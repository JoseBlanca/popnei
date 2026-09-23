/**
 * The association study: which of the variants of a dataset are associated
 * with a trait that was measured on its individuals.
 *
 * A trait is one number per individual, a measurement such as the height of
 * a plant. Covariates are other numbers per individual whose effect on the
 * trait is not of interest but has to be taken out, such as the field the
 * plant grew in. The model is fitted once with no variant in it, which is
 * the null model, and then every variant is tested against what that model
 * left unexplained: `beta` is the effect of one more copy of a non major
 * allele, in the units of the trait, `se` is how uncertain that effect is,
 * and `pValue` is the probability of seeing an effect that far from 0 when
 * the variant has none.
 *
 * What `calcGwas` fits is the continuous half of the study, the two models
 * of a trait that is a measurement. Without a kinship it is the linear
 * model, which is what plink2's `--glm` computes, and its test is the t test
 * of the effect: the effect divided by its standard error, which under the
 * hypothesis that the variant has none follows a Student t distribution
 * with as many degrees of freedom as there are individuals left once the
 * covariates and the variant have been fitted, and `pValue` is the chance
 * that such a t falls further from 0 than this one did, either way.
 *
 * With a kinship it is the linear mixed model, which a panel with families
 * in it needs: the trait carries a random effect whose covariance is the
 * kinship times a variance, so that a variant which only marks the ancestry
 * of the panel does not look associated. Its two tests are rrBLUP's Wald
 * test and GMMAT's score test.
 *
 * The two logistic models of a binomial trait are being written, and so is
 * the GRAMMAR-Gamma approximation a mixed model can take instead of the
 * exact denominator of its test; asking for one is an `Error` that says so.
 * `docs/specs/gwas.md` has the four models.
 */

import {
  default_transform_to_biallelic as defaultTransformToBiallelic,
} from "../wasm/popnei.js";

import { aBoolean, aString, whatWasGiven } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Kinship } from "./kinship.js";
import { theValuesOf } from "./pca.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * The name of the effect of the column of ones that every design has, which
 * is what the trait is worth when every covariate is 0.
 *
 * It is the key of `covariateEffects` that the intercept comes back under,
 * and it is pyNei's name for the same row of its series.
 */
const INTERCEPT = "intercept";

/**
 * What `float` takes in Python and `Number` does not, and the other way
 * round: a string is read as the number it holds only when it is written the
 * way `float` writes one.
 *
 * `Number` gives a value where `float` raises, so a string that holds no
 * number would become a phenotype of 0 or a silently missing individual, and
 * `float` takes what `Number` refuses, the underscores of `1_000`. What is
 * matched is an optional sign and then digits with an optional point and an
 * optional exponent, with `_` allowed between digits, which is `float`'s
 * grammar for a decimal literal. The words `nan`, `inf` and `infinity`,
 * which `float` also takes, are matched apart, since `Number` reads only the
 * second of the three.
 */
const A_DIGIT_RUN = "\\d+(?:_\\d+)*";
const WRITTEN_AS_FLOAT_WRITES_ONE = new RegExp(
  `^[+-]?(?:${A_DIGIT_RUN}(?:\\.(?:${A_DIGIT_RUN})?)?|\\.${A_DIGIT_RUN})` +
    `(?:[eE][+-]?${A_DIGIT_RUN})?$`,
);
const WRITTEN_AS_A_WORD_FLOAT_TAKES = /^([+-]?)(nan|inf|infinity)$/i;

/** What was measured on each individual. */
export type TraitType = "continuous" | "binomial";

/**
 * Which test is made of every variant: `wald` fits the model again with the
 * variant in it and measures its effect in its own standard errors away
 * from 0, and `score` never fits it again and measures at the null model
 * how steeply the fit would improve if the variant's effect were let off 0.
 *
 * Under the hypothesis that the variant has no effect the two have the same
 * distribution in large samples; they differ in what they cost. The linear
 * model has the Wald test alone, which for it is the t test of the effect,
 * and asking it for the score test is an `Error`.
 */
export type TestType = "wald" | "score";

/**
 * Which of the four models a study fits, which the trait and the kinship
 * together decide: `lm` a linear model, `lmm` a linear mixed model, `glm` a
 * logistic regression and `glmm` a logistic mixed model.
 */
export type GwasModel = "lm" | "lmm" | "glm" | "glmm";

/**
 * The answers of a study, one row for each variant it was given, each column
 * as its own array.
 *
 * The rows are in the order the variants came. A variant whose dosages are
 * all the same among the tested individuals has no variance and cannot be
 * tested: its row is here with its `alleleFreq`, and its `beta`, its `se`
 * and its `pValue` are NaN.
 */
export interface GwasStats {
  /**
   * The name of the chromosome of each variant, and `undefined` when the
   * source has no such column.
   */
  readonly chrom?: readonly string[];
  /**
   * The position of each variant, 1 based as in a VCF, and `undefined` when
   * the source has no such column.
   */
  readonly pos?: Float64Array;
  /**
   * The id of each variant, and `undefined` when the source has no such
   * column.
   */
  readonly id?: readonly string[];
  /**
   * The frequency of the alleles that are not the major one, over the tested
   * individuals: the mean dosage over the ploidy, which for a biallelic
   * variant is the minor allele frequency and which is what plink2 reports
   * as `A1_FREQ`.
   */
  readonly alleleFreq: Float64Array;
  /**
   * The effect of one more copy of a non major allele, in the units of the
   * trait, and NaN for a variant that has no answer.
   */
  readonly beta: Float64Array;
  /** The standard error of that effect, and NaN where `beta` is NaN. */
  readonly se: Float64Array;
  /**
   * The p-value of the test that the effect is 0, and NaN where `beta` is
   * NaN.
   */
  readonly pValue: Float64Array;
}

/**
 * The model a study fitted with no variant in it, which every variant was
 * then tested against.
 */
export interface GwasNullModel {
  /** Which of the four models it is. */
  readonly model: GwasModel;
  /**
   * The effect of the intercept, under the name `intercept`, and of every
   * covariate, under the name it was given.
   */
  readonly covariateEffects: Readonly<Record<string, number>>;
  /**
   * What the model left unexplained, and `undefined` for a binomial trait,
   * whose variance is decided by its mean.
   */
  readonly residualVariance?: number;
  /**
   * The variance of the random effect of the kinship, and `undefined`
   * without a kinship.
   */
  readonly geneticVariance?: number;
  /**
   * The genetic variance over the sum of the two, only for the linear mixed
   * model.
   */
  readonly heritability?: number;
  /** How many individuals the study tested. */
  readonly numIndividuals: number;
}

/** What an association study gives back. */
export interface GwasResult {
  /** One row for each variant the study was given, column by column. */
  readonly stats: GwasStats;
  /** The model fitted with no variant in it. */
  readonly nullModel: GwasNullModel;
  /** What was measured. */
  readonly trait: TraitType;
  /** Which test was made of every variant. */
  readonly test: TestType;
  /**
   * The names of the individuals that were tested, in the order the source
   * has them, which is the order their phenotype and their design were read
   * in.
   */
  readonly individuals: readonly string[];
  /** Whether the GRAMMAR-Gamma approximation was used. */
  readonly usedGrammarGammaApprox: boolean;
  /**
   * How many variants the pass gave, used or not, and what each filter of
   * the `Variants` was given and kept.
   */
  readonly passStats: PassStats;
}

/** What a study is asked for. */
export interface CalcGwasOptions {
  /**
   * What was measured on each individual, under its name: the individuals
   * that are tested are those that have a value here and that the
   * `Variants` has.
   *
   * A value that is not a number is read as one, as `float` reads it in the
   * Python package and in pyNei: the string `"1.7"`, which is how a trait
   * read from a file arrives, and `true` and `false`. A value whose `float`
   * would be NaN is an individual with no phenotype, which is left untested:
   * NaN itself, the string `"nan"`, and a key that is not in the object at
   * all. What holds no number is an `Error` naming the individual, `null`,
   * `undefined`, the empty string and any other string that is not a number
   * among them, and so is an infinity.
   */
  phenotype: Readonly<Record<string, number | string | boolean>>;
  /**
   * What the trait is: `"continuous"`, a measurement, or `"binomial"`, 0 for
   * an individual that has not a condition and 1 for one that has. A
   * binomial trait is a logistic regression, which is being written.
   */
  trait: TraitType;
  /**
   * The other numbers per individual whose effect on the trait is taken out,
   * each of them under its own name and holding one value for each tested
   * individual.
   *
   * The top principal components of the panel, from
   * `kinship.principalComponents` or `doPcaFromVariants`, go in here, which
   * is how the ancestry of individuals that are not close relatives is
   * accounted for. A covariate whose values are names and not numbers, the
   * field a plant grew in, is given as one covariate for each of its values,
   * 1 for the individuals of that value and 0 for the others.
   *
   * A value that is not a number is read as one here as it is in the
   * phenotype, the string `"1.7"` and `true` and `false` among them, and
   * one that holds no number is an `Error` naming the covariate and the
   * individual.
   */
  covariates?: Readonly<
    Record<string, Readonly<Record<string, number | string | boolean>>>
  >;
  /**
   * Which test is made of every variant, and the default of the model when
   * it is not given. What `calcGwas` fits is the linear model, whose only
   * test is `wald`, so `score` is an `Error` that says so.
   */
  test?: TestType;
  /**
   * Whether every allele that is not the major one counts the same, which is
   * what gives a variant of more than two alleles a dosage. False when it is
   * not given, and such a variant is then an `Error`.
   */
  transformToBiallelic?: boolean;
  /**
   * The relatedness of every pair of individuals, which the linear mixed
   * model takes as the covariance of a random effect, so that a variant
   * which only marks the ancestry of a panel does not look associated.
   *
   * It is what `calcKinship` gives, or a `Kinship` built over the matrix
   * another program wrote. It has to hold every individual that is tested,
   * and a tested individual it has not is an `Error` naming them; the ones
   * it holds over are left out, as `Kinship` of some of a panel would be.
   * Without it the structure of a panel goes in as the top principal
   * components among the covariates, which is enough for individuals that
   * are not close relatives.
   */
  kinship?: Kinship;
  /**
   * Whether the GRAMMAR-Gamma approximation is made, which stands in for the
   * denominator of a mixed model's test and which only a mixed model has.
   * It is being written, and asking for it is an `Error`: with no kinship
   * because there is no such denominator to approximate, and with one
   * because popnei cannot approximate it yet.
   */
  useGrammarGammaApprox?: boolean;
}

/**
 * Which of the variants of `variants`, after its steps, are associated with
 * the trait of `phenotype`.
 *
 * Each variant becomes one number per individual, its dosage: how many
 * alleles of the genotype are not the major allele of that variant, which is
 * the most frequent among its called alleles. A genotype with any allele
 * missing takes the mean of the dosages of its variant. The dosages, the
 * mean and the frequency are over the tested individuals and not over the
 * whole panel, so a phenotype that leaves individuals out gives a variant
 * another frequency than the panel's.
 *
 * The individuals that are tested are those that have a phenotype and that
 * the `Variants` has, in the order the `Variants` has them, whatever order
 * the phenotype was written in. The design is the column of ones of the
 * intercept and one column for each covariate, and its columns have to be
 * independent: a covariate that is constant or a copy of another is an
 * `Error`.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, and nothing of
 * the size of the variants x the individuals is held.
 *
 * It is pyNei's `calc_gwas`, whose `samples` is `individuals` here.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when the
 * options are not given at all; when `phenotype` is not an object of a name
 * to a value that holds a number, which `null`, `undefined` and the empty
 * string do not; when a name of it is of nobody the pass gives; when `trait`
 * is not one of the two names and when `test` is not one of the two; when a
 * covariate is not an object of a name to such a value, does not cover a
 * tested individual, or holds a value that is missing or holds no number;
 * when a covariate is named `intercept`, which is the name the effect of the
 * column of ones comes back under; when `kinship` is not a `Kinship` or has
 * not an individual that is tested; when `useGrammarGammaApprox` is asked
 * for, which is being written; when the score test is asked of a linear
 * model, which has it not; when no individual is tested or they are fewer
 * than the columns of the design plus two; when a phenotype or a covariate
 * is not a finite number once it is read as one; when the trait is the same
 * in every tested individual, which leaves nothing for a variant to be
 * associated with; when the columns of the design are not independent; when the trait is binomial, which is a
 * logistic model and is being written; when the source cannot be read, a
 * wrong line of a VCF among the causes; when a variant has more than two
 * alleles among its called genotypes and `transformToBiallelic` is false;
 * when the pass gives no variant; when the linear algebra of the fit or of a
 * test could not be done; and when `init` has not been awaited.
 */
export function calcGwas(
  variants: Variants,
  options: CalcGwasOptions,
): GwasResult {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  theOptions(options);
  const trait = aString("trait", options.trait);
  const test =
    options.test === undefined ? undefined : aString("test", options.test);
  const transformToBiallelic =
    options.transformToBiallelic === undefined
      ? defaultTransformToBiallelic()
      : aBoolean("transformToBiallelic", options.transformToBiallelic);
  // An option written and left `undefined` is one that was not given, which
  // is what spreading an object of options over a call leaves behind.
  const useGrammarGammaApprox =
    options.useGrammarGammaApprox === undefined
      ? false
      : aBoolean("useGrammarGammaApprox", options.useGrammarGammaApprox);
  const kinship = theKinship(options.kinship);
  const individuals = theTestedIndividuals(
    options.phenotype,
    variants.individuals,
  );
  const phenotype = Float64Array.from(
    individuals.map(({ phenotype: value }) => value),
  );
  const names = individuals.map(({ name }) => name);
  const covariates = theCovariates(options.covariates, names);
  const design = theDesign(covariates, names);
  const positions = Uint32Array.from(
    individuals.map(({ position }) => position),
  );
  const numCoefs = covariates.length + 1;
  // The steps of the pass are a copy of the list, made after the arguments
  // were checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = source.calc_gwas(
    positions,
    phenotype,
    design,
    numCoefs,
    trait,
    test,
    kinship === undefined ? undefined : theKinshipOfTheTested(kinship, names),
    useGrammarGammaApprox,
    transformToBiallelic,
    steps.of_a_pass(),
  );
  try {
    // Every array is moved out of the result as it is read, and not cloned:
    // the generated code copies the values into an array of the JavaScript
    // heap and frees the memory of wasm after it.
    const effects = theValuesOf(
      calculated.covariate_effects(),
      "covariateEffects",
    );
    const stats: GwasStats = {
      chrom: calculated.chroms(),
      pos: calculated.poss(),
      id: calculated.ids(),
      alleleFreq: theValuesOf(calculated.allele_freq(), "alleleFreq"),
      beta: theValuesOf(calculated.beta(), "beta"),
      se: theValuesOf(calculated.se(), "se"),
      pValue: theValuesOf(calculated.p_value(), "pValue"),
    };
    return {
      stats,
      nullModel: {
        model: calculated.model() as GwasModel,
        covariateEffects: theEffects(effects, covariates),
        residualVariance: calculated.residual_variance(),
        geneticVariance: calculated.genetic_variance(),
        heritability: calculated.heritability(),
        numIndividuals: calculated.num_individuals(),
      },
      trait: trait as TraitType,
      test: calculated.test() as TestType,
      individuals: Object.freeze(names),
      usedGrammarGammaApprox: calculated.used_grammar_gamma_approx(),
      passStats: passStatsOf(calculated.pass_stats()),
    };
  } finally {
    calculated.free();
  }
}

/** One tested individual: its name, where it is in the pass and its trait. */
interface TestedIndividual {
  name: string;
  position: number;
  phenotype: number;
}

/**
 * The individuals that are tested, in the order `ofThePass` has them: those
 * that have a phenotype there.
 *
 * The order is the source's and not the phenotype's, because the phenotype,
 * the rows of the design and the dosages of a variant are read together row
 * by row: a study that took them in the order the phenotype was written in
 * would measure one individual's trait against another's genotypes.
 *
 * @throws {Error} When `phenotype` is not an object of a name to a value
 * that holds a number, when a name of it is of nobody the pass gives, and
 * when a value of it holds no number.
 */
function theTestedIndividuals(
  phenotype: Readonly<Record<string, number | string | boolean>>,
  ofThePass: readonly string[],
): TestedIndividual[] {
  if (
    typeof phenotype !== "object" ||
    phenotype === null ||
    Array.isArray(phenotype)
  ) {
    throw new Error(
      "popnei: `phenotype` is an object of the name of an individual to what " +
        `was measured on it, {ind00: 1.7}, and ${whatWasGiven(phenotype)} ` +
        "was given",
    );
  }
  const ofTheSource = new Set(ofThePass);
  for (const name of Object.keys(phenotype)) {
    if (!ofTheSource.has(name)) {
      throw new Error(
        `popnei: \`${name}\` has a phenotype and is not one of the ` +
          `${ofThePass.length} individuals these variants give`,
      );
    }
  }
  const tested: TestedIndividual[] = [];
  for (const [position, name] of ofThePass.entries()) {
    // An individual the object has no key for has no phenotype and is not
    // tested, which is what a name the series has not says in Python.
    if (!Object.hasOwn(phenotype, name)) {
      continue;
    }
    const value: unknown = phenotype[name];
    if (value === null || value === undefined) {
      throw new Error(
        `popnei: the phenotype of \`${name}\` is ${whatWasGiven(value)}, and ` +
          "a trait is a number: an individual with no phenotype is left out " +
          "of `phenotype` altogether, which is what leaves it untested",
      );
    }
    const number = theNumberOf(value);
    if (number === undefined) {
      throw new Error(
        `popnei: the phenotype of \`${name}\` is ${whatWasGiven(value)}, and ` +
          "a trait is a number: a trait whose values are names is not a " +
          "trait of a study, and one that is 0 and 1 is written as those " +
          'numbers with trait: "binomial"',
      );
    }
    // A value whose `float` is NaN is an individual with no phenotype, and
    // is left untested rather than refused: that is the rule the spec
    // settled on 23 September 2026, by the oracle, and it is what the
    // Python package does with the NaN of a table read from a file. NaN
    // itself and the string `nan` are what arrive here.
    if (Number.isNaN(number)) {
      continue;
    }
    // An infinity is refused here, where the individual has a name: it
    // would carry through the null model into the effect of every variant.
    // The core refuses it as well, by the place of the individual, which is
    // what a caller of the wasm module reads. pyNei accepts it and then
    // gives NaN for every variant.
    if (!Number.isFinite(number)) {
      throw new Error(
        `popnei: the phenotype of \`${name}\` is ${number}, and a study is ` +
          "fitted on numbers: an individual with no phenotype is left out " +
          "of `phenotype` altogether, which is what leaves it untested",
      );
    }
    tested.push({ name, position, phenotype: number });
  }
  return tested;
}

/**
 * `value` as `float` reads it in Python, and `undefined` where `float`
 * raises.
 *
 * `float` is the rule and not `Number`, because `float` is what pyNei and
 * the Python package of popnei read a phenotype with and the two layers
 * answer the same thing or they are two libraries: a value means an
 * individual with no phenotype exactly where this gives NaN, and it is
 * refused exactly where this gives `undefined`. `Number` agrees with `float`
 * nowhere that matters, which is why nothing here is handed to it but a
 * string already matched: `Number("abc")` gives NaN where `float` raises, so
 * a typo would vanish instead of being reported; `Number(null)` and
 * `Number("")` give 0, so a blank cell would become a phenotype of zero; and
 * `Number("0x10")` gives 16 where `float` raises.
 *
 * A number, a whole number and a boolean are read as `float` reads them, and
 * a string is read when it is written the way `float` writes a number, which
 * is how a trait read from a file arrives. `"nan"` comes back as NaN, so it
 * is an individual with no phenotype and not a refusal, which is the case
 * nobody guesses.
 */
function theNumberOf(value: unknown): number | undefined {
  if (
    typeof value === "number" ||
    typeof value === "boolean" ||
    typeof value === "bigint"
  ) {
    return Number(value);
  }
  if (typeof value !== "string") {
    return undefined;
  }
  // `float` skips the whitespace around the number and nothing else.
  const written = value.trim();
  const word = WRITTEN_AS_A_WORD_FLOAT_TAKES.exec(written);
  if (word !== null) {
    const sign = word[1] === "-" ? -1 : 1;
    return (word[2] as string).toLowerCase() === "nan"
      ? Number.NaN
      : sign * Number.POSITIVE_INFINITY;
  }
  if (!WRITTEN_AS_FLOAT_WRITES_ONE.test(written)) {
    return undefined;
  }
  return Number(written.replaceAll("_", ""));
}

/** One covariate: its name and its value for each tested individual. */
interface Covariate {
  name: string;
  values: number[];
}

/**
 * The covariates, each with one value for each of the individuals of
 * `tested`, in their order.
 *
 * @throws {Error} When `covariates` is not an object of a name to an object
 * of a name to a value that holds a number, when one of them is named
 * `intercept`, when one does not cover a tested individual, and when a value
 * of one is missing or holds no number.
 */
function theCovariates(
  covariates:
    | Readonly<
        Record<string, Readonly<Record<string, number | string | boolean>>>
      >
    | undefined,
  tested: readonly string[],
): Covariate[] {
  if (covariates === undefined) {
    return [];
  }
  if (
    typeof covariates !== "object" ||
    covariates === null ||
    Array.isArray(covariates)
  ) {
    throw new Error(
      "popnei: `covariates` is an object of the name of a covariate to its " +
        `value for each individual, {cov1: {ind00: 1.7}}, and ` +
        `${whatWasGiven(covariates)} was given`,
    );
  }
  const asked = [];
  for (const [name, values] of Object.entries(covariates)) {
    if (name === INTERCEPT) {
      throw new Error(
        "popnei: a covariate is named `intercept`, which is the name the " +
          "effect of the column of ones every design has comes back under in " +
          "`nullModel.covariateEffects`: the two would be one entry of it, " +
          "so give the covariate another name",
      );
    }
    if (
      typeof values !== "object" ||
      values === null ||
      Array.isArray(values)
    ) {
      throw new Error(
        `popnei: the covariate \`${name}\` is an object of the name of an ` +
          "individual to its value, {ind00: 1.7}, and " +
          `${whatWasGiven(values)} was given`,
      );
    }
    asked.push({ name, values: theValuesOfTheCovariate(name, values, tested) });
  }
  return asked;
}

/**
 * The value of the covariate `name` for each of the individuals of `tested`,
 * in their order.
 *
 * @throws {Error} When an individual has no value, and when a value is
 * missing or holds no number. A value that is not a number is read as one
 * first, as the phenotype's is. A covariate whose values are names is
 * refused here, and what a user does with one is to give one covariate for
 * each of its values, 1 for the individuals of that value and 0 for the
 * others.
 */
function theValuesOfTheCovariate(
  name: string,
  values: Readonly<Record<string, number | string | boolean>>,
  tested: readonly string[],
): number[] {
  return tested.map((individual) => {
    const value: unknown = values[individual];
    if (value === undefined || value === null) {
      throw new Error(
        `popnei: the covariate \`${name}\` has no value for ` +
          `\`${individual}\`, and a covariate holds one for every individual ` +
          "that is tested",
      );
    }
    const number = theNumberOf(value);
    if (number === undefined) {
      throw new Error(
        `popnei: the value of the covariate \`${name}\` at \`${individual}\` ` +
          `is ${whatWasGiven(value)}, and a covariate is a number: one ` +
          "whose values are names is given as one covariate for each of " +
          "them, 1 for the individuals of that value and 0 for the others",
      );
    }
    // A NaN and an infinity are numbers to JavaScript and not to a fit, and
    // this is the layer that has the name of the covariate and of the
    // individual: the core refuses them as well, by their places among the
    // columns and the rows, which is what a caller of the core crate reads.
    if (!Number.isFinite(number)) {
      throw new Error(
        `popnei: the value of the covariate \`${name}\` at \`${individual}\` ` +
          `is ${number}, and a study is fitted on numbers: it would carry ` +
          "through the null model into the effect of every variant",
      );
    }
    return number;
  });
}

/**
 * The design of the study: one row for each individual of `tested`, the 1 of
 * the intercept and then the value of each covariate, row after row.
 */
function theDesign(
  covariates: readonly Covariate[],
  tested: readonly string[],
): Float64Array {
  const numCoefs = covariates.length + 1;
  const design = new Float64Array(tested.length * numCoefs);
  for (let individual = 0; individual < tested.length; individual += 1) {
    const row = individual * numCoefs;
    design[row] = 1;
    for (const [covariate, { values }] of covariates.entries()) {
      design[row + covariate + 1] = values[individual] as number;
    }
  }
  return design;
}

/**
 * The effect of each column of the design under its name: the intercept
 * first, which is the column of ones, and then the covariates in the order
 * they were given.
 *
 * @throws {Error} When the model answered another number of effects than the
 * design has columns, which is a defect of popnei.
 */
function theEffects(
  effects: Float64Array,
  covariates: readonly Covariate[],
): Readonly<Record<string, number>> {
  const names = [INTERCEPT, ...covariates.map(({ name }) => name)];
  if (effects.length !== names.length) {
    throw new Error(
      `popnei: the null model gave ${effects.length} effects and its design ` +
        `has ${names.length} columns, which is a defect of popnei; please ` +
        "report it",
    );
  }
  const named: Record<string, number> = {};
  for (const [column, name] of names.entries()) {
    named[name] = effects[column] as number;
  }
  return Object.freeze(named);
}

/**
 * Refuses a call with no options at all.
 *
 * `calcGwas(variants)` is what a user writes who has read the signature of
 * `calcKinship`, and what it gave was the `TypeError` of `Object.keys` on
 * `undefined`, which names neither popnei nor what is missing.
 *
 * @throws {Error} When `options` is not an object.
 */
function theOptions(options: CalcGwasOptions): void {
  if (
    typeof options !== "object" ||
    options === null ||
    Array.isArray(options)
  ) {
    throw new Error(
      "popnei: a study is asked for with the trait of the individuals and " +
        "what was measured, calcGwas(variants, {phenotype, trait: " +
        `"continuous"}), and ${whatWasGiven(options)} was given`,
    );
  }
}

/**
 * What was written in `kinship`, given back unless it is neither a `Kinship`
 * nor missing.
 *
 * An option written and left `undefined` is one that was not given, which is
 * what spreading an object of options over a call leaves behind, so a study
 * with `kinship: undefined` is a study with no kinship and no error.
 *
 * @throws {Error} When it holds something that is not a `Kinship`. What a
 * user gives instead is usually the matrix itself, and what that gave was
 * the error of a `Float64Array` read off `undefined`.
 */
function theKinship(kinship: Kinship | undefined): Kinship | undefined {
  if (kinship === undefined || kinship instanceof Kinship) {
    return kinship;
  }
  throw new Error(
    `popnei: \`kinship\` is ${whatWasGiven(kinship)}, and the kinship a ` +
      "mixed model takes is a `Kinship`: give it what `calcKinship` gives, " +
      "or build one over the matrix another program wrote, new " +
      "Kinship(matrix, individuals, numVars)",
  );
}

/**
 * The kinship cut to the individuals of `tested`, in their order, as the
 * core reads it: one row and one column for each of them, row after row.
 *
 * The core is given numbers and holds no name to cut a matrix by, so the
 * cutting is here, where the names are. An individual the kinship holds over
 * is left out, as a `Kinship` of some of a panel leaves it out: the kinship
 * a user has is of their panel, and a phenotype that leaves individuals out
 * does not make it another matrix.
 *
 * @throws {Error} When a tested individual is not one of the kinship's: the
 * random effect of a mixed model is the relatedness of every pair that is
 * tested, and there is nothing to put in the row of an individual the matrix
 * has not.
 */
function theKinshipOfTheTested(
  kinship: Kinship,
  tested: readonly string[],
): Float64Array {
  const ofTheMatrix = new Map(
    kinship.individuals.map((name, row) => [name, row]),
  );
  const rows = tested.map((name) => {
    const row = ofTheMatrix.get(name);
    if (row === undefined) {
      throw new Error(
        `popnei: \`${name}\` is tested and is not one of the ` +
          `${kinship.individuals.length} individuals of the \`kinship\`, ` +
          "which holds the relatedness of every pair that is tested: give a " +
          "kinship of them, `calcKinship(variants)` for instance, or leave " +
          "that individual out of the phenotype",
      );
    }
    return row;
  });
  const ofThePanel = kinship.individuals.length;
  const cut = new Float64Array(rows.length * rows.length);
  for (const [row, ofTheRow] of rows.entries()) {
    for (const [column, ofTheColumn] of rows.entries()) {
      cut[row * rows.length + column] = kinship.matrix[
        ofTheRow * ofThePanel + ofTheColumn
      ] as number;
    }
  }
  return cut;
}
