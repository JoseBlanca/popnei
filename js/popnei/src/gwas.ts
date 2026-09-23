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
 * What `calcGwas` fits is the linear model, a continuous trait with no
 * kinship, which is what plink2's `--glm` computes, and its test is the t
 * test of the effect. The linear mixed model, which takes a kinship so that
 * a variant that only marks the ancestry of a panel does not look
 * associated, and the two logistic models of a binomial trait are being
 * written; asking for one is an `Error` that says so. `docs/specs/gwas.md`
 * has the four of them.
 */

import {
  default_transform_to_biallelic as defaultTransformToBiallelic,
} from "../wasm/popnei.js";

import { aBoolean, aString, whatWasGiven } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
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
 * The options of `calcGwas` that the linear mixed model brings and that this
 * build has not: a kinship, which model with it, and the approximation that
 * a mixed model can take instead of a fit per variant.
 *
 * They are refused by name rather than ignored: a call that gave a kinship
 * and got a study without one would be a linear model reported as a mixed
 * one, with nothing to show it.
 */
const OF_THE_MIXED_MODEL = ["kinship", "test", "useGrammarGammaApprox"];

/** What was measured on each individual. */
export type TraitType = "continuous" | "binomial";

/** Which test is made of every variant. */
export type GwasTestType = "wald" | "score";

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
  readonly test: GwasTestType;
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
   * that are tested are those that have a number here and that the
   * `Variants` has. An individual whose value is NaN, `null` or `undefined`
   * has no phenotype and is not tested.
   */
  phenotype: Readonly<Record<string, number | null | undefined>>;
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
   */
  covariates?: Readonly<Record<string, Readonly<Record<string, number>>>>;
  /**
   * Whether every allele that is not the major one counts the same, which is
   * what gives a variant of more than two alleles a dosage. False when it is
   * not given, and such a variant is then an `Error`.
   */
  transformToBiallelic?: boolean;
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
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * `phenotype` is not an object of a name to a number; when a name of it is
 * of nobody the pass gives; when `trait` is not one of the two names; when a
 * covariate is not an object of a name to a number, does not cover a tested
 * individual, or holds a value that is missing or is not a number; when a
 * covariate is named `intercept`, which is the name the effect of the column
 * of ones comes back under; when `kinship`, `test` or `useGrammarGammaApprox`
 * is given, which the linear mixed model brings; when no individual is
 * tested or they are fewer than the columns of the design plus two; when a
 * phenotype or a covariate is not a finite number; when the columns of the
 * design are not independent; when the trait is binomial, which is a
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
  theOptionsOfTheMixedModel(options);
  const trait = aString("trait", options.trait);
  const transformToBiallelic =
    options.transformToBiallelic === undefined
      ? defaultTransformToBiallelic()
      : aBoolean("transformToBiallelic", options.transformToBiallelic);
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
      test: calculated.test() as GwasTestType,
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
 * @throws {Error} When `phenotype` is not an object of a name to a number,
 * when a name of it is of nobody the pass gives, and when a value of it is
 * neither a number nor missing.
 */
function theTestedIndividuals(
  phenotype: Readonly<Record<string, number | null | undefined>>,
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
    const value = phenotype[name];
    // A missing phenotype is an individual that is not tested, which is what
    // a name that is not in the object is: `null` is what a phenotype read
    // from JSON holds where a NaN was written.
    if (value === undefined || value === null) {
      continue;
    }
    if (typeof value !== "number") {
      throw new Error(
        `popnei: the phenotype of \`${name}\` is ` +
          `${whatWasGiven(value)}, and a trait is a number`,
      );
    }
    if (Number.isNaN(value)) {
      continue;
    }
    tested.push({ name, position, phenotype: value });
  }
  return tested;
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
 * of a name to a number, when one of them is named `intercept`, when one
 * does not cover a tested individual, and when a value of one is missing or
 * is not a number.
 */
function theCovariates(
  covariates:
    | Readonly<Record<string, Readonly<Record<string, number>>>>
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
        "popnei: a covariate is not named `intercept`, which is the name the " +
          "effect of the column of ones every design has comes back under in " +
          "`nullModel.covariateEffects`",
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
 * missing or is not a number. A covariate whose values are names is refused
 * here, and what a user does with one is to give one covariate for each of
 * its values, 1 for the individuals of that value and 0 for the others.
 */
function theValuesOfTheCovariate(
  name: string,
  values: Readonly<Record<string, number>>,
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
    if (typeof value !== "number" || Number.isNaN(value)) {
      throw new Error(
        `popnei: the value of the covariate \`${name}\` at \`${individual}\` ` +
          `is ${whatWasGiven(value)}, and a covariate is a number: one ` +
          "whose values are names is given as one covariate for each of " +
          "them, 1 for the individuals of that value and 0 for the others",
      );
    }
    return value;
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
 * Refuses the options that the linear mixed model brings, which this build
 * has not.
 *
 * @throws {Error} When `kinship`, `test` or `useGrammarGammaApprox` is
 * given.
 */
function theOptionsOfTheMixedModel(options: CalcGwasOptions): void {
  const given = Object.keys(options).filter((option) =>
    OF_THE_MIXED_MODEL.includes(option),
  );
  if (given.length > 0) {
    throw new Error(
      `popnei: \`${given.join("`, `")}\` belongs to the linear mixed model, ` +
        "which accounts for the relatedness of a panel and is being written; " +
        "what calcGwas fits is the linear model, a continuous trait with no " +
        "kinship",
    );
  }
}
