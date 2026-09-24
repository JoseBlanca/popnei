/**
 * How much variety each population of a dataset holds.
 *
 * `calcPopDiversity` makes one pass over the variants and gives, for each
 * population, how many alleles it called, how many of those no other
 * population called, how many of the variants vary in it, how its variants
 * are spread over the count of their rarer allele, and how far its genotypes
 * are from the proportions its allele frequencies would give. The first
 * three come both as they stand and standardized to a common number of
 * called alleles, so that a population of 20 individuals and one of 200 can
 * be compared.
 *
 * Standardized means drawn: take `numCalledAlleles` of the called alleles a
 * population has at a variant, without replacement, and take the expectation
 * over every such draw. Applied to a count of alleles that operation is
 * called rarefaction and applied to the spectrum it is called projection,
 * and both come from that one argument.
 *
 * A population is a named set of individuals, and `pops` names them as
 * `calcPerVarDistribs` does: an object of population name to the names of
 * its individuals, with no `pops` meaning one population, `pop`, of every
 * individual.
 *
 * `docs/specs/diversity.md` has what each of the five statistics is, the
 * program each of its numbers was checked against and the numbers the tests
 * assert.
 */

import {
  default_min_num_individuals as defaultMinNumIndividuals,
  diversity_stats_without_a_draw as statsWithoutADraw,
} from "../wasm/popnei.js";

import {
  LARGEST_WHOLE_NUMBER,
  namesOf,
  popsOfTheObject,
  whatWasGiven,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * The five statistics `calcPopDiversity` calculates, each named as the field
 * of the result that holds it is named in Python.
 */
const THE_STATISTICS = [
  "num_alleles",
  "private_alleles",
  "variable_vars_ratio",
  "folded_sfs",
  "fis",
] as const;

/**
 * The name of one of the five statistics of a population: how many alleles
 * it called, how many of those no other population of the call called at the
 * same variant, how many of the variants it called more than one allele at,
 * how its variants are spread over the count of their rarer allele in a
 * draw, and F_IS, how far its genotypes are from the proportions its allele
 * frequencies would give if its individuals paired at random.
 */
export type PopDiversityStat = (typeof THE_STATISTICS)[number];

/** What `calcPopDiversity` calculates, for which populations and how. */
export interface CalcPopDiversityOptions {
  /**
   * Which of the five statistics to calculate. The default is the four that
   * need no draw: the alleles a population called, the private ones among
   * them, the variants that vary in it and F_IS. Asking for fewer is a
   * saving of work and changes no value, and a result holds `null` for one
   * nobody asked for.
   *
   * `folded_sfs` is the fifth and is not in the default, its bins being
   * counts of the rarer allele in a draw of `numCalledAlleles`: a user who
   * wants it names it here and gives that draw, and asking for it without one
   * is an `Error`.
   */
  stats?: readonly PopDiversityStat[];

  /**
   * The populations: an object of population name to the names of its
   * individuals, which are looked up among the individuals the pass gives.
   * With no `pops` there is one population, `pop`, of every individual, and
   * then every allele it called is private.
   *
   * The populations of every result are in the order the keys of this object
   * iterate in. An individual may be in more than one population and is
   * counted in each of them.
   */
  pops?: Record<string, readonly string[]>;

  /**
   * How many called alleles every population is brought down to, so that a
   * population of 20 individuals and one of 200 can be compared: draw that
   * many of the alleles a population called at a variant, without
   * replacement, and take the expectation over every such draw.
   *
   * A number that is no whole number of 2 or more is refused as the wrong
   * argument it is, a draw of one allele finding one allele whatever the
   * population holds, and so is a draw of more alleles than the dataset holds
   * gene copies, its individuals times their ploidy, which no variant of any
   * population could reach. A draw at or below that number which this
   * dataset's missing genotypes leave no population able to fill is no error:
   * every `inDraw` value is then NaN, every bin of the spectrum is 0 and
   * `numVars.inDraw` says why.
   *
   * With no `numCalledAlleles`, which is the default, every `inDraw` value is
   * NaN, `numVars.inDraw` is 0 and the spectrum cannot be asked for.
   *
   * It is one number and not a sequence of them, so a user who wants the
   * curve of allelic richness against the number of alleles drawn, which
   * shows whether a population has been sampled enough, calls this once per
   * point of it.
   */
  numCalledAlleles?: number;

  /**
   * How many called genotypes a population needs at a variant for the
   * variant to count for it, 20 when it is not given, measured as the called
   * alleles of the population over the ploidy, so a half called genotype
   * counts as half an individual. Strictly fewer and the variant does not
   * count for that population.
   *
   * A population that called nothing at a variant does not count it whatever
   * `minNumIndividuals` is, so a 0 does not put variants with no data into
   * the totals.
   */
  minNumIndividuals?: number;
}

/**
 * A count of alleles of every population: as it stands, as a mean over the
 * variants, and standardized to a draw of `numCalledAlleles`.
 *
 * Every array holds one value per population, in the order of `pops` of the
 * result.
 */
export interface PopAlleleCounts {
  /** The count summed over the variants that it is counted over. */
  total: Uint32Array;

  /**
   * That total over those variants, and NaN for a population no variant
   * counted for.
   */
  mean: Float64Array;

  /**
   * What a draw of `numCalledAlleles` is expected to show, averaged over the
   * variants in the draw for the population, and NaN where there was no draw
   * or no variant reached it.
   */
  inDraw: Float64Array;
}

/**
 * How many of the variants each population called more than one allele at:
 * that count, that count over the variants that counted for the population,
 * and the same in a draw.
 *
 * Every array holds one value per population, in the order of `pops` of the
 * result.
 */
export interface VariableVarsRatio {
  /** The variants the population called more than one allele at. */
  total: Uint32Array;

  /**
   * That count over `numVars.withData`, and NaN for a population no variant
   * counted for.
   */
  ratio: Float64Array;

  /**
   * The chance that a draw of `numCalledAlleles` is not all of one allele,
   * averaged over the variants in the draw, and NaN where there was no draw or
   * no variant reached it.
   */
  inDraw: Float64Array;
}

/**
 * How many variants counted for each population, one value per population in
 * the order of `pops` of the result.
 */
export interface DiversityNumVars {
  /**
   * The variants the population called something at and had
   * `minNumIndividuals` called genotypes in, which is the divisor of every
   * `mean` and of the `ratio`.
   */
  withData: Uint32Array;

  /**
   * Of those, the ones whose called alleles reached `numCalledAlleles`, which
   * is the divisor of every `inDraw` value. It is 0 for a call that gave no
   * draw and for a draw the population reached at no variant.
   */
  inDraw: Uint32Array;
}

/**
 * What `calcPopDiversity` gives back. A statistic that was not asked for is
 * `null`.
 *
 * A population for which no variant counted has 0 in every count, NaN in
 * every mean and ratio and NaN in `fis`, and a column of zeros in the
 * spectrum. It is not an error: a user whose filters left one population
 * with nothing still wants the others.
 */
export interface PopDiversity {
  /**
   * The names of the populations, in the order the keys of `pops` iterate
   * in, which is the order of every array of the result.
   */
  pops: readonly string[];

  /**
   * How many alleles each population called. An allele numbered 1 at one
   * variant is not the allele numbered 1 at the next, so `total` is a sum of
   * per variant counts and never a count of distinct things across the
   * dataset, and `mean` is the allelic richness a user compares between
   * populations.
   */
  numAlleles: PopAlleleCounts | null;

  /**
   * How many alleles each population called that no other population of the
   * call called at the same variant.
   *
   * The divisor of `mean` is `numVarsEveryPop` and not `numVars.withData`: a
   * variant where one population has too little called is out of the private
   * alleles of every population, since there every allele of every other
   * population would be private and the count would measure the missing
   * data. `inDraw` is over `numVarsEveryPopInDraw` for the same reason.
   *
   * With one population every allele it called is private, since there is no
   * other population to hold it, so `total` is then the `total` of
   * `numAlleles`.
   */
  privateAlleles: PopAlleleCounts | null;

  /** How many of the variants each population called more than one allele
   * at. */
  variableVarsRatio: VariableVarsRatio | null;

  /**
   * How many variants of each population show each count of their rarer
   * allele in a draw of `numCalledAlleles`, under the name of the
   * population: one value per count of the rarer allele, 0 to
   * `numCalledAlleles / 2` rounded down.
   *
   * It is folded because nothing in a VCF says which allele is the ancestral
   * one, so the counts `j` and `numCalledAlleles - j` are one bin.
   *
   * The values are not whole numbers: each variant in the draw for the
   * population gives every bin the chance that a draw shows that many copies of
   * the rarer allele there, so a spectrum sums to `numVars.inDraw` of that
   * population. It is `null` when nobody asked for it, which a call that names
   * no statistic does not.
   */
  foldedSfs: Record<string, Float64Array> | null;

  /**
   * One minus the mean observed heterozygosity of each population over its
   * mean unbiased expected one.
   *
   * It is 0 when the genotypes of the population are in the proportions its
   * allele frequencies would give if its individuals paired at random,
   * positive when it holds fewer heterozygous genotypes than that, which
   * inbreeding, selfing and a population split into unmixed groups all
   * produce, and negative when it holds more. It is Nei's F_IS, read on one
   * population on its own, and not Weir and Cockerham's, which comes out of
   * a decomposition of the variance across populations.
   *
   * It is NaN for a population that has no F_IS, in four cases. One that no
   * variant counted for. One no counted variant of which carries both
   * heterozygosities, which a population whose counted variants hold no
   * whole called genotype reaches: individuals whose genotypes are all half
   * called count their variants and have no observed heterozygosity at any
   * of them. One whose mean unbiased expected heterozygosity is 0, every
   * variant it counted having held one allele. And every population of a
   * haploid dataset, where no genotype can be heterozygous and the ratio
   * would be 1 wherever the population has any diversity. The draw does not
   * touch it: the observed heterozygosity is a property of whole genotypes
   * and not of a sample of alleles.
   */
  fis: Float64Array | null;

  /** How many variants counted for each population. */
  numVars: DiversityNumVars;

  /**
   * The variants that counted for every population, which is the divisor of
   * the `mean` of `privateAlleles`.
   */
  numVarsEveryPop: number;

  /**
   * Of those, the ones every population reached `numCalledAlleles` at, which
   * is the divisor of the `inDraw` of `privateAlleles`.
   */
  numVarsEveryPopInDraw: number;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  passStats: PassStats;
}

/**
 * How much variety each population of `variants` holds, in one pass over
 * them.
 *
 * The five statistics are how many alleles a population called, how many of
 * those no other population called, how many of the variants vary in it, how
 * its variants are spread over the count of their rarer allele, and F_IS,
 * how far its genotypes are from the proportions its allele frequencies
 * would give. All five are built from the same thing, how often each
 * population called each allele at each variant, so asking for several of
 * them costs one reading of the variants and not five.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, and the `Variants`
 * is as it was afterwards. The populations are resolved against the
 * individuals that pass gives, which are the ones a `filterIndividuals` kept
 * when the `Variants` carries one.
 *
 * The default `stats` is the four that need no draw, so a call that gives no
 * options runs and gives the alleles each population called, the private ones
 * among them, the variants that vary in it and F_IS, each over the called
 * alleles the population has. The fifth, the folded spectrum, has bins that
 * are counts of the rarer allele in a draw of `numCalledAlleles`, so a user
 * who wants it names it in `stats` and gives that draw.
 *
 * pyNei has none of the five, so no result here is compared with it.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * `stats` is not an array of names, when a name of it is of no statistic and
 * when it names none at all; when `pops` is not an object of names to arrays
 * of names, when a population names an individual the pass does not give,
 * names one twice or names none, and when it holds no population; when
 * `numCalledAlleles` is not a whole number of 2 or more, when it is above the
 * individuals of the dataset times their ploidy and when
 * `minNumIndividuals` is not a whole number of 0 or more; when `folded_sfs` is
 * among the statistics and no `numCalledAlleles` was given; when the
 * source cannot be read, a wrong line of a VCF among the causes; when the
 * pass gives no variant, whether the source holds none or the steps kept
 * none; and when `init` has not been awaited.
 */
export function calcPopDiversity(
  variants: Variants,
  options: CalcPopDiversityOptions = {},
): PopDiversity {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  const stats = theStats(options.stats);
  const pops = thePops(options.pops);
  const numCalledAlleles = theNumCalledAlleles(options.numCalledAlleles);
  const minNumIndividuals =
    options.minNumIndividuals === undefined
      ? defaultMinNumIndividuals()
      : wholeNumberOfZeroOrMore("minNumIndividuals", options.minNumIndividuals);
  // The steps of the pass are a copy of the list, made after every argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const diversity = source.calc_pop_diversity(
    steps.of_a_pass(),
    stats,
    pops.names,
    pops.individuals,
    pops.numIndividualsPerPop,
    numCalledAlleles,
    minNumIndividuals,
  );
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the
  // user gets are the copies.
  try {
    // The name of a population and its counts are read at the same place of
    // these arrays, so a result with no names would give a user numbers with
    // nothing to say which population each of them is of.
    const names = diversity.pop_names();
    if (names === undefined) {
      throw new Error(
        "popnei: the pass named none of the populations it ran over",
      );
    }
    const numVars = {
      withData: arrayOf(diversity.num_vars_with_data(), "variants with data"),
      inDraw: arrayOf(diversity.num_vars_in_draw(), "variants in the draw"),
    };
    return {
      pops: Object.freeze(names),
      numAlleles: meanOf(
        countsOf(
          diversity.num_alleles_total(),
          diversity.num_alleles_in_draw(),
          numVars.withData,
          "alleles",
        ),
      ),
      privateAlleles: meanOf(
        countsOf(
          diversity.private_alleles_total(),
          diversity.private_alleles_in_draw(),
          everyPop(diversity.num_vars_every_pop(), numVars.withData.length),
          "private alleles",
        ),
      ),
      variableVarsRatio: ratioOf(
        countsOf(
          diversity.variable_vars_total(),
          diversity.variable_vars_in_draw(),
          numVars.withData,
          "variable variants",
        ),
      ),
      foldedSfs: spectraOf(
        diversity.folded_sfs(),
        diversity.num_sfs_bins(),
        names,
      ),
      fis: diversity.fis() ?? null,
      numVars,
      numVarsEveryPop: diversity.num_vars_every_pop(),
      numVarsEveryPopInDraw: diversity.num_vars_every_pop_in_draw(),
      passStats: passStatsOf(diversity.pass_stats()),
    };
  } finally {
    diversity.free();
  }
}

/**
 * The names of the statistics a user asked for, each once and in the order
 * they named them, and the four that need no draw when they named none.
 *
 * Those four are the set the core holds, so a statistic that needs no draw is
 * added there and is in the default of this package and of the Python one
 * with nothing written in either.
 *
 * What is refused here is what is no array of names at all. Which names there
 * are is the core's rule, and so is a `stats` that names none, which the core
 * refuses with the five names and with what such a pass would do, read every
 * variant of the source and compute nothing of them: a copy of that sentence
 * here would be the same rule written in a third language.
 */
function theStats(stats: readonly PopDiversityStat[] | undefined): string[] {
  if (stats === undefined) {
    return statsWithoutADraw();
  }
  const asked = namesOf("stats", stats, {
    oneOfThem: "statistic",
    anExample: "fis",
  });
  return [...new Set(asked)];
}

/**
 * The populations a user named, as the flat arrays the binding crate takes:
 * their names in the order the keys iterate in, the names of the individuals
 * of every one of them one after another, and how many individuals each of
 * them holds. The names are `undefined` when the user named no population,
 * which is one population of every individual.
 *
 * The names of the individuals are not looked up here: they are resolved
 * against the individuals the pass gives, which are those of the source
 * after a filter of individuals when the `Variants` has one, and only the
 * pass knows them.
 */
function thePops(pops: Record<string, readonly string[]> | undefined): {
  names: string[] | undefined;
  individuals: string[];
  numIndividualsPerPop: Uint32Array;
} {
  if (pops === undefined) {
    return {
      names: undefined,
      individuals: [],
      numIndividualsPerPop: new Uint32Array(0),
    };
  }
  return popsOfTheObject(pops);
}

/**
 * How many called alleles every population is brought down to, and
 * `undefined` for a pass that takes no draw.
 *
 * The rule is written here and not taken from the checks that every other
 * whole number goes through, because those say "0 or more" and 0 and 1 are
 * no draw: a draw of one allele finds one allele whatever the population
 * holds. A user who reads a bound and writes a number inside it has to be
 * taken.
 *
 * @throws {Error} When `numCalledAlleles` is not a whole number of 2 or more
 * that the core holds. What the generated code would do with the numbers
 * this refuses is to throw the fraction of 20.5 away and turn -1 into a
 * count of about four thousand million.
 */
function theNumCalledAlleles(
  numCalledAlleles: number | undefined,
): number | undefined {
  if (numCalledAlleles === undefined) {
    return undefined;
  }
  if (
    typeof numCalledAlleles !== "number" ||
    !Number.isSafeInteger(numCalledAlleles) ||
    numCalledAlleles < 2 ||
    numCalledAlleles > LARGEST_WHOLE_NUMBER
  ) {
    throw new Error(
      "popnei: `numCalledAlleles` is how many called alleles every " +
        "population is brought down to, a whole number of 2 or more and at " +
        `most ${LARGEST_WHOLE_NUMBER}, or left out for no draw at all, and ` +
        `${whatWasGiven(numCalledAlleles)} was given`,
    );
  }
  return numCalledAlleles;
}

/**
 * The folded spectrum of each population under its name, or `null` when
 * nobody asked for it.
 *
 * The bins of the populations cross in one array, the first population's
 * first, because an array of arrays is not one of the types that cross from
 * wasm, and `numBins` is what cuts them. Each population's run is copied into
 * an array of its own, so that a user who holds one spectrum does not hold the
 * memory of the others.
 *
 * @throws {Error} When the array does not hold `numBins` values for each of
 * the populations, which is a defect of popnei.
 */
function spectraOf(
  binsOfEveryPop: Float64Array | undefined,
  numBins: number,
  names: readonly string[],
): Record<string, Float64Array> | null {
  if (binsOfEveryPop === undefined) {
    return null;
  }
  if (binsOfEveryPop.length !== numBins * names.length) {
    throw new Error(
      `popnei: the pass gave ${binsOfEveryPop.length} bins of the folded ` +
        `spectrum and not the ${numBins} of each of ${names.length} ` +
        "populations",
    );
  }
  const spectra: Record<string, Float64Array> = {};
  for (const [pop, name] of names.entries()) {
    spectra[name] = binsOfEveryPop.slice(pop * numBins, (pop + 1) * numBins);
  }
  return spectra;
}

/**
 * One count of the result: its total for every population, that total over
 * the variants it is counted over, and its standardized value.
 *
 * The middle one is `mean` for the two counts of alleles and `ratio` for the
 * variable variants, and it is the one arithmetic of this layer.
 */
interface CountOfEveryPop {
  total: Uint32Array;
  over: Float64Array;
  inDraw: Float64Array;
}

/**
 * The count `total` of every population with its standardized value, and
 * `null` when nobody asked for the statistic.
 *
 * @throws {Error} When the pass gave the total of the statistic or its
 * standardized value and not both, and when it gave one value per population
 * of neither of the two, which are both a defect of popnei.
 */
function countsOf(
  total: Uint32Array | undefined,
  inDraw: Float64Array | undefined,
  divisor: Uint32Array,
  statistic: string,
): CountOfEveryPop | null {
  if (total === undefined && inDraw === undefined) {
    return null;
  }
  if (total === undefined || inDraw === undefined) {
    throw new Error(
      `popnei: the pass counted the ${statistic} and gave the total or the ` +
        "value in a draw and not both",
    );
  }
  if (total.length !== divisor.length || inDraw.length !== divisor.length) {
    throw new Error(
      `popnei: the pass counted the variants of ${divisor.length} ` +
        `populations and the ${statistic} of ${total.length}`,
    );
  }
  return { total, over: overTheVariants(total, divisor), inDraw };
}

/** A count of alleles as a user reads it, its middle column the mean over
 * the variants. */
function meanOf(counts: CountOfEveryPop | null): PopAlleleCounts | null {
  return counts === null
    ? null
    : { total: counts.total, mean: counts.over, inDraw: counts.inDraw };
}

/** The variable variants as a user reads them, their middle column the
 * ratio over the variants. */
function ratioOf(counts: CountOfEveryPop | null): VariableVarsRatio | null {
  return counts === null
    ? null
    : { total: counts.total, ratio: counts.over, inDraw: counts.inDraw };
}

/**
 * `total` over `divisor`, one value per population, and NaN where the
 * divisor is 0: a population no variant counted for has no mean and no
 * ratio.
 */
function overTheVariants(
  total: Uint32Array,
  divisor: Uint32Array,
): Float64Array {
  return Float64Array.from(total, (count, pop) => {
    const variants = divisor[pop] as number;
    return variants === 0 ? Number.NaN : count / variants;
  });
}

/** The same count for every one of `numPops` populations, which is what the
 * private alleles are divided by. */
function everyPop(count: number, numPops: number): Uint32Array {
  return new Uint32Array(numPops).fill(count);
}

/**
 * The counts the pass gave for each of its populations.
 *
 * @throws {Error} When the pass gave none of them, which is a defect of
 * popnei: every pass counts the variants of each of its populations,
 * whatever statistics it was asked for.
 */
function arrayOf(counts: Uint32Array | undefined, what: string): Uint32Array {
  if (counts === undefined) {
    throw new Error(`popnei: the pass counted the ${what} of no population`);
  }
  return counts;
}
