/**
 * The statistics of the variants, per population, and of the individuals.
 *
 * A population is a named set of individuals that a calculation treats as a
 * group, and `pops` is how a user names them: an object of population name
 * to the names of its individuals. Every statistic here is calculated for
 * each population over its individuals alone, and every result holds one
 * value per population, in the order the keys of `pops` iterate in. With no
 * `pops` there is one population, named `pop`, of every individual.
 *
 * `calcPerVarDistribs` makes one pass over the variants and gives, for each
 * statistic and population, the mean over the variants that had a value and
 * a histogram of them. The per variant values are not kept: a million of
 * them for each population do not fit a browser tab, and a user who wants
 * them takes the genotypes with `iterBlocks`.
 *
 * `calcPerIndividualStats` makes a pass of its own and gives two numbers for
 * each individual instead: the share of the variants at which its genotype
 * is missing and the share of its called genotypes at which it is
 * heterozygous. It takes no `pops`, since each of its values is of one
 * individual.
 */

import {
  default_bin_type as defaultBinType,
  default_hist_range as defaultHistRange,
  default_min_num_individuals as defaultMinNumIndividuals,
  default_num_bins as defaultNumBins,
  default_poly_threshold as defaultPolyThreshold,
} from "../wasm/popnei.js";

import {
  aNumber,
  aString,
  namesOf,
  whatWasGiven,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * The five statistics `calcPerVarDistribs` calculates, each named as the
 * field of the result that holds it is named in Python.
 */
const THE_STATISTICS = [
  "obs_het",
  "maf",
  "exp_het",
  "unbiased_exp_het",
  "poly_vars_ratio",
] as const;

/**
 * The name of one of the five statistics of a variant: the observed
 * heterozygosity, the heterozygous genotypes of a population over its called
 * ones; the major allele frequency, the count of its commonest allele over
 * its called alleles; the expected heterozygosity, the chance that gene
 * copies taken at random from the population are not all of the same allele,
 * plain and corrected for the frequencies being estimated from the copies
 * the statistic is computed over; and the polymorphism ratio, how many of
 * the variants vary in the population, which is a count and not a
 * distribution.
 */
export type PerVarStat = (typeof THE_STATISTICS)[number];

/** The two kinds of bins a histogram is made of. */
export type BinType = "linear" | "logarithmic";

/** The histogram the values of every statistic are counted in. */
export interface HistKwargs {
  /**
   * The two ends of the histogram, `[0, 1]` when it is not given, which is
   * where these statistics live.
   */
  range?: readonly [number, number];

  /** How many bins it holds, 40 when it is not given. */
  numBins?: number;

  /**
   * `"linear"` for bins of equal width or `"logarithmic"` for bins of equal
   * ratio, whose `range` has to start above 0. Equal widths when it is not
   * given. pyNei spells the first one `lineal`, the Spanish word.
   */
  binType?: BinType;
}

/** What `calcPerVarDistribs` calculates, for which populations and how. */
export interface PerVarDistribsOptions {
  /**
   * Which of the five statistics to calculate, all of them when it is not
   * given. Asking for fewer is a saving of work and changes no value, and a
   * result holds `null` for one nobody asked for.
   */
  stats?: readonly PerVarStat[];

  /**
   * The populations: an object of population name to the names of its
   * individuals, which are looked up among the individuals the pass gives.
   * With no `pops` there is one population, `pop`, of every individual.
   *
   * The populations of every result are in the order the keys of this object
   * iterate in, which is the order they were written in, and the numeric
   * order for names that are whole numbers, as JavaScript has it.
   */
  pops?: Record<string, readonly string[]>;

  /**
   * How many called genotypes a population needs at a variant for the
   * variant to have a value there, 20 when it is not given. The test is on
   * the called data counted in genotypes, the called alleles of the
   * population over the ploidy, which is a half when a genotype is half
   * called, and the variant has no value when that number is strictly less
   * than the threshold. A variant with no value in a population is out of
   * the mean and in no bin of the histogram of that population.
   */
  minNumIndividuals?: number;

  /** The histogram every statistic is counted in. */
  histKwargs?: HistKwargs;

  /**
   * The two expected heterozygosities' alone: the number the allele
   * frequencies are raised to, and how many copies the unbiased one draws,
   * which is the ploidy of the variants when it is not given.
   */
  ploidy?: number;

  /**
   * The polymorphism ratio's alone: below that major allele frequency a
   * variant is polymorphic in a population, 0.95 when it is not given. A
   * number from 0 to 1, both included.
   */
  polyThreshold?: number;
}

/**
 * The distribution of one statistic over the variants of a pass, per
 * population.
 *
 * A variant has no value of a statistic in a population when the population
 * has too little data at it, and such a variant is out of the mean and in no
 * bin, so the histograms of two populations can count different numbers of
 * variants. A value outside the range of the bins is in the mean and in no
 * bin.
 */
export interface StatsDistrib {
  /**
   * The mean over the variants that had a value, one number per population,
   * in the order of `pops` of the result, and NaN for a population in which
   * no variant had one.
   */
  mean: Float64Array;

  /**
   * The edges of the bins, one more number than there are bins.
   *
   * The four distributions of one result share this array, as pyNei's do,
   * so it is read only: a number written into the edges of one statistic
   * would be in the edges of the other three.
   */
  histBinEdges: Readonly<Float64Array>;

  /**
   * How many variants fell in each bin, the bins of one population after
   * the bins of the one before it: the count of the bin `b` of the
   * population `p` is at `p * numBins + b`.
   */
  histCounts: Uint32Array;
}

/**
 * How many of the variants of a pass vary in each population.
 *
 * A variant is polymorphic in a population when its major allele frequency
 * there is below `polyThreshold`, strictly, and variable when that frequency
 * is below 1. Both are counted among the variants that have a major allele
 * frequency in the population. Every array holds one value per population,
 * in the order of `pops` of the result.
 */
export interface PolyVarsStats {
  /** The polymorphic variants of each population. */
  numPoly: Uint32Array;

  /** `numPoly` over `totNumVariantsWithData`, NaN when the latter is 0. */
  polyRatio: Float64Array;

  /** `numPoly` over `numVariable`, NaN when the latter is 0. */
  polyRatioOverVariables: Float64Array;

  /** The variable variants of each population. */
  numVariable: Uint32Array;

  /**
   * The variants that have a major allele frequency in each population,
   * which the other two counts are among.
   */
  totNumVariantsWithData: Uint32Array;
}

/**
 * What `calcPerVarDistribs` gives back. A statistic that was not asked for
 * is `null`.
 */
export interface PerVarDistribs {
  /**
   * The name of each population, in the order the keys of `pops` iterate
   * in, which is the order of every array of the result.
   */
  pops: readonly string[];

  /** The distribution of the observed heterozygosity. */
  obsHet: StatsDistrib | null;

  /** The distribution of the major allele frequency. */
  maf: StatsDistrib | null;

  /** The distribution of the plain expected heterozygosity. */
  expHet: StatsDistrib | null;

  /** The distribution of the unbiased expected heterozygosity. */
  unbiasedExpHet: StatsDistrib | null;

  /** The counts of the polymorphism ratio. */
  polyVarsRatio: PolyVarsStats | null;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  passStats: PassStats;
}

/**
 * Up to five statistics of every variant and every population of `variants`,
 * in one pass over them, as a mean and a histogram each.
 *
 * The statistics are the observed heterozygosity, the heterozygous genotypes
 * of a population over its called ones; the major allele frequency, the
 * count of the commonest allele over the called alleles; the expected
 * heterozygosity, the chance that gene copies taken at random from the
 * population are not all of the same allele, plain and corrected for the
 * frequencies being estimated from the copies the statistic is computed
 * over; and the polymorphism ratio, how many of the variants vary in the
 * population, which is a count and not a distribution.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, and the `Variants`
 * is as it was afterwards. The populations are resolved against the
 * individuals that pass gives, which are the ones a `filterIndividuals` kept
 * when the `Variants` carries one.
 *
 * A value falls in the bin whose left edge is at most the value and whose
 * right edge is above it, and the last bin takes its right edge too, as
 * `numpy.histogram` does; a value outside the range of the bins is in no bin
 * and in the mean.
 *
 * It is pyNei's `calc_per_var_distribs` under the names of this package,
 * with these differences: `expHet` of the result is the plain expected
 * heterozygosity and `unbiasedExpHet` the unbiased one, where pyNei's
 * `exp_het` holds whichever its `unbiased_exp_het` argument chose, the
 * unbiased one by default, so a user who reads `exp_het` of both libraries
 * reads two numbers; there is no `num_threads`, since wasm has one thread;
 * the observed heterozygosity is held to `minNumIndividuals` too, which
 * pyNei exempts; a block in which nothing is called gives no value, where
 * pyNei gives the expected heterozygosity of such a block a 1; the unbiased
 * correction is the one of the ploidy in hand, where pyNei applies the
 * diploid one at every ploidy; `ploidy` is the exponent alone, where pyNei
 * also counts with it the alleles the individuals are expected to hold; a
 * duplicated name in a population, an empty population and an empty `pops`
 * are refused; and the result has `passStats`.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * `stats` is not an array of names, when a name of it is of no statistic and
 * when it names none at all; when `pops` is not an object of names to arrays
 * of names, when a population names an individual the pass does not give,
 * names one twice or names none, and when it holds no population; when
 * `minNumIndividuals` or `ploidy` is not a whole number of 0 or more, and
 * when the ploidy is 0 or above 255; when `histKwargs` holds a key that is
 * none of the three, a range that does not run from a number up to a larger
 * one, no bin, a kind of bins that is neither of the two, or a logarithmic
 * range that starts at 0 or below; when `polyThreshold` is not a number from
 * 0 to 1; when the source cannot be read, a wrong line of a VCF among the
 * causes; when the pass gives no variant, whether the source holds none or
 * the steps kept none; and when `init` has not been awaited.
 */
export function calcPerVarDistribs(
  variants: Variants,
  options: PerVarDistribsOptions = {},
): PerVarDistribs {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  const stats = theStats(options.stats);
  const pops = thePops(options.pops);
  const minNumIndividuals =
    options.minNumIndividuals === undefined
      ? defaultMinNumIndividuals()
      : wholeNumberOfZeroOrMore("minNumIndividuals", options.minNumIndividuals);
  const histogram = theHistogram(options.histKwargs);
  const ploidy =
    options.ploidy === undefined
      ? undefined
      : wholeNumberOfZeroOrMore("ploidy", options.ploidy);
  const polyThreshold =
    options.polyThreshold === undefined
      ? defaultPolyThreshold()
      : aNumber("polyThreshold", options.polyThreshold);
  // The steps of the pass are a copy of the list, made after every argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const distribs = source.calc_per_var_distribs(
    steps.of_a_pass(),
    stats,
    pops.names,
    pops.individuals,
    pops.numIndividualsPerPop,
    minNumIndividuals,
    histogram.start,
    histogram.end,
    histogram.numBins,
    histogram.binType,
    ploidy,
    polyThreshold,
  );
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the
  // user gets are the copies.
  try {
    const edges = distribs.hist_bin_edges();
    return {
      pops: Object.freeze(distribs.pop_names()),
      obsHet: distribOf(
        edges,
        distribs.obs_het_mean(),
        distribs.obs_het_hist_counts(),
        "obs_het",
      ),
      maf: distribOf(
        edges,
        distribs.maf_mean(),
        distribs.maf_hist_counts(),
        "maf",
      ),
      expHet: distribOf(
        edges,
        distribs.exp_het_mean(),
        distribs.exp_het_hist_counts(),
        "exp_het",
      ),
      unbiasedExpHet: distribOf(
        edges,
        distribs.unbiased_exp_het_mean(),
        distribs.unbiased_exp_het_hist_counts(),
        "unbiased_exp_het",
      ),
      polyVarsRatio: polyVarsStatsOf(
        distribs.num_poly(),
        distribs.poly_ratio(),
        distribs.poly_ratio_over_variables(),
        distribs.num_variable(),
        distribs.num_vars_with_data(),
      ),
      passStats: passStatsOf(distribs.pass_stats()),
    };
  } finally {
    distribs.free();
  }
}

/**
 * The names of the statistics a user asked for, each once and in the order
 * they named them, and the five of them when they named none.
 *
 * Which names there are is the binding crate's rule, as the two kinds of
 * bins are: what is refused here is what is no array of names and an array
 * of none, which would make a pass over the whole source for nothing.
 */
function theStats(stats: readonly PerVarStat[] | undefined): string[] {
  if (stats === undefined) {
    return [...THE_STATISTICS];
  }
  const asked = namesOf("stats", stats, {
    oneOfThem: "statistic",
    anExample: "maf",
  });
  if (asked.length === 0) {
    throw new Error(
      "popnei: `stats` names no statistic, and a result holds the ones that " +
        "were asked for: leave `stats` out for the five of them",
    );
  }
  return [...new Set(asked)];
}

/**
 * The populations a user named, as the flat arrays the binding crate takes:
 * their names in the order the keys iterate in, the names of the individuals
 * of every one of them one after another, and how many individuals each of
 * them holds. The names are `undefined` when the user named no population.
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
  if (typeof pops !== "object" || pops === null || Array.isArray(pops)) {
    throw new Error(
      "popnei: `pops` is an object of the name of a population to the names " +
        `of its individuals, {pop1: ["ind00", "ind01"]}, and ` +
        `${whatWasGiven(pops)} was given`,
    );
  }
  const names = Object.keys(pops);
  const individuals: string[] = [];
  const numIndividualsPerPop = new Uint32Array(names.length);
  for (const [which, pop] of names.entries()) {
    const ofThePop = namesOf(`pops.${pop}`, pops[pop], {
      oneOfThem: "individual",
      anExample: "ind00",
    });
    individuals.push(...ofThePop);
    numIndividualsPerPop[which] = ofThePop.length;
  }
  return { names, individuals, numIndividualsPerPop };
}

/** The three keys the histogram is given under. */
const HIST_KEYS = ["range", "numBins", "binType"];

/**
 * The two ends of the histogram, how many bins it holds and of which kind,
 * out of the object a user gave, which is read and not changed.
 *
 * A key that is none of the three is refused: pyNei ignores such a key, so a
 * user who writes `nbins` gets the 40 bins of the default with nothing said,
 * which is a result that is not the one they asked for and that says so
 * nowhere.
 */
function theHistogram(histKwargs: HistKwargs | undefined): {
  start: number;
  end: number;
  numBins: number;
  binType: string;
} {
  for (const key of Object.keys(histKwargs ?? {})) {
    if (!HIST_KEYS.includes(key)) {
      throw new Error(
        `popnei: \`${key}\` is not a key of \`histKwargs\`, whose keys are ` +
          "`range`, the two ends of the histogram, `numBins` and `binType`",
      );
    }
  }
  const range = histKwargs?.range;
  let start: number;
  let end: number;
  if (range === undefined) {
    const theDefault = defaultHistRange();
    const [defaultStart, defaultEnd] = theDefault;
    if (defaultStart === undefined || defaultEnd === undefined) {
      throw new Error(
        "popnei: the default range of the histogram is not two ends but " +
          `${theDefault.length} numbers`,
      );
    }
    start = defaultStart;
    end = defaultEnd;
  } else if (!Array.isArray(range) || range.length !== 2) {
    throw new Error(
      "popnei: `histKwargs.range` is the two ends of the histogram, [0, 1], " +
        `and ${whatWasGiven(range)} was given`,
    );
  } else {
    start = aNumber("histKwargs.range[0]", range[0]);
    end = aNumber("histKwargs.range[1]", range[1]);
  }
  const numBins =
    histKwargs?.numBins === undefined
      ? defaultNumBins()
      : wholeNumberOfZeroOrMore("histKwargs.numBins", histKwargs.numBins);
  const binType =
    histKwargs?.binType === undefined
      ? defaultBinType()
      : aString("histKwargs.binType", histKwargs.binType);
  return { start, end, numBins, binType };
}

/**
 * The distribution of one statistic, or `null` when nobody asked for it: the
 * mean of each population and the counts of its bins, which the binding
 * crate gives, under the edges of the bins of the pass.
 */
function distribOf(
  histBinEdges: Readonly<Float64Array>,
  mean: Float64Array | undefined,
  histCounts: Uint32Array | undefined,
  statistic: string,
): StatsDistrib | null {
  if (mean === undefined && histCounts === undefined) {
    return null;
  }
  if (mean === undefined || histCounts === undefined) {
    throw new Error(
      `popnei: the pass calculated the ${statistic} and gave its mean or its ` +
        "histogram and not both",
    );
  }
  return { mean, histBinEdges, histCounts };
}

/**
 * The counts of the polymorphism ratio, or `null` when nobody asked for
 * them.
 */
function polyVarsStatsOf(
  numPoly: Uint32Array | undefined,
  polyRatio: Float64Array | undefined,
  polyRatioOverVariables: Float64Array | undefined,
  numVariable: Uint32Array | undefined,
  totNumVariantsWithData: Uint32Array | undefined,
): PolyVarsStats | null {
  const counts = [
    numPoly,
    polyRatio,
    polyRatioOverVariables,
    numVariable,
    totNumVariantsWithData,
  ];
  if (counts.every((count) => count === undefined)) {
    return null;
  }
  if (
    numPoly === undefined ||
    polyRatio === undefined ||
    polyRatioOverVariables === undefined ||
    numVariable === undefined ||
    totNumVariantsWithData === undefined
  ) {
    throw new Error(
      "popnei: the pass counted the polymorphism ratio and gave some of its " +
        "five arrays and not every one of them",
    );
  }
  return {
    numPoly,
    polyRatio,
    polyRatioOverVariables,
    numVariable,
    totNumVariantsWithData,
  };
}

/** What `calcPerIndividualStats` gives back. */
export interface PerIndividualStats {
  /**
   * The name of each individual the pass gave, in its order, which is the
   * order of the source unless a `filterIndividuals` named them in another
   * one, and the order of the two arrays below.
   */
  individuals: readonly string[];

  /**
   * The variants at which the genotype of the individual is missing, a half
   * called genotype among them, over the variants of the pass.
   */
  missingGtRate: Float64Array;

  /**
   * The variants at which the genotype of the individual is called and its
   * alleles are not all the same, over its called genotypes, and NaN for an
   * individual that called none of them.
   */
  obsHetRate: Float64Array;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  passStats: PassStats;
}

/**
 * The missing rate and the heterozygosity rate of every individual of
 * `variants`, in one pass over them.
 *
 * The missing rate is the share of the variants at which the individual has
 * no genotype, and a half called genotype is missing and not heterozygous.
 * The heterozygosity rate is the share of its called genotypes at which its
 * alleles are not all the same. The first tells a user which individuals
 * were badly genotyped, and the second which ones are more heterozygous than
 * the rest, a sign of a mixed sample or of an outcrossed individual among
 * inbred ones. An individual that called no genotype has a missing rate of 1
 * and no heterozygosity rate, NaN.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, and the `Variants`
 * is as it was afterwards. The individuals are the ones that pass gives,
 * which a `filterIndividuals` kept, in the order the user named them there.
 *
 * It is pyNei's `calc_per_sample_stats` under the names of this package,
 * with these differences: the heterozygosity rate divides by the called
 * genotypes of the individual, where pyNei divides by every variant, so an
 * individual with more missing data looks less heterozygous there, and
 * popnei's number is what plink2's `--het` gives, with the missing rate
 * beside it saying what pyNei's one number said; pyNei's sample is popnei's
 * individual; there is no `num_threads`, since wasm has one thread; and the
 * result has `passStats`.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when the
 * source cannot be read, a wrong line of a VCF among the causes; when the
 * pass gives no variant, whether the source holds none or the steps kept
 * none; and when `init` has not been awaited.
 */
export function calcPerIndividualStats(
  variants: Variants,
): PerIndividualStats {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  // The steps of the pass are a copy of the list: the call takes it over and
  // frees it.
  const stats = source.calc_per_individual_stats(steps.of_a_pass());
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the user
  // gets are the copies.
  try {
    return {
      individuals: Object.freeze(stats.individuals()),
      missingGtRate: stats.missing_gt_rate(),
      obsHetRate: stats.obs_het_rate(),
      passStats: passStatsOf(stats.pass_stats()),
    };
  } finally {
    stats.free();
  }
}
