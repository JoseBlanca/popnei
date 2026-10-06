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
 *
 * `calcVarDensity` takes no `pops` either: it counts the variants in windows
 * along each chromosome, and reads no genotype.
 *
 * `calcVariantsSummary` gives the three in one pass, where they take three,
 * each the same to the bit as its own call gives it.
 *
 * The three build their results at the end of the pass from totals they add
 * up block by block, so each of them, and `calcVariantsSummary`, can give,
 * while the pass runs, the result over the variants read so far, to a
 * function of the application, `onSoFar`, which a page draws the histograms
 * from as they fill. Python has no such option, and no
 * `calcVariantsSummary`: it has no page to draw on.
 */

import type {
  PerIndividualStats as PerIndividualStatsOfTheCore,
  PerVarDistribs as PerVarDistribsOfTheCore,
  VarDensityOfAPass,
  VariantsSummaryOfAPass,
} from "../wasm/popnei.js";
import {
  ArgumentsOfTheDensity,
  ArgumentsOfThePass,
  default_bin_type as defaultBinType,
  default_hist_range as defaultHistRange,
  default_min_num_individuals as defaultMinNumIndividuals,
  default_num_bins as defaultNumBins,
  default_poly_threshold as defaultPolyThreshold,
} from "../wasm/popnei.js";

import {
  anObjectOfOptions,
  aNumber,
  aString,
  distanceInBasePairs,
  namesOf,
  popsOfTheObject,
  whatWasGiven,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * The six statistics `calcPerVarDistribs` calculates, each named as the
 * field of the result that holds it is named in Python.
 */
const THE_STATISTICS = [
  "obs_het",
  "maf",
  "exp_het",
  "unbiased_exp_het",
  "poly_vars_ratio",
  "missing_rate",
] as const;

/**
 * The name of one of the six statistics of a variant: the observed
 * heterozygosity, the heterozygous genotypes of a population over its called
 * ones; the major allele frequency, the count of its commonest allele over
 * its called alleles; the expected heterozygosity, the chance that gene
 * copies taken at random from the population are not all of the same allele,
 * plain and corrected for the frequencies being estimated from the copies
 * the statistic is computed over; the polymorphism ratio, how many of the
 * variants vary in the population, which is a count and not a distribution;
 * and the missing rate, the missing genotypes of the population over its
 * individuals, called or not, a half called genotype being missing.
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

/**
 * The two options of the result so far, which `calcPerVarDistribs`,
 * `calcPerIndividualStats`, `calcVarDensity` and `calcVariantsSummary` take,
 * `T` being the result the calculation returns.
 */
export interface SoFarOptions<T> {
  /**
   * A function that is given, while the pass runs, the result over the
   * variants read so far, of the type the calculation returns, with its
   * `passStats` as they stand: `passStats.numVars` is how many variants it
   * covers. It is called after a block of the pass, the last one too, when
   * `soFarEvery` seconds have gone by since the pass started or since the
   * last call, so the last call, when there is one, is given what the
   * calculation then returns.
   *
   * The result so far is the one the calculation would return over those
   * variants alone. The bins of the histograms are those of the last result
   * from the first call; the density with no length for a chromosome has
   * its windows up to the last variant read so far. Each call builds the
   * result as the final one is built, a few hundred numbers for the
   * histograms, two for each individual and four for each window of the
   * density, so a density of millions of windows asks for a longer
   * `soFarEvery`.
   *
   * A value it throws ends the pass and is what the calculation throws, as
   * a value thrown by the function of `Variants.onProgress` is, and that
   * function is told nothing more of the pass. From inside it, the `free()`
   * of the variants being read is refused as it is from inside the function
   * of `onProgress`.
   *
   * It is called and not awaited: a function that returns a promise is
   * called, the pass goes on without waiting for it, and what the promise
   * rejects with does not stop the pass.
   */
  onSoFar?: (soFar: T) => void;

  /**
   * How many seconds go by between two calls of `onSoFar`, a finite number
   * of 0 or more, 2 when it is not given; 0 calls it after every block. It is
   * an `Error` without `onSoFar`, since alone it does nothing.
   */
  soFarEvery?: number;
}

/** How many seconds go by between two calls of `onSoFar` when not given. */
const DEFAULT_SO_FAR_EVERY = 2;

/** The keys of `PerVarDistribsOptions` but the two of the result so far. */
const PER_VAR_KEYS = [
  "stats",
  "pops",
  "minNumIndividuals",
  "histKwargs",
  "ploidy",
  "polyThreshold",
];

/** The two keys of the options of the result so far. */
const SO_FAR_KEYS = ["onSoFar", "soFarEvery"];

/** What `calcPerVarDistribs` calculates, for which populations and how. */
export interface PerVarDistribsOptions extends SoFarOptions<PerVarDistribs> {
  /**
   * Which of the six statistics to calculate, all of them when it is not
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
   * than the threshold, for every statistic but the missing rate, which
   * every variant has. A variant with no value in a population is out of
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
   * The five distributions of one result share this array, as pyNei's do,
   * so it is read only: a number written into the edges of one statistic
   * would be in the edges of the other four.
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
   * The distribution of the missing rate, which every variant has in every
   * population.
   */
  missingRate: StatsDistrib | null;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  passStats: PassStats;
}

/**
 * Up to six statistics of every variant and every population of `variants`,
 * in one pass over them, as a mean and a histogram each.
 *
 * The statistics are the observed heterozygosity, the heterozygous genotypes
 * of a population over its called ones; the major allele frequency, the
 * count of the commonest allele over the called alleles; the expected
 * heterozygosity, the chance that gene copies taken at random from the
 * population are not all of the same allele, plain and corrected for the
 * frequencies being estimated from the copies the statistic is computed
 * over; the polymorphism ratio, how many of the variants vary in the
 * population, which is a count and not a distribution; and the missing
 * rate, the missing genotypes of the population over its individuals, called
 * or not, a half called genotype being missing. The missing rate has a value
 * at every variant, whatever `minNumIndividuals` is, and a variant with
 * nothing called in a population has a rate of 1 there.
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
 * `onSoFar` is given the distributions over the variants read so far while
 * the pass runs, every `soFarEvery` seconds, as `SoFarOptions` says.
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
 * are refused; the result has `passStats`; pyNei has no missing rate; and
 * neither pyNei nor the Python package has `onSoFar` or `soFarEvery`.
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
 * 0 to 1; when `onSoFar` is not a function, when `soFarEvery` is not a
 * finite number of 0 or more and when it is given without `onSoFar`; when
 * the source cannot be read, a wrong line of a VCF among the causes; when the
 * pass gives no variant, whether the source holds none or the steps kept
 * none; and when `init` has not been awaited. It throws what `onSoFar`
 * threw.
 */
export function calcPerVarDistribs(
  variants: Variants,
  options: PerVarDistribsOptions = {},
): PerVarDistribs {
  theWasmHasToBeLoaded();
  anObjectOfOptions("calcPerVarDistribs", options, [
    ...PER_VAR_KEYS,
    ...SO_FAR_KEYS,
  ]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const asked = thePerVarArguments(options);
  const soFar = theResultSoFar(options, distribsOf);
  // The steps of the pass and its arguments are objects of the binding
  // crate, made after every argument was checked so that nothing refused
  // here leaves one behind: the call takes them over and frees them.
  const distribs = whileTheRunReads(() =>
    source.calc_per_var_distribs(
      steps.of_a_pass(),
      argumentsOfThePassOf(asked),
      soFar.told,
      soFar.every,
    ),
  );
  return distribsOf(distribs);
}

/**
 * The arguments of the distributions of the statistics of each variant, out
 * of the options a user gave `calcPerVarDistribs` or the `perVar` of
 * `calcVariantsSummary`, each checked and with its default where it was not
 * given.
 */
interface PerVarArguments {
  stats: string[];
  pops: ReturnType<typeof thePops>;
  minNumIndividuals: number;
  histogram: ReturnType<typeof theHistogram>;
  ploidy: number | undefined;
  polyThreshold: number;
}

/**
 * The arguments of the distributions out of `options`, whose keys the caller
 * has checked.
 *
 * @throws {Error} What the checks of `calcPerVarDistribs` refuse, each
 * named in its doc comment.
 */
function thePerVarArguments(
  options: Omit<PerVarDistribsOptions, keyof SoFarOptions<PerVarDistribs>>,
): PerVarArguments {
  return {
    stats: theStats(options.stats),
    pops: thePops(options.pops),
    minNumIndividuals:
      options.minNumIndividuals === undefined
        ? defaultMinNumIndividuals()
        : wholeNumberOfZeroOrMore(
            "minNumIndividuals",
            options.minNumIndividuals,
          ),
    histogram: theHistogram(options.histKwargs),
    ploidy:
      options.ploidy === undefined
        ? undefined
        : wholeNumberOfZeroOrMore("ploidy", options.ploidy),
    polyThreshold:
      options.polyThreshold === undefined
        ? defaultPolyThreshold()
        : aNumber("polyThreshold", options.polyThreshold),
  };
}

/**
 * `asked` as the object of the binding crate that carries it, which the call
 * it is given to takes over and frees.
 */
function argumentsOfThePassOf(asked: PerVarArguments): ArgumentsOfThePass {
  return new ArgumentsOfThePass(
    asked.stats,
    asked.pops.names,
    asked.pops.individuals,
    asked.pops.numIndividualsPerPop,
    asked.minNumIndividuals,
    asked.histogram.start,
    asked.histogram.end,
    asked.histogram.numBins,
    asked.histogram.binType,
    asked.ploidy,
    asked.polyThreshold,
  );
}

/**
 * The distributions of a pass, or of its first blocks, out of what the
 * binding crate gives, which is freed here.
 */
function distribsOf(distribs: PerVarDistribsOfTheCore): PerVarDistribs {
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
      missingRate: distribOf(
        edges,
        distribs.missing_rate_mean(),
        distribs.missing_rate_hist_counts(),
        "missing_rate",
      ),
      passStats: passStatsOf(distribs.pass_stats()),
    };
  } finally {
    distribs.free();
  }
}

/**
 * The function the binding crate calls with the result so far, which builds
 * the result of this package out of what crosses with `built` and gives it
 * to `onSoFar`, and the seconds between two calls; no function when the
 * application gave no `onSoFar`.
 *
 * @throws {Error} When `onSoFar` is not a function, when `soFarEvery` is not
 * a finite number of 0 or more, and when `soFarEvery` is given without
 * `onSoFar`, which would be a value that does nothing.
 */
function theResultSoFar<OfTheCore, T>(
  options: SoFarOptions<T>,
  built: (ofTheCore: OfTheCore) => T,
): { told: ((ofTheCore: OfTheCore) => void) | undefined; every: number } {
  const { onSoFar, soFarEvery } = options;
  if (onSoFar === undefined) {
    if (soFarEvery !== undefined) {
      throw new Error(
        "popnei: `soFarEvery` was given with no `onSoFar`, and it is how often " +
          "`onSoFar` is called, so alone it does nothing",
      );
    }
    return { told: undefined, every: DEFAULT_SO_FAR_EVERY };
  }
  if (typeof onSoFar !== "function") {
    throw new Error(
      "popnei: `onSoFar` is a function that is given the result so far, and " +
        `${whatWasGiven(onSoFar)} was given`,
    );
  }
  let every = DEFAULT_SO_FAR_EVERY;
  if (soFarEvery !== undefined) {
    if (
      typeof soFarEvery !== "number" ||
      !Number.isFinite(soFarEvery) ||
      soFarEvery < 0
    ) {
      throw new Error(
        "popnei: `soFarEvery` is a number of seconds, finite and 0 or more, " +
          `and ${whatWasGiven(soFarEvery)} was given`,
      );
    }
    every = soFarEvery;
  }
  return { told: (ofTheCore) => onSoFar(built(ofTheCore)), every };
}

/**
 * The names of the statistics a user asked for, each once and in the order
 * they named them, and the six of them when they named none.
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
        "were asked for: leave `stats` out for the six of them",
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
  return popsOfTheObject(pops);
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

/** The options of `calcPerIndividualStats`, the two of the result so far. */
export type PerIndividualStatsOptions = SoFarOptions<PerIndividualStats>;

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
 * `onSoFar` is given the rates over the variants read so far while the pass
 * runs, every `soFarEvery` seconds, as `SoFarOptions` says.
 *
 * It is pyNei's `calc_per_sample_stats` under the names of this package,
 * with these differences: the heterozygosity rate divides by the called
 * genotypes of the individual, where pyNei divides by every variant, so an
 * individual with more missing data looks less heterozygous there, and
 * popnei's number is what plink2's `--het` gives, with the missing rate
 * beside it saying what pyNei's one number said; pyNei's sample is popnei's
 * individual; there is no `num_threads`, since wasm has one thread; the
 * result has `passStats`; and neither pyNei nor the Python package has
 * `onSoFar` or `soFarEvery`.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when the
 * options are not an object or hold a key that is not one of the two; when
 * `onSoFar` is not a function, when `soFarEvery` is not a finite number of 0
 * or more and when it is given without `onSoFar`; when the source cannot be
 * read, a wrong line of a VCF among the causes; when the pass gives no
 * variant, whether the source holds none or the steps kept none; and when
 * `init` has not been awaited. It throws what `onSoFar` threw.
 */
export function calcPerIndividualStats(
  variants: Variants,
  options: PerIndividualStatsOptions = {},
): PerIndividualStats {
  theWasmHasToBeLoaded();
  anObjectOfOptions("calcPerIndividualStats", options, [
    "onSoFar",
    "soFarEvery",
  ]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const soFar = theResultSoFar(options, ratesOf);
  // The steps of the pass are a copy of the list, made after every argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const stats = whileTheRunReads(() =>
    source.calc_per_individual_stats(
      steps.of_a_pass(),
      soFar.told,
      soFar.every,
    ),
  );
  return ratesOf(stats);
}

/**
 * The rates of a pass, or of its first blocks, out of what the binding crate
 * gives, which is freed here.
 */
function ratesOf(stats: PerIndividualStatsOfTheCore): PerIndividualStats {
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

/** The options of `calcVarDensity`. */
export interface VarDensityOptions extends SoFarOptions<VarDensity> {
  /**
   * The length of each chromosome, an object of chromosome name to length,
   * which replaces the lengths of the source for every chromosome: one it
   * does not name has no length. When it is not given the lengths are those
   * of the source, the `##contig` lines of a VCF that have a `length` and
   * what a vars file keeps of them.
   */
  chromLengths?: Record<string, number>;
}

/** What `calcVarDensity` gives back, one entry of each array per window. */
export interface VarDensity {
  /** The name of the chromosome of each window. */
  chroms: readonly string[];

  /** The first position of each window, counted from 1. */
  start: Float64Array;

  /** The last position of each window, included. */
  end: Float64Array;

  /** How many variants of the pass are at a position from start to end. */
  numVars: Uint32Array;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  passStats: PassStats;
}

/**
 * How many variants fall in each window of `windowSize` base pairs along
 * each chromosome, in one pass over `variants`.
 *
 * A user sees with it where the variants are crowded, where there are none,
 * a centromere or a region that did not map, and how evenly a filter took
 * variants out. The windows of a chromosome are laid end to end from the
 * position 1 and do not overlap: window k, counted from 0, holds the
 * positions from k x `windowSize` + 1 to (k + 1) x `windowSize`. A window
 * with no variant is in the result with a count of 0.
 *
 * With the length of a chromosome the windows cover it to its end, and the
 * last one ends at the length. Without one, the windows go up to the one
 * that holds the last variant of the chromosome, and that one ends at its
 * full width. A chromosome with a length is in the result whether or not it
 * has a variant. The chromosomes are in the order of the lengths, and after
 * them those with variants and no length, in the order their first variant
 * came; the windows of each in the order of their positions. The order of
 * `chromLengths` is the order JavaScript gives the keys of an object, which
 * puts the keys that are whole numbers first, in ascending order, and then
 * the others in the order they were written: `{X: 1, "10": 1, "2": 1}` gives
 * 2, 10 and X, where the same dict in Python gives X, 10 and 2.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, reading only the
 * chromosome and the position of each variant, and the `Variants` is as it
 * was afterwards. The variants need not be sorted. pyNei has no density of
 * the variants; the result is the columns of the frame of Python, one array
 * each.
 *
 * `onSoFar` is given the density over the variants read so far while the
 * pass runs, every `soFarEvery` seconds, as `SoFarOptions` says: a
 * chromosome with a length has all its windows from the first call, and one
 * with none has them up to its last variant read so far. Python has neither
 * option.
 *
 * @throws {Error} When `windowSize` or a length of `chromLengths` is not a
 * whole number from 1 to 2^53 - 1; when `chromLengths` is not a plain
 * object, a `Map` among the rest; when a variant is past the length of its
 * chromosome, which the message says came from `chromLengths` or from the
 * source, or at the position 0; when the density would have more than 10
 * million windows; when a window ends past 2^53, which a number of
 * JavaScript would round; when the arrays of the windows do not fit in the
 * memory the page has left; when `onSoFar` is not a function, when
 * `soFarEvery` is not a finite number of 0 or more and when it is given
 * without `onSoFar`; when the source cannot be read; when the pass gives no
 * variant; and when `init` has not been awaited. It throws what `onSoFar`
 * threw.
 */
export function calcVarDensity(
  variants: Variants,
  windowSize: number,
  options: VarDensityOptions = {},
): VarDensity {
  theWasmHasToBeLoaded();
  anObjectOfOptions("calcVarDensity", options, [
    "chromLengths",
    ...SO_FAR_KEYS,
  ]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const asked = theDensityArguments(
    "windowSize",
    windowSize,
    options.chromLengths,
  );
  const soFar = theResultSoFar(options, windowsOf);
  // The steps of the pass and its arguments are objects of the binding
  // crate, made after every argument was checked so that nothing refused
  // here leaves one behind: the call takes them over and frees them.
  const density = whileTheRunReads(() =>
    source.calc_var_density(
      steps.of_a_pass(),
      argumentsOfTheDensityOf(asked),
      soFar.told,
      soFar.every,
    ),
  );
  return windowsOf(density);
}

/**
 * The width of the windows of the density and the lengths of the
 * chromosomes, out of what a user gave `calcVarDensity` or the `density` of
 * `calcVariantsSummary`, checked: the lengths are `undefined` when they were
 * not given. `windowSizeName` is what the errors call the width,
 * `windowSize` or `density.windowSize`.
 *
 * @throws {Error} When `windowSize` or a length is not a whole number from 1
 * to 2^53 - 1, and when `chromLengths` is not a plain object.
 */
function theDensityArguments(
  windowSizeName: string,
  windowSize: unknown,
  chromLengths: unknown,
): {
  windowSize: number;
  lengths: ReturnType<typeof theChromLengths>;
} {
  return {
    windowSize: distanceInBasePairs(windowSizeName, windowSize, 1),
    lengths: theChromLengths(chromLengths),
  };
}

/**
 * `asked` as the object of the binding crate that carries it, which the call
 * it is given to takes over and frees.
 */
function argumentsOfTheDensityOf(
  asked: ReturnType<typeof theDensityArguments>,
): ArgumentsOfTheDensity {
  return new ArgumentsOfTheDensity(
    asked.windowSize,
    asked.lengths?.names,
    asked.lengths?.lengths ?? new Float64Array(0),
  );
}

/**
 * The windows of a pass, or of its first blocks, out of what the binding
 * crate gives, which is freed here.
 */
function windowsOf(density: VarDensityOfAPass): VarDensity {
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the user
  // gets are the copies.
  try {
    const names = density.chroms();
    const windowsPerChrom = density.windows_per_chrom();
    const chroms: string[] = [];
    for (const [which, name] of names.entries()) {
      for (let window = 0; window < (windowsPerChrom[which] ?? 0); window++) {
        chroms.push(name);
      }
    }
    return {
      chroms: Object.freeze(chroms),
      start: density.starts(),
      end: density.ends(),
      numVars: density.num_vars(),
      passStats: passStatsOf(density.pass_stats()),
    };
  } finally {
    density.free();
  }
}

/**
 * The names and the lengths of `chromLengths`, in the order its keys
 * iterate in, or `undefined` when it was not given.
 *
 * @throws {Error} When it is not a plain object, one whose prototype is
 * `Object.prototype` or `null`, a `Map` among the rest, or a length is not a
 * whole number from 1 to 2^53 - 1.
 */
function theChromLengths(
  value: unknown,
): { names: string[]; lengths: Float64Array } | undefined {
  if (value === undefined) {
    return undefined;
  }
  // A `Map`, and any object of a class of its own, keeps its entries where
  // `Object.keys` does not look, so it would be read as no lengths and say
  // nothing of it: the windows up to each length would not be there, and a
  // variant past a length would not be refused.
  const prototype =
    typeof value === "object" && value !== null
      ? Object.getPrototypeOf(value)
      : undefined;
  if (prototype !== Object.prototype && prototype !== null) {
    throw new Error(
      "popnei: `chromLengths` is a plain object of chromosome name to length, " +
        `{chr1: 248956422}, and ${whatWasGiven(value)} was given`,
    );
  }
  const given = value as Record<string, unknown>;
  const names = Object.keys(given);
  const lengths = new Float64Array(names.length);
  for (const [which, name] of names.entries()) {
    lengths[which] = distanceInBasePairs(
      `chromLengths.${name}`,
      given[name],
      1,
    );
  }
  return { names, lengths };
}

/**
 * The options of `calcPerVarDistribs` but its `onSoFar` and `soFarEvery`,
 * which `perVar` of `VariantsSummaryOptions` takes: the result so far of a
 * summary is given by the `onSoFar` of the summary.
 */
export type PerVarOptionsOfASummary = Omit<
  PerVarDistribsOptions,
  keyof SoFarOptions<PerVarDistribs>
>;

/**
 * The options of `calcVariantsSummary`: which of the three statistics it
 * gives, each with the options of its own call, and the two of the result
 * so far. A statistic whose key is not there, or is `undefined`, is not
 * given.
 */
export interface VariantsSummaryOptions extends SoFarOptions<VariantsSummary> {
  /**
   * The distributions of the statistics of each variant, with the options
   * of `calcPerVarDistribs` but its `onSoFar` and `soFarEvery`; `{}` for
   * the six statistics of one population of every individual.
   */
  perVar?: PerVarOptionsOfASummary;

  /**
   * `{}` for the rates of each individual of `calcPerIndividualStats`,
   * which takes no option but the two of the result so far.
   */
  perIndividual?: Record<string, never>;

  /**
   * The density of the variants of `calcVarDensity`: the width of a window
   * in base pairs, and the lengths of the chromosomes as its `chromLengths`
   * has them.
   */
  density?: {
    windowSize: number;
    chromLengths?: Record<string, number>;
  };
}

/**
 * What `calcVariantsSummary` gives back: each of the three statistics of
 * the type its own call returns, or `null` when it was not asked for, and
 * the counts of the one pass.
 */
export interface VariantsSummary {
  /** The distributions of the statistics of each variant. */
  perVar: PerVarDistribs | null;

  /** The missing rate and the heterozygosity rate of each individual. */
  perIndividual: PerIndividualStats | null;

  /** The number of variants in each window along each chromosome. */
  density: VarDensity | null;

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept, which each of the three
   * carries too.
   */
  passStats: PassStats;
}

/**
 * What `calcPerVarDistribs`, `calcPerIndividualStats` and `calcVarDensity`
 * give, in one pass over `variants` where the three calls take three.
 *
 * A page that shows the three when a file is opened reads the file once.
 * Each of the three is what its own call gives with the same options, to
 * the bit: the three are added up from the same blocks by the same code,
 * and only the pass is shared, so the code of a page that draws one from
 * its own call draws it from here unchanged. The pass reads the genotypes
 * when `perVar` or `perIndividual` is asked for, and the chromosome and the
 * position when `density` is. What it saves depends on the file: over a
 * plain VCF of 403 MB, 100000 variants of 1000 diploid individuals, one
 * pass took 44% less than the three, 1.247 s against 2.240 s, since the
 * genotypes are parsed once; over the vars file of the same variants, whose
 * blocks take almost no time to read, 2% less, 0.391 s against 0.400 s.
 * Both were measured natively on an Apple M5 Pro on 7 October 2026, as
 * "The three statistics of a file in one pass" of
 * `docs/specs/js_sources.md` has it.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, and the
 * `Variants` is as it was afterwards.
 *
 * `onSoFar` is given the three over the variants read so far while the pass
 * runs, every `soFarEvery` seconds, as `SoFarOptions` says, each the result
 * so far of its own call.
 *
 * An error of any of the three ends the pass, and none of them is given: a
 * page that wants the other two when the density refuses a variant past the
 * length of its chromosome calls them on their own. pyNei has no such
 * function, and neither has the Python package of popnei.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * the options are not an object or hold a key that is none of the five;
 * when none of `perVar`, `perIndividual` and `density` is given; when
 * `perVar` is not an object, holds `onSoFar`, `soFarEvery` or a key that is
 * no option of `calcPerVarDistribs`, or holds what `calcPerVarDistribs`
 * refuses; when `perIndividual` is not an object or
 * holds a key; when `density` is not an object, holds a key that is neither
 * `windowSize` nor `chromLengths`, or holds what `calcVarDensity` refuses;
 * when `onSoFar` is not a function, when `soFarEvery` is not a finite
 * number of 0 or more and when it is given without `onSoFar`; what the pass
 * of any of the three refuses, which `calcPerVarDistribs`,
 * `calcPerIndividualStats` and `calcVarDensity` name; when the pass gives
 * no variant; and when `init` has not been awaited. It throws what
 * `onSoFar` threw.
 */
export function calcVariantsSummary(
  variants: Variants,
  options: VariantsSummaryOptions = {},
): VariantsSummary {
  theWasmHasToBeLoaded();
  anObjectOfOptions("calcVariantsSummary", options, [
    "perVar",
    "perIndividual",
    "density",
    ...SO_FAR_KEYS,
  ]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const { perVar, perIndividual, density } = options;
  if (
    perVar === undefined &&
    perIndividual === undefined &&
    density === undefined
  ) {
    throw new Error(
      "popnei: `calcVariantsSummary` was asked for none of its three " +
        "statistics: give `perVar: {}` for the distributions of the " +
        "statistics of each variant, `perIndividual: {}` for the rates of " +
        "each individual, or `density: {windowSize}` for the density of the " +
        "variants",
    );
  }
  let perVarAsked: PerVarArguments | undefined;
  if (perVar !== undefined) {
    anObjectOfOptions("calcVariantsSummary.perVar", perVar, PER_VAR_KEYS);
    perVarAsked = thePerVarArguments(perVar);
  }
  if (perIndividual !== undefined) {
    anObjectOfOptions("calcVariantsSummary.perIndividual", perIndividual, []);
  }
  let densityAsked: ReturnType<typeof theDensityArguments> | undefined;
  if (density !== undefined) {
    anObjectOfOptions("calcVariantsSummary.density", density, [
      "windowSize",
      "chromLengths",
    ]);
    densityAsked = theDensityArguments(
      "density.windowSize",
      density.windowSize,
      density.chromLengths,
    );
  }
  const soFar = theResultSoFar(options, summaryOf);
  // The steps of the pass and the arguments of each statistic are objects
  // of the binding crate, made after every argument was checked so that
  // nothing refused here leaves one behind: the call takes them over and
  // frees them. The arguments are made before the steps, and the first is
  // freed when the second cannot be made, so that neither leaves the other
  // behind.
  const summary = whileTheRunReads(() => {
    const perVarOfThePass =
      perVarAsked === undefined ? undefined : argumentsOfThePassOf(perVarAsked);
    let densityOfThePass: ArgumentsOfTheDensity | undefined;
    try {
      densityOfThePass =
        densityAsked === undefined
          ? undefined
          : argumentsOfTheDensityOf(densityAsked);
    } catch (thrown: unknown) {
      perVarOfThePass?.free();
      throw thrown;
    }
    return source.calc_variants_summary(
      steps.of_a_pass(),
      perVarOfThePass,
      perIndividual !== undefined,
      densityOfThePass,
      soFar.told,
      soFar.every,
    );
  });
  return summaryOf(summary);
}

/**
 * The three statistics of a pass, or of its first blocks, out of what the
 * binding crate gives, which is freed here with each of the three.
 */
function summaryOf(summary: VariantsSummaryOfAPass): VariantsSummary {
  // Each of the three is taken out and built as its own call builds it,
  // which frees it; one that a throw leaves inside is freed with the rest.
  try {
    const perVar = summary.take_per_var();
    const perVarBuilt = perVar === undefined ? null : distribsOf(perVar);
    const perIndividual = summary.take_per_individual();
    const perIndividualBuilt =
      perIndividual === undefined ? null : ratesOf(perIndividual);
    const density = summary.take_density();
    const densityBuilt = density === undefined ? null : windowsOf(density);
    return {
      perVar: perVarBuilt,
      perIndividual: perIndividualBuilt,
      density: densityBuilt,
      passStats: passStatsOf(summary.pass_stats()),
    };
  } finally {
    summary.free();
  }
}
