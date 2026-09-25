/**
 * How strongly two variants go together: the matrix of r² of every pair of
 * the variants of a dataset, and how that r² falls off as the two variants
 * of a pair move apart along a chromosome.
 *
 * Two variants are in linkage disequilibrium when the genotype of one tells
 * something about the genotype of the other, which happens when they sit
 * close enough on a chromosome that few recombinations have separated them.
 * r² measures it: each variant becomes one number per individual, its
 * dosage, how many alleles of the genotype are not the major allele of that
 * variant, and r² is the square of the correlation between the two dosage
 * vectors, 0 when one variant says nothing about the other and 1 when it
 * says everything. `docs/specs/ld.md` has what is computed.
 *
 * `calcRogersHuffR2Matrix` gives that r² for every pair of the variants of
 * a dataset. `calcLdAndDistPerPop` gives, for each population on its own,
 * the mean r² of the pairs of each bin of distance, which is the fall-off:
 * two variants that sit close together have had fewer recombinations
 * between them than two that sit far apart, and how fast r² falls is a
 * property of the population. Beside those bins it gives the curve fitted
 * to the pairs of each population, which carries the distance at which r²
 * has fallen to half.
 */

import {
  default_max_allowed_maf as defaultMaxAllowedMaf,
  default_max_dist as defaultMaxDist,
  default_max_num_vars as defaultMaxNumVars,
  default_min_dist as defaultMinDist,
  default_num_dist_bins as defaultNumDistBins,
} from "../wasm/popnei.js";

import {
  aNumber,
  distanceInBasePairs,
  popsOfTheObject,
  varsOfTheMatrixOfEveryPair,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/** How the matrix of r² is calculated: the variants it is taken of. */
export interface CalcRogersHuffR2MatrixOptions {
  /**
   * How many variants the calculation takes before it refuses, a whole
   * number from 1 to 65535. 5000 when it is not given, which is the default
   * of the core crate and 200 MB of matrix.
   *
   * The matrix holds one r² for each pair of the variants of the pass, 8
   * bytes each, so it grows with the square of them: 5000 variants are 200
   * MB, 10000 are 800 MB and 23170 are the 4 GB a page holds of everything
   * that is open in it at a time. It is the one calculation of popnei whose
   * result grows with the square of its input, so a pass of more variants
   * than this is an `Error` and not a matrix the tab is asked for the
   * memory of. A user who wants the matrix of more variants, and has the
   * memory, raises it; a user who has more variants than memory puts a
   * filter on the `Variants` first.
   *
   * 65535 is where it stops in a browser, and a cap above it is an `Error`
   * at the call: the values of the matrix are counted in a whole number of
   * 32 bits there, and 65536 variants hold more of them than it counts. The
   * memory of the tab is reached at 23170, so the number a user can write
   * here is not one a browser runs. Python, whose whole numbers are 64
   * bits wide, takes a `max_num_vars` of up to 4294967295.
   */
  maxNumVars?: number;
}

/**
 * The r² of every pair of the variants of a dataset, with the chromosome
 * and the position of each of them and the counts of the pass.
 *
 * It is the `R2Matrix` of the Python package with the names of TypeScript:
 * a `Float64Array` where Python has a numpy array, and an array of strings
 * where Python has a tuple.
 */
export interface R2Matrix {
  /** How many variants the pass gave, which is the side of the matrix. */
  readonly numVars: number;
  /**
   * The r² of every pair, `numVars` x `numVars` row after row: the r² of
   * the variants `i` and `j` is the value at `i * numVars + j`, and the
   * value at `j * numVars + i` is the same one.
   *
   * A pair that has no r² is NaN, which is a pair holding a variant whose
   * called genotypes have one dosage only, and the diagonal of a variant
   * that has two dosages at least is 1. The individuals whose genotype is
   * missing at either variant of a pair are left out of that pair, so each
   * pair has its own number of individuals.
   */
  readonly r2: Float64Array;
  /**
   * The name of the chromosome of each variant, in the order the pass gave
   * them, which is the order of the rows and of the columns of the matrix.
   *
   * The array is frozen, as the names of a `Distances` are: writing into
   * it throws, and a user who wants one of their own copies it.
   */
  readonly chroms: readonly string[];
  /**
   * The position of each variant, 1 based as in a VCF. The distance of a
   * pair is the difference of two of them, and a pair whose variants are on
   * two chromosomes has no distance and an r² like any other pair.
   */
  readonly poss: Float64Array;
  /**
   * How many variants the calculation took, after the steps, and how many
   * variants each filter of the pass was given and kept.
   */
  readonly passStats: PassStats;
}

/**
 * The r² of every pair of the variants of `variants`.
 *
 * Each variant becomes one number per individual, its dosage: how many
 * alleles of the genotype are not the major allele of that variant, which
 * is the most frequent among its called alleles and the lowest numbered of
 * two that are equally frequent. r² is the square of the correlation
 * between the dosages of the two variants of a pair, over the individuals
 * whose genotype is missing at neither of them, so each pair has its own
 * number of individuals; a variant whose called genotypes hold one dosage
 * only has no correlation with anything and is NaN in its whole row, in its
 * whole column and in its diagonal cell.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, so the matrix
 * is of the variants the filters kept, and the `variants` are as they were
 * afterwards.
 *
 * It is `calc_rogers_huff_r2_matrix` of the Python package, which it
 * mirrors, and pyNei's function of that name with two differences a caller
 * sees: pyNei gives r and this gives its square, which loses the sign and
 * nothing else, and pyNei leaves an individual with a missing genotype in
 * the pair with a dosage of -1 where this leaves it out, which is what
 * plink2 does and what makes plink2 the reference program of this
 * calculation.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `maxNumVars` is not a whole number from 1 to 65535, when the pass gives
 * more variants than `maxNumVars`, whose message has both numbers and the
 * memory the matrix would have needed, when the pass gives no variant, whose
 * message says whether the source held none or the steps kept none and how
 * many variants each filter was given and kept, when the source cannot be
 * read, a wrong line of a VCF among the causes, when a position of the
 * source is above 2^53, which a number of JavaScript rounds, when the
 * memory of the tab does not take the matrix, 8 bytes a pair, and when
 * `init` has not been awaited.
 */
export function calcRogersHuffR2Matrix(
  variants: Variants,
  options: CalcRogersHuffR2MatrixOptions = {},
): R2Matrix {
  theWasmHasToBeLoaded();
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  // The default is the core's, as the ploidy of `openVcf` is, so that
  // Python and TypeScript cannot drift apart on how many variants the
  // matrix is taken of when the caller says nothing.
  const maxNumVars =
    options.maxNumVars === undefined
      ? defaultMaxNumVars()
      : varsOfTheMatrixOfEveryPair("maxNumVars", options.maxNumVars);
  // The steps of the pass are a copy of the list, made after the arguments
  // were checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = whileTheRunReads(() =>
    source.calc_rogers_huff_r2_matrix(maxNumVars, steps.of_a_pass()),
  );
  try {
    // Each array is moved out of the result as it is read, and not cloned.
    // The generated code still copies the values into an array of the
    // JavaScript heap and frees the memory of wasm after it, so the matrix
    // is held twice while that copy is made, 400 MB at 5000 variants, and
    // once afterwards.
    return {
      numVars: calculated.num_vars(),
      r2: thePartOfTheResult(calculated.r2(), "r2"),
      // The names are frozen, as the names of a `Distances` and the
      // individuals of a `Variants` are: what Python gives here is a tuple,
      // and a user who wrote into the array would be writing into the
      // result of a pass that is over.
      chroms: Object.freeze(thePartOfTheResult(calculated.chroms(), "chroms")),
      poss: thePartOfTheResult(calculated.poss(), "poss"),
      passStats: passStatsOf(calculated.pass_stats()),
    };
  } finally {
    calculated.free();
  }
}

/**
 * One array of a result of this module on its way out of WebAssembly, the
 * matrix of r² or the bins of the fall-off.
 *
 * Each of them leaves the memory of wasm the first time it is asked for, so
 * the call after that gives nothing; the two calls read each of them once,
 * and the `Error` here is a defect of this package and not something a
 * caller can do. It is `theValuesOf` of `pca.ts` for the arrays of these
 * results, one of which is an array of names.
 *
 * @throws {Error} When the array had been read already.
 */
function thePartOfTheResult<Part>(part: Part | undefined, name: string): Part {
  if (part === undefined) {
    throw new Error(
      `popnei: \`${name}\` was read twice out of the memory of WebAssembly, ` +
        "which is a defect of popnei; please report it",
    );
  }
  return part;
}

/** How the fall-off of r² with distance is calculated, and over whom. */
export interface CalcLdAndDistPerPopOptions {
  /**
   * The populations the fall-off is read for, an object of the name of a
   * population to the names of its individuals, which are looked up among
   * the individuals the pass gives.
   *
   * With no `pops` there is one population, `pop`, of every individual,
   * which is what pyNei names it. An individual may be in more than one
   * population, and one in none is read by none of them.
   *
   * The curve of each population is read on its own because the
   * recombination it has had and the number of individuals it has been
   * through shape it: a population that went through few individuals keeps
   * linkage disequilibrium over longer stretches.
   */
  pops?: Record<string, readonly string[]>;

  /**
   * The smallest distance in base pairs at which a pair of variants is
   * counted, that distance included, 1 when it is not given.
   *
   * At 1 the only pairs left out are those of two variants at one position,
   * a SNP and an indel at the same base, whose distance is 0.
   */
  minDist?: number;

  /**
   * The largest distance in base pairs at which a pair is counted, that
   * distance included, 1000000 when it is not given, which is plink2's own
   * `--ld-window-kb 1000`.
   *
   * It is also how far back the pass holds the variants it has read, so it
   * is what the memory of the call grows with.
   */
  maxDist?: number;

  /**
   * How many bins of equal width the distances from `minDist` to `maxDist`
   * are cut into, 50 when it is not given.
   *
   * The width is (`maxDist` − `minDist` + 1) / `numBins`, and a pair at the
   * distance d falls in the bin that floor((d − `minDist`) / width) gives,
   * the last bin taking anything the rounding would put past it.
   */
  numBins?: number;

  /**
   * The largest major allele frequency a variant has in a population and is
   * still counted there, both ends included, 0.95 when it is not given.
   *
   * The major allele frequency is the count of the commonest allele over
   * the called alleles, worked out over the individuals of that population
   * alone, so two populations of one pass count different variants; a
   * variant that no individual of a population called has none there and is
   * left out of it, and can still be counted in another population. Those
   * variants are left out because the r² of a variant that hardly varies
   * rests on the one or two individuals that carry the rare allele, and
   * keeping them raises the curve everywhere.
   */
  maxAllowedMaf?: number;
}

/**
 * The bins of distance of one population: one value in each of the five
 * arrays for each bin, in the order of the distances.
 *
 * It is the pandas frame the Python package gives for a population, with a
 * `Float64Array` for each of its columns. The distances and the counts are
 * `Float64Array` too, as the positions of a block are: a number of
 * JavaScript is a float64 and holds every whole number up to 2^53.
 */
export interface LdBins {
  /** The smallest distance of each bin in base pairs, that one included. */
  readonly smallestDist: Float64Array;
  /** The largest distance of each bin, that one included. */
  readonly largestDist: Float64Array;
  /** How many pairs of variants fell in each bin. */
  readonly numPairs: Float64Array;
  /**
   * The mean of the r² of the pairs of each bin, and NaN for a bin with no
   * pair. A bin of one pair has that pair's r².
   */
  readonly meanR2: Float64Array;
  /**
   * The standard deviation of those r², with the pairs of the bin as the
   * divisor, and NaN for a bin with no pair. A bin of one pair has 0.
   */
  readonly sdR2: Float64Array;
}

/**
 * The curve of r² against distance fitted to the pairs of one population,
 * and the distance at which it has fallen to half.
 *
 * The curve is the r² that two variants of a population are expected to
 * have at a given recombination between them, under drift and
 * recombination: Hill and Weir (1988), with the correction of Weir and Hill
 * (1986) for r² being measured on a sample of individuals and not on the
 * whole population, which is what holds the curve up at long distances. It
 * has one number to fit, `rhoPerBp`, and the number of individuals of the
 * population enters it as it is given.
 *
 * The curve is fitted to every pair the pass counted, each at its own
 * distance, so the bins do not move it and two runs over one dataset give
 * the same three numbers.
 *
 * A population whose pairs no curve was fitted to has NaN in all three, and
 * that is no error. It happens to a population whose pairs fall at fewer
 * than two distances, one distance saying nothing about a fall-off, which
 * includes a population with no pair at all and one left with two variants;
 * and to a population whose best fit falls at an end of the range of the ρ
 * per base pair that is searched, 1e-12 to 100, a curve flat across
 * `maxDist` or fallen before the second base pair being no fall-off that
 * these pairs pin down.
 *
 * It is the `LdDecay` of the Python package with the names of TypeScript.
 */
export interface LdDecay {
  /**
   * The fitted 4Nr: four times the effective size of the population times
   * the recombination per base pair, which is by how much ρ, the scaled
   * recombination of the curve, grows with each base pair between the two
   * variants of a pair.
   *
   * The effective size and the recombination rate enter the curve only as
   * that product, and one pass over one dataset does not separate them, so
   * neither is given on its own.
   */
  readonly rhoPerBp: number;

  /**
   * The fitted curve at a distance of 0, which is its own ceiling: two
   * variants that never recombine still do not reach an r² of 1, because
   * their allele frequencies drift apart.
   *
   * How many individuals the population has fixes it on its own,
   * 0.46198347107438015 at 100 of them, and no pair of the dataset moves
   * it.
   */
  readonly r2AtZero: number;

  /**
   * The distance in base pairs at which the fitted curve has fallen to half
   * of `r2AtZero`.
   *
   * It is read off the curve and not off the pairs, so it can fall beyond
   * every distance the pass counted a pair at, and a plot that marks it on
   * the bins may have to run past them. Three individuals of the dataset
   * `docs/specs/ld.md` verifies against, taken as one population at the
   * defaults, give 1868334.8 base pairs where `maxDist` is 1000000 and
   * their furthest pair is under 260000.
   *
   * What is halved is the curve at a distance of 0 and not the mean r² of
   * the shortest bin, so `numBins` does not move this distance; `minDist`
   * and `maxDist` do, through `rhoPerBp`, since they choose which pairs the
   * curve is fitted to.
   *
   * It is NaN when the other two are, and NaN on its own for a curve that
   * never falls to half of its value at 0, which happens below three
   * individuals: what the curve falls towards as ρ grows is 1 over the
   * individuals, which is above that half at one and at two of them. No
   * pass of `calcLdAndDistPerPop` reaches such a population, since one
   * individual has no variant with variance and so no pair, and two give
   * every pair an r² of 1, which the three NaN above cover.
   */
  readonly halfDist: number;
}

/**
 * How the r² of a pair of variants falls off with the distance between
 * them, for each population: in bins of distance, and as a curve fitted to
 * its pairs with the distance at which that curve has fallen to half.
 *
 * It is the `LdAndDistPerPop` of the Python package with the names of
 * TypeScript: five `Float64Array` where Python has the columns of a pandas
 * frame.
 */
export interface LdAndDistPerPop {
  /**
   * The bins of each population, under its name and in the order the keys
   * of the `pops` that was given iterate in.
   */
  readonly perPop: Record<string, LdBins>;

  /**
   * How many variants each population kept at its major allele frequency,
   * under its name.
   *
   * It is worked out over the individuals of that population alone, so two
   * populations of one pass count different variants, and a variant that no
   * individual of a population has called is out of it.
   */
  readonly numVarsPerPop: Record<string, number>;

  /**
   * The curve fitted to the pairs of each population, under its name and in
   * the same order.
   *
   * The `LdDecay` of a population holds the fitted ρ per base pair, the
   * curve at a distance of 0 and the distance at which it has fallen to
   * half of that, as three numbers and not as arrays: one curve is fitted
   * to every pair of the population and not one to each bin. The three are
   * NaN for a population whose pairs no curve was fitted to.
   */
  readonly decayPerPop: Record<string, LdDecay>;

  /**
   * How many variants the pass gave, before the major allele frequency of
   * any population, and how many variants each filter of it was given and
   * kept.
   */
  readonly passStats: PassStats;
}

/**
 * How the r² of a pair of variants falls off as the two move apart along a
 * chromosome, for each population of `options.pops`.
 *
 * A pair of variants is counted in a population when both of its variants
 * passed the major allele frequency of that population, when the two are on
 * one chromosome, and when their distance, the difference of their
 * positions, is from `minDist` to `maxDist`, both included. The pairs are
 * put into `numBins` bins of equal width across that range, and each bin
 * gets how many pairs it holds, the mean of their r² and its standard
 * deviation. A pair that has no r², one whose individuals called at both
 * variants hold one dosage at one of them, is in no bin, and so is a pair
 * whose two variants are on two chromosomes.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, which serves
 * every population, and the `variants` are as they were afterwards.
 *
 * Beside the bins each population gets a curve fitted to its pairs, in the
 * `decayPerPop` of the result: an `LdDecay` with the fitted ρ per base
 * pair, the curve at a distance of 0 and the distance at which it has
 * fallen to half of that. The curve is fitted to every pair and at the
 * distance of each pair, so `numBins` does not move it, and a population
 * whose pairs no curve was fitted to has NaN in all three, which `LdDecay`
 * says when and which is no error.
 *
 * A population in which every variant was left out, and a dataset whose
 * variants are all further apart than `maxDist` or each on a chromosome of
 * their own, give every bin empty, which is no error: those bins hold no
 * pair and have NaN for their mean and their standard deviation.
 *
 * It is `calc_ld_and_dist_per_pop` of the Python package, which it mirrors,
 * and pyNei's function of that name with these differences: it gives the
 * bins over every pair, where pyNei gives a sample of at most
 * `max_num_measures_to_keep` pairs drawn with no seed, so that two runs
 * over one dataset give the same numbers here and different points there;
 * it fits the curve and gives the half distance, where pyNei hands its
 * sample of pairs over and leaves the fitting to the user;
 * it gives r² where pyNei gives r, and a missing genotype takes its
 * individual out of that pair where pyNei leaves it in with a dosage of -1;
 * `minDist` counts the pair at that distance, where pyNei keeps the pairs
 * strictly beyond it; `maxDist` has a default where pyNei's is `None`; the
 * bins and the standard deviation are new; there is no `method`, pyNei's
 * two counting each unordered pair once and twice; and it makes one pass
 * that serves every population, where pyNei makes one per population.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * `minDist` or `maxDist` is not a whole number of base pairs of 0 or more,
 * and when `minDist` is above `maxDist`; when `numBins` is not a whole
 * number of 0 or more, and when it is 0; when `maxAllowedMaf` is not a
 * number, and when it is not one from 0 to 1; when `pops` is not an object
 * of names to arrays of names, when a population names an individual the
 * pass does not give, names one twice or names none, and when it holds no
 * population; when the source cannot be read, a wrong line of a VCF among
 * the causes; when the pass gives no variant, whether the source holds none
 * or the steps kept none; when the memory of the tab does not take the bins
 * or the variants the pass holds within `maxDist`; and when `init` has not
 * been awaited.
 */
export function calcLdAndDistPerPop(
  variants: Variants,
  options: CalcLdAndDistPerPopOptions = {},
): LdAndDistPerPop {
  theWasmHasToBeLoaded();
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  // The populations cross flat, the names of the individuals of every one
  // of them in one array, and no name at all is one population of every
  // individual, which the core names `pop`.
  const pops =
    options.pops === undefined
      ? {
          names: undefined,
          individuals: [] as string[],
          numIndividualsPerPop: new Uint32Array(0),
        }
      : popsOfTheObject(options.pops);
  // The four defaults are the core's, as the cap of the matrix is, so that
  // Python and TypeScript cannot drift apart on what a call that says
  // nothing counts.
  const minDist =
    options.minDist === undefined
      ? defaultMinDist()
      : distanceInBasePairs("minDist", options.minDist, 0);
  const maxDist =
    options.maxDist === undefined
      ? defaultMaxDist()
      : distanceInBasePairs("maxDist", options.maxDist, 0);
  const numBins =
    options.numBins === undefined
      ? defaultNumDistBins()
      : wholeNumberOfZeroOrMore("numBins", options.numBins);
  // Whether the frequency is one from 0 to 1 is the core's rule, which it
  // holds for every pass: what is refused here is what is no number at all,
  // which would reach the core as a NaN and take every variant out.
  const maxAllowedMaf =
    options.maxAllowedMaf === undefined
      ? defaultMaxAllowedMaf()
      : aNumber("maxAllowedMaf", options.maxAllowedMaf);
  // The steps of the pass are a copy of the list, made after every argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = whileTheRunReads(() =>
    source.calc_ld_and_dist_per_pop(
      steps.of_a_pass(),
      pops.names,
      pops.individuals,
      pops.numIndividualsPerPop,
      minDist,
      maxDist,
      numBins,
      maxAllowedMaf,
    ),
  );
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the
  // user gets are the copies.
  try {
    const popNames = thePartOfTheResult(calculated.pop_names(), "popNames");
    const binsOfEveryPop = calculated.num_bins();
    const smallestDist = thePartOfTheResult(
      calculated.smallest_dist(),
      "smallestDist",
    );
    const largestDist = thePartOfTheResult(
      calculated.largest_dist(),
      "largestDist",
    );
    const numPairs = thePartOfTheResult(calculated.num_pairs(), "numPairs");
    const meanR2 = thePartOfTheResult(calculated.mean_r2(), "meanR2");
    const sdR2 = thePartOfTheResult(calculated.sd_r2(), "sdR2");
    const varsOfEachPop = thePartOfTheResult(
      calculated.num_vars_per_pop(),
      "numVarsPerPop",
    );
    // The curve of a population is three numbers and not three arrays, so
    // each of these three arrays holds one value for each population and is
    // read at the index of the population, where the arrays of the bins are
    // cut into `binsOfEveryPop` values each.
    const rhoPerBp = thePartOfTheResult(calculated.rho_per_bp(), "rhoPerBp");
    const r2AtZero = thePartOfTheResult(calculated.r2_at_zero(), "r2AtZero");
    const halfDist = thePartOfTheResult(calculated.half_dist(), "halfDist");
    const perPop: Record<string, LdBins> = {};
    const numVarsPerPop: Record<string, number> = {};
    const decayPerPop: Record<string, LdDecay> = {};
    const numPops = popNames.length;
    for (const [which, pop] of popNames.entries()) {
      // The bins of one population lie together in each array, the bins of
      // the population before it first, as the values of one measure of
      // `calcPopDists` do.
      const first = which * binsOfEveryPop;
      const pastTheLast = first + binsOfEveryPop;
      perPop[pop] = {
        smallestDist: smallestDist.subarray(first, pastTheLast),
        largestDist: largestDist.subarray(first, pastTheLast),
        numPairs: numPairs.subarray(first, pastTheLast),
        meanR2: meanR2.subarray(first, pastTheLast),
        sdR2: sdR2.subarray(first, pastTheLast),
      };
      numVarsPerPop[pop] = theValueOfThePop(
        varsOfEachPop,
        which,
        numPops,
        "how many variants it kept",
      );
      decayPerPop[pop] = {
        rhoPerBp: theValueOfThePop(
          rhoPerBp,
          which,
          numPops,
          "the ρ per base pair",
        ),
        r2AtZero: theValueOfThePop(
          r2AtZero,
          which,
          numPops,
          "the r² at a distance of 0",
        ),
        halfDist: theValueOfThePop(
          halfDist,
          which,
          numPops,
          "the half distance",
        ),
      };
    }
    return {
      perPop,
      numVarsPerPop,
      decayPerPop,
      passStats: passStatsOf(calculated.pass_stats()),
    };
  } finally {
    calculated.free();
  }
}

/**
 * The value of the population `which` of an array that holds one for each
 * population of the pass, `what` naming what it is.
 *
 * The variants a population kept and the three values of its curve cross in
 * one array each, one value for each population, and this reads the one of a
 * population out of such an array. Every per-population value of the result
 * that is not an array of bins is read through it, so a value added later
 * is checked the same way as the four there are.
 *
 * @throws {Error} When the array holds no value for that population, which
 * is a defect of popnei and not something a caller can do: the binding
 * crate fills each of the four arrays with one value for each population it
 * counted.
 */
function theValueOfThePop(
  values: Float64Array,
  which: number,
  numPops: number,
  what: string,
): number {
  const value = values[which];
  if (value === undefined) {
    throw new Error(
      `popnei: the pass counted ${numPops} populations and gave ${what} of ` +
        `${values.length} of them, which is a defect of popnei; please ` +
        "report it",
    );
  }
  return value;
}
