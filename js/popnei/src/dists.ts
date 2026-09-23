/**
 * The distances between individuals: the Kosman distance of every pair, and
 * the `Distances` that every distance calculation of popnei gives.
 */

import { wholeNumberOfZeroOrMore } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/** How the Kosman distances are calculated: the variants a pair needs. */
export interface CalcPairwiseKosmanDistsOptions {
  /**
   * How many variants a pair of individuals needs before it gets a
   * distance: a pair that was called at fewer than this many variants has
   * none, and one called at exactly this many keeps its distance. 0 when it
   * is not given, which is every pair that was called at one variant at
   * least.
   *
   * It keeps the name pyNei gives it, although the variants need not be
   * SNPs.
   */
  minNumSnps?: number;
}

/**
 * The distance of every pair of individuals of a dataset, with the names of
 * those individuals and the counts of the pass the distances were
 * calculated over.
 *
 * `distVector` holds one value for each pair in the order (0, 1), (0, 2),
 * ..., (0, N-1), (1, 2), ..., the upper triangle of the square matrix of
 * the distances row by row, and NaN for a pair that has no distance.
 * `squareDists` gives that square matrix.
 *
 * The pairs of the measures of `calcPopDists` are pairs of populations, and
 * everything here reads the same way with the populations in the place of
 * the individuals: their names are the `names`.
 *
 * It is the result of every distance calculation of popnei, as
 * `docs/specs/dists.md` has it, and it is the `Distances` of the Python
 * package with the names of TypeScript: a `Float64Array` where Python has a
 * numpy array, and a square matrix of numbers where Python has a pandas
 * frame indexed by the names.
 */
export class Distances {
  /**
   * The distance of every pair, NaN for a pair with no distance, in the
   * order (0, 1), (0, 2), ..., (1, 2), ....
   */
  readonly distVector: Float64Array;

  /**
   * The names of the individuals, in the order the source has them, or of
   * the populations of `calcPopDists`, in the order they were named in.
   */
  readonly names: readonly string[];

  /**
   * How many variants the calculation took, after the steps, and how many
   * variants each filter of the pass was given and kept.
   */
  readonly passStats: PassStats;

  /**
   * How far each distance would move if the variants it was calculated over
   * were drawn again, one value for each pair in the order of `distVector`,
   * and `null` when the calculation gave none.
   *
   * The Kosman distances between individuals never give one, and
   * `calcPopDists` gives one for each measure only when it was asked to cut
   * the variants into resampling groups. A pair that has a distance and no
   * standard error is NaN.
   */
  readonly standardErrors: Float64Array | null;

  /**
   * The distances of the pairs of `names`, which is what
   * `calcPairwiseKosmanDists` and `calcPopDists` build.
   *
   * @throws {Error} When `distVector`, or `standardErrors` where there are
   * any, does not hold one value for each pair of `names`, which is
   * `names.length * (names.length - 1) / 2` of them.
   */
  constructor(
    distVector: Float64Array,
    names: readonly string[],
    passStats: PassStats,
    standardErrors: Float64Array | null = null,
  ) {
    const numPairs = (names.length * (names.length - 1)) / 2;
    if (distVector.length !== numPairs) {
      throw new Error(
        `popnei: the distances of ${names.length} individuals are ` +
          `${numPairs} pairs, and ${distVector.length} values were given`,
      );
    }
    if (standardErrors !== null && standardErrors.length !== numPairs) {
      throw new Error(
        `popnei: the standard error of each pair goes beside the distance ` +
          `of that pair, and ${numPairs} distances were given with ` +
          `${standardErrors.length} standard errors`,
      );
    }
    this.distVector = distVector;
    this.names = Object.freeze([...names]);
    this.passStats = passStats;
    this.standardErrors = standardErrors;
  }

  /**
   * The distances as the square matrix, N x N values row by row: the
   * distance of the individuals `i` and `j` is the value at `i * N + j`,
   * and the value at `j * N + i` is the same one.
   *
   * The diagonal is 0, also for an individual that has no called genotype,
   * and a pair with no distance is NaN in both of its cells. The array is
   * the caller's own: it is built at every call.
   */
  squareDists(): Float64Array {
    const numIndividuals = this.names.length;
    const square = new Float64Array(numIndividuals * numIndividuals);
    let pair = 0;
    for (let first = 0; first < numIndividuals; first += 1) {
      for (
        let second = first + 1;
        second < numIndividuals;
        second += 1, pair += 1
      ) {
        const dist = this.distVector[pair] as number;
        square[first * numIndividuals + second] = dist;
        square[second * numIndividuals + first] = dist;
      }
    }
    return square;
  }

  /**
   * The standard errors as the square matrix, N x N values row by row, and
   * `null` where `standardErrors` is `null`, so that a caller who asked for
   * none reads the same answer from the method and from the field.
   *
   * The standard error of a pair is in both of its cells and the diagonal is
   * NaN, where `squareDists` has 0: a distance of an individual or a
   * population with itself is 0 and known, and how far that 0 would move is
   * nothing the calculation gives. The array is the caller's own: it is
   * built at every call.
   */
  squareStandardErrors(): Float64Array | null {
    if (this.standardErrors === null) {
      return null;
    }
    const numNames = this.names.length;
    const square = new Float64Array(numNames * numNames).fill(Number.NaN);
    let pair = 0;
    for (let first = 0; first < numNames; first += 1) {
      for (let second = first + 1; second < numNames; second += 1, pair += 1) {
        const error = this.standardErrors[pair] as number;
        square[first * numNames + second] = error;
        square[second * numNames + first] = error;
      }
    }
    return square;
  }
}

/**
 * The Kosman distance of every pair of individuals of `variants`.
 *
 * At one variant, the two called genotypes are laid side by side, each
 * allele of one paired with an allele of the other in the pairing that
 * leaves the fewest pairs of different alleles, and d is that number of
 * pairs divided by the ploidy: for two diploids, 0 when they hold the same
 * alleles, 0.5 when they share one, and 1 when they share none. The
 * distance of the pair is the mean of d over the variants at which both
 * genotypes were called, so it runs from 0 for two individuals with the
 * same genotype everywhere to 1 for two that share no allele anywhere. It
 * is the distance of Kosman and Leonard (2005, Molecular Ecology 14: 415)
 * for codominant markers, at any ploidy, and it is what a tree or a
 * principal coordinate analysis of individuals is built from.
 *
 * A genotype is missing for the pair when one of its alleles at least was
 * not called, `0/.` among them, and such a variant counts for neither of
 * the two individuals at it. Each pair has its own missing genotypes, so
 * each has its own number of variants; `minNumSnps` is how many a pair
 * needs before it gets a distance, and a pair below it, or one with no
 * variant at all, is NaN in the vector.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, so the
 * distances are over the variants the filters kept, and the `variants` are
 * as they were afterwards.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `minNumSnps` is not a whole number of 0 or more, when the pass gives no
 * variant, whose message says whether the source held none or the steps
 * kept none and how many variants each filter was given and kept, when the
 * source cannot be read, a wrong line of a VCF among the
 * causes, when the memory of the tab does not take the two counts popnei
 * keeps for every pair, 8 bytes a pair, and when `init` has not been
 * awaited.
 */
export function calcPairwiseKosmanDists(
  variants: Variants,
  options: CalcPairwiseKosmanDistsOptions = {},
): Distances {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  const minNumSnps =
    options.minNumSnps === undefined
      ? 0
      : wholeNumberOfZeroOrMore("minNumSnps", options.minNumSnps);
  // The steps of the pass are a copy of the list, made after the argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = source.calc_pairwise_kosman_dists(
    minNumSnps,
    steps.of_a_pass(),
  );
  try {
    // The vector and the names are moved out of the result as they are
    // read, and not cloned. The generated code still copies the values into
    // an array of the JavaScript heap and frees the memory of wasm after
    // it, so the vector is held twice while that copy is made, 800 MB at
    // 10000 individuals, and once afterwards.
    const distVector = calculated.dist_vector();
    const names = calculated.names();
    if (distVector === undefined || names === undefined) {
      throw new Error(
        "popnei: the distances of this pass gave no vector or no names",
      );
    }
    return new Distances(
      distVector,
      names,
      passStatsOf(calculated.pass_stats()),
    );
  } finally {
    calculated.free();
  }
}
