/**
 * The principal coordinates of distances, and Lingoes' correction of
 * distances that are not Euclidean.
 *
 * A principal coordinate analysis, PCoA, places the individuals on
 * components as a principal component analysis does, from the distance of
 * every pair of them instead of a table of their values: the individuals are
 * put in a space where the straight line between each two of them is as long
 * as their distance, and the components are the directions of that space
 * along which they vary most. It is Gower's method (1966, Biometrika 53:
 * 325), which R's `cmdscale` and ape's `pcoa` compute, and
 * `docs/specs/pca.md` has what it computes.
 *
 * Distances are Euclidean when some space has points whose straight line
 * distances are those distances. Kosman distances with missing genotypes
 * often are not, and a PCoA of them would draw directions that no space has,
 * so `doPcoa` refuses them. `correctDistsByLingoes` makes them Euclidean by
 * adding one constant to every squared distance, and says how large that
 * constant is and how much of the variance the directions that no space has
 * held, which is what an application warns its user from.
 */

import {
  correct_dists_by_lingoes as correctDistsByLingoesOfTheCore,
  pcoa as pcoaOfTheCore,
  room_for_the_principal_coordinates_of as roomForThePrincipalCoordinatesOf,
} from "../wasm/popnei.js";

import { whatWasGiven } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Distances } from "./dists.js";
import { theValuesOf } from "./pca.js";
import type { PassStats } from "./variant.js";

/**
 * What the principal coordinates of a `Distances` give.
 *
 * Only the components of the positive eigenvalues of the matrix of the
 * squared distances are here: the centering of that matrix always leaves an
 * eigenvalue of 0, so 5 individuals give 4 components at most. In each
 * component the projection of the largest absolute value is positive, which
 * is the rule that makes the numbers the same whichever library did the
 * decomposition, in TypeScript as in Python.
 */
export interface PcoaResult {
  /**
   * The names of the individuals, or of the populations, of the
   * `Distances`, in its order, which is the order of the rows of
   * `projections`.
   */
  readonly names: readonly string[];
  /** How many components there are. */
  readonly numComps: number;
  /**
   * Where each individual falls along each component, the individuals x
   * `numComps`, row after row.
   */
  readonly projections: Float64Array;
  /**
   * The variance each component holds as a percentage of the variance of
   * the individuals placed at their distances, which is the sum of every
   * eigenvalue, one number per component. They add up to 100.
   */
  readonly explainedVariancePercent: Float64Array;
  /**
   * The constant of Lingoes' correction, which is 0: `doPcoa` corrects
   * nothing, and the distances it was given may be corrected ones.
   */
  readonly lingoesConstant: number;
  /**
   * The share of the negative eigenvalues of the distances before a
   * correction, which is 0: `doPcoa` refuses distances that have any.
   */
  readonly negativeEigenvaluesPercent: number;
  /** The counts of the pass of the `Distances`, as it holds them. */
  readonly passStats: PassStats;
}

/** What Lingoes' correction of a `Distances` gives. */
export interface LingoesCorrection {
  /**
   * The corrected distances, sqrt(d² + 2c) for every distance d of the ones
   * given, with their names and their `passStats`, and no standard errors:
   * those of the distances given are not those of the corrected ones. They
   * are the distances given when those were Euclidean.
   */
  readonly distances: Distances;
  /**
   * c, the absolute value of the most negative eigenvalue of the matrix of
   * the squared distances, in the units of a squared distance. 0 when the
   * distances were Euclidean.
   */
  readonly constant: number;
  /**
   * 100 times the sum of the absolute values of the negative eigenvalues
   * over the sum of every eigenvalue, of the distances given: how much of
   * their variance was in directions that no space has. 0 when they were
   * Euclidean.
   */
  readonly negativeEigenvaluesPercent: number;
}

/**
 * The principal coordinates of `distances`: where each individual falls on
 * the components of the space in which the straight line between each two
 * of them is their distance.
 *
 * With d the distance of two individuals, the matrix that is decomposed is
 * that of the -d²/2, centered by its rows and by its columns, and the
 * projections of a component are its eigenvector times the square root of
 * its eigenvalue. It is ape's `pcoa` of R with no correction, and pyNei's
 * `do_pcoa`, which gives a component for every eigenvalue, the negative ones
 * and the 0 of the centering included.
 *
 * Distances that are not Euclidean are refused, and the message names
 * `correctDistsByLingoes`, whose distances this then takes. It corrects
 * nothing itself, so a `lingoesConstant` and a `negativeEigenvaluesPercent`
 * of 0 are in its result, and what a correction added and took is in the
 * result of the correction.
 *
 * A page holds 4 GB at a time, and the matrix of the individuals, its
 * eigenvectors, the workspace of the eigendecomposition and the projections
 * are about 56.8 bytes per pair of individuals, so distances of more than
 * 8695 individuals are an `Error` here and are analysed by a program outside
 * the browser, popnei in Python among them.
 *
 * @throws {Error} When `distances` is not a `Distances`; when it is of more
 * than 8695 individuals; when it is of fewer than 2; when a pair has no
 * distance, whose message says how many do not, names the first of them and
 * the individual that is in the most of them, and says to give each of those
 * pairs a distance or take one of its two individuals out; when a distance
 * is negative or infinite, a negative F_ST or f_2 of two populations the
 * dataset cannot tell apart among them; when every distance is 0; when the
 * distances are not Euclidean, whose message says how many eigenvalues are
 * negative and their share, and names `correctDistsByLingoes`; when the
 * linear algebra of the analysis could not be done; and when `init` has not
 * been awaited.
 */
export function doPcoa(distances: Distances): PcoaResult {
  theWasmHasToBeLoaded();
  const given = aDistances(distances);
  // The page is asked before the vector is copied into the memory of wasm,
  // which is the first thing a vector the page cannot hold the analysis of
  // would not fit in.
  roomForThePrincipalCoordinatesOf(given.names.length);
  const result = pcoaOfTheCore(given.distVector, [...given.names]);
  try {
    return {
      names: given.names,
      numComps: result.num_comps(),
      projections: theValuesOf(result.projections(), "projections"),
      explainedVariancePercent: theValuesOf(
        result.explained_variance_percent(),
        "explainedVariancePercent",
      ),
      lingoesConstant: result.lingoes_constant(),
      negativeEigenvaluesPercent: result.negative_eigenvalues_percent(),
      passStats: given.passStats,
    };
  } finally {
    result.free();
  }
}

/**
 * Lingoes' correction of `distances`, which makes them Euclidean.
 *
 * With c the absolute value of the most negative eigenvalue of the matrix
 * that `doPcoa` decomposes, every distance d becomes sqrt(d² + 2c): the
 * matrix of the corrected distances has the eigenvectors of that one and
 * every eigenvalue but the 0 of the centering c larger, so the most negative
 * becomes 0 and none is below it. It is ape's `pcoa` of R with `correction =
 * "lingoes"` (Lingoes 1971). Distances that are Euclidean give a constant of
 * 0 and themselves back.
 *
 * The correction pushes every pair apart by the same amount of squared
 * distance, which moves the nearest pairs the most, so what it did is in
 * the result: `constant`, and the `negativeEigenvaluesPercent` of the
 * distances given. Which values a user should be warned at is the
 * application's to decide.
 *
 * It holds less of the memory of a page than `doPcoa`, since it writes no
 * projections, and takes the same limit of 8695 individuals.
 *
 * @throws {Error} When `distances` is not a `Distances`; when it is of more
 * than 8695 individuals; when it is of fewer than 2; when a pair has no
 * distance, whose message says how many do not, names the first of them and
 * the individual that is in the most of them, and says to give each of those
 * pairs a distance or take one of its two individuals out; when a distance
 * is negative or infinite, a negative F_ST or f_2 of two populations the
 * dataset cannot tell apart among them; when every distance is 0; when the
 * constant, which is in the units of a squared distance, is beyond what a 64
 * bit float holds, which distances above 1.3e154 or all below 1e-150 reach and
 * whose message says to divide them by a number near the largest first; when
 * the linear algebra of the correction could not be done; and when `init`
 * has not been awaited.
 */
export function correctDistsByLingoes(distances: Distances): LingoesCorrection {
  theWasmHasToBeLoaded();
  const given = aDistances(distances);
  roomForThePrincipalCoordinatesOf(given.names.length);
  const correction = correctDistsByLingoesOfTheCore(given.distVector, [
    ...given.names,
  ]);
  try {
    const distVector = correction.dist_vector();
    if (distVector === undefined) {
      throw new Error(
        "popnei: `distVector` was read twice out of the memory of " +
          "WebAssembly, which is a defect of popnei; please report it",
      );
    }
    return {
      distances: new Distances(distVector, given.names, given.passStats),
      constant: correction.constant(),
      negativeEigenvaluesPercent: correction.negative_eigenvalues_percent(),
    };
  } finally {
    correction.free();
  }
}

/**
 * `value`, when it is a `Distances`.
 *
 * @throws {Error} When it is not, which names what was given: what a user
 * gives instead is usually the vector itself.
 */
function aDistances(value: unknown): Distances {
  if (value instanceof Distances) {
    return value;
  }
  throw new Error(
    `popnei: \`distances\` is a \`Distances\`, and ${whatWasGiven(value)} ` +
      "was given: give it what calcPairwiseKosmanDists or calcPopDists " +
      "gives, or build one, new Distances(distVector, names, passStats)",
  );
}
