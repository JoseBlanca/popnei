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
 *
 * `doPcoaFromVariants` is the Kosman distances of `calcPairwiseKosmanDists`
 * followed by `doPcoa`, in one pass over the variants of a `Variants` and
 * with no `Distances` built between the two, and it applies Lingoes'
 * correction inside when it is asked to.
 */

import {
  correct_dists_by_lingoes as correctDistsByLingoesOfTheCore,
  pcoa as pcoaOfTheCore,
  room_for_the_principal_coordinates_of as roomForThePrincipalCoordinatesOf,
} from "../wasm/popnei.js";

import {
  aBoolean,
  anObjectOfOptions,
  whatWasGiven,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Distances } from "./dists.js";
import { theValuesOf } from "./pca.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

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
 * are about 56.8 bytes per cell of the individuals x individuals matrix, so
 * distances of more than 8695 individuals are an `Error` here and are
 * analysed by a program outside the browser, popnei in Python among them.
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
 * constant, which is in the units of a squared distance, is beyond the range
 * of a 64 bit float, above 1.8e308 or below the smallest normal one,
 * 2.2e-308, which distances of about 1.3e154 and above or all below about
 * 1e-154 reach and whose message says to divide them by a number near the largest first; when
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
 * What the principal coordinates of the variants of a dataset give.
 *
 * The components are those of the positive eigenvalues, as in
 * `PcoaResult`, and the fields that the result of `doPcaFromVariants` has
 * for drawing the individuals have the same names and shapes here, so that
 * one piece of code draws either. A PCoA has no weights of variants, so
 * there is no `numPrinComps`, `princomps` or `usedVars`.
 */
export interface VariantsPcoaResult {
  /**
   * The names of the individuals, in the order the pass gives them, which is
   * the order of the rows of `projections`: those of the source, or the ones
   * `filterIndividuals` kept, in the order it named them.
   */
  readonly individuals: readonly string[];
  /** How many components there are. */
  readonly numComps: number;
  /**
   * Where each individual falls along each component, the individuals x
   * `numComps`, row after row.
   */
  readonly projections: Float64Array;
  /**
   * The variance each component holds as a percentage of the sum of every
   * eigenvalue, of the corrected distances when there was a correction, one
   * number per component. They add up to 100.
   */
  readonly explainedVariancePercent: Float64Array;
  /**
   * c, the constant of Lingoes' correction, in the units of a squared
   * distance: 0 when `correctByLingoes` was false or the distances were
   * Euclidean.
   */
  readonly lingoesConstant: number;
  /**
   * 100 times the sum of the absolute values of the negative eigenvalues
   * over the sum of every eigenvalue, of the distances before the
   * correction: how much of their variance was in directions that no space
   * has. 0 when `correctByLingoes` was false, since distances with a negative
   * eigenvalue are then refused.
   */
  readonly negativeEigenvaluesPercent: number;
  /**
   * How many variants the pass gave, called in a pair or not, and what each
   * filter of the `Variants` was given and kept.
   */
  readonly passStats: PassStats;
}

/** How the principal coordinates of the variants of a dataset are taken. */
export interface DoPcoaFromVariantsOptions {
  /**
   * How many variants a pair of individuals needs to be called together at
   * before it gets a distance, a whole number of 0 or more: a pair called
   * together at fewer, or at none, has no distance, which is an `Error`
   * here. 0 when it is not given, which is every pair called together at
   * one variant at least. It is the `minNumSnps` of
   * `calcPairwiseKosmanDists`.
   */
  minNumSnps?: number;
  /**
   * Whether distances that are not Euclidean are corrected by Lingoes'
   * correction inside the analysis, as `correctDistsByLingoes` corrects a
   * `Distances`. False when it is not given, and such distances are then an
   * `Error`.
   */
  correctByLingoes?: boolean;
}

/**
 * The principal coordinates of the Kosman distances of the individuals of
 * `variants`, after its steps.
 *
 * It is `doPcoa` of the distances of `calcPairwiseKosmanDists`, made in one
 * pass over the source, through the steps the `Variants` has when it is
 * called, and with no vector of distances built: the matrix that is
 * decomposed is built from what the pass counted for each pair. Kosman
 * distances with missing genotypes are often not Euclidean, and then the
 * analysis is refused, unless `correctByLingoes` is true: the correction is
 * then applied inside, as `correctDistsByLingoes` would apply it to the
 * distances, and what it added and took is in `lingoesConstant` and
 * `negativeEigenvaluesPercent`. It is `doPcoa` of the corrected distances
 * within the rounding of their square roots.
 *
 * It is pyNei's `do_pcoa_from_variants`, which has no correction and gives a
 * component for every eigenvalue, the negative ones included.
 *
 * A page holds 4 GB at a time, and the analysis holds about 56.8 bytes per
 * cell of the individuals x individuals matrix, so a pass of more than 8695
 * individuals is an `Error` here, before the source is read, and is
 * analysed by a program outside the browser, popnei in Python among them.
 * The individuals counted are those the pass gives, so a `filterIndividuals`
 * that keeps fewer is analysed.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when the
 * options are not an object or hold a key that is neither `minNumSnps` nor
 * `correctByLingoes`; when `minNumSnps` is not a whole number from 0 to
 * 4294967295; when `correctByLingoes` is not a boolean; when the pass gives
 * more than 8695 individuals; when it gives fewer than 2; when the pass
 * gives no variant, whose message says whether the source held none or the
 * steps kept none and how many variants each filter was given and kept;
 * when the source cannot be read, a wrong line of a VCF among the causes;
 * when the counts of a pair go above 4294967295; when the memory the page
 * has left does not take the two counts popnei keeps for every pair, 8
 * bytes a pair, or the matrix of the analysis; when a pair has no
 * distance, whose message says how many do not, names the first of them and
 * the individual that is in the most of them, and says to take that
 * individual out with `filterIndividuals`, lower `minNumSnps` or run the PCA
 * of the variants; when every distance is 0; when the distances are not
 * Euclidean and `correctByLingoes` is false, whose message says how many
 * eigenvalues are negative and their share, and names `correctByLingoes`;
 * when the linear algebra of the analysis could not be done; when the
 * eigenvalue of 0 that the centering of the matrix always gives is not
 * found, which is a defect of popnei; and when `init` has not been awaited.
 */
export function doPcoaFromVariants(
  variants: Variants,
  options: DoPcoaFromVariantsOptions = {},
): VariantsPcoaResult {
  theWasmHasToBeLoaded();
  anObjectOfOptions("doPcoaFromVariants", options, [
    "minNumSnps",
    "correctByLingoes",
  ]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const minNumSnps =
    options.minNumSnps === undefined
      ? 0
      : wholeNumberOfZeroOrMore("minNumSnps", options.minNumSnps);
  const correctByLingoes =
    options.correctByLingoes === undefined
      ? false
      : aBoolean("correctByLingoes", options.correctByLingoes);
  // The names are read before the call: the pass gives the individuals the
  // steps say, and nothing of the run changes them.
  const individuals = variants.individuals;
  // The steps of the pass are a copy of the list, made after the arguments
  // were checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const result = whileTheRunReads(() =>
    source.pcoa_of_variants(minNumSnps, correctByLingoes, steps.of_a_pass()),
  );
  try {
    return {
      individuals,
      numComps: result.num_comps(),
      projections: theValuesOf(result.projections(), "projections"),
      explainedVariancePercent: theValuesOf(
        result.explained_variance_percent(),
        "explainedVariancePercent",
      ),
      lingoesConstant: result.lingoes_constant(),
      negativeEigenvaluesPercent: result.negative_eigenvalues_percent(),
      passStats: passStatsOf(result.pass_stats()),
    };
  } finally {
    result.free();
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
