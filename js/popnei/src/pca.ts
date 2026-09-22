/**
 * The principal component analysis of a table that the page holds, and of
 * the variants of a dataset.
 *
 * A principal component analysis places the rows of a table, the
 * individuals, on a few axes that hold as much of the variation between them
 * as that many axes can. Each column of the table, a trait, is centered, its
 * mean taken from it, and standardized, divided by its standard deviation,
 * and the components are the directions in the space of the traits along
 * which the rows vary most: the first has the largest variance that any
 * direction has, the second the largest among the directions at a right
 * angle to the first, and so on.
 *
 * `doPca` takes the table the user brings, and `doPcaFromVariants` makes the
 * table out of the variants of a dataset: each variant becomes one number
 * per individual, its dosage, how many alleles of its genotype are not the
 * major allele of that variant. `docs/specs/pca.md` has what each of the two
 * computes. `PcaResult` is what both give, and `VariantsPcaResult` adds the
 * names of the individuals, the variants that were used and the counts of
 * the pass to it.
 */

import {
  default_center_data as defaultCenterData,
  default_num_prin_comps as defaultNumPrinComps,
  default_standardize_data as defaultStandardizeData,
  default_transform_to_biallelic as defaultTransformToBiallelic,
  pca as pcaOfTheCore,
} from "../wasm/popnei.js";

import {
  aBoolean,
  valuesOfATable,
  wholeNumberOfOneOrMore,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * What a principal component analysis gives.
 *
 * Only the components that have variance are here: centering takes one
 * dimension out of the data, so a table of 8 rows and 30 traits has 7
 * components and not 8. In each component the projection of the largest
 * absolute value is positive, which is the rule that makes the numbers the
 * same whichever library did the decomposition, in TypeScript as in Python,
 * and the weights of that component have the sign that rule gave it.
 *
 * The names of the rows and of the traits stay with the application: the
 * rows of `projections` are the rows of the table in the order they were
 * given, and the columns of `princomps` are its traits in their order.
 */
export interface PcaResult {
  /** How many components have variance, which is how many are given. */
  readonly numComps: number;
  /**
   * Where each row falls along each component, the rows of the table x
   * `numComps`, row after row.
   */
  readonly projections: Float64Array;
  /**
   * The variance each component holds as a percentage of the variance of
   * every component of the data, the ones with no variance counted, one
   * number per component.
   */
  readonly explainedVariancePercent: Float64Array;
  /**
   * How many components the weights are given for, which for a table is
   * `numComps`.
   */
  readonly numPrinComps: number;
  /**
   * The weight of each trait in each component, `numPrinComps` x the traits
   * of the table, row after row.
   */
  readonly princomps: Float64Array;
}

/** The two steps on the traits of a table before its components are taken. */
export interface DoPcaOptions {
  /**
   * Whether the mean of each trait is taken from it. True when it is not
   * given. Without it the first component mostly points at the mean of the
   * data.
   */
  centerData?: boolean;
  /**
   * Whether each trait is then divided by its standard deviation, which puts
   * traits measured in different units on one scale; without it the traits
   * with the largest numbers dominate. True when it is not given, and true
   * with `centerData` false is an `Error`: the standard deviation it divides
   * by is the one the trait has once it is centered.
   *
   * The divisor of that standard deviation is the number of rows and not the
   * number of rows less one, which is pyNei's. R's `prcomp` divides by the
   * number of rows less one, so the projections it gives for a standardized
   * table are these times the square root of (n - 1) / n, 0.9975 at 200
   * rows; the percentages and the weights are the same in both, and so are
   * the projections of a table that is not standardized.
   */
  standardizeData?: boolean;
}

/**
 * The principal components of `data`, a table of `numRows` rows of `numCols`
 * values each, row after row, the rows being the individuals and the columns
 * the traits.
 *
 * Every trait is a column of the weights, the one with no variance included,
 * which gets a weight of 0 when the table is not standardized. No value may
 * be missing: a table comes whole.
 *
 * It is pyNei's `do_pca` and popnei's `do_pca` in Python, which take the
 * table as a pandas frame and put the names of its index and of its columns
 * on the result; here the two sides of the table are given as numbers and
 * the names stay with the application.
 *
 * A table whose smaller side is more than 9381 is an `Error` here, and in
 * Python the same table is analysed. What is decomposed is the square of
 * that side, and it, its eigenvectors and the workspace of the
 * eigendecomposition are about 6 times 8 bytes per pair of it, which at 9382
 * is more than the 4 GB a page holds at a time. It is a limit of the browser
 * and not of popnei: WebAssembly addresses 4 GB in one page, and an
 * allocation that does not fit ends the module where an `Error` belongs.
 *
 * @throws {Error} When `data` is not a `Float64Array` of `numRows` times
 * `numCols` values, when a side of the table is not a whole number of 1 or
 * more, when an option is not a boolean, when the smaller side of the table
 * is more than 9381, when a value of the table is not finite, when the table
 * is to be standardized and not centered, when it has fewer than 2 rows or
 * no traits, when it is standardized and a trait has no variance, when no
 * trait of it has any, when the values of a trait are so large or so small
 * that its mean or its standard deviation is not a number the analysis can
 * use, when the linear algebra of the analysis could not be done, and when
 * `init` has not been awaited.
 */
export function doPca(
  data: Float64Array,
  numRows: number,
  numCols: number,
  options: DoPcaOptions = {},
): PcaResult {
  theWasmHasToBeLoaded();
  const rows = wholeNumberOfOneOrMore("numRows", numRows);
  const cols = wholeNumberOfOneOrMore("numCols", numCols);
  const values = valuesOfATable("data", data, rows, cols);
  // The two defaults are the core's, as the ploidy of `openVcf` is, so that
  // Python and TypeScript cannot drift apart on what a table is analysed as
  // when the caller says nothing.
  const centerData =
    options.centerData === undefined
      ? defaultCenterData()
      : aBoolean("centerData", options.centerData);
  const standardizeData =
    options.standardizeData === undefined
      ? defaultStandardizeData()
      : aBoolean("standardizeData", options.standardizeData);
  const result = pcaOfTheCore(values, rows, cols, centerData, standardizeData);
  try {
    return {
      numComps: result.num_comps(),
      projections: theValuesOf(result.projections(), "projections"),
      explainedVariancePercent: theValuesOf(
        result.explained_variance_percent(),
        "explainedVariancePercent",
      ),
      numPrinComps: result.num_prin_comps(),
      princomps: theValuesOf(result.princomps(), "princomps"),
    };
  } finally {
    result.free();
  }
}

/**
 * What the principal components of the variants of a dataset give, which is
 * what a table's analysis gives with the names of the individuals, the
 * variants that were used and the counts of the pass beside it.
 */
export interface VariantsPcaResult extends PcaResult {
  /**
   * The names of the individuals, in the order the source has them, which
   * is the order of the rows of `projections`. They are the names of the
   * rows of the table that a user of `doPca` keeps for themselves: here the
   * source carries them.
   */
  readonly individuals: readonly string[];
  /**
   * Where each individual falls along each component, the individuals x
   * `numComps`, row after row. Every component that has variance is here,
   * whatever `numPrinComps` was.
   */
  readonly projections: Float64Array;
  /**
   * How many components the weights are given for: `numPrinComps` of the
   * call, or the components that have variance when fewer were found, or 0
   * when no weights were asked for.
   */
  readonly numPrinComps: number;
  /**
   * The weight of each variant that was used in each component,
   * `numPrinComps` x `usedVars.length`, row after row. It has no row when no
   * weights were asked for, and `usedVars` names its columns then as well.
   */
  readonly princomps: Float64Array;
  /**
   * The position of each variant that was used, from 0, among the variants
   * the pass gave, the ones with no variance included. They are the columns
   * of `princomps`, and they are given whether or not any weights were
   * asked for.
   */
  readonly usedVars: Uint32Array;
  /**
   * How many variants the pass gave, used or not, and what each filter of
   * the `Variants` was given and kept. The analysis reads the source twice
   * when it is asked for weights, both passes count the same, and these are
   * the counts of the first.
   */
  readonly passStats: PassStats;
}

/** How the principal components of the variants of a dataset are taken. */
export interface DoPcaFromVariantsOptions {
  /**
   * Whether every allele that is not the major one counts the same, which
   * is what gives a variant of more than two alleles a dosage. False when
   * it is not given, and such a variant is then an `Error`.
   */
  transformToBiallelic?: boolean;
  /**
   * How many components the weight of each variant is given for, a whole
   * number of 0 or more. 10 when it is not given, more than the components
   * there are gives those there are, and 0 gives no weight and reads the
   * source once instead of twice.
   *
   * It cuts the weights and nothing else: `projections` and
   * `explainedVariancePercent` are of every component that has variance
   * whatever it is. What it is for is the weights of a million variants,
   * which are 80 MB for 10 components and 8 GB for every component of a
   * dataset of 1000 individuals.
   */
  numPrinComps?: number;
}

/**
 * The principal components of the variants of `variants`, after its steps.
 *
 * Each variant becomes one number per individual, its dosage: how many
 * alleles of the genotype are not the major allele of that variant, which is
 * the most frequent among its called alleles. A genotype with an allele
 * missing takes the mean of the dosages of its variant, so that after
 * centering it pulls its individual nowhere. A variant whose called
 * genotypes all have one dosage has no variance and is left out, a variant
 * with one allele and one where every individual is heterozygous among them;
 * `usedVars` is the ones that were used.
 *
 * The source is read once for the components and a second time for the
 * weights, since a weight needs the eigenvectors and those are known when
 * the first reading ends. Nothing of the size of the variants x the
 * individuals is held: what stays between two blocks is the individuals x
 * individuals matrix. Both readings go through the steps the `Variants` has
 * when the call starts.
 *
 * It is pyNei's `do_pca_from_variants`, whose filters are steps of the
 * `Variants` here, and which gives the weights of every variant where this
 * gives those of the first `numPrinComps` components.
 *
 * A page holds 4 GB at a time, and the matrix of the individuals, its
 * eigenvectors and the workspace of the eigendecomposition are about 6 times
 * 8 bytes per pair of individuals, so a dataset of more than 9381
 * individuals is an `Error` here and is analysed by a program outside the
 * browser, popnei in Python among them.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `transformToBiallelic` is not a boolean, when `numPrinComps` is not a
 * whole number of 0 or more and at most 4294967295, when the analysis of
 * these individuals does not fit in the memory of a page, when the source
 * cannot be read, a wrong line of a VCF among the causes, when a variant has
 * more than two alleles among its called genotypes and `transformToBiallelic`
 * is false, when the pass gives no variant or no variant with variance, when
 * the source has no individual, when the ploidy is above 254, the
 * individuals are more than 46340 or the variants or the weights are more
 * than a whole number of WebAssembly counts, when the linear algebra could
 * not be done, and when `init` has not been awaited.
 */
export function doPcaFromVariants(
  variants: Variants,
  options: DoPcaFromVariantsOptions = {},
): VariantsPcaResult {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  const transformToBiallelic =
    options.transformToBiallelic === undefined
      ? defaultTransformToBiallelic()
      : aBoolean("transformToBiallelic", options.transformToBiallelic);
  const numPrinComps =
    options.numPrinComps === undefined
      ? defaultNumPrinComps()
      : wholeNumberOfZeroOrMore("numPrinComps", options.numPrinComps);
  // The steps of the two passes are a copy of the list, made after the
  // arguments were checked so that nothing refused here leaves one behind:
  // the call takes it over and frees it.
  const result = source.pca_of_variants(
    transformToBiallelic,
    numPrinComps,
    steps.of_a_pass(),
  );
  try {
    return {
      individuals: variants.individuals,
      numComps: result.num_comps(),
      projections: theValuesOf(result.projections(), "projections"),
      explainedVariancePercent: theValuesOf(
        result.explained_variance_percent(),
        "explainedVariancePercent",
      ),
      numPrinComps: result.num_prin_comps(),
      princomps: theValuesOf(result.princomps(), "princomps"),
      usedVars: thePositionsOf(result.used_vars(), "usedVars"),
      passStats: passStatsOf(result.pass_stats()),
    };
  } finally {
    result.free();
  }
}

/**
 * The positions of one array of the result.
 *
 * It is `theValuesOf` for the array of positions, which is a
 * `Uint32Array`: each array of a result leaves the memory of wasm the first
 * time it is asked for.
 *
 * @throws {Error} When the array had been read already.
 */
function thePositionsOf(
  positions: Uint32Array | undefined,
  name: string,
): Uint32Array {
  if (positions === undefined) {
    throw new Error(
      `popnei: \`${name}\` was read twice out of the memory of WebAssembly, ` +
        "which is a defect of popnei; please report it",
    );
  }
  return positions;
}

/**
 * The values of one array of the result.
 *
 * Each array leaves the memory of wasm the first time it is asked for, so
 * the call after that gives nothing; `doPca` and `doPcaFromVariants` read
 * each of them once, and the `Error` here is a defect of this package and
 * not something a caller can do.
 *
 * @throws {Error} When the array had been read already.
 */
export function theValuesOf(
  values: Float64Array | undefined,
  name: string,
): Float64Array {
  if (values === undefined) {
    throw new Error(
      `popnei: \`${name}\` was read twice out of the memory of WebAssembly, ` +
        "which is a defect of popnei; please report it",
    );
  }
  return values;
}
