/**
 * The principal component analysis of a table that the page holds.
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
 * `docs/specs/pca.md` has what it computes. The principal components of the
 * variants of a dataset are the other analysis of that spec and are not
 * written yet; `PcaResult` is what both give, and the one over variants adds
 * the names of the individuals, the variants it used and the counts of its
 * pass to it.
 */

import { do_pca as doPcaOfTheCore } from "../wasm/popnei.js";

import {
  aBoolean,
  valuesOfATable,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";

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
 * @throws {Error} When `data` is not a `Float64Array` of `numRows` times
 * `numCols` values, when a side of the table is not a whole number of 1 or
 * more, when an option is not a boolean, when a value of the table is not
 * finite, when the table is to be standardized and not centered, when it has
 * fewer than 2 rows or no traits, when it is standardized and a trait has no
 * variance, and when `init` has not been awaited.
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
  // The two defaults are pyNei's, which `docs/specs/pca.md` keeps: a table
  // is centered and standardized when the caller says nothing.
  const centerData =
    options.centerData === undefined
      ? true
      : aBoolean("centerData", options.centerData);
  const standardizeData =
    options.standardizeData === undefined
      ? true
      : aBoolean("standardizeData", options.standardizeData);
  const result = doPcaOfTheCore(values, rows, cols, centerData, standardizeData);
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
 * The values of one array of the result.
 *
 * Each array leaves the memory of wasm the first time it is asked for, so
 * the call after that gives nothing; `doPca` reads each of them once, and
 * the `Error` here is a defect of this package and not something a caller
 * can do.
 *
 * @throws {Error} When the array had been read already.
 */
function theValuesOf(
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
