/**
 * The principal components of a table, at `doPca`.
 *
 * The table is iris, 150 rows x 4 traits, the one of pyNei's
 * `test/datasets.py` that `tests/reference/pca/make_reference.py` writes as
 * `tests/reference/pca/iris.tsv`. The numbers asserted are the literals of
 * "How it is verified" of "The PCA of a table" of `docs/specs/pca.md`, which
 * R's `prcomp` gave: where the first row falls on each of the four
 * components, and the percentage of the variance each component holds, with
 * the traits standardized and with them centered only. The Python tests
 * assert the same numbers.
 *
 * The analysis itself is tested in the core crate, over every row of the
 * reference files. What these tests say is that the table reaches the core
 * as it was written, row after row and not trait after trait, that the
 * arrays come back with the shape of the result, and that an error of the
 * core is thrown as an `Error` with the message it has in Rust.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { doPca, init } from "popnei";

import { referenceTable } from "./reference.ts";

await init();

/**
 * The tolerance of "How it is verified" of `docs/specs/pca.md`: the
 * reference files write 12 significant digits, and the two routes to the
 * components differ by 3e-15 in a projection.
 */
const TOLERANCE = 1e-9;

/** Iris, 150 rows x 4 traits, read once for every test here. */
const IRIS = await referenceTable("iris.tsv");

/** Where the first row of iris falls on each component, standardized. */
const STANDARDIZED_ROW_0 = [
  -2.26470280881, 0.480026596521, -0.1277060223, -0.0241682038555,
];

/** The variance each component of iris holds, standardized, as a percentage. */
const STANDARDIZED_PERCENT = [
  72.9624454133, 22.8507617867, 3.66892188928, 0.517870910715,
];

/** The same two, with the traits centered and not standardized. */
const CENTERED_ROW_0 = [
  -2.68412562597, 0.319397246585, -0.0279148275894, -0.00226243707132,
];

const CENTERED_PERCENT = [
  92.4618723202, 5.30664831171, 1.71026098079, 0.521218387328,
];

/**
 * The 3 rows x 5 traits of "How it is verified" of the same part, which has
 * fewer rows than traits: the matrix that is decomposed is then the 3 x 3
 * product of the rows and not the 5 x 5 product of the traits, and a table
 * of 150 rows x 4 traits never reaches that half of the analysis.
 */
const FIVE_TRAITS = Float64Array.from([
  1, 2, 3, 4, 5, 2, 4, 1, 3, 2, 5, 1, 4, 2, 6,
]);

/** Its first row and its percentages, neither centered nor standardized. */
const FIVE_TRAITS_ROW_0 = [7.07789159873, 1.12190772301, 1.90912901021];

const FIVE_TRAITS_PERCENT = [87.0392022767, 9.31343669027, 3.64736103304];

/** Each value against the one of the reference, within the tolerance. */
function assertClose(
  got: ArrayLike<number>,
  expected: readonly number[],
  what: string,
): void {
  assert.equal(got.length, expected.length, `${what}: the count of values`);
  for (const [position, reference] of expected.entries()) {
    const value = got[position] as number;
    assert.ok(
      Math.abs(value - reference) <= TOLERANCE,
      `${what}: the value at ${position} is ${value} and the reference is ${reference}`,
    );
  }
}

test("iris standardized gives the projections and the percentages of R", () => {
  const result = doPca(IRIS.values, IRIS.numRows, IRIS.numCols);
  assert.equal(result.numComps, 4);
  assertClose(
    result.projections.subarray(0, 4),
    STANDARDIZED_ROW_0,
    "the projections of the first row",
  );
  assertClose(
    result.explainedVariancePercent,
    STANDARDIZED_PERCENT,
    "the percentages",
  );
});

test("iris centered and not standardized gives the other numbers of R", () => {
  const result = doPca(IRIS.values, IRIS.numRows, IRIS.numCols, {
    standardizeData: false,
  });
  assertClose(
    result.projections.subarray(0, 4),
    CENTERED_ROW_0,
    "the projections of the first row",
  );
  assertClose(
    result.explainedVariancePercent,
    CENTERED_PERCENT,
    "the percentages",
  );
});

test("the result of iris has a row per row and a weight per trait", () => {
  const result = doPca(IRIS.values, IRIS.numRows, IRIS.numCols);
  assert.equal(result.numComps, 4);
  assert.equal(result.numPrinComps, 4);
  assert.equal(result.projections.length, 150 * 4);
  assert.equal(result.princomps.length, 4 * 4);
  assert.ok(result.projections instanceof Float64Array);
  assert.ok(result.explainedVariancePercent instanceof Float64Array);
  assert.ok(result.princomps instanceof Float64Array);
});

test("a table of fewer rows than traits, neither centered nor standardized", () => {
  const result = doPca(FIVE_TRAITS, 3, 5, {
    centerData: false,
    standardizeData: false,
  });
  assert.equal(result.numComps, 3);
  assertClose(
    result.projections.subarray(0, 3),
    FIVE_TRAITS_ROW_0,
    "the projections of the first row",
  );
  assertClose(
    result.explainedVariancePercent,
    FIVE_TRAITS_PERCENT,
    "the percentages",
  );
  assert.equal(result.princomps.length, 3 * 5);
});

test("a value that is not finite is an Error that says where it is", () => {
  const values = Float64Array.from(IRIS.values);
  values[4 + 2] = Number.POSITIVE_INFINITY;
  assert.throws(
    () => doPca(values, IRIS.numRows, IRIS.numCols),
    /the value at row 1, trait 2 of the table is inf/,
  );
});

test("standardizing without centering is an Error", () => {
  assert.throws(
    () =>
      doPca(IRIS.values, IRIS.numRows, IRIS.numCols, {
        centerData: false,
        standardizeData: true,
      }),
    /standardized and not centered/,
  );
});

test("a table of fewer than two rows is an Error", () => {
  assert.throws(
    () => doPca(IRIS.values.subarray(0, 4), 1, 4),
    /the table is 1 x 4/,
  );
});

test("a trait with no variance is an Error that names the trait", () => {
  // The three traits of `test_pca_refuses_traits_with_no_variance` of
  // pyNei, the second of which is 5 in every row.
  const fixed = Float64Array.from([1, 5, 3, 2, 5, 1, 3, 5, 2]);
  assert.throws(() => doPca(fixed, 3, 3), /1 of the 3 traits have no variance/);
});

test("a table whose values are not as many as its two sides is an Error", () => {
  assert.throws(
    () => doPca(IRIS.values, IRIS.numRows, 3),
    /popnei: the table was given as 150 x 3/,
  );
});

test("a table that is not a Float64Array is an Error that says so", () => {
  assert.throws(
    () => doPca([1, 2, 3, 4] as unknown as Float64Array, 2, 2),
    /popnei: `data` is a Float64Array/,
  );
});

test("a side of the table that is not a whole number is an Error", () => {
  assert.throws(
    () => doPca(IRIS.values, 150.5, 4),
    /popnei: `numRows` is a whole number/,
  );
});

test("an option that is not a boolean is an Error", () => {
  assert.throws(
    () =>
      doPca(IRIS.values, IRIS.numRows, IRIS.numCols, {
        centerData: "yes" as unknown as boolean,
      }),
    /popnei: `centerData` is true or false/,
  );
});
