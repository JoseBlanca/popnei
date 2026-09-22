/**
 * The principal components of a table, at `doPca`, and of the variants of a
 * dataset, at `doPcaFromVariants`.
 *
 * The table is iris, 150 rows x 4 traits, the one of pyNei's
 * `test/datasets.py` that `tests/reference/pca/make_reference.py` writes as
 * `tests/reference/pca/iris.tsv`. The numbers asserted are the literals of
 * "How it is verified" of "The PCA of a table" of `docs/specs/pca.md`, which
 * R's `prcomp` gave: where the first row falls on each of the four
 * components, and the percentage of the variance each component holds, with
 * the traits standardized and with them centered only. The tests of the core
 * crate assert the same numbers; the Python tests hold no literal of R and
 * compare with pyNei instead.
 *
 * The dataset of `doPcaFromVariants` is the worked example of "How it is
 * verified" of "The PCA of the variants", `tests/reference/pca/worked.vcf`,
 * 5 individuals and 5 variants of which the third and the fourth have no
 * variance and are left out, and `worked3.vcf`, which adds a variant of
 * three alleles.
 *
 * The analysis itself is tested in the core crate, over every row of the
 * reference files. What these tests say is that the table reaches the core
 * as it was written, row after row and not trait after trait, that the two
 * readers of the variants are opened over the source and the steps of the
 * `Variants`, that the arrays come back with the shape of the result, and
 * that an error of the core is thrown as an `Error` with the message it has
 * in Rust.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { doPca, doPcaFromVariants, init, openVcf } from "popnei";

import { room_for_the_analysis as roomForTheAnalysis } from "../wasm/popnei.js";
import { theValuesOf } from "../dist/pca.js";
import { referencePcaVcf, referenceTable } from "./reference.ts";

await init();

/**
 * The tolerance of "How it is verified" of `docs/specs/pca.md`, which the
 * numbers here are far inside.
 *
 * The reference files write 12 significant digits, and what is compared with
 * them is what popnei computes, so what the tolerance has to cover is the
 * distance between popnei and R. Measured here on 22 September 2026, over
 * every number of the three files of the worked example: 5.73e-15 in a
 * projection, 4.93e-12 in a percentage and 4.43e-14 in a weight. The
 * percentages are the widest apart because they are numbers of about 70 and
 * the others are of about 1. For a table, the spec measures 3e-15 in a
 * projection between the two routes to the components, the
 * eigendecomposition of the product of the rows and that of the product of
 * the traits.
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

test("a table in which no trait has variance is an Error", () => {
  // Centered and not standardized, a table whose traits are each one number
  // repeated is no error at the traits and has no component left: every
  // value of it is 0 once the means are taken out.
  const flat = Float64Array.from([4, 7, 4, 7, 4, 7]);
  assert.throws(
    () => doPca(flat, 3, 2, { standardizeData: false }),
    /no trait has variance, there is nothing to do a PCA with/,
  );
});

test("a trait whose values are too large to standardize is an Error", () => {
  // The squares of the deviations of the first trait sum above the largest
  // float64, so its standard deviation is an infinity and the trait would
  // become a column of zeros, which is what a trait with no variance gives.
  const huge = Float64Array.from([1e154, 1, -1e154, 2, 0, 3]);
  assert.throws(
    () => doPca(huge, 3, 2),
    /the trait at the position 0 cannot be centered or standardized: the squares of its deviations sum above the largest f64/,
  );
});

test("the values of an array of the result are read once", () => {
  // Each array leaves the memory of wasm the first time it is asked for, so
  // a second read gives nothing; `doPca` reads each of them once and what
  // this asserts is that a defect that read one twice would be said and not
  // handed out as an empty array.
  assert.deepEqual(
    theValuesOf(Float64Array.from([1, 2]), "projections"),
    Float64Array.from([1, 2]),
  );
  assert.throws(
    () => theValuesOf(undefined, "projections"),
    /popnei: `projections` was read twice out of the memory of WebAssembly/,
  );
});

test("a table whose buffer was transferred away is an Error", () => {
  const values = Float64Array.from([1, 2, 3, 4, 5, 6]);
  // What a page does when it sends the values to a web worker: the buffer
  // moves and the array that is left has nothing behind it, and the
  // generated code would read it as a table of its own.
  structuredClone(values.buffer, { transfer: [values.buffer] });
  // `detached` is of ES2024, which is later than the library this package
  // is compiled against.
  assert.equal((values.buffer as { detached?: boolean }).detached, true);
  assert.throws(
    () => doPca(values, 3, 2),
    /popnei: the buffer of `data` was transferred/,
  );
});

/** The bytes of the two VCFs of the worked example, read once. */
const WORKED_VCF = await referencePcaVcf("worked.vcf");
const WORKED3_VCF = await referencePcaVcf("worked3.vcf");

/** Where each of the 5 individuals of the worked example falls. */
const WORKED_PROJECTIONS = [
  -0.76228776254, 0.994550847306, -0.320852228237, -0.863187556071,
  -0.87503028832, -0.007193635382, 3.025289193432, -0.097439191351,
  -0.02164630288, -0.536626318751, 0.852948920685, 0.35688580188,
  -0.863187556071, -0.87503028832, -0.007193635382,
];

const WORKED_PERCENT = [76.74407104469, 21.71669104093, 1.53923791438];

/**
 * Where the first individual of `worked3.vcf` falls, and the variance each
 * of its components holds, from `worked3.r.projections.tsv` and
 * `worked3.r.percent.tsv`: the dataset of the worked example with a sixth
 * variant of three alleles, read with every allele that is not the major one
 * counting the same.
 */
const WORKED3_ROW_0 = [
  -0.540318123982, 1.0305897036, 0.564890574265,
];

const WORKED3_PERCENT = [70.7485198149, 26.495544613, 2.75593557214];

/** The weight of each of the three variants used, component after component. */
const WORKED_PRINCOMPS = [
  0.5019679573733, 0.6484646269251, 0.5722952012706, -0.7978606574051,
  0.0917847782057, 0.5958136670595, -0.3338360992098, 0.7556911949445,
  -0.5634574311763,
];

/**
 * The analysis of `bytes` read as a VCF, which frees the handle afterwards.
 *
 * Every call of `doPcaFromVariants` reads the source again, so one handle
 * would serve every test here; each builds its own so that a test that
 * leaves a step on it cannot reach the next.
 */
function theVariantsPca(
  bytes: Uint8Array,
  options?: Parameters<typeof doPcaFromVariants>[1],
): ReturnType<typeof doPcaFromVariants> {
  const variants = openVcf(bytes);
  try {
    return doPcaFromVariants(variants, options);
  } finally {
    variants.free();
  }
}

test("the worked example gives the projections and the weights of R", () => {
  const result = theVariantsPca(WORKED_VCF, { numPrinComps: 3 });
  assert.deepEqual(result.individuals, ["i0", "i1", "i2", "i3", "i4"]);
  assert.equal(result.numComps, 3);
  assertClose(result.projections, WORKED_PROJECTIONS, "the projections");
  assertClose(
    result.explainedVariancePercent,
    WORKED_PERCENT,
    "the percentages",
  );
  assert.deepEqual(result.usedVars, Uint32Array.from([0, 1, 4]));
  assert.equal(result.numPrinComps, 3);
  assertClose(result.princomps, WORKED_PRINCOMPS, "the weights");
});

test("the worked example with no option gives its counts and 10 weights", () => {
  const result = theVariantsPca(WORKED_VCF);
  // The five variants of the file, the two with no variance included: the
  // pass gave them and the analysis left them out, and no filter ran.
  assert.deepEqual(result.passStats, { numVars: 5, filtering: {} });
  // The default of `numPrinComps` is 10, which is more components than this
  // dataset has, so the weights are of the 3 there are: with a default of 0
  // there would be none.
  assert.equal(result.numPrinComps, 3);
  assert.equal(result.princomps.length, 9);
  assertClose(result.princomps, WORKED_PRINCOMPS, "the weights");
});

test("no weights are asked for and princomps has no row", () => {
  const result = theVariantsPca(WORKED_VCF, { numPrinComps: 0 });
  assert.equal(result.numComps, 3);
  assert.equal(result.numPrinComps, 0);
  assert.equal(result.princomps.length, 0);
  // The variants that were used are the columns of the weights there would
  // have been, and they are given whether or not any were asked for.
  assert.deepEqual(result.usedVars, Uint32Array.from([0, 1, 4]));
  assertClose(result.projections, WORKED_PROJECTIONS, "the projections");
});

test("more components than there are gives the weights of those there are", () => {
  const result = theVariantsPca(WORKED_VCF, { numPrinComps: 10 });
  assert.equal(result.numPrinComps, 3);
  assertClose(result.princomps, WORKED_PRINCOMPS, "the weights");
});

test("a variant of three alleles is an Error that names it", () => {
  assert.throws(
    () => theVariantsPca(WORKED3_VCF),
    /the variant at the position 5 among those given has 3 different alleles among its called genotypes/,
  );
});

test("every allele that is not the major one counts the same", () => {
  const result = theVariantsPca(WORKED3_VCF, {
    transformToBiallelic: true,
    numPrinComps: 3,
  });
  // R gives four components for this dataset and popnei the three that have
  // variance, the fourth holding 1.09e-31 of the variance, so the three
  // percentages are the first three of R's.
  assert.equal(result.numComps, 3);
  assert.deepEqual(result.usedVars, Uint32Array.from([0, 1, 4, 5]));
  assertClose(
    result.projections.subarray(0, 3),
    WORKED3_ROW_0,
    "the projections of the first individual",
  );
  assertClose(
    result.explainedVariancePercent,
    WORKED3_PERCENT,
    "the percentages",
  );
  assert.equal(result.princomps.length, 3 * 4);
});

test("a numPrinComps that is not a whole number of 0 or more is an Error", () => {
  assert.throws(
    () => theVariantsPca(WORKED_VCF, { numPrinComps: -1 }),
    /popnei: `numPrinComps` is a whole number of 0 or more/,
  );
});

test("variants that were freed cannot be analysed", () => {
  const variants = openVcf(WORKED_VCF);
  variants.free();
  assert.throws(() => doPcaFromVariants(variants), /were freed/);
});

/**
 * A VCF of `numIndividuals` individuals and 2 variants, for the datasets
 * that are too many individuals for a page.
 */
function vcfOfManyIndividuals(numIndividuals: number): Uint8Array {
  const names = Array.from(
    { length: numIndividuals },
    (_unused, individual) => `i${individual}`,
  );
  const genotypes = (variant: number) =>
    Array.from(
      { length: numIndividuals },
      (_unused, individual) => ["0/0", "0/1", "1/1"][(individual + variant) % 3],
    ).join("\t");
  return new TextEncoder().encode(
    [
      "##fileformat=VCFv4.2",
      `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
      `1\t1\tv0\tA\tC\t.\t.\t.\tGT\t${genotypes(0)}`,
      `1\t2\tv1\tA\tC\t.\t.\t.\tGT\t${genotypes(1)}`,
      "",
    ].join("\n"),
  );
}

test("more individuals than a page holds is an Error and not a trap", () => {
  // 10000 individuals, which the objectives of popnei name and the core
  // takes: before this was checked, the eigendecomposition allocated its
  // workspace, the allocation failed, and an allocation that fails in wasm
  // aborts, which ended the module with `RuntimeError: unreachable` and left
  // every later call to popnei broken. The dataset is refused before any
  // variant is read, so nothing of it is allocated.
  assert.throws(
    () => theVariantsPca(vcfOfManyIndividuals(10000), { numPrinComps: 1 }),
    /the principal components of 10000 individuals hold about 5 GB/,
  );
});

test("the individuals a page holds are the ones measured under node", () => {
  // 9410 individuals ran under node and 9415 ended the module, so the
  // largest dataset popnei takes is under both. Neither of the two is run
  // here: the one that works takes five minutes, since the time of the
  // eigendecomposition goes with the cube of the individuals.
  roomForTheAnalysis(9381);
  assert.throws(
    () => roomForTheAnalysis(9382),
    /the principal components of 9382 individuals hold about 5 GB/,
  );
  // The worked example, and every dataset a page really holds, passes.
  roomForTheAnalysis(5);
  roomForTheAnalysis(0);
});
