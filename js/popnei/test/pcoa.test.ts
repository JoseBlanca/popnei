/**
 * The principal coordinates of a `Distances`, at `doPcoa`, Lingoes'
 * correction of one, at `correctDistsByLingoes`, and the principal
 * coordinates of the Kosman distances of a `Variants`, at
 * `doPcoaFromVariants`.
 *
 * The distances are the worked example of "How it is verified" of "The
 * principal coordinates of distances" of `docs/specs/pca.md`: the ten
 * distances of pyNei's `test_pcoa`, of five individuals `i1` to `i5`, which
 * are not Euclidean. The numbers asserted are the literals of that part,
 * which R's `ape::pcoa` with `correction = "lingoes"` gave, and which
 * `tests/reference/pca/small.lingoes.r.*.tsv` hold. The variants are the
 * panel of that part, `tests/reference/dists/panel.vcf.gz`, 200 individuals
 * and 1200 variants, whose corrected principal coordinates R gave in
 * `tests/reference/pca/panel.lingoes.r.*.tsv`.
 *
 * The analysis and the correction are tested in the core crate, the twin and
 * the Kosman distances of `four_alleles.vcf.gz` among them. What these tests
 * say is that the vector reaches the core in its order, that the arrays come
 * back with the shape of the result and the names and the counts of the
 * `Distances` or of the pass, that the corrected distances are a
 * `Distances` of the same individuals, that an error of the core is thrown
 * with the names a TypeScript user writes, and that a page refuses
 * distances, or variants, of more individuals than it holds the analysis
 * of.
 */

import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

import type { PassStats, Progress, Variants } from "popnei";
import {
  calcPairwiseKosmanDists,
  correctDistsByLingoes,
  Distances,
  doPcoa,
  doPcoaFromVariants,
  init,
  openVcf,
} from "popnei";

import { room_for_the_principal_coordinates_of as roomForThePrincipalCoordinates } from "../wasm/popnei.js";

import { referenceDists } from "./reference.ts";

await init();

/**
 * The tolerance of "How it is verified" of `docs/specs/pca.md`: the
 * reference files write 15 significant digits of R, and the core crate finds
 * popnei within 3e-14 of them on these distances.
 */
const TOLERANCE = 1e-9;

/** The individuals of the worked example. */
const NAMES = ["i1", "i2", "i3", "i4", "i5"];

/**
 * The ten distances of pyNei's `test_pcoa`, in the order (i1, i2), (i1, i3),
 * ..., (i4, i5).
 */
const TEN_DISTANCES = [0.2, 0.3, 0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2];

/**
 * Counts of a pass that no calculation gives, so that a test can tell them
 * from the empty ones: a `Distances` a user builds carries what the user put
 * in it, and the result carries them on.
 */
const PASS_STATS: PassStats = {
  numVars: 17,
  filtering: { maf: { varsProcessed: 20, varsKept: 17 } },
};

/** The ten distances as the `Distances` a user builds. */
function theTenDistances(): Distances {
  return new Distances(Float64Array.from(TEN_DISTANCES), NAMES, PASS_STATS);
}

/** Lingoes' constant of the ten distances, c, the most negative eigenvalue. */
const CONSTANT = 0.0640069399611263;

/** The share of the negative eigenvalues of the ten distances. */
const NEGATIVE_PERCENT = 7.88262807403034;

/** The distance of i1 and i2 corrected, sqrt(0.2² + 2c). */
const FIRST_CORRECTED = 0.409894962060102;

/** Where each individual falls on the 3 components of the corrected ones. */
const CORRECTED_PROJECTIONS = [
  [-0.431869046368213, -0.0415086068851603, 0.224881556185394],
  [-0.283479006142767, -0.158680403878776, -0.138801060767122],
  [-0.269028184151739, 0.2124683967908, -0.131478395434635],
  [0.492920681079785, 0.195146288470815, 0.0612591236372794],
  [0.491455555582935, -0.207425674497679, -0.015861223620916],
];

/** The share of the variance of each of the 3 components. */
const CORRECTED_PERCENT = [77.1278402980914, 14.3397713765961, 8.53238832531248];

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

test("the ten distances are refused, and the message names the correction of TypeScript", () => {
  assert.throws(
    () => doPcoa(theTenDistances()),
    (error: unknown) => {
      assert.ok(error instanceof Error);
      assert.match(
        error.message,
        /^1 of the 5 eigenvalues of the matrix of the squared distances are negative, 7\.88 percent of the sum of all of them/,
      );
      assert.match(
        error.message,
        /`correctDistsByLingoes` makes them Euclidean by adding the same amount to every squared distance$/,
      );
      assert.doesNotMatch(error.message, /correct_dists_by_lingoes/);
      return true;
    },
  );
});

test("Lingoes' correction of the ten distances gives R's constant and a Distances of the same individuals", () => {
  const given = theTenDistances();
  const correction = correctDistsByLingoes(given);
  assert.ok(Math.abs(correction.constant - CONSTANT) <= TOLERANCE);
  assert.ok(
    Math.abs(correction.negativeEigenvaluesPercent - NEGATIVE_PERCENT) <=
      TOLERANCE,
  );
  const corrected = correction.distances;
  assert.ok(corrected instanceof Distances);
  assert.equal(corrected.distVector.length, 10);
  assert.ok(
    Math.abs((corrected.distVector[0] as number) - FIRST_CORRECTED) <=
      TOLERANCE,
  );
  // Every distance is sqrt(d² + 2c) of its own, which the order of the vector
  // is what puts beside the right d: the last pair, i4 and i5, is at 0.2 as
  // the first is, and the ninth at 0.8.
  assert.ok(
    Math.abs((corrected.distVector[9] as number) - FIRST_CORRECTED) <=
      TOLERANCE,
  );
  assert.ok(
    Math.abs(
      (corrected.distVector[8] as number) - Math.sqrt(0.8 * 0.8 + 2 * CONSTANT),
    ) <= TOLERANCE,
  );
  assert.deepEqual(corrected.names, NAMES);
  assert.deepEqual(corrected.passStats, PASS_STATS);
  assert.equal(corrected.standardErrors, null);
  // The distances given are the user's and are left as they were.
  assert.deepEqual([...given.distVector], TEN_DISTANCES);
});

test("the principal coordinates of the corrected distances are R's", () => {
  const corrected = correctDistsByLingoes(theTenDistances()).distances;
  const result = doPcoa(corrected);
  assert.equal(result.numComps, 3);
  assert.ok(result.projections instanceof Float64Array);
  assert.ok(result.explainedVariancePercent instanceof Float64Array);
  assertClose(
    result.projections,
    CORRECTED_PROJECTIONS.flat(),
    "the projections, individual after individual",
  );
  assertClose(result.explainedVariancePercent, CORRECTED_PERCENT, "the percentages");
  // doPcoa corrects nothing, so both numbers of a correction are 0 in its
  // result, whatever the distances went through before.
  assert.equal(result.lingoesConstant, 0);
  assert.equal(result.negativeEigenvaluesPercent, 0);
  assert.deepEqual(result.names, NAMES);
  assert.deepEqual(result.passStats, PASS_STATS);
});

test("distances that are Euclidean are not corrected", () => {
  // Three individuals on a line, at 0, 1 and 3: the straight line distances
  // of points are Euclidean, so the constant and the share are 0 and the
  // distances come back as they were.
  const onALine = new Distances(Float64Array.from([1, 3, 2]), ["a", "b", "c"], PASS_STATS);
  const correction = correctDistsByLingoes(onALine);
  assert.equal(correction.constant, 0);
  assert.equal(correction.negativeEigenvaluesPercent, 0);
  assert.deepEqual([...correction.distances.distVector], [1, 3, 2]);
  const result = doPcoa(onALine);
  assert.equal(result.numComps, 1);
  assert.equal(result.explainedVariancePercent.length, 1);
  assert.ok(
    Math.abs((result.explainedVariancePercent[0] as number) - 100) <= TOLERANCE,
  );
});

test("pairs with no distance are refused with the names of the individuals", () => {
  // The pairs (i2, i4) and (i2, i5) have none, so i2 is in both of them and
  // the first in the order of the vector is i2 with i4.
  const withNoDistance = Float64Array.from(TEN_DISTANCES);
  withNoDistance[5] = Number.NaN;
  withNoDistance[6] = Number.NaN;
  const distances = new Distances(withNoDistance, NAMES, PASS_STATS);
  for (const call of [
    () => doPcoa(distances),
    () => correctDistsByLingoes(distances),
  ]) {
    assert.throws(call, (error: unknown) => {
      assert.ok(error instanceof Error);
      assert.equal(
        error.message,
        "2 of the 10 pairs of individuals have no distance, the first of " +
          "them `i2` and `i4`, and `i2` is in 2 of them; a principal " +
          "coordinate analysis places every individual by its distance to " +
          "every other, so each of those pairs has to be given a distance or " +
          "one of its two individuals taken out of the distances",
      );
      return true;
    });
  }
});

test("a negative or an infinite distance is refused with the names of its pair", () => {
  // The second distance of the vector is of i1 and i3.
  const withANegative = Float64Array.from(TEN_DISTANCES);
  withANegative[1] = -0.3;
  const negative = new Distances(withANegative, NAMES, PASS_STATS);
  // The ninth is of i3 and i5.
  const withAnInfinity = Float64Array.from(TEN_DISTANCES);
  withAnInfinity[8] = Number.POSITIVE_INFINITY;
  const infinite = new Distances(withAnInfinity, NAMES, PASS_STATS);
  // A distance below 1e-6 is written as JavaScript writes it, with an
  // exponent, and not as the 300 zeros of the number in full.
  const withATinyNegative = Float64Array.from(TEN_DISTANCES);
  withATinyNegative[1] = -1e-300;
  const tinyNegative = new Distances(withATinyNegative, NAMES, PASS_STATS);
  for (const [distances, pair, value] of [
    [negative, "`i1` and `i3`", "-0.3"],
    [infinite, "`i3` and `i5`", "Infinity"],
    [tinyNegative, "`i1` and `i3`", "-1e-300"],
  ] as const) {
    for (const call of [
      () => doPcoa(distances),
      () => correctDistsByLingoes(distances),
    ]) {
      assert.throws(call, (error: unknown) => {
        assert.ok(error instanceof Error);
        assert.equal(
          error.message,
          `the distance of ${pair} is ${value}, and a principal coordinate ` +
            "analysis needs every distance finite and 0 or above; a negative " +
            "F_ST or f_2 is of two populations the dataset cannot tell apart",
        );
        return true;
      });
    }
  }
});

test("distances that are all 0, and fewer than two individuals, are refused", () => {
  const atOnePoint = new Distances(new Float64Array(10), NAMES, PASS_STATS);
  for (const call of [
    () => doPcoa(atOnePoint),
    () => correctDistsByLingoes(atOnePoint),
  ]) {
    assert.throws(
      call,
      /every distance is 0, so the individuals are all at one point and there is nothing to do a PCoA with$/,
    );
  }
  const alone = new Distances(new Float64Array(0), ["i1"], PASS_STATS);
  for (const call of [() => doPcoa(alone), () => correctDistsByLingoes(alone)]) {
    assert.throws(
      call,
      /there is 1 individual, and a principal coordinate analysis places 2 at least/,
    );
  }
});

test("something that is not a Distances is refused and named", () => {
  const vector = Float64Array.from(TEN_DISTANCES);
  assert.throws(
    () => doPcoa(vector as unknown as Distances),
    /`distances` is a `Distances`, and an object of the type `Float64Array` was given/,
  );
  assert.throws(
    () => correctDistsByLingoes(vector as unknown as Distances),
    /`distances` is a `Distances`, and an object of the type `Float64Array` was given/,
  );
});

test("the individuals a page holds the principal coordinates of are 8695", () => {
  // 56.8 bytes for each cell of the individuals x individuals matrix, which
  // "How it runs" of the spec counts, fit 8695 individuals in the 4 GiB of a
  // page and not 8696. Neither is run here: the one that works takes minutes,
  // since the time of the eigendecomposition goes with the cube of the
  // individuals.
  roomForThePrincipalCoordinates(8695);
  assert.throws(
    () => roomForThePrincipalCoordinates(8696),
    /the principal coordinates of 8696 individuals hold about 5 GB, the individuals x individuals matrix of the analysis/,
  );
  roomForThePrincipalCoordinates(5);
  roomForThePrincipalCoordinates(0);
});

test("a Distances too large for a page is refused by both functions", () => {
  // The 37805860 distances of 8696 individuals are 302 MB, which node holds.
  // The package asks the page before it copies them into the memory of wasm
  // and the binding crate asks again, so what this says is that neither
  // function reaches the analysis; which of the two refused it is not seen
  // from here, since both give the same message.
  const names = Array.from({ length: 8696 }, (_, position) => `s${position}`);
  const distances = new Distances(
    new Float64Array((8696 * 8695) / 2),
    names,
    PASS_STATS,
  );
  for (const call of [
    () => doPcoa(distances),
    () => correctDistsByLingoes(distances),
  ]) {
    assert.throws(
      call,
      /the principal coordinates of 8696 individuals hold about 5 GB/,
    );
  }
});

/** The panel of "How it is verified", 200 individuals and 1200 variants. */
const PANEL_VCF = await referenceDists("panel.vcf.gz");

/** How many variants the panel holds, every one of which a pass gives. */
const PANEL_NUM_VARS = 1200;

/**
 * A table of R's for the panel corrected by Lingoes, one of
 * `tests/reference/pca/panel.lingoes.r.*.tsv`: the names of its rows, when
 * its first line names the columns and the first field of every other line
 * names the row, and its values row after row.
 *
 * @throws {Error} When a value is not a finite number: a table read wrong
 * would be the check of something else.
 */
async function rsTableOfThePanel(
  name: string,
  withNames: boolean,
): Promise<{ rows: string[]; values: number[] }> {
  const text = await readFile(
    new URL(`../../../tests/reference/pca/${name}`, import.meta.url),
    "utf8",
  );
  const lines = text.split("\n").filter((line) => line.trim() !== "");
  const body = withNames ? lines.slice(1) : lines;
  const rows: string[] = [];
  const values: number[] = [];
  for (const line of body) {
    const fields = line.split("\t");
    if (withNames) {
      rows.push(fields.shift() as string);
    }
    for (const field of fields) {
      const value = Number(field.trim());
      if (!Number.isFinite(value)) {
        throw new Error(`${name}: the field \`${field}\` is not a number`);
      }
      values.push(value);
    }
  }
  return { rows, values };
}

/** R's projections of the panel corrected, 200 individuals x 198 components. */
const PANEL_PROJECTIONS = await rsTableOfThePanel(
  "panel.lingoes.r.projections.tsv",
  true,
);

/** R's share of the variance of each of the 198 components. */
const PANEL_PERCENT = (
  await rsTableOfThePanel("panel.lingoes.r.percent.tsv", false)
).values;

/** R's constant of Lingoes' correction of the panel, and its percent. */
const [PANEL_CONSTANT, PANEL_NEGATIVE_PERCENT] = (
  await rsTableOfThePanel("panel.lingoes.r.constant.tsv", false)
).values as [number, number];

/**
 * The principal coordinates of the panel, with the steps `withSteps` puts
 * on its `Variants` before the call.
 */
function pcoaOfThePanel(
  options: Parameters<typeof doPcoaFromVariants>[1],
  withSteps?: (variants: Variants) => void,
): ReturnType<typeof doPcoaFromVariants> {
  const variants = openVcf(PANEL_VCF, { onlyPassed: false });
  try {
    withSteps?.(variants);
    return doPcoaFromVariants(variants, options);
  } finally {
    variants.free();
  }
}

test("the panel is refused without correctByLingoes, with its 44 of 200 and its percent", () => {
  for (const options of [{}, { correctByLingoes: false }]) {
    assert.throws(
      () => pcoaOfThePanel(options),
      (error: unknown) => {
        assert.ok(error instanceof Error);
        assert.match(
          error.message,
          /^44 of the 200 eigenvalues of the matrix of the squared distances are negative, 2\.98 percent of the sum of all of them/,
        );
        assert.match(
          error.message,
          /; `correctByLingoes` makes them Euclidean by adding the same amount to every squared distance$/,
        );
        assert.doesNotMatch(error.message, /correct_by_lingoes|correctDistsByLingoes/);
        return true;
      },
    );
  }
});

test("the principal coordinates of the panel with correctByLingoes are R's", () => {
  const result = pcoaOfThePanel({ correctByLingoes: true });
  assert.equal(result.numComps, 198);
  assert.ok(result.projections instanceof Float64Array);
  assert.ok(result.explainedVariancePercent instanceof Float64Array);
  assert.ok(Math.abs(result.lingoesConstant - 0.014182298472042) <= TOLERANCE);
  assert.ok(
    Math.abs(result.negativeEigenvaluesPercent - 2.98343616373556) <= TOLERANCE,
  );
  // The literals of the spec, which the files of R hold as well.
  assert.ok(Math.abs(result.lingoesConstant - PANEL_CONSTANT) <= TOLERANCE);
  assert.ok(
    Math.abs(result.negativeEigenvaluesPercent - PANEL_NEGATIVE_PERCENT) <=
      TOLERANCE,
  );
  const row = (individual: number) =>
    [...result.projections.subarray(individual * 198, individual * 198 + 3)];
  assertClose(
    row(0),
    [0.0131009923566961, 0.103593034190948, -0.0461610570616341],
    "the first three projections of s000",
  );
  assertClose(
    row(1),
    [0.0188069490905941, 0.102449570787996, -0.0506243303501242],
    "the first three projections of s001",
  );
  assertClose(
    row(199),
    [-0.0728375733863349, -0.0163081366424616, 0.0108102147666687],
    "the first three projections of s199",
  );
  assertClose(
    result.explainedVariancePercent.subarray(0, 3),
    [9.62407114041929, 6.65101405201371, 1.73580899943624],
    "the first three percentages",
  );
  // And every number of R: the rows of its file are the individuals in the
  // order of the VCF, which is the order of `individuals`.
  assert.deepEqual(result.individuals, PANEL_PROJECTIONS.rows);
  assertClose(result.projections, PANEL_PROJECTIONS.values, "the projections");
  assertClose(result.explainedVariancePercent, PANEL_PERCENT, "the percentages");
  const expected: PassStats = { numVars: PANEL_NUM_VARS, filtering: {} };
  assert.deepEqual(result.passStats, expected);
  // What a user of the result of the PCA reads has the same names here, and
  // what a PCoA has not is not here.
  assert.ok(!("numPrinComps" in result));
  assert.ok(!("princomps" in result));
  assert.ok(!("usedVars" in result));
});

test("minNumSnps 1105 leaves 35 pairs of the panel with no distance, named in TypeScript", () => {
  // The pairs are looked at before the eigenvalues, so the message is the
  // same whether or not the correction was asked for.
  for (const correctByLingoes of [true, false]) {
    assert.throws(
      () => pcoaOfThePanel({ minNumSnps: 1105, correctByLingoes }),
      {
        name: "Error",
        message:
          "35 of the 19900 pairs of individuals have no distance, the first " +
          "of them `s001` and `s082`, and `s082` is in 17 of them; those " +
          "pairs were called together at fewer variants than `minNumSnps`, " +
          "or at none; take that individual out with `filterIndividuals`, " +
          "lower `minNumSnps`, or run the PCA of the variants, which gives " +
          "every individual a projection",
      },
    );
  }
});

test("the counts of a filter of the panel are those the Kosman distances give over it", () => {
  const withTheFilter = (variants: Variants) => {
    variants.filterByMaf(0.9);
  };
  const result = pcoaOfThePanel({ correctByLingoes: true }, withTheFilter);
  const variants = openVcf(PANEL_VCF, { onlyPassed: false });
  let distances: Distances;
  try {
    withTheFilter(variants);
    distances = calcPairwiseKosmanDists(variants);
  } finally {
    variants.free();
  }
  assert.deepEqual(result.passStats, distances.passStats);
  const maf = result.passStats.filtering.maf;
  assert.ok(maf !== undefined, "the counts of the filter are there");
  assert.equal(maf.varsProcessed, PANEL_NUM_VARS);
  // The filter took some variants out, or the counts would say nothing a
  // pass without it does not.
  assert.ok(maf.varsKept < PANEL_NUM_VARS, `it kept ${maf.varsKept}`);
  assert.equal(result.passStats.numVars, maf.varsKept);
});

test("filterIndividuals gives the kept individuals, in their order, at the numbers of the corrected distances", () => {
  // The analysis of the variants with the correction is `doPcoa` of the
  // corrected Kosman distances, within the rounding of the square roots, as
  // the spec says: the two routes put each individual on its own row only
  // if the names and the rows are in the same order.
  const kept = ["s150", "s003", "s042", "s199", "s000", "s077", "s121", "s010"];
  const keep = (variants: Variants) => {
    variants.filterIndividuals(kept);
  };
  const result = pcoaOfThePanel({ correctByLingoes: true }, keep);
  assert.deepEqual(result.individuals, kept);
  const variants = openVcf(PANEL_VCF, { onlyPassed: false });
  let distances: Distances;
  try {
    keep(variants);
    distances = calcPairwiseKosmanDists(variants);
  } finally {
    variants.free();
  }
  const correction = correctDistsByLingoes(distances);
  const fromTheDistances = doPcoa(correction.distances);
  assert.deepEqual(fromTheDistances.names, kept);
  assert.equal(result.numComps, fromTheDistances.numComps);
  assertClose(
    result.projections,
    [...fromTheDistances.projections],
    "the projections",
  );
  assertClose(
    result.explainedVariancePercent,
    [...fromTheDistances.explainedVariancePercent],
    "the percentages",
  );
  assert.ok(
    Math.abs(result.lingoesConstant - correction.constant) <= TOLERANCE,
  );
  assert.ok(
    Math.abs(
      result.negativeEigenvaluesPercent - correction.negativeEigenvaluesPercent,
    ) <= TOLERANCE,
  );
});

/**
 * A VCF of `numIndividuals` diploid individuals and two variants, whose
 * genotypes cycle through the three of one alternative allele.
 */
function vcfOfManyIndividuals(numIndividuals: number): Uint8Array {
  const names = Array.from(
    { length: numIndividuals },
    (_unused, individual) => `i${individual}`,
  );
  const genotypes = (variant: number) =>
    Array.from(
      { length: numIndividuals },
      (_unused, individual) =>
        ["0/0", "0/1", "1/1"][(individual + variant) % 3],
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

test("variants of more individuals than a page holds are refused before the source is read", () => {
  const variants = openVcf(vcfOfManyIndividuals(8696));
  const calls: Progress[] = [];
  variants.onProgress((progress) => {
    calls.push(progress);
  });
  try {
    assert.throws(
      () => doPcoaFromVariants(variants, { correctByLingoes: true }),
      /the principal coordinates of 8696 individuals hold about 5 GB/,
    );
    // No pass started: its first read would have told the page.
    assert.deepEqual(calls, []);
    // The matrix is of the individuals the pass gives, so the same file
    // with three of them kept is analysed.
    variants.filterIndividuals(["i0", "i1", "i2"]);
    const result = doPcoaFromVariants(variants, { correctByLingoes: true });
    assert.deepEqual(result.individuals, ["i0", "i1", "i2"]);
    assert.ok(calls.length > 0);
  } finally {
    variants.free();
  }
});

test("the options of doPcoaFromVariants are checked before the source is read", () => {
  const variants = openVcf(PANEL_VCF, { onlyPassed: false });
  try {
    assert.throws(
      () =>
        doPcoaFromVariants(variants, {
          correctByLingoe: true,
        } as unknown as Parameters<typeof doPcoaFromVariants>[1]),
      /`correctByLingoe`/,
    );
    assert.throws(
      () => doPcoaFromVariants(variants, { minNumSnps: -1 }),
      /`minNumSnps` is a whole number of 0 or more and at most 4294967295/,
    );
    assert.throws(
      () => doPcoaFromVariants(variants, { minNumSnps: 4294967296 }),
      /`minNumSnps` is a whole number of 0 or more and at most 4294967295/,
    );
    assert.throws(
      () =>
        doPcoaFromVariants(variants, {
          correctByLingoes: "yes",
        } as unknown as Parameters<typeof doPcoaFromVariants>[1]),
      /`correctByLingoes` is true or false/,
    );
    assert.throws(
      () => doPcoaFromVariants({} as unknown as Variants),
      /`variants`/,
    );
  } finally {
    variants.free();
  }
});
