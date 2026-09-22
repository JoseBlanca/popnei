/**
 * The Kosman distances between individuals from TypeScript:
 * `calcPairwiseKosmanDists` and the `Distances` it gives.
 *
 * "How it is verified" of `docs/specs/dists.md` has the numbers. The first
 * dataset is the panel, `tests/reference/dists/panel.vcf.gz`: 200 diploid
 * individuals, `s000` to `s199`, 1200 biallelic variants, 3 in 100
 * genotypes missing. Its distances are those that `gd.kosman` of the R
 * package PopGenReport 3.1.3 gives, which are in
 * `tests/reference/dists/panel.gdkosman.tsv`, one line for each pair in the
 * order of the distance vector; five of them are in the table of the spec
 * and are written here as literals. The second dataset is the worked
 * example of the spec, 4 variants of 3 individuals, written as a VCF here:
 * its three distances are 1/4, 5/6 and 2/6, two integers divided once, so
 * they are asserted exactly.
 *
 * Nothing here computes an expected value with popnei: every number comes
 * from the spec or from the file that R wrote.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Variants } from "popnei";
import { Distances, calcPairwiseKosmanDists, init, openVcf } from "popnei";

import { referenceDists, vcfOf } from "./reference.ts";

await init();

/** How many variants the panel holds, every one of them given. */
const PANEL_NUM_VARS = 1200;

/** How many individuals the panel holds, `s000` to `s199`. */
const PANEL_NUM_INDIVIDUALS = 200;

/** The tolerance of a distance, the digits of the shortest of the spec. */
const TOLERANCE = 1e-9;

/** The bytes of the panel, read once for every test that runs on it. */
const PANEL_VCF = await referenceDists("panel.vcf.gz");

/**
 * Five pairs of the panel with the distance `gd.kosman` gives for each, the
 * table of "How it is verified" of the spec: the first two pairs, the last
 * pair, the largest distance and the smallest one.
 */
const PANEL_LITERALS: { first: string; second: string; dist: number }[] = [
  { first: "s000", second: "s001", dist: 0.1657754010695187 },
  { first: "s000", second: "s002", dist: 0.16666666666666666 },
  { first: "s198", second: "s199", dist: 0.15476190476190477 },
  { first: "s010", second: "s033", dist: 0.3579697239536955 },
  { first: "s116", second: "s119", dist: 0.13680494263018536 },
];

/**
 * The worked example of the spec: 4 variants of 3 diploid individuals, with
 * a half called genotype at the third variant and a missing one at the
 * fourth. Every variant declares the three alleles `A`, `C` and `G`, so
 * that the allele 2 of the second and the third variants is one the file
 * names.
 */
const WORKED_EXAMPLE = vcfOfThreeIndividuals(
  ["s0", "s1", "s2"],
  ["0/0\t0/1\t1/1", "0/1\t0/1\t1/2", "0/0\t0/.\t2/2", "./.\t1/1\t1/1"],
);

/**
 * The three distances of the worked example, the sums of d 1, 5 and 2 over
 * the numbers of variants 2, 3 and 3, each of them the two integers of the
 * spec divided once.
 */
const WORKED_EXAMPLE_DISTS = [1 / 4, 5 / 6, 2 / 6];

/**
 * The haploid worked example of the spec: 4 variants of 3 haploid
 * individuals, the third of them called in `h1` and `h2` alone. Its
 * distances are 1/3, 2/3 and 1/2, the sums of d 1, 2 and 2 over the numbers
 * of variants 3, 3 and 4.
 */
const HAPLOID_WORKED_EXAMPLE = vcfOfThreeIndividuals(
  ["h0", "h1", "h2"],
  ["0\t0\t1", "0\t1\t2", ".\t1\t1", "0\t0\t0"],
);

/** The tetraploid dataset of the spec, 12 individuals and 200 variants. */
const TETRAPLOID_VCF = await referenceDists("tetraploid.vcf.gz");
const TETRAPLOID_NUM_INDIVIDUALS = 12;
const TETRAPLOID_NUM_VARS = 200;

/**
 * Three pairs of that dataset with the distance `gd.kosman` gives for each,
 * the tetraploid rows of the table of "How it is verified" of the spec.
 */
const TETRAPLOID_LITERALS: { first: string; second: string; dist: number }[] = [
  { first: "t00", second: "t01", dist: 0.375 },
  { first: "t00", second: "t02", dist: 0.38797814207650272 },
  { first: "t00", second: "t03", dist: 0.40163934426229508 },
];

/**
 * A VCF of the three individuals `names`, one line for each text of
 * `genotypes`, which holds the genotype of each of them.
 *
 * Every variant declares the three alleles `A`, `C` and `G`, so that the
 * allele 2 of a genotype is one the file names, and no genotype is written
 * with a ploidy of its own: the ploidy is what `openVcf` is told.
 */
function vcfOfThreeIndividuals(
  names: readonly string[],
  genotypes: readonly string[],
): Uint8Array {
  const lines = genotypes.map(
    (variant, index) =>
      `chr1\t${(index + 1) * 100}\t.\tA\tC,G\t.\tPASS\t.\tGT\t${variant}`,
  );
  const header = [
    "##fileformat=VCFv4.4",
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
  ];
  return new TextEncoder().encode([...header, ...lines, ""].join("\n"));
}

/**
 * Where the pair of the individuals at `first` and at `second` sits in the
 * distance vector of `numIndividuals` individuals, which runs (0, 1),
 * (0, 2), ..., (1, 2), ...
 *
 * It counts the pairs of the rows before `first` one row at a time, so that
 * it does not repeat a formula of the code it checks.
 */
function indexOfThePair(
  first: number,
  second: number,
  numIndividuals: number,
): number {
  let index = 0;
  for (let row = 0; row < first; row += 1) {
    index += numIndividuals - 1 - row;
  }
  return index + (second - first - 1);
}

/** The distances of the panel, and the counts of the pass that gave them. */
function distancesOfThePanel(
  options: { minNumSnps?: number } = {},
  filter?: (variants: Variants) => void,
): Distances {
  const variants = openVcf(PANEL_VCF, { onlyPassed: false });
  try {
    filter?.(variants);
    return calcPairwiseKosmanDists(variants, options);
  } finally {
    variants.free();
  }
}

/** The distances of the worked example. */
function distancesOfTheWorkedExample(
  options: { minNumSnps?: number } = {},
): Distances {
  const variants = openVcf(WORKED_EXAMPLE, { onlyPassed: false });
  try {
    return calcPairwiseKosmanDists(variants, options);
  } finally {
    variants.free();
  }
}

/**
 * The distance of every pair that `gd.kosman` gave, in the order of the
 * distance vector, out of the reference file `name`.
 *
 * Its three columns are the distance, the variants both individuals were
 * called at and the ploidy times the sum of d; the first is the one popnei
 * gives a TypeScript user, and the header is checked so that a file written
 * with its columns in another order is not read as this one.
 */
async function distancesOfTheReference(name: string): Promise<number[]> {
  const text = new TextDecoder().decode(await referenceDists(name));
  const lines = text.trimEnd().split("\n");
  const header = lines.shift();
  assert.equal(header, "dist\tn\tk_sum");
  return lines.map((line) => Number(line.split("\t")[0]));
}

test("the five distances of the panel that the spec gives", () => {
  const distances = distancesOfThePanel();

  assert.equal(distances.names.length, PANEL_NUM_INDIVIDUALS);
  assert.equal(distances.names[0], "s000");
  assert.equal(
    distances.distVector.length,
    (PANEL_NUM_INDIVIDUALS * (PANEL_NUM_INDIVIDUALS - 1)) / 2,
  );
  for (const { first, second, dist } of PANEL_LITERALS) {
    const pair = indexOfThePair(
      distances.names.indexOf(first),
      distances.names.indexOf(second),
      PANEL_NUM_INDIVIDUALS,
    );
    assert.ok(
      Math.abs((distances.distVector[pair] as number) - dist) < TOLERANCE,
      `${first}, ${second}: ${String(distances.distVector[pair])} and not ${dist}`,
    );
  }
});

test("the counts of the pass over the panel are its 1200 variants", () => {
  const distances = distancesOfThePanel();

  const expected: PassStats = { numVars: PANEL_NUM_VARS, filtering: {} };
  assert.deepEqual(distances.passStats, expected);
});

test("a filter that keeps every variant of the panel counts them and changes no distance", () => {
  // A missing rate of at most 1 is every variant, whatever its missing
  // genotypes, so the pass is the one of the test above with the counts of
  // one filter beside it.
  const distances = distancesOfThePanel({}, (variants) => {
    variants.filterByMissingData(1);
  });

  const expected: PassStats = {
    numVars: PANEL_NUM_VARS,
    filtering: {
      missing_data: {
        varsProcessed: PANEL_NUM_VARS,
        varsKept: PANEL_NUM_VARS,
      },
    },
  };
  assert.deepEqual(distances.passStats, expected);
  for (const { first, second, dist } of PANEL_LITERALS) {
    const pair = indexOfThePair(
      distances.names.indexOf(first),
      distances.names.indexOf(second),
      PANEL_NUM_INDIVIDUALS,
    );
    assert.ok(
      Math.abs((distances.distVector[pair] as number) - dist) < TOLERANCE,
    );
  }
});

test("every pair of the panel is the distance R gave for it", async () => {
  const reference = await distancesOfTheReference("panel.gdkosman.tsv");
  const distances = distancesOfThePanel();

  assert.equal(reference.length, distances.distVector.length);
  let largestDifference = 0;
  for (const [pair, dist] of reference.entries()) {
    const difference = Math.abs((distances.distVector[pair] as number) - dist);
    largestDifference = Math.max(largestDifference, difference);
  }
  assert.ok(
    largestDifference < TOLERANCE,
    `the largest difference with gd.kosman is ${largestDifference}`,
  );
});

test("the three distances of the worked example", () => {
  const distances = distancesOfTheWorkedExample();

  assert.deepEqual(distances.names, ["s0", "s1", "s2"]);
  // The core divides two integers once, so these are the bits of 1/4, 5/6
  // and 2/6 and not a number near them.
  assert.deepEqual(Array.from(distances.distVector), WORKED_EXAMPLE_DISTS);
  const expected: PassStats = { numVars: 4, filtering: {} };
  assert.deepEqual(distances.passStats, expected);
});

test("minNumSnps leaves the pair below it without a distance", () => {
  // The pair s0, s1 was called at 2 variants and the other two at 3.
  const distances = distancesOfTheWorkedExample({ minNumSnps: 3 });

  assert.deepEqual(Array.from(distances.distVector), [
    Number.NaN,
    WORKED_EXAMPLE_DISTS[1],
    WORKED_EXAMPLE_DISTS[2],
  ]);
});

test("a minNumSnps above every pair leaves no distance", () => {
  const distances = distancesOfTheWorkedExample({ minNumSnps: 4 });

  assert.deepEqual(Array.from(distances.distVector), [
    Number.NaN,
    Number.NaN,
    Number.NaN,
  ]);
});

test("a minNumSnps of 0 is what no minNumSnps is", () => {
  const distances = distancesOfTheWorkedExample({ minNumSnps: 0 });

  assert.deepEqual(Array.from(distances.distVector), WORKED_EXAMPLE_DISTS);
});

test("a negative minNumSnps is refused and names what was given", () => {
  const variants = openVcf(WORKED_EXAMPLE, { onlyPassed: false });
  try {
    assert.throws(
      () => calcPairwiseKosmanDists(variants, { minNumSnps: -1 }),
      { message: /minNumSnps.*-1/ },
    );
  } finally {
    variants.free();
  }
});

test("squareDists has the pair with no distance in both of its cells", () => {
  const distances = distancesOfTheWorkedExample({ minNumSnps: 3 });

  const square = distances.squareDists();
  assert.equal(square.length, 9);
  assert.deepEqual(Array.from(square), [
    0,
    Number.NaN,
    WORKED_EXAMPLE_DISTS[1],
    Number.NaN,
    0,
    WORKED_EXAMPLE_DISTS[2],
    WORKED_EXAMPLE_DISTS[1],
    WORKED_EXAMPLE_DISTS[2],
    0,
  ]);
});

test("the three tetraploid distances that the spec gives", () => {
  // Every pair of a tetraploid dataset: the sums of d are divided by 4
  // times the variants of the pair, so a layer that took the ploidy for 2
  // would give twice these distances.
  const variants = openVcf(TETRAPLOID_VCF, { ploidy: 4, onlyPassed: false });
  let distances: Distances;
  try {
    distances = calcPairwiseKosmanDists(variants);
  } finally {
    variants.free();
  }

  assert.equal(distances.names.length, TETRAPLOID_NUM_INDIVIDUALS);
  const expected: PassStats = {
    numVars: TETRAPLOID_NUM_VARS,
    filtering: {},
  };
  assert.deepEqual(distances.passStats, expected);
  for (const { first, second, dist } of TETRAPLOID_LITERALS) {
    const pair = indexOfThePair(
      distances.names.indexOf(first),
      distances.names.indexOf(second),
      TETRAPLOID_NUM_INDIVIDUALS,
    );
    assert.ok(
      Math.abs((distances.distVector[pair] as number) - dist) < TOLERANCE,
      `${first}, ${second}: ${String(distances.distVector[pair])} and not ${dist}`,
    );
  }
});

test("the three distances of the haploid worked example", () => {
  // One allele in each genotype, so d is 0 for the same allele and 1 for
  // two different ones, and the third variant is called in h1 and h2 alone.
  const variants = openVcf(HAPLOID_WORKED_EXAMPLE, {
    ploidy: 1,
    onlyPassed: false,
  });
  let distances: Distances;
  try {
    distances = calcPairwiseKosmanDists(variants);
  } finally {
    variants.free();
  }

  assert.deepEqual(distances.names, ["h0", "h1", "h2"]);
  assert.deepEqual(Array.from(distances.distVector), [1 / 3, 2 / 3, 1 / 2]);
});

test("a pair that was never called in both has no distance", () => {
  // s0 and s1 are called at no variant in common, so their n is 0 and they
  // have no distance whatever minNumSnps is. The other two pairs keep
  // theirs: 0/0 against 0/1 is a d of 0.5 over one variant, and 1/1 against
  // 0/0 a d of 1 over one variant.
  const variants = openVcf(
    vcfOfThreeIndividuals(["s0", "s1", "s2"], ["0/0\t./.\t0/1", "./.\t1/1\t0/0"]),
    { onlyPassed: false },
  );
  let distances: Distances;
  try {
    distances = calcPairwiseKosmanDists(variants);
  } finally {
    variants.free();
  }

  assert.deepEqual(Array.from(distances.distVector), [Number.NaN, 0.5, 1]);
  assert.deepEqual(Array.from(distances.squareDists()), [
    0,
    Number.NaN,
    0.5,
    Number.NaN,
    0,
    1,
    0.5,
    1,
    0,
  ]);
});

test("the variants are as they were after the calculation", () => {
  // The call is a consumer: it makes one pass over the source and puts
  // nothing on the `Variants`, so the steps are the one that was there and
  // a second call reads the source again and gives the same vector.
  const variants = openVcf(WORKED_EXAMPLE, { onlyPassed: false });
  try {
    variants.filterByMissingData(1);
    const stepsBefore = variants.steps;

    const first = calcPairwiseKosmanDists(variants);
    assert.deepEqual(variants.steps, stepsBefore);
    const second = calcPairwiseKosmanDists(variants);

    assert.deepEqual(
      Array.from(second.distVector),
      Array.from(first.distVector),
    );
    assert.deepEqual(second.passStats, first.passStats);
  } finally {
    variants.free();
  }
});

test("distances of a vector that is not one value for each pair are refused", () => {
  // The vector of three individuals holds three values. Four is the length
  // that pyNei takes and puts the fourth value in no cell of the square
  // matrix, which "Its Python function" of the spec has among the
  // differences: popnei refuses it in both languages.
  assert.throws(
    () =>
      new Distances(Float64Array.from([0.1, 0.2, 0.3, 0.4]), ["a", "b", "c"], {
        numVars: 4,
        filtering: {},
      }),
    { message: /3 pairs, and 4 values/ },
  );
});

test("a source with no variant is refused and says that the source has none", () => {
  // A VCF whose header names three individuals and that has no data line.
  // The whole message is asserted, because it is the one the Python
  // function gives for the same source: the two languages say the same
  // thing, and Python writes the path of the file before it, which the bytes
  // a TypeScript user gives have not.
  const variants = openVcf(vcfOf([]));
  try {
    assert.throws(() => calcPairwiseKosmanDists(variants), {
      message:
        "the source has no variant, and a calculation needs 1 variant at least",
    });
  } finally {
    variants.free();
  }
});

test("a source with no variant names the counts of the filters that were on it", () => {
  // Every filter of the pass ran over the nothing the source gave, and the
  // counts say so: they are what a failed pass otherwise loses.
  const variants = openVcf(vcfOf([]));
  try {
    variants.filterByMaf(0.95);
    assert.throws(() => calcPairwiseKosmanDists(variants), {
      message:
        "the source has no variant, and a calculation needs 1 variant at " +
        "least: the filter `maf` was given 0 variants and kept 0",
    });
  } finally {
    variants.free();
  }
});

test("steps that keep no variant are refused with the counts of each filter", () => {
  // The major allele frequency of a variant is at least one over its
  // alleles, so a threshold of 0 keeps none of the 4 variants, and a missing
  // rate of at most 1 keeps every one of them: the two filters are named in
  // the order of the steps.
  const variants = openVcf(WORKED_EXAMPLE, { onlyPassed: false });
  try {
    variants.filterByMissingData(1);
    variants.filterByMaf(0);
    assert.throws(() => calcPairwiseKosmanDists(variants), {
      message:
        "the steps kept no variant of the 4 the source gave, and a " +
        "calculation needs 1 variant at least: the filter `missing_data` was " +
        "given 4 variants and kept 4, the filter `maf` was given 4 variants " +
        "and kept 0",
    });
  } finally {
    variants.free();
  }
});
