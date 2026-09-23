/**
 * The kinship from TypeScript: `calcKinship` and the `Kinship` it gives.
 *
 * "How it is verified" of `docs/specs/kinship.md` has the numbers. The first
 * dataset is `tests/reference/kinship/panel_called.vcf.gz`, 200 diploid
 * individuals, `s000` to `s199`, and 1200 biallelic variants with every
 * genotype called; the five entries of `s000` that the table of the spec
 * gives are the literals here, and plink2 v2.0.0-a.7.7 wrote them. The
 * second is the worked example of the same section, 4 variants of 4 diploid
 * individuals, written as a VCF here: two of its variants have no variance
 * and are left out, one genotype is missing, and the matrix of the two that
 * are left is whole numbers.
 *
 * The calculation is tested in the core crate, over every entry of both
 * reference panels. What these tests say is that the matrix reaches
 * TypeScript with the shape and the order it has in the core, that the
 * individuals a user names are the ones the frequencies are of, and that an
 * error of the core is thrown as an `Error` with the message it has in Rust.
 *
 * Nothing here computes an expected value with popnei: every number comes
 * from the spec, which has them from plink2 and from pyNei.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Variants } from "popnei";
import { calcKinship, init, Kinship, openVcf } from "popnei";

import { referenceKinship } from "./reference.ts";

await init();

/** How many individuals the panel holds, `s000` to `s199`. */
const PANEL_NUM_INDIVIDUALS = 200;

/** How many variants of the panel have variance and are used. */
const PANEL_NUM_VARS = 1200;

/**
 * The tolerance of an entry of the panel: plink2 writes six significant
 * digits, so an entry near 1 is rounded by up to 5e-6, and "How it is
 * verified" of the spec compares within 1e-5 absolute, one unit of the last
 * digit it prints for an entry of that size.
 */
const OF_PLINK2 = 1e-5;

/**
 * The tolerance of the worked example, whose entries are whole numbers: run
 * through pyNei, every one of them is a whole number within 4.4e-16, and the
 * spec asserts them within 1e-12 absolute.
 */
const OF_THE_WORKED_EXAMPLE = 1e-12;

/** The bytes of the panel, read once for every test that runs on it. */
const PANEL_VCF = await referenceKinship("panel_called.vcf.gz");

/**
 * The five entries of `s000` that plink2 gave, the rows of the table of
 * "How it is verified" that name `s000` on `panel_called`: its diagonal, its
 * two full sibs `s001` and `s002`, the unrelated `s004` and the last
 * individual of the panel, `s199`.
 */
const OF_S000: { other: string; ofPlink2: number }[] = [
  { other: "s000", ofPlink2: 1.09309 },
  { other: "s001", ofPlink2: 0.648081 },
  { other: "s002", ofPlink2: 0.615611 },
  { other: "s004", ofPlink2: -0.0945533 },
  { other: "s199", ofPlink2: -0.0760273 },
];

/**
 * The worked example of the spec: 4 diploid individuals, `i0` to `i3`, and 4
 * variants. `v2`, where every individual is heterozygous, and `v3`, which
 * has one allele, have no variance and are left out; the genotype `./.` of
 * `i2` at `v1` takes the mean dosage of its variant and is in the
 * denominator of no pair.
 */
const WORKED_EXAMPLE = vcfOf(
  ["i0", "i1", "i2", "i3"],
  [
    "0/0\t0/1\t1/1\t0/1",
    "0/0\t0/1\t./.\t1/1",
    "0/1\t0/1\t0/1\t0/1",
    "0/0\t0/0\t0/0\t0/0",
  ],
);

/**
 * The matrix of the worked example, row after row: the denominator of every
 * pair with `i2` is 1 where the others have 2, which is why the -2 of `i0`
 * and `i2` stays -2 and the -2 of `i0` and `i3` becomes -1.
 */
const OF_THE_WORKED_EXAMPLE_MATRIX = [
  2, 0, -2, -1, 0, 0, 0, 0, -2, 0, 2, 0, -1, 0, 0, 1,
];

/**
 * The matrix of the same four variants over `i0` and `i3` alone, which is
 * what `individuals` asks for: the frequencies are those two individuals',
 * so `v0`, whose dosages there are 0 and 1, is divided by
 * `sqrt(2 * 0.25 * 0.75)` where over the four it is divided by
 * `sqrt(2 * 0.5 * 0.5)`. pyNei at commit ef0ca6e gave it.
 */
const OF_I0_AND_I3 = [4 / 3, -4 / 3, -4 / 3, 4 / 3];

/**
 * A VCF of the individuals `names`, one line for each text of `genotypes`,
 * which holds the genotype of each of them.
 *
 * Every variant declares the two alleles `A` and `C`, and no genotype is
 * written with a ploidy of its own: the ploidy is what `openVcf` is told.
 */
function vcfOf(
  names: readonly string[],
  genotypes: readonly string[],
): Uint8Array {
  const lines = genotypes.map(
    (variant, index) =>
      `chr1\t${(index + 1) * 100}\t.\tA\tC\t.\tPASS\t.\tGT\t${variant}`,
  );
  const header = [
    "##fileformat=VCFv4.4",
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
  ];
  return new TextEncoder().encode([...header, ...lines, ""].join("\n"));
}

/** The kinship of `bytes`, with `options` as a user would give them. */
function kinshipOf(
  bytes: Uint8Array,
  options: {
    individuals?: readonly string[];
    transformToBiallelic?: boolean;
  } = {},
): Kinship {
  const variants = openVcf(bytes, { onlyPassed: false });
  try {
    return calcKinship(variants, options);
  } finally {
    variants.free();
  }
}

/**
 * The kinship of `bytes` over the variants a step leaves, with `filter` put
 * on the `Variants` before the call.
 */
function kinshipAfter(
  bytes: Uint8Array,
  filter: (variants: Variants) => void,
): Kinship {
  const variants = openVcf(bytes, { onlyPassed: false });
  try {
    filter(variants);
    return calcKinship(variants);
  } finally {
    variants.free();
  }
}

/** The entry of the individuals `one` and `other` of `kinship`. */
function entryOf(kinship: Kinship, one: string, other: string): number {
  const row = kinship.individuals.indexOf(one);
  const column = kinship.individuals.indexOf(other);
  assert.ok(row >= 0 && column >= 0, `${one} or ${other} is not in the matrix`);
  return kinship.matrix[row * kinship.individuals.length + column] as number;
}

test("the five entries of s000 of the panel that plink2 gives", () => {
  const kinship = kinshipOf(PANEL_VCF);

  assert.equal(kinship.individuals.length, PANEL_NUM_INDIVIDUALS);
  assert.equal(kinship.individuals[0], "s000");
  assert.equal(kinship.numVars, PANEL_NUM_VARS);
  assert.equal(
    kinship.matrix.length,
    PANEL_NUM_INDIVIDUALS * PANEL_NUM_INDIVIDUALS,
  );
  for (const { other, ofPlink2 } of OF_S000) {
    const entry = entryOf(kinship, "s000", other);
    assert.ok(
      Math.abs(entry - ofPlink2) < OF_PLINK2,
      `the entry of s000 and ${other} is ${entry} and plink2 gives ${ofPlink2}`,
    );
  }
});

test("the counts of the pass of the panel are its 1200 variants", () => {
  const kinship = kinshipOf(PANEL_VCF);

  assert.equal(kinship.passStats?.numVars, PANEL_NUM_VARS);
  assert.deepEqual(kinship.passStats?.filtering, {});
});

test("the worked example gives the whole numbers of the spec", () => {
  const kinship = kinshipOf(WORKED_EXAMPLE);

  assert.deepEqual(kinship.individuals, ["i0", "i1", "i2", "i3"]);
  assert.equal(kinship.numVars, 2);
  for (const [at, expected] of OF_THE_WORKED_EXAMPLE_MATRIX.entries()) {
    const entry = kinship.matrix[at] as number;
    assert.ok(
      Math.abs(entry - expected) < OF_THE_WORKED_EXAMPLE,
      `the entry ${at} of the matrix is ${entry} and the spec gives ${expected}`,
    );
  }
});

test("the kinship of two individuals takes their own frequencies", () => {
  const kinship = kinshipOf(WORKED_EXAMPLE, { individuals: ["i0", "i3"] });

  assert.deepEqual(kinship.individuals, ["i0", "i3"]);
  assert.equal(kinship.numVars, 2);
  for (const [at, expected] of OF_I0_AND_I3.entries()) {
    const entry = kinship.matrix[at] as number;
    assert.ok(
      Math.abs(entry - expected) < OF_THE_WORKED_EXAMPLE,
      `the entry ${at} of the matrix is ${entry} and pyNei gives ${expected}`,
    );
  }
});

test("an individual that is not in the dataset is refused by its name", () => {
  assert.throws(
    () => kinshipOf(WORKED_EXAMPLE, { individuals: ["i0", "i9"] }),
    {
      message: /i9/,
    },
  );
});

test("filterIndividuals takes the rows and the columns and calculates nothing", () => {
  const kinship = kinshipOf(WORKED_EXAMPLE);

  const some = kinship.filterIndividuals(["i3", "i0"]);

  assert.deepEqual(some.individuals, ["i3", "i0"]);
  assert.equal(some.numVars, kinship.numVars);
  assert.equal(some.passStats?.numVars, kinship.passStats?.numVars);
  // The rows and the columns of the matrix of the four, turned around, and
  // not the kinship of the two, which is `OF_I0_AND_I3`.
  for (const [at, expected] of [1, -1, -1, 2].entries()) {
    const entry = some.matrix[at] as number;
    assert.ok(
      Math.abs(entry - expected) < OF_THE_WORKED_EXAMPLE,
      `the entry ${at} of the matrix is ${entry} and the four rows give ${expected}`,
    );
  }
});

test("filterIndividuals refuses a name that is of nobody in the matrix", () => {
  const kinship = kinshipOf(WORKED_EXAMPLE);

  assert.throws(() => kinship.filterIndividuals(["i0", "i9"]), {
    message: /i9/,
  });
});

test("a kinship a user builds by hand has no counts of a pass", () => {
  const built = new Kinship(
    Float64Array.from(OF_THE_WORKED_EXAMPLE_MATRIX),
    ["i0", "i1", "i2", "i3"],
    2,
  );

  assert.equal(built.passStats, undefined);
  assert.equal(built.numVars, 2);
  assert.equal(built.individuals.length, 4);
});

test("a matrix that is not symmetric is refused, naming the two individuals", () => {
  const notSymmetric = Float64Array.from([1, 0.5, 0.4, 1]);

  assert.throws(() => new Kinship(notSymmetric, ["i0", "i1"], 2), {
    message: /a kinship is symmetric, and the pair of `i0` and `i1` is 0.5/,
  });
});

test("a matrix that is not one value for each pair is refused", () => {
  assert.throws(
    () => new Kinship(Float64Array.from([1, 0.5, 0.5]), ["i0", "i1"], 2),
    { message: /4 values, and 3 were given/ },
  );
});

test("a variant of more than two alleles is refused unless it is read as biallelic", () => {
  const threeAlleles = new TextEncoder().encode(
    [
      "##fileformat=VCFv4.4",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3",
      "chr1\t100\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/0\t0/1\t1/2\t0/1",
      "chr1\t200\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1",
      "",
    ].join("\n"),
  );

  assert.throws(() => kinshipOf(threeAlleles), { message: /alleles/ });

  const read = kinshipOf(threeAlleles, { transformToBiallelic: true });

  assert.equal(read.numVars, 2);
});

test("the counts of the pass are the variants it gave and not the ones used", () => {
  // `filterByMissingData(0)` drops `v1`, the one with the missing genotype,
  // so the pass gives `v0`, `v2` and `v3`, of which `v2`, where every
  // individual is heterozygous, and `v3`, which has one allele, have no
  // variance. The two counts are 3 and 1: a count that read the variants
  // that were used would give 1 here, and one that read the source 4, and
  // on a pass with no filter all three are the same number.
  const kinship = kinshipAfter(WORKED_EXAMPLE, (variants) => {
    variants.filterByMissingData(0);
  });

  assert.equal(kinship.passStats?.numVars, 3);
  assert.equal(kinship.numVars, 1);
  assert.deepEqual(kinship.passStats?.filtering, {
    missing_data: { varsProcessed: 4, varsKept: 3 },
  });
});

test("a pair with no variant called in both is refused by their names", () => {
  // Three individuals and three variants of "Missing genotypes, variants
  // with no variance, and what pyNei asserts" of the spec: `i0` and `i2`
  // are never called at the same variant, so their entry of the kinship
  // would be divided by 0. The core names the positions the two have among
  // the individuals of the kinship, which with `individuals` on the call
  // are not even the ones of the file, and a user drops a name.
  const neverTogether = vcfOf(
    ["i0", "i1", "i2"],
    ["0/0\t0/1\t./.", "./.\t0/1\t1/1", "0/0\t1/1\t./."],
  );

  let thrown = "";
  try {
    kinshipOf(neverTogether);
  } catch (error) {
    thrown = (error as Error).message;
  }

  assert.match(
    thrown,
    /the individuals `i0` and `i2` have no variant called in both of them/,
  );
  assert.match(
    thrown,
    /2 variants are called in the first and 1 in the second/,
  );
  assert.ok(
    !thrown.includes("at the positions"),
    `the message writes a position: ${thrown}`,
  );
});

test("an entry that is not a number is refused, naming its cell", () => {
  // On the diagonal, where the check of the symmetry compares no pair of
  // cells and would let it through.
  const withANaN = Float64Array.from([1, 0, 0, NaN]);

  assert.throws(() => new Kinship(withANaN, ["i0", "i1"], 2), {
    message: /entry of the kinship of `i1` and `i1` is NaN/,
  });
});

test("an infinity is refused as a value and not as an asymmetry", () => {
  const withAnInfinity = Float64Array.from([1, Infinity, Infinity, 1]);

  assert.throws(() => new Kinship(withAnInfinity, ["i0", "i1"], 2), {
    message: /entry of the kinship of `i0` and `i1` is Infinity/,
  });
});

test("a kinship of no individual is refused", () => {
  assert.throws(() => new Kinship(Float64Array.from([]), [], 2), {
    message: /no individual was named/,
  });
});

test("a numVars that is not a whole number of 0 or more is refused", () => {
  assert.throws(() => new Kinship(Float64Array.from([1]), ["i0"], -1), {
    message: /`numVars` is a whole number of 0 or more/,
  });
});

test("a matrix that is not a Float64Array says what was given", () => {
  assert.throws(
    () => new Kinship([1, 0, 0, 1] as unknown as Float64Array, ["i0", "i1"], 2),
    { message: /an object of the type `Array` was given/ },
  );
});
