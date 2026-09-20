/**
 * `openVcf` and `iterBlocks`: what a TypeScript user reads from a VCF.
 *
 * The cases are the ones the two specs give to the test under node: the four
 * variants of `cases.vcf` and the two of `differences.vcf` of
 * `docs/specs/io_vcf.md`, with the default and with every variant given, and
 * the blocks of three variants of `docs/specs/block.md`. The numbers are the
 * literals of those tables, which come from bcftools 1.24, and the files are
 * the ones the Python tests read.
 *
 * It imports `popnei`, the name of this package, which node resolves to the
 * built `dist/node.js`, so the entry point a user of node gets is the one
 * that is tested. `dist/variant.js`, which no entry point re-exports, is
 * imported by its path for the count of the passes that hold memory of wasm.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Block, IterBlocksOptions, Variants } from "popnei";
import { init, openVcf } from "popnei";

import { numberOfOpenPasses } from "../dist/variant.js";
import { referenceVcf, vcfOf } from "./reference.ts";

await init();

/** Every field a block can carry besides the genotypes. */
const ALL_FIELDS = ["chrom", "pos", "id", "alleles", "qual"];

/** One variant of a block, with the fields the tables of the specs give. */
interface Row {
  chrom: string;
  pos: number;
  id: string | null;
  alleles: string[];
  qual: number | null;
  gts: number[];
}

/** The four variants of `cases.vcf`, the table of `docs/specs/io_vcf.md`. */
const CASES: Row[] = [
  {
    chrom: "chr1",
    pos: 100,
    id: "rs1",
    alleles: ["A", "T"],
    qual: 29.5,
    gts: [0, 0, 0, 1, 1, 1],
  },
  {
    chrom: "chr1",
    pos: 200,
    id: null,
    alleles: ["A", "T"],
    qual: null,
    gts: [-1, -1, 0, 1, -1, 0],
  },
  {
    chrom: "chr1",
    pos: 300,
    id: null,
    alleles: ["A", "G", "T"],
    qual: 67,
    gts: [1, 2, 2, 1, 2, 2],
  },
  {
    chrom: "chr1",
    pos: 400,
    id: null,
    alleles: ["T"],
    qual: 47,
    gts: [0, 0, 0, 0, 0, 0],
  },
];

/** The FILTER of the second variant of `cases.vcf` is `q10`, and no other. */
const CASES_THAT_PASSED = CASES.filter((_row, variant) => variant !== 1);

/** The two variants of `differences.vcf`, both with `PASS`. */
const DIFFERENCES: Row[] = [
  {
    chrom: "chr2",
    pos: 50,
    id: "ms1",
    alleles: ["GTC", "G", "GTCT"],
    qual: 50,
    gts: [0, 1, 0, 2, -1, -1],
  },
  {
    chrom: "chr2",
    pos: 60,
    id: null,
    alleles: ["A", "<DEL>", "*"],
    qual: null,
    gts: [0, 1, 2, 2, 0, 0],
  },
];

/** The variants of every block of `variants`, one after another. */
function rowsOf(
  variants: Variants,
  options: IterBlocksOptions = { fields: ALL_FIELDS },
): Row[] {
  const rows: Row[] = [];
  for (const block of variants.iterBlocks(options)) {
    rows.push(...rowsOfTheBlock(block, variants.numIndividuals * variants.ploidy));
  }
  return rows;
}

/** The variants of one block, `allelesPerVariant` genotype numbers each. */
function rowsOfTheBlock(block: Block, allelesPerVariant: number): Row[] {
  const { chrom, pos, id, alleles, qual } = block;
  assert.ok(chrom !== null && pos !== null && id !== null);
  assert.ok(alleles !== null && qual !== null);
  return Array.from({ length: block.numVars }, (_unused, variant) => ({
    chrom: chrom[variant] as string,
    pos: pos[variant] as number,
    id: id[variant] ?? null,
    alleles: alleles[variant] as string[],
    // The core gives NaN to a variant with no quality, which the table of
    // the spec has as none.
    qual: Number.isNaN(qual[variant]) ? null : (qual[variant] as number),
    gts: [
      ...block.gts.subarray(
        variant * allelesPerVariant,
        (variant + 1) * allelesPerVariant,
      ),
    ],
  }));
}

test("the individuals and the ploidy are known when the VCF is opened", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  assert.equal(variants.numIndividuals, 3);
  assert.equal(variants.ploidy, 2);
  variants.free();
});

test("every variant of cases.vcf is given with onlyPassed false", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"), {
    onlyPassed: false,
  });
  assert.deepEqual(rowsOf(variants), CASES);
  variants.free();
});

test("by default the variant of cases.vcf that failed a filter is left out", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.deepEqual(rowsOf(variants), CASES_THAT_PASSED);
  variants.free();
});

test("the gzipped cases.vcf gives what the plain one gives", async () => {
  const plain = openVcf(await referenceVcf("cases.vcf"));
  const gzipped = openVcf(await referenceVcf("cases.vcf.gz"));
  assert.deepEqual(rowsOf(gzipped), rowsOf(plain));
  plain.free();
  gzipped.free();
});

test("the two variants of differences.vcf are read as bcftools reads them", async () => {
  const bytes = await referenceVcf("differences.vcf");
  const byDefault = openVcf(bytes);
  const everyVariant = openVcf(bytes, { onlyPassed: false });
  // Both variants have `PASS`, so the default and every variant give the
  // same two rows, one with a leading separator in two of its genotypes.
  assert.deepEqual(rowsOf(byDefault), DIFFERENCES);
  assert.deepEqual(rowsOf(everyVariant), DIFFERENCES);
  byDefault.free();
  everyVariant.free();
});

test("the blocks are cut by the count of the variants", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"), {
    onlyPassed: false,
  });
  const blocks = [
    ...variants.iterBlocks({ fields: ALL_FIELDS, numVarsPerBlock: 3 }),
  ];
  assert.deepEqual(
    blocks.map((block) => block.numVars),
    [3, 1],
  );
  assert.deepEqual(
    blocks.flatMap((block) => rowsOfTheBlock(block, 6)),
    CASES,
  );
  variants.free();
});

test("the genotypes alone leave every other column out", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  const [block] = [...variants.iterBlocks({ fields: [] })];
  assert.ok(block !== undefined);
  assert.equal(block.numVars, 3);
  assert.equal(block.chrom, null);
  assert.equal(block.pos, null);
  assert.equal(block.id, null);
  assert.equal(block.alleles, null);
  assert.equal(block.qual, null);
  assert.deepEqual([...block.gts], [0, 0, 0, 1, 1, 1, 1, 2, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0]);
  variants.free();
});

test("the chromosome and the position travel together", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  const [block] = [...variants.iterBlocks({ fields: ["pos"] })];
  assert.ok(block !== undefined);
  assert.deepEqual(block.chrom, ["chr1", "chr1", "chr1"]);
  assert.deepEqual([...(block.pos ?? [])], [100, 300, 400]);
  variants.free();
});

test("every pass over the variants reads the source again", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  const first = rowsOf(variants);
  const second = rowsOf(variants, { fields: ALL_FIELDS, numVarsPerBlock: 2 });
  assert.deepEqual(second, first);
  variants.free();
});

test("a field that is not one of the five is refused at the call", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.throws(() => variants.iterBlocks({ fields: ["chrom", "depth"] }), {
    name: "Error",
    message: /depth/,
  });
  variants.free();
});

test("blocks of no variant are refused at the call", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.throws(() => variants.iterBlocks({ numVarsPerBlock: 0 }), {
    name: "Error",
    message: /0 variants/,
  });
  variants.free();
});

test("bytes that are not a VCF are refused when they are opened", () => {
  assert.throws(() => openVcf(new TextEncoder().encode("chr1\t100\n")), {
    name: "Error",
    message: /not a VCF/,
  });
});

test("a genotype of another ploidy is an error at the block that holds it", () => {
  const tetraploid = vcfOf([
    "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1/1/1\t0/0/0/0",
  ]);
  const variants = openVcf(tetraploid);
  assert.throws(() => [...variants.iterBlocks()], {
    name: "Error",
    message: /ind1/,
  });
  variants.free();
});

test("the ploidy of the reader is the one that was asked for", () => {
  const tetraploid = vcfOf([
    "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1/1/1\t0/0/0/0",
  ]);
  const variants = openVcf(tetraploid, { ploidy: 4 });
  assert.equal(variants.ploidy, 4);
  const [block] = [...variants.iterBlocks()];
  assert.ok(block !== undefined);
  assert.deepEqual([...block.gts], [0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0]);
  variants.free();
});

test("a VCF with no variant gives no block", () => {
  const variants = openVcf(vcfOf([]));
  assert.deepEqual([...variants.iterBlocks()], []);
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  variants.free();
});

test("an iteration that ends gives back the memory of wasm of its pass", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.equal([...variants.iterBlocks()].length, 1);
  assert.equal(numberOfOpenPasses(), 0);
  variants.free();
});

test("an iteration that is left with a break gives back its pass too", async () => {
  const variants = openVcf(await referenceVcf("many.vcf"));
  for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
    assert.equal(block.numVars, 100);
    break;
  }
  assert.equal(numberOfOpenPasses(), 0);
  variants.free();
});

test("an iteration that throws gives back its pass too", () => {
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1\t1/1",
    ]),
  );
  assert.throws(() => [...variants.iterBlocks({ numVarsPerBlock: 1 })], {
    name: "Error",
  });
  assert.equal(numberOfOpenPasses(), 0);
  variants.free();
});

test("a position above 2^53 is refused, because a float64 rounds it", () => {
  // 9007199254740993 is 2^53 + 1, the first whole number a float64 does not
  // hold: it would reach a user as 9007199254740992, where a Python user of
  // the same file reads the number the source has.
  const variants = openVcf(
    vcfOf(["chr1\t9007199254740993\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"]),
  );
  assert.throws(() => [...variants.iterBlocks()], {
    name: "Error",
    message: /9007199254740993/,
  });
  variants.free();
});

test("a position of 2^53 is read, the largest a float64 holds", () => {
  const variants = openVcf(
    vcfOf(["chr1\t9007199254740992\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"]),
  );
  const [block] = [...variants.iterBlocks()];
  assert.ok(block !== undefined);
  assert.deepEqual([...(block.pos ?? [])], [9007199254740992]);
  variants.free();
});

test("variants that were freed cannot be read again", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  variants.free();
  // The names and the ploidy are in JavaScript, so they still answer.
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  assert.equal(variants.ploidy, 2);
  assert.throws(() => variants.iterBlocks(), {
    name: "Error",
    message: /freed/,
  });
  // A second free is not an error: it has nothing left to give back.
  variants.free();
});
