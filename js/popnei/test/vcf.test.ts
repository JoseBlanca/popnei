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

import type { Block, Field, IterBlocksOptions, Variants } from "popnei";
import { init, openVcf } from "popnei";

import loadTheWasm from "../wasm/popnei.js";
import { numberOfOpenPasses } from "../dist/variant.js";
import { manyVariantsVcf, referenceVcf, vcfOf } from "./reference.ts";

await init();

/** The WebAssembly of the core, for the tests that watch its memory. */
const wasm = await loadTheWasm();

/** How many bytes the memory of wasm holds. */
function memoryOfWasm(): number {
  return wasm.memory.buffer.byteLength;
}

/** Every field a block can carry besides the genotypes. */
const ALL_FIELDS: Field[] = ["chrom", "pos", "id", "alleles", "qual"];

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

test("a block carries the individuals and the ploidy its genotypes are read with", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  const [block] = [...variants.iterBlocks()];
  assert.ok(block !== undefined);
  assert.equal(block.numIndividuals, 3);
  assert.equal(block.ploidy, 2);
  // The alleles of the individual 1 of the variant 1, the `2|1` of chr1
  // 300, found with the two numbers of the block alone.
  const first = (1 * block.numIndividuals + 1) * block.ploidy;
  assert.deepEqual([...block.gts.subarray(first, first + block.ploidy)], [2, 1]);
  variants.free();
});

test("the default fields are the chromosomes and the positions", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  const [block] = [...variants.iterBlocks()];
  assert.ok(block !== undefined);
  assert.deepEqual(block.chrom, ["chr1", "chr1", "chr1"]);
  assert.deepEqual([...(block.pos ?? [])], [100, 300, 400]);
  assert.equal(block.id, null);
  assert.equal(block.alleles, null);
  assert.equal(block.qual, null);
  variants.free();
});

test("the columns of a block are copies, and the memory of wasm may grow", () => {
  const variants = openVcf(manyVariantsVcf(20));
  const [kept] = [...variants.iterBlocks({ fields: ["chrom", "pos"] })];
  assert.ok(kept !== undefined);
  const gts = [...kept.gts];
  const pos = [...(kept.pos ?? [])];
  const chrom = kept.chrom;
  variants.free();
  // A VCF of some megabytes, which the memory of wasm has to grow to hold.
  // A column that were a view into that memory would be detached by the
  // growth, and its length would be 0.
  const memory = memoryOfWasm();
  const big = openVcf(manyVariantsVcf(100000));
  for (const block of big.iterBlocks()) {
    assert.ok(block.numVars > 0);
  }
  big.free();
  assert.ok(
    memoryOfWasm() > memory,
    `the memory of wasm did not grow: ${memoryOfWasm()} bytes`,
  );
  assert.deepEqual([...kept.gts], gts);
  assert.deepEqual([...(kept.pos ?? [])], pos);
  assert.deepEqual(kept.chrom, chrom);
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
  // `depth` does not compile, which is what the type of `fields` is for,
  // and a user of JavaScript, who has no types, gets this error.
  const fields = ["chrom", "depth"] as Field[];
  assert.throws(() => variants.iterBlocks({ fields }), {
    name: "Error",
    message: /depth/,
  });
  variants.free();
});

test("blocks of no variant are refused at the call", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.throws(() => variants.iterBlocks({ numVarsPerBlock: 0 }), {
    name: "Error",
    message: /`numVarsPerBlock` is a whole number of 1 or more/,
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
  const before = numberOfOpenPasses();
  const variants = openVcf(await referenceVcf("cases.vcf"));
  assert.equal([...variants.iterBlocks()].length, 1);
  assert.equal(numberOfOpenPasses(), before);
  variants.free();
});

test("an iteration that is left with a break gives back its pass too", async () => {
  const before = numberOfOpenPasses();
  const variants = openVcf(await referenceVcf("many.vcf"));
  for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
    assert.equal(block.numVars, 100);
    assert.equal(numberOfOpenPasses(), before + 1);
    break;
  }
  assert.equal(numberOfOpenPasses(), before);
  variants.free();
});

test("an iteration that throws gives back its pass, and throws its error", () => {
  const before = numberOfOpenPasses();
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1\t1/1",
    ]),
  );
  // The error a user reads is the one of the core, the genotype of another
  // ploidy of the second line, and not one of the free of the pass that the
  // `finally` of the iteration does after it.
  assert.throws(() => [...variants.iterBlocks({ numVarsPerBlock: 1 })], {
    name: "Error",
    message: /line 4 of the VCF, the column of ind1/,
  });
  assert.equal(numberOfOpenPasses(), before);
  variants.free();
});

test("an iterator that is never started holds its pass", async () => {
  const before = numberOfOpenPasses();
  const variants = openVcf(await referenceVcf("cases.vcf"));
  variants.iterBlocks();
  // The pass is open in the memory of wasm and no `finally` will free it:
  // a generator that never ran its first line never runs its last. What
  // gives it back is the `FinalizationRegistry` of wasm-bindgen, when the
  // garbage collector reaches the iterator.
  assert.equal(numberOfOpenPasses(), before);
  variants.free();
});

test("a source that is not a Uint8Array is refused", async () => {
  const text = new TextDecoder().decode(await referenceVcf("cases.vcf"));
  // The bytes of the file as text, which is what readFileSync(path, "utf8")
  // gives: the core would read the memory of the string and say that the
  // source starts with bytes the user never had.
  assert.throws(() => openVcf(text as unknown as Uint8Array), {
    name: "Error",
    message: /Uint8Array/,
  });
  assert.throws(() => openVcf(null as unknown as Uint8Array), {
    name: "Error",
    message: /null/,
  });
});

test("a ploidy that is not a whole number of 1 or more is refused", () => {
  const bytes = vcfOf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"]);
  // 2.5 and 2^32 + 2 both reached the core as a ploidy of 2, and -1 as
  // 4294967295.
  for (const ploidy of [2.5, 4294967298, -1, 0, Number.NaN]) {
    assert.throws(() => openVcf(bytes, { ploidy }), {
      name: "Error",
      message: /`ploidy` is a whole number/,
    });
  }
});

test("an onlyPassed that is not a boolean is refused", () => {
  const bytes = vcfOf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"]);
  assert.throws(
    () => openVcf(bytes, { onlyPassed: 0 as unknown as boolean }),
    { name: "Error", message: /`onlyPassed` is true or false/ },
  );
});

test("a numVarsPerBlock that is not a whole number of 1 or more is refused", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  // 2^32 + 1 gave blocks of one variant, and NaN blocks of none.
  for (const numVarsPerBlock of [4294967297, -1, 0, 1.5, Number.NaN]) {
    assert.throws(() => variants.iterBlocks({ numVarsPerBlock }), {
      name: "Error",
      message: /`numVarsPerBlock` is a whole number/,
    });
  }
  variants.free();
});

test("one field written where the array of fields goes is refused", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  // A string spread into an array is its letters, and popnei would look for
  // a field called `a`.
  assert.throws(
    () => variants.iterBlocks({ fields: "alleles" as unknown as Field[] }),
    { name: "Error", message: /array of names/ },
  );
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
  // The steps live in the memory of wasm beside the source, so they go
  // with it.
  assert.throws(() => variants.steps, {
    name: "Error",
    message: /freed/,
  });
  // A second free is not an error: it has nothing left to give back.
  variants.free();
});

test("a Variants is freed by the using of a block too", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  // What `using variants = openVcf(...)` calls at the end of its block.
  variants[Symbol.dispose]();
  assert.throws(() => variants.iterBlocks(), {
    name: "Error",
    message: /freed/,
  });
});
