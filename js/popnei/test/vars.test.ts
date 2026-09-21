/**
 * `writeVars` and `openVars`: the vars file a TypeScript user writes and
 * reads back.
 *
 * The round trip is the one `docs/specs/io_vars.md` gives to the test under
 * node: `cases.vcf` read with every variant, written with batches of 3
 * variants and read back, compared with the table of the four variants of
 * `docs/specs/io_vcf.md`, whose numbers come from bcftools 1.24 and are the
 * literals below. The file that popnei cannot write, compressed with zstd,
 * is `tests/reference/vars/zstd.vars`, which the Python tests read too.
 *
 * It imports `popnei`, the name of this package, which node resolves to the
 * built `dist/node.js`. `dist/variant.js`, which no entry point re-exports,
 * is imported by its path for the count of the passes that hold memory of
 * wasm.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Block, Field, Variants } from "popnei";
import { init, openVars, openVcf, writeVars } from "popnei";

import loadTheWasm, { room_for_bytes as roomForBytes } from "../wasm/popnei.js";
import { numberOfOpenPasses } from "../dist/variant.js";
import {
  manyVariantsVcf,
  referenceVars,
  referenceVcf,
  vcfOf,
} from "./reference.ts";

await init();

/** The WebAssembly of the core, for the tests that watch its memory. */
const wasm = await loadTheWasm();

/** How many bytes the memory of wasm holds. */
function memoryOfWasm(): number {
  return wasm.memory.buffer.byteLength;
}

/** Every field a block can carry besides the genotypes. */
const ALL_FIELDS: Field[] = ["chrom", "pos", "id", "alleles", "qual"];

/** The six bytes an arrow IPC file starts and ends with. */
const ARROW1 = "ARROW1";

/** One variant of a block, with the fields the table of the spec gives. */
interface Row {
  chrom: string;
  pos: number;
  id: string | null;
  alleles: string[];
  qual: number | null;
  gts: number[];
}

/**
 * The four variants of `cases.vcf`, the table of `docs/specs/io_vcf.md`,
 * which is what the file written from that VCF holds.
 */
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

/** The variants of every block of `variants`, one after another. */
function rowsOf(variants: Variants, numVarsPerBlock?: number): Row[] {
  const rows: Row[] = [];
  for (const block of variants.iterBlocks({
    fields: ALL_FIELDS,
    numVarsPerBlock,
  })) {
    rows.push(
      ...rowsOfTheBlock(block, variants.numIndividuals * variants.ploidy),
    );
  }
  return rows;
}

/**
 * The bytes of the vars file of the four variants of `cases.vcf`, batches
 * of `numVarsPerBlock` variants.
 */
async function varsFileOfCases(numVarsPerBlock: number): Promise<Uint8Array> {
  const variants = openVcf(await referenceVcf("cases.vcf"), {
    onlyPassed: false,
  });
  try {
    return writeVars(variants, { numVarsPerBlock }).bytes;
  } finally {
    variants.free();
  }
}

/** What `bytes` says at `first`, as text: the mark of an arrow file. */
function textAt(bytes: Uint8Array, first: number, length: number): string {
  return new TextDecoder().decode(bytes.subarray(first, first + length));
}

/**
 * How many variants a batch of the vars file in `bytes` holds, read from
 * the `popnei` key of its schema.
 *
 * The keys of the schema of an arrow file are text that nothing compresses,
 * so the json is in the file as it was written and this is what says at
 * which size the batches were written, which the blocks that come out of a
 * pass do not: they are cut to the size the pass asks for. Each byte is one
 * character here, which is what `latin1` gives, so the bytes of the
 * genotypes do not move what is found.
 */
function numVarsPerBlockOf(bytes: Uint8Array): number {
  const text = new TextDecoder("latin1").decode(bytes);
  const written = /"num_vars_per_block":(\d+)/.exec(text)?.[1];
  assert.ok(written !== undefined, "the file has no `popnei` key");
  return Number(written);
}

test("the four variants of cases.vcf come back from the vars file", async () => {
  const bytes = await varsFileOfCases(3);
  const variants = openVars(bytes);
  // The batches of the file hold 3 variants and 1, and a pass that asks for
  // blocks of 3 gets them as they are: the two are cut nowhere else.
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

test("the blocks of a vars file are joined to the size that is asked for", async () => {
  const variants = openVars(await varsFileOfCases(3));
  // With no size the blocks are the ones popnei chooses for three
  // individuals, which holds the two batches of the file together.
  const blocks = [...variants.iterBlocks({ fields: ALL_FIELDS })];
  assert.deepEqual(
    blocks.map((block) => block.numVars),
    [4],
  );
  assert.deepEqual(
    blocks.flatMap((block) => rowsOfTheBlock(block, 6)),
    CASES,
  );
  variants.free();
});

test("the bytes of a vars file are those of an arrow file", async () => {
  const bytes = await varsFileOfCases(3);
  assert.equal(textAt(bytes, 0, ARROW1.length), ARROW1);
  assert.equal(textAt(bytes, bytes.length - ARROW1.length, ARROW1.length), ARROW1);
});

test("the batches of the file hold the variants that were asked for", async () => {
  // What a pass gives says nothing about this: it cuts the blocks to the
  // size it was asked for whatever the file holds. The `popnei` key of the
  // schema is what the file itself says.
  assert.equal(numVarsPerBlockOf(await varsFileOfCases(3)), 3);
  assert.equal(numVarsPerBlockOf(await varsFileOfCases(1)), 1);
  const variants = openVcf(await referenceVcf("cases.vcf"), {
    onlyPassed: false,
  });
  // With no size it is the one popnei chooses for three diploid
  // individuals, which holds the four variants in one batch.
  const chosen = numVarsPerBlockOf(writeVars(variants).bytes);
  assert.ok(chosen > 4, `popnei chose batches of ${chosen} variants`);
  variants.free();
});

test("the bytes of a vars file are a copy, and the memory of wasm may grow", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"), {
    onlyPassed: false,
  });
  const bytes = writeVars(variants, { numVarsPerBlock: 3 }).bytes;
  const kept = [...bytes];
  variants.free();
  assert.ok(bytes.length > 0);
  // A VCF of 13 MB, more than what the tests before this one left free in
  // the memory of wasm, so that memory has to grow to hold it. Bytes that
  // were a view into it would be detached by the growth: their length would
  // be 0 and reading them would throw.
  const memory = memoryOfWasm();
  const big = openVcf(manyVariantsVcf(300000));
  for (const block of big.iterBlocks()) {
    assert.ok(block.numVars > 0);
  }
  big.free();
  assert.ok(
    memoryOfWasm() > memory,
    `the memory of wasm did not grow: ${memoryOfWasm()} bytes`,
  );
  assert.deepEqual([...bytes], kept);
});

test("the individuals and the ploidy are known when the vars file is opened", async () => {
  const variants = openVars(await varsFileOfCases(3));
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  assert.equal(variants.numIndividuals, 3);
  assert.equal(variants.ploidy, 2);
  variants.free();
});

test("bytes that are not a vars file are refused when they are opened", async () => {
  // A VCF, which is the file a user has beside the vars file.
  const vcf = await referenceVcf("cases.vcf");
  assert.throws(() => openVars(vcf), {
    name: "Error",
    message: /not a vars file/,
  });
});

test("a vars file compressed with zstd opens and throws at its first block", async () => {
  // popnei cannot write it: no build of it carries the zstd crate, and
  // pyarrow wrote `tests/reference/vars/zstd.vars` for this test. Arrow
  // decompresses a batch when it reads it, so the file opens and the error
  // comes with the first block.
  const variants = openVars(await referenceVars("zstd.vars"));
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  assert.throws(() => [...variants.iterBlocks()], {
    name: "Error",
    message: /zstd/,
  });
  variants.free();
});

test("every pass over a vars file reads the same bytes again", async () => {
  const variants = openVars(await varsFileOfCases(3));
  const first = rowsOf(variants);
  const second = rowsOf(variants, 2);
  assert.deepEqual(second, first);
  assert.deepEqual(first, CASES);
  variants.free();
});

test("a vars file is written again from the variants of one", async () => {
  const read = openVars(await varsFileOfCases(3));
  const written = openVars(writeVars(read, { numVarsPerBlock: 2 }).bytes);
  const blocks = [
    ...written.iterBlocks({ fields: ALL_FIELDS, numVarsPerBlock: 2 }),
  ];
  assert.deepEqual(
    blocks.map((block) => block.numVars),
    [2, 2],
  );
  assert.deepEqual(
    blocks.flatMap((block) => rowsOfTheBlock(block, 6)),
    CASES,
  );
  assert.deepEqual(written.individuals, read.individuals);
  read.free();
  written.free();
});

test("a tetraploid VCF gives a vars file of tetraploid genotypes", () => {
  const tetraploid = vcfOf([
    "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1/1/1\t0/0/0/0",
    "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t1/1/1/1\t0/0/0/1\t./././.",
  ]);
  const vcf = openVcf(tetraploid, { ploidy: 4 });
  const bytes = writeVars(vcf, { numVarsPerBlock: 2 }).bytes;
  vcf.free();
  const variants = openVars(bytes);
  // The ploidy comes from the `popnei` key of the file, and the width of
  // its `gts` column is the individuals times that ploidy.
  assert.equal(variants.ploidy, 4);
  assert.equal(variants.numIndividuals, 3);
  const [block] = [...variants.iterBlocks()];
  assert.ok(block !== undefined);
  assert.equal(block.ploidy, 4);
  assert.deepEqual(
    [...block.gts],
    [0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, -1, -1, -1, -1],
  );
  variants.free();
});

test("a source that fails half way is thrown and leaves the handle usable", () => {
  const before = numberOfOpenPasses();
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\tten\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    ]),
  );
  // The error a user reads is the one of the VCF that was being read, and
  // no file of the variants that were written before it comes back.
  assert.throws(() => writeVars(variants, { numVarsPerBlock: 1 }), {
    name: "Error",
    message: /line 4 of the VCF, the column POS/,
  });
  assert.equal(numberOfOpenPasses(), before);
  // The handle is the one it was: the bytes of the VCF are still there and
  // the variants before the wrong line are still read.
  assert.deepEqual(variants.individuals, ["ind1", "ind2", "ind3"]);
  let numBlocks = 0;
  for (const block of variants.iterBlocks({ numVarsPerBlock: 1 })) {
    assert.deepEqual([...(block.pos ?? [])], [10]);
    numBlocks += 1;
    // The wrong line is the next block, and this pass is left before it.
    break;
  }
  assert.equal(numBlocks, 1);
  variants.free();
});

test("a pass over a vars file gives itself back however it ends", async () => {
  const before = numberOfOpenPasses();
  const variants = openVars(await varsFileOfCases(1));
  assert.equal([...variants.iterBlocks({ numVarsPerBlock: 1 })].length, 4);
  assert.equal(numberOfOpenPasses(), before);
  for (const block of variants.iterBlocks({ numVarsPerBlock: 1 })) {
    assert.equal(block.numVars, 1);
    assert.equal(numberOfOpenPasses(), before + 1);
    break;
  }
  assert.equal(numberOfOpenPasses(), before);
  variants.free();
  const zstd = openVars(await referenceVars("zstd.vars"));
  assert.throws(() => [...zstd.iterBlocks()], { name: "Error" });
  assert.equal(numberOfOpenPasses(), before);
  zstd.free();
});

test("a source of openVars that is not a Uint8Array is refused", async () => {
  const text = new TextDecoder().decode(await referenceVcf("cases.vcf"));
  assert.throws(() => openVars(text as unknown as Uint8Array), {
    name: "Error",
    message: /Uint8Array/,
  });
  assert.throws(() => openVars(null as unknown as Uint8Array), {
    name: "Error",
    message: /and null was given/,
  });
  // `undefined` said "the undefined undefined", and an array "a Array".
  assert.throws(() => openVars(undefined as unknown as Uint8Array), {
    name: "Error",
    message: /and undefined was given/,
  });
  assert.throws(() => openVars([65] as unknown as Uint8Array), {
    name: "Error",
    message: /an object of the type `Array` was given/,
  });
});

test("bytes that do not fit in the memory of wasm are refused", () => {
  // What the package asks before it lets the generated code copy an array
  // into the memory of wasm, which allocates the whole length first and
  // traps when that fails, leaving the module unusable. It is called here
  // with a length alone: an array of 4 GB in node, to reach it through
  // `openVars`, is 4 GB of the machine this test runs on.
  //
  // A wasm module addresses 4 GB and some megabytes of them are already
  // open here, so neither of these fits. The first is above what a whole
  // number of wasm counts, 2^32 - 1, and the second is under it.
  for (const numBytes of [4294967296, 4294000000]) {
    assert.throws(() => roomForBytes(numBytes), {
      name: "Error",
      message: /do not fit in the memory popnei has left/,
    });
  }
  // What fits is not an error, and what it grew stays for the copy.
  roomForBytes(1024);
});

test("a source whose buffer was transferred away is refused", async () => {
  const bytes = await varsFileOfCases(3);
  // What a page does when it sends the bytes of a file to a web worker:
  // the buffer moves and the array that is left has nothing behind it. The
  // generated code read it as bytes of its own and threw a `TypeError`
  // that named neither the argument nor what had happened.
  structuredClone(bytes.buffer, { transfer: [bytes.buffer] });
  // `detached` is of ES2024, which is later than the library this package
  // is compiled against.
  assert.equal((bytes.buffer as { detached?: boolean }).detached, true);
  assert.throws(() => openVars(bytes), {
    name: "Error",
    message: /the buffer of `source` was transferred/,
  });
});

test("variants that writeVars was not given by popnei are refused", async () => {
  const bytes = await varsFileOfCases(3);
  // The bytes of the file where the handle over them goes is the mistake
  // that is easiest to make.
  assert.throws(() => writeVars(bytes as unknown as Variants), {
    name: "Error",
    message:
      /`variants` is what openVcf or openVars gives, and an object of the type `Uint8Array`/,
  });
  assert.throws(() => writeVars(null as unknown as Variants), {
    name: "Error",
    message: /null/,
  });
});

test("variants that were freed cannot be written", async () => {
  const variants = openVars(await varsFileOfCases(3));
  variants.free();
  assert.throws(() => writeVars(variants), {
    name: "Error",
    message: /freed/,
  });
});

test("a numVarsPerBlock that writeVars cannot use is refused", async () => {
  const variants = openVcf(await referenceVcf("cases.vcf"));
  // 2^32 + 1 reached the core as batches of one variant, and NaN as none.
  for (const numVarsPerBlock of [4294967297, -1, 0, 1.5, Number.NaN]) {
    assert.throws(() => writeVars(variants, { numVarsPerBlock }), {
      name: "Error",
      message: /`numVarsPerBlock` is a whole number/,
    });
  }
  variants.free();
});
