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

import loadTheWasm from "../wasm/popnei.js";
import { numberOfOpenPasses } from "../dist/variant.js";
import { referenceVars, referenceVcf } from "./reference.ts";

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
 * The bytes of a VCF of `numIndividuals` diploid individuals and `numVars`
 * variants, for the test that watches the memory of wasm, which needs a
 * file of some megabytes.
 *
 * The genotypes are drawn with a generator of its own, so that the file
 * does not compress to nothing: a column of one repeated genotype would be
 * a few kilobytes of lz4 whatever its number of variants.
 */
function vcfOfDrawnGenotypes(
  numVars: number,
  numIndividuals: number,
): Uint8Array {
  const names = Array.from(
    { length: numIndividuals },
    (_unused, individual) => `ind${individual + 1}`,
  );
  const lines = [
    "##fileformat=VCFv4.4",
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
  ];
  // A linear congruential generator, the one of numerical recipes: the
  // numbers are the same at every run, so the size of the file is too.
  let drawn = 1;
  const nextAllele = (): number => {
    drawn = (Math.imul(drawn, 1664525) + 1013904223) >>> 0;
    return drawn >>> 30;
  };
  for (let variant = 0; variant < numVars; variant += 1) {
    const genotypes = Array.from({ length: numIndividuals }, () => {
      const first = nextAllele();
      const second = nextAllele();
      return `${first > 2 ? "." : first}/${second > 2 ? "." : second}`;
    });
    lines.push(
      `chr1\t${variant + 1}\t.\tA\tC,G\t.\tPASS\t.\tGT\t${genotypes.join("\t")}`,
    );
  }
  return new TextEncoder().encode([...lines, ""].join("\n"));
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

test("writing a vars file holds no pass when it returns", async () => {
  const before = numberOfOpenPasses();
  const variants = openVcf(await referenceVcf("cases.vcf"));
  // The pass over the source is the core's, inside the one call, so nothing
  // of it is left in the memory of wasm afterwards.
  assert.ok(writeVars(variants).bytes.length > 0);
  assert.equal(numberOfOpenPasses(), before);
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

test("the passes of a vars file share the bytes it was opened with", () => {
  const vcf = openVcf(vcfOfDrawnGenotypes(4000, 300), { onlyPassed: false });
  const bytes = writeVars(vcf, { numVarsPerBlock: 100 }).bytes;
  vcf.free();
  const variants = openVars(bytes);
  const before = memoryOfWasm();
  // Twelve passes at once, each with its reader and the batch it decoded in
  // the memory of wasm. A reader that copied the bytes of the file, which
  // is what the reader of a VCF did before the passes shared them, would
  // hold twelve copies of the file: that was measured at 7.0 MB of growth
  // for the 1.75 MB file this writes, where the readers that share the bytes
  // grow the memory by nothing at all. One pass alone shows neither: a copy
  // of 1.8 MB fits in what the writing of the file left free.
  const passes = Array.from({ length: 12 }, () => {
    const pass = variants.iterBlocks({ numVarsPerBlock: 100 });
    assert.equal(pass.next().value?.numVars, 100);
    return pass;
  });
  const grew = memoryOfWasm() - before;
  for (const pass of passes) {
    pass.return?.();
  }
  assert.ok(
    grew < bytes.length,
    `twelve passes over a vars file of ${bytes.length} bytes grew the memory ` +
      `of wasm by ${grew} bytes`,
  );
  variants.free();
});

test("a source of openVars that is not a Uint8Array is refused", async () => {
  const text = new TextDecoder().decode(await referenceVcf("cases.vcf"));
  assert.throws(() => openVars(text as unknown as Uint8Array), {
    name: "Error",
    message: /Uint8Array/,
  });
  assert.throws(() => openVars(null as unknown as Uint8Array), {
    name: "Error",
    message: /null/,
  });
});

test("variants that writeVars was not given by popnei are refused", async () => {
  const bytes = await varsFileOfCases(3);
  // The bytes of the file where the handle over them goes is the mistake
  // that is easiest to make.
  assert.throws(() => writeVars(bytes as unknown as Variants), {
    name: "Error",
    message: /`variants` is what openVcf or openVars gives, and a Uint8Array/,
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
