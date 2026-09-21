/**
 * That every pass over a vars file reads the bytes `openVars` was given,
 * and none copies them.
 *
 * It is alone in its file because the memory of wasm never shrinks: what
 * the tests of another file freed would stay there as room a copy of the
 * file fits into, and the measurement would then say nothing. node's test
 * runner gives each test file its own process, as `before_init.test.ts`
 * uses too.
 *
 * The measurement is twelve passes open at once. One pass alone does not
 * tell the two apart: its copy fits in what writing the file left free,
 * and the memory of wasm does not grow at all.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVars, openVcf, writeVars } from "popnei";

import loadTheWasm from "../wasm/popnei.js";

await init();

/** The WebAssembly of the core, whose memory this file watches. */
const wasm = await loadTheWasm();

/** How many bytes the memory of wasm holds. */
function memoryOfWasm(): number {
  return wasm.memory.buffer.byteLength;
}

/** How many passes are open at once while the memory is read. */
const NUM_PASSES = 12;

/**
 * The bytes of a VCF of `numIndividuals` diploid individuals and `numVars`
 * variants, whose vars file is some megabytes.
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

test("the passes of a vars file share the bytes it was opened with", () => {
  const vcf = openVcf(vcfOfDrawnGenotypes(14000, 600), { onlyPassed: false });
  const bytes = writeVars(vcf, { numVarsPerBlock: 100 });
  vcf.free();
  const variants = openVars(bytes);
  const before = memoryOfWasm();
  const passes = Array.from({ length: NUM_PASSES }, () => {
    const pass = variants.iterBlocks({ numVarsPerBlock: 100 });
    assert.equal(pass.next().value?.numVars, 100);
    return pass;
  });
  const grew = memoryOfWasm() - before;
  for (const pass of passes) {
    pass.return?.();
  }
  variants.free();
  // Each pass holds its reader and the batch it decoded, 100 variants of
  // 600 individuals, and the twelve together grow the memory of wasm by
  // nothing at all. A reader that copied the bytes of the file, which is
  // what the reader of a VCF did before the passes shared them, holds
  // twelve copies of it: measured at 36765696 bytes of growth for the
  // 12231602 byte file this writes, three times what is asked here.
  assert.ok(
    grew < bytes.length,
    `${NUM_PASSES} passes over a vars file of ${bytes.length} bytes grew ` +
      `the memory of wasm by ${grew} bytes`,
  );
});
