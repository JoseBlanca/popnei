/**
 * What writing a file with `onBytes` costs in the memory of wasm, against
 * writing it whole, as "The file of a writer in pieces" of
 * `docs/specs/js_sources.md` has it.
 *
 * What a write in pieces holds is what its pass holds for a block, the
 * reader's block and what the writer makes of it, and not the file, so the
 * file is made large enough that a block is well under half of it: 30000
 * variants of 600 individuals, 72.9 MB of text.
 *
 * The tests are alone in their file, and in this order, because the memory
 * of wasm never shrinks: what a test frees stays there as room the next one
 * fits into. The two writes in pieces come first, and the write of the
 * whole file last, which fits in what they left only for the room of a
 * block. node's test runner gives each test file its own process.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVcf, writeVars, writeVcf } from "popnei";

import loadTheWasm from "../wasm/popnei.js";
import { vcfOfDrawnGenotypes } from "./reference.ts";

await init();

/** The WebAssembly of the core, whose memory this file watches. */
const wasm = await loadTheWasm();

/** How many bytes the memory of wasm holds. */
function memoryOfWasm(): number {
  return wasm.memory.buffer.byteLength;
}

const variants = openVcf(vcfOfDrawnGenotypes(30000, 600), { onlyPassed: false });

/** How many bytes the plain VCF holds, which the write in pieces counts. */
let numBytesOfTheVcf = 0;

test("a vars file written in pieces grows the memory of wasm by less than half the file", () => {
  let numBytes = 0;
  const before = memoryOfWasm();
  writeVars(variants, {
    numVarsPerBlock: 100,
    onBytes: (piece) => {
      numBytes += piece.length;
    },
  });
  const grew = memoryOfWasm() - before;
  assert.ok(numBytes > 20_000_000, `the file is ${numBytes} bytes`);
  assert.ok(
    grew < numBytes / 2,
    `writing a vars file of ${numBytes} bytes in pieces grew the memory of wasm by ${grew} bytes`,
  );
});

test("a plain VCF written in pieces grows the memory of wasm by less than half the file", () => {
  const before = memoryOfWasm();
  writeVcf(variants, {
    bgzip: false,
    onBytes: (piece) => {
      numBytesOfTheVcf += piece.length;
    },
  });
  const grew = memoryOfWasm() - before;
  assert.ok(numBytesOfTheVcf > 70_000_000, `the file is ${numBytesOfTheVcf} bytes`);
  assert.ok(
    grew < numBytesOfTheVcf / 2,
    `writing a VCF of ${numBytesOfTheVcf} bytes in pieces grew the memory of wasm by ${grew} bytes`,
  );
});

test("the same VCF written whole grows the memory of wasm by more than half the file", () => {
  const before = memoryOfWasm();
  const written = writeVcf(variants, { bgzip: false });
  const grew = memoryOfWasm() - before;
  variants.free();
  assert.equal(written.bytes.length, numBytesOfTheVcf);
  assert.ok(
    grew > numBytesOfTheVcf / 2,
    `writing a VCF of ${numBytesOfTheVcf} bytes whole grew the memory of wasm by ${grew} bytes`,
  );
});
