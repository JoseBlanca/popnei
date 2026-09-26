/**
 * What a BED of some megabytes costs in the memory of wasm, which never
 * shrinks: what `filterByRegions` holds while it reads one is memory the
 * page keeps.
 *
 * The bytes of the BED are copied into that memory by the code wasm-bindgen
 * generates, and the core reads them there, with no copy of its own when
 * the BED is not gzipped, into 16 bytes for each region. A BED of 20 MB of
 * regions that do not overlap grew the memory by 74.5 MB when the core
 * copied the bytes and joined the regions into a second vector, on 26
 * September 2026, and by 54.5 MB when it read them where they were and
 * joined them in place; 89.8 MB grew it by 315 MB and 225 MB.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVcf } from "popnei";

import loadTheWasm from "../wasm/popnei.js";
import { referenceVcf } from "./reference.ts";

await init();

const wasm = await loadTheWasm();
const MANY_VCF = await referenceVcf("many.vcf");

/** A BED of about `bytes` bytes, of regions of 5 bases 10 apart on three chromosomes. */
function aBedOf(bytes: number): Uint8Array {
  const lines: string[] = [];
  let length = 0;
  for (let line = 0; length < bytes; line += 1) {
    const text = `chr${1 + (line % 3)}\t${line * 10}\t${line * 10 + 5}\n`;
    lines.push(text);
    length += text.length;
  }
  return new TextEncoder().encode(lines.join(""));
}

test("a BED of 20 MB grows the memory of wasm by less than three times its size", () => {
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  const bed = aBedOf(20_000_000);
  const before = wasm.memory.buffer.byteLength;
  variants.filterByRegions(bed);
  const grew = wasm.memory.buffer.byteLength - before;
  assert.ok(
    grew < 3 * bed.length,
    `a BED of ${bed.length} bytes grew the memory of wasm by ${grew} bytes`,
  );
  assert.equal(variants.steps[0]?.args["numRegions"], 962_963);
});
