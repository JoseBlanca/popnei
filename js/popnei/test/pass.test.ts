/**
 * That a pass over a VCF ends at its own error, which the class the binding
 * crate exports is what shows.
 *
 * `docs/specs/block.md` asks of every reader that it give no block after an
 * error: a reader that went on would hand out the variants that follow the
 * wrong one as if nothing had happened. The errors of the pass itself, a
 * block that is not of its own size and a position that a number of
 * JavaScript does not hold, happen after the block was taken from the
 * reader, so the reader knows nothing of them and only the pass can keep
 * that promise.
 *
 * The package closes its generator when a block throws, so what a user of
 * the package sees is already the end of the iteration. This test goes to
 * `Blocks`, the class of wasm that the package is written on and that
 * anything importing the generated JavaScript reaches, and asks it for a
 * block after its error.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init } from "popnei";

import { Steps, open_vcf } from "../wasm/popnei.js";
import { vcfOf } from "./reference.ts";

await init();

/**
 * The first position a number of JavaScript does not hold, 2^53 + 1, which
 * the binding crate refuses once per block: rounding it would give a
 * TypeScript user another position than a Python user.
 */
const POSITION_ABOVE_THE_LARGEST = "9007199254740993";

/** Six variants, the third of them at a position JavaScript cannot hold. */
const SIX_VARIANTS = vcfOf([
  "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  "chr1\t200\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  `chr1\t${POSITION_ABOVE_THE_LARGEST}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1`,
  "chr1\t400\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  "chr1\t500\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  "chr1\t600\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
]);

test("a pass gives no block after an error of its own", () => {
  const source = open_vcf(SIX_VARIANTS, 2, false);
  // A pass runs the steps it is given, and a pass with none is asked for
  // with an empty list of them: the argument cannot be left out. The steps
  // are built with the individuals of the source, which is what a filter of
  // individuals resolves its names against, and this one has no step.
  const pass = source.blocks(["chrom", "pos"], 2, new Steps([]));
  try {
    const first = pass.next_block();
    assert.ok(first !== undefined);
    assert.deepEqual(Array.from(first.pos() ?? []), [100, 200]);
    first.free();
    // The block of the third and the fourth variants holds the position
    // that JavaScript cannot hold, so it is lost with the error.
    assert.throws(() => pass.next_block(), {
      message: new RegExp(POSITION_ABOVE_THE_LARGEST),
    });
    // The variant at 400 was in that block, and the two that follow it
    // would come out of a pass that went on as if nothing had happened.
    assert.equal(pass.next_block(), undefined);
    assert.equal(pass.next_block(), undefined);
  } finally {
    pass.free();
    source.free();
  }
});
