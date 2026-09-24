/**
 * What an application gets back when the function it is told the progress
 * with throws: the value it threw, and not an error of popnei.
 *
 * `onProgress` sets a function that a pass calls while it reads, and a
 * function that throws ends the pass where it was reading, as "What the
 * source tells the page" of `docs/specs/js_sources.md` says. The read fails,
 * the error travels out through the readers of the core, and the consumer
 * throws the value the function threw, which the application recognises with
 * `===` and without reading a message. A pass is one reading of a source
 * from its start, a run is one call of one consumer with the passes it
 * makes, and a consumer is a function of the package that reads a
 * `Variants`, the iteration of `iterBlocks` among them, as
 * `docs/glossary.md` has the three words.
 *
 * Which error the core gives for the failed read depends on the reader that
 * was reading: a VCF in the middle of a line and a gzipped VCF at the end of
 * its last member fail in different words, and the file that was stopped at
 * its end would otherwise be the error of a file that was cut short. So the
 * tests here stop a pass at three places, and every one of them gives the
 * value the function threw.
 *
 * The files are `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
 * individuals in 117346 bytes, its bgzipped `many.vcf.gz` of 21904, and a
 * VCF of 150000 variants of 3 individuals that the test writes, 5.2 MB,
 * which is more than the 4 MiB a pass reads between two calls and is
 * therefore told of three times.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Progress, Variants } from "popnei";
import {
  calcPerIndividualStats,
  calcPerVarDistribs,
  doPcaFromVariants,
  init,
  openVcf,
} from "popnei";

import { manyVariantsVcf, referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The bytes of `many.vcf.gz`, which bgzip wrote as four gzip members. */
const MANY_VCF_GZ = await referenceVcf("many.vcf.gz");

/**
 * How many variants `many.vcf` holds, which the tests here read with
 * `onlyPassed` false.
 *
 * 25 of the 500 failed a filter of the file, as `docs/specs/io_vcf.md` says,
 * so a pass that leaves those out gives 475: what a test of the stop asserts
 * is the whole file, read again after the run that was stopped.
 */
const VARIANTS_OF_MANY_VCF = 500;

/** `many.vcf` read with every variant of it, the 25 that failed a filter
 * among them. */
const EVERY_VARIANT = { onlyPassed: false } as const;

/**
 * What the function of these tests throws, which is not an `Error`: an
 * application cancels a run with a value of its own and tells it from a file
 * that could not be read by what it is, and not by a message.
 */
const THE_CANCEL = { whyThePassEnded: "the user pressed the button" };

/**
 * What `run` threw, and an assertion that fails when it threw nothing.
 *
 * `assert.throws` is not what these tests use: what is thrown here is not an
 * `Error`, and what they assert of it is that it is one value and not
 * another one that reads the same.
 */
function whatWasThrownBy(run: () => void): unknown {
  try {
    run();
  } catch (thrown: unknown) {
    return thrown;
  }
  return assert.fail("the run gave its result instead of being stopped");
}

/**
 * Sets a function on `variants` that throws [`THE_CANCEL`] at the call
 * `stopAt`, the first being 0, and does nothing at the others, and gives
 * back the calls it was made with.
 *
 * It throws once and not at every call after that: a `Variants` whose run
 * was stopped is the one it was, and the tests read its variants afterwards
 * with the same function still set on it.
 */
function stoppedAt(variants: Variants, stopAt: number): Progress[] {
  const calls: Progress[] = [];
  variants.onProgress((progress) => {
    calls.push(progress);
    if (calls.length === stopAt + 1) {
      throw THE_CANCEL;
    }
  });
  return calls;
}

test("a function that throws at the first read stops the run with its value", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  const calls = stoppedAt(variants, 0);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerVarDistribs(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.ok(
      !(thrown instanceof Error),
      "the value of the application arrived as an error of popnei",
    );
    // The first read of the pass is the header of the VCF, read when the
    // reader is built, so the pass was stopped before it gave a variant.
    assert.equal(calls.length, 1);
    assert.equal(calls.at(0)?.bytesRead, 0);
    // Nothing of popnei is left in a state a later call notices: the same
    // `Variants` reads the file again from its start, with the function that
    // threw still set on it.
    let numVars = 0;
    for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
      numVars += block.numVars;
    }
    assert.equal(numVars, VARIANTS_OF_MANY_VCF);
  } finally {
    variants.free();
  }
});

test("a function that throws in the middle of a VCF stops the run with its value", () => {
  // 150000 variants of 3 individuals, 5.2 MB, which a pass is told of at its
  // first read, at the read that finds a range of 4 MiB read since then, and
  // at the read that finds no more bytes. The second of the three is inside
  // the file, where the reader is in the middle of its lines.
  const variants = openVcf(manyVariantsVcf(150000));
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerIndividualStats(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 2);
    assert.ok(
      (calls.at(1)?.bytesRead ?? 0) >= 4 * 1024 * 1024,
      `the call that threw says ${calls.at(1)?.bytesRead} bytes read`,
    );
  } finally {
    variants.free();
  }
});

test("a function that throws at the end of a gzipped VCF is not a file cut short", () => {
  const variants = openVcf(MANY_VCF_GZ);
  // The last call of a pass over a VCF is the read that finds no more bytes,
  // which for a bgzipped file is the decoder looking for another member
  // after the last one. What that read fails with in the core is the error
  // of a file that ends where a member should start, and what the
  // application is given is its own value all the same.
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerVarDistribs(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 2);
    assert.equal(calls.at(1)?.bytesRead, calls.at(1)?.numBytes);
  } finally {
    variants.free();
  }
});

test("a function that throws at the first read of the second pass stops the pca", () => {
  const variants = openVcf(MANY_VCF);
  // The principal components of the variants build both of their readers
  // before they ask either for a block, and a reader reads when it is built,
  // so the calls of a run are pass 1 at 0 bytes, pass 2 at 0, and then the
  // two reads that find the end of the file. The second call is the first
  // read of the second pass.
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      doPcaFromVariants(variants, {
        numPrinComps: 10,
        transformToBiallelic: true,
      });
    });
    assert.equal(thrown, THE_CANCEL);
    assert.deepEqual(
      calls.map((call) => ({ pass: call.pass, bytesRead: call.bytesRead })),
      [
        { pass: 1, bytesRead: 0 },
        { pass: 2, bytesRead: 0 },
      ],
    );
  } finally {
    variants.free();
  }
});

test("an iteration of blocks that is stopped throws the value of the application", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  // The first read of the pass is the header, and the read that finds the
  // end of the file is the second call: this one is thrown at the end of the
  // iteration, when every block has been given.
  const calls = stoppedAt(variants, 1);
  try {
    let numVars = 0;
    const thrown = whatWasThrownBy(() => {
      for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
        numVars += block.numVars;
      }
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 2);
    assert.equal(numVars, VARIANTS_OF_MANY_VCF);
  } finally {
    variants.free();
  }
});

test("a function that frees the variants while an iteration reads them does not trap", () => {
  // The VCF of 150000 variants, whose pass is told at its first read, in the
  // middle of the file and at its end: the source is freed at the call in
  // the middle, so the iteration has 1 MB of the file left to read after it.
  const variants = openVcf(manyVariantsVcf(150000));
  let calls = 0;
  let freedAt = 0;
  variants.onProgress(() => {
    calls += 1;
    // What the source keeps in JavaScript, the function itself, stays until
    // the pass that is reading is done with it, so the call at the end of
    // the file is made through it after the source was freed.
    if (calls === 2) {
      freedAt = calls;
      variants.free();
    }
  });
  let numVars = 0;
  for (const block of variants.iterBlocks({ numVarsPerBlock: 1000 })) {
    numVars += block.numVars;
  }
  assert.equal(freedAt, 2);
  assert.equal(calls, 3);
  assert.equal(numVars, 150000);
  // A second `free` is not an error: it has nothing left to give back.
  variants.free();
});

test("a function that runs a consumer over the same variants does not trap", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  const insideTheFunction: number[] = [];
  let alreadyRan = false;
  variants.onProgress(() => {
    // The function is called with no table of the binding crate borrowed, so
    // a consumer started from inside it runs as any other call does. It is
    // started once: the pass of that consumer tells the page through the
    // same function, and a function that started one at every call would
    // start them without end.
    if (alreadyRan) {
      return;
    }
    alreadyRan = true;
    insideTheFunction.push(calcPerIndividualStats(variants).passStats.numVars);
  });
  try {
    const distribs = calcPerVarDistribs(variants);
    // The run that was reading when the other one started gives the variants
    // of the whole file, and so does the one that ran inside it.
    assert.equal(distribs.passStats.numVars, VARIANTS_OF_MANY_VCF);
    assert.deepEqual(insideTheFunction, [VARIANTS_OF_MANY_VCF]);
  } finally {
    variants.free();
  }
});

test("free from inside the function of a run that is reading is refused", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  // While a consumer runs, wasm-bindgen holds the source for the length of
  // that call, and the free of a value it holds throws inside wasm. What
  // the generated `free` does before that call is to zero the pointer of
  // the handle and take it out of the `FinalizationRegistry`, so a free that
  // was let through left the source in the memory of wasm with nothing left
  // to free it and the next call of these variants reading a null pointer.
  const refused: unknown[] = [];
  variants.onProgress(() => {
    refused.push(
      whatWasThrownBy(() => {
        variants.free();
      }),
    );
  });
  try {
    const stats = calcPerIndividualStats(variants);
    assert.equal(stats.passStats.numVars, VARIANTS_OF_MANY_VCF);
    assert.ok(refused.length > 0, "the run told the page nothing");
    for (const thrown of refused) {
      assert.ok(thrown instanceof Error, `the free threw ${String(thrown)}`);
      assert.match(thrown.message, /a run is reading these variants/);
    }
    // The handle is the one it was: the source was not freed, and the
    // variants read the file again from its start. The function is taken
    // off first, because between two blocks of an iteration no call holds
    // the source and a free from inside it would be taken.
    variants.onProgress();
    let numVars = 0;
    for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
      numVars += block.numVars;
    }
    assert.equal(numVars, VARIANTS_OF_MANY_VCF);
  } finally {
    variants.free();
  }
  // The free after the run goes through, and what it freed is gone: the
  // message is popnei's for variants that were freed, and not the
  // `null pointer passed to rust` of a handle that was half freed.
  assert.throws(() => variants.iterBlocks(), {
    name: "Error",
    message: /these variants were freed/,
  });
});

test("the error of a free that was refused stops the run that was reading", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  // The function does not catch what `free` threw, so it leaves the
  // function as any other value an application throws does: the pass ends
  // there and the consumer gives that error back.
  variants.onProgress(() => {
    variants.free();
  });
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerVarDistribs(variants);
    });
    assert.ok(thrown instanceof Error, `the run threw ${String(thrown)}`);
    assert.match(thrown.message, /a run is reading these variants/);
  } finally {
    variants.free();
  }
});
