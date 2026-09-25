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
 * Each of the twelve consumers is stopped in two ways, in the two loops
 * over `THE_CONSUMERS` of `consumers.ts`: with a value of the application,
 * which the consumer has to give back as it is, and with the `free()` of the
 * variants that are being read, which popnei refuses and whose refusal then
 * stops the pass as any other thrown value does. Both need every consumer to
 * hold its source while its run reads: a consumer that hands the source to
 * the core outside that hold lets a `free()` from the function through, and
 * what the generated `free` does before the call it fails in is to zero the
 * pointer of the handle and take it out of the `FinalizationRegistry`, so
 * the source stays in the memory of wasm with no handle left to free it.
 *
 * Which error the core gives for the failed read depends on the reader that
 * was reading: a plain VCF in the middle of a line, a gzipped one inside its
 * decompressor and a vars file inside the read of a batch, which the reader
 * of that file turns into the error of a file that was cut short. So the
 * tests here stop a pass over each of the three, and every one of them gives
 * the value the function threw.
 *
 * The calls that end a run are the ones that cannot stop anything: they are
 * made when the run is over, and a value thrown in one of them is dropped,
 * as "What the source tells the page" of `docs/specs/js_sources.md` says.
 *
 * The files are `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
 * individuals in 117346 bytes, its bgzipped `many.vcf.gz` of 21904, and a
 * VCF of 150000 variants of 3 individuals that the test writes, 6188977
 * bytes, which is more than the 4 MiB a pass reads between two calls and is
 * therefore told of three times.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Progress, Variants } from "popnei";
import {
  calcPerIndividualStats,
  calcPerVarDistribs,
  calcPopDiversity,
  doPcaFromVariants,
  init,
  openVars,
  openVcf,
  writeVars,
} from "popnei";

import { THE_CONSUMERS, theConsumersTheCrateNames } from "./consumers.ts";
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

/**
 * How many variants of `many.vcf` passed the filter the file itself records,
 * which is what a pass over it gives with the `onlyPassed` the two loops
 * over the consumers open it with.
 */
const VARIANTS_THAT_PASSED = 475;

/** How many bytes `many.vcf` holds, which `docs/specs/io_vcf.md` gives. */
const BYTES_OF_MANY_VCF = 117346;

/** How many bytes `many.vcf.gz` holds. */
const BYTES_OF_MANY_VCF_GZ = 21904;

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

/**
 * How many variants `variants` gives when they are read again, with no
 * function set on them.
 *
 * It is what says that the `Variants` is the one it was: a source that was
 * left half freed answers a read with the error of variants that were freed,
 * and one whose handle was broken with a null pointer passed to Rust.
 */
function numVarsReadAgain(variants: Variants): number {
  variants.onProgress();
  let numVars = 0;
  for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
    numVars += block.numVars;
  }
  return numVars;
}

test("the twelve consumers the tests stop are the twelve the crate names", () => {
  // The loops below are worth what their list holds: a consumer left out of
  // it is never stopped and nothing says so. The crate's own list is what
  // the message of a name that is of no consumer gives, and it is written in
  // no file of TypeScript.
  assert.deepEqual(
    [...THE_CONSUMERS.map((consumer) => consumer.name)].sort(),
    [...theConsumersTheCrateNames()].sort(),
  );
});

for (const consumer of THE_CONSUMERS) {
  test(`${consumer.name} stopped at its first read throws the value of the application`, () => {
    const variants = openVcf(MANY_VCF);
    const calls = stoppedAt(variants, 0);
    try {
      const thrown = whatWasThrownBy(() => {
        consumer.run(variants);
      });
      assert.equal(thrown, THE_CANCEL);
      assert.ok(
        !(thrown instanceof Error),
        `${consumer.name} gave an error of popnei: ${String(thrown)}`,
      );
      // The first read of the pass is the header of the VCF, read when the
      // reader is built, so every one of the twelve was stopped before it
      // was given a variant.
      assert.equal(calls.length, 1);
      assert.equal(calls.at(0)?.bytesRead, 0);
      assert.equal(numVarsReadAgain(variants), VARIANTS_THAT_PASSED);
    } finally {
      variants.free();
    }
  });
}

for (const consumer of THE_CONSUMERS) {
  test(`free from inside the function while ${consumer.name} reads is refused`, () => {
    const variants = openVcf(MANY_VCF);
    // While a consumer runs, wasm-bindgen holds the source for the length of
    // that call, and the free of a value it holds throws inside wasm after
    // the generated code has zeroed the pointer of the handle and taken it
    // out of the `FinalizationRegistry`. So the package counts the run and
    // refuses the free itself, with an `Error` of its own, and that error
    // leaves the function as any other thrown value does and stops the pass.
    // A consumer that hands the source to the core without counting its run
    // gives wasm-bindgen's own sentence here and leaves a `Variants` that
    // says it was freed over a source nothing can free.
    let calls = 0;
    variants.onProgress(() => {
      calls += 1;
      if (calls > 1) {
        return;
      }
      variants.free();
    });
    try {
      const thrown = whatWasThrownBy(() => {
        consumer.run(variants);
      });
      assert.ok(
        thrown instanceof Error,
        `${consumer.name} threw ${String(thrown)}`,
      );
      assert.match(thrown.message, /a run is reading these variants/);
      assert.equal(numVarsReadAgain(variants), VARIANTS_THAT_PASSED);
    } finally {
      variants.free();
    }
  });
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
  // 150000 variants of 3 individuals, 6188977 bytes, which a pass is told of
  // at its first read, at the read that finds a range of 4 MiB read since
  // then, and at the end of the run. The second of the three is inside the
  // file, where the reader is in the middle of its lines.
  const variants = openVcf(manyVariantsVcf(150000));
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerIndividualStats(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 2);
    const stoppedAtBytes = calls.at(1);
    assert.ok(
      (stoppedAtBytes?.bytesRead ?? 0) >= 4 * 1024 * 1024,
      `the call that threw says ${stoppedAtBytes?.bytesRead} bytes read`,
    );
    // The call it stopped at is inside the file and not the one at its end:
    // a pass stopped where no read follows would end the same way whether
    // the stop worked or not. It says 4194313 bytes read of 6188977, with
    // 1994664 left.
    assert.ok(
      (stoppedAtBytes?.bytesRead ?? 0) < (stoppedAtBytes?.numBytes ?? 0),
      `the call that threw says the whole file, ${JSON.stringify(stoppedAtBytes)}`,
    );
  } finally {
    variants.free();
  }
});

test("a function that throws in the middle of a VCF stops the diversity with its value", () => {
  // The same file and the same call of the two above, with the diversity of
  // the populations as the consumer: a consumer that reads the source
  // without opening a run of its own is told the progress of its first read
  // and of nothing after it, and gives popnei's error for the failed read
  // where the application is owed the value it threw. The populations are
  // left out, which is one population of the three individuals of the file.
  const variants = openVcf(manyVariantsVcf(150000));
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPopDiversity(variants, { minNumIndividuals: 1 });
    });
    assert.equal(thrown, THE_CANCEL);
    assert.ok(
      !(thrown instanceof Error),
      "the value of the application arrived as an error of popnei",
    );
    assert.equal(calls.length, 2);
    const stoppedAtBytes = calls.at(1);
    assert.ok(
      (stoppedAtBytes?.bytesRead ?? 0) >= 4 * 1024 * 1024,
      `the call that threw says ${stoppedAtBytes?.bytesRead} bytes read`,
    );
    assert.ok(
      (stoppedAtBytes?.bytesRead ?? 0) < (stoppedAtBytes?.numBytes ?? 0),
      `the call that threw says the whole file, ${JSON.stringify(stoppedAtBytes)}`,
    );
  } finally {
    variants.free();
  }
});

test("a function that throws over a gzipped VCF gives its value and not a broken stream", () => {
  const variants = openVcf(MANY_VCF_GZ);
  // The read that fails is inside the decompressor, which is reading the
  // first member of the file for the header of the VCF: what it makes of a
  // read that failed is the error of a gzip stream that ends in the middle,
  // and what the application is given is its own value all the same.
  const calls = stoppedAt(variants, 0);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerVarDistribs(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 1);
    assert.equal(calls.at(0)?.numBytes, BYTES_OF_MANY_VCF_GZ);
  } finally {
    variants.free();
  }
});

test("a function that throws over a vars file gives its value and not a file cut short", () => {
  const vcf = openVcf(MANY_VCF, EVERY_VARIANT);
  const file = writeVars(vcf).bytes;
  vcf.free();
  const variants = openVars(file);
  // The first read of a pass over a vars file is the ten last bytes of the
  // file, which say how long its footer is, and it is made with the
  // `read_exact` of `bytes_at`: what that turns a read of fewer bytes into
  // is the error of a vars file that was cut short, and what the application
  // is given is its own value.
  const calls = stoppedAt(variants, 0);
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerIndividualStats(variants);
    });
    assert.equal(thrown, THE_CANCEL);
    assert.ok(
      !(thrown instanceof Error),
      "the value of the application arrived as an error of popnei",
    );
    assert.equal(calls.length, 1);
  } finally {
    variants.free();
  }
});

test("a function that throws only in the calls that end a run does not stop it", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  // The calls of this run are the first read of the pass, at 0 bytes, and
  // the one the end of the run makes, at the 117346 of the file. The second
  // is made when the run is over and there is no read left for it to end,
  // so what it throws is dropped and the consumer gives its result.
  const calls = stoppedAt(variants, 1);
  try {
    const distribs = calcPerVarDistribs(variants);
    assert.equal(distribs.passStats.numVars, VARIANTS_OF_MANY_VCF);
    assert.deepEqual(
      calls.map((call) => call.bytesRead),
      [0, BYTES_OF_MANY_VCF],
    );
  } finally {
    variants.free();
  }
});

test("a function that throws at every call is called once in the pass it stopped", () => {
  // The VCF of 150000 variants, whose pass is told at its first read, in the
  // middle of the file and at the end of the run. The function throws a
  // value of its own at every call, so a pass that read on after the first
  // one, or a call that ended a pass that was stopped, would be a second
  // call and a second value.
  const variants = openVcf(manyVariantsVcf(150000));
  const thrownByTheFunction: unknown[] = [];
  variants.onProgress(() => {
    const cancel = { whichCall: thrownByTheFunction.length };
    thrownByTheFunction.push(cancel);
    throw cancel;
  });
  try {
    const thrown = whatWasThrownBy(() => {
      calcPerIndividualStats(variants);
    });
    assert.equal(thrownByTheFunction.length, 1);
    assert.equal(thrown, thrownByTheFunction.at(0));
  } finally {
    variants.free();
  }
});

test("a function that throws at the first read of the second pass stops the pca", () => {
  const variants = openVcf(MANY_VCF);
  // The principal components of the variants build both of their readers
  // before they ask either for a block, and a reader reads when it is built,
  // so the calls of a run are pass 1 at 0 bytes, pass 2 at 0, and then the
  // two the end of the run makes. The second call is the first read of the
  // second pass.
  const calls = stoppedAt(variants, 1);
  try {
    const thrown = whatWasThrownBy(() => {
      doPcaFromVariants(variants, {
        numPrinComps: 10,
        transformToBiallelic: true,
      });
    });
    assert.equal(thrown, THE_CANCEL);
    // The third call is the one the end of the run makes for its first pass,
    // which the application did not stop: that pass read the 617 bytes of
    // the header of `many.vcf` when its reader was built and no more, since
    // the run ended before either reader was asked for a block. The pass
    // that was stopped, the second, is the one no call ends.
    assert.deepEqual(
      calls.map((call) => ({ pass: call.pass, bytesRead: call.bytesRead })),
      [
        { pass: 1, bytesRead: 0 },
        { pass: 2, bytesRead: 0 },
        { pass: 1, bytesRead: 617 },
      ],
    );
  } finally {
    variants.free();
  }
});

test("an iteration of blocks that is stopped throws the value of the application", () => {
  const variants = openVcf(MANY_VCF, EVERY_VARIANT);
  // The first read of the pass is the header of the VCF, read when the
  // iteration is opened, so the pass is stopped before it gave a block.
  const calls = stoppedAt(variants, 0);
  try {
    let numVars = 0;
    const thrown = whatWasThrownBy(() => {
      for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
        numVars += block.numVars;
      }
    });
    assert.equal(thrown, THE_CANCEL);
    assert.equal(calls.length, 1);
    assert.equal(numVars, 0);
  } finally {
    variants.free();
  }
});

test("a function that frees the variants while an iteration reads them does not trap", () => {
  // The VCF of 150000 variants, 6188977 bytes, whose pass is told at its
  // first read, in the middle of the file and at the end of the run: the
  // source is freed at the call in the middle, which says 4194313 bytes
  // read, so the iteration has 1994664 bytes of the file left to read after
  // it.
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

test("the counts of a pass read from inside the function while it reads a block are refused", () => {
  // The VCF of 150000 variants, 6188977 bytes, whose pass is told at its
  // first read, in the middle of the file and at the end of the run. The
  // call in the middle is made from inside the read of a block, where
  // wasm-bindgen holds the pass for the length of that read: a read of the
  // counts there fails that borrow and leaves behind a reference to the
  // pass that nothing drops, so the pass and the bytes of the file it holds
  // stay in the memory of wasm for as long as the page lives.
  const variants = openVcf(manyVariantsVcf(150000));
  const blocks = variants.iterBlocks({ numVarsPerBlock: 1000 });
  // The function is set after the pass was built, because the call at its
  // first read is made inside `iterBlocks`, where there is no pass yet to
  // read the counts of.
  const whatTheCountsGave: unknown[] = [];
  variants.onProgress(() => {
    try {
      whatTheCountsGave.push(blocks.passStats.numVars);
    } catch (thrown: unknown) {
      whatTheCountsGave.push(thrown);
    }
  });
  let numVars = 0;
  try {
    for (const block of blocks) {
      numVars += block.numVars;
    }
  } finally {
    variants.free();
  }
  // The pass gave every variant of the file: the refusal is the counts', and
  // it changed nothing of the reading.
  assert.equal(numVars, 150000);
  assert.equal(whatTheCountsGave.length, 2);
  const [whileItRead, whenItEnded] = whatTheCountsGave;
  assert.ok(
    whileItRead instanceof Error,
    `the counts gave ${String(whileItRead)} while the pass read a block`,
  );
  assert.match(whileItRead.message, /while it is reading a block/);
  // The call that ends the run is made when the pass is over, and there the
  // counts are the ones it was freed with.
  assert.equal(whenItEnded, 150000);
});
