/**
 * That `numPassesOf` says how many times a consumer will read the source
 * before it is started.
 *
 * A consumer is the function that runs a `Variants`, `calcKinship` or the
 * iteration of `iterBlocks`, and a pass is one reading of the source from
 * its start. Every consumer of popnei makes one pass, except two. The
 * principal components of the variants read the source a second time for
 * the weights of the variants: those need the eigenvectors, which are known
 * when the first pass ends. The association study reads it a second time
 * when it is asked for the GRAMMAR-Gamma approximation: the factor of that
 * approximation is estimated from the first block of that second pass. A
 * page that draws one bar for a whole run asks this function how many
 * passes the bar covers, as `docs/specs/js_sources.md` has it under "How
 * many passes a consumer makes".
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { ConsumerName } from "popnei";
import { init, numPassesOf } from "popnei";

import { theConsumersTheCrateNames } from "./consumers.ts";

/**
 * How many passes each consumer makes when it is asked with no options.
 *
 * Eleven of the twelve read the source once, the principal components of the
 * variants being the one that reads it twice for the weights its ten
 * components ask for. The association study is among the eleven: it reads
 * the source twice only when it is asked for the GRAMMAR-Gamma
 * approximation, which the default does not ask for, and the two tests
 * below are of that option.
 *
 * It is a `Record` over `ConsumerName` and not an array of names, because
 * TypeScript makes a `Record` hold every member of the union and lets an
 * array hold as few as it likes: a name added to `ConsumerName` and to no
 * consumer of the crate compiled and type-checked while this was an array,
 * and now it is a `TS2741` here and the comparison with the crate's own
 * names below is what refuses it at run time.
 */
const THE_PASSES_OF_EACH_CONSUMER: Record<ConsumerName, number> = {
  calcPerVarDistribs: 1,
  calcPerIndividualStats: 1,
  calcPairwiseKosmanDists: 1,
  calcPopDists: 1,
  calcPopDiversity: 1,
  calcRogersHuffR2Matrix: 1,
  calcLdAndDistPerPop: 1,
  calcKinship: 1,
  calcGwas: 1,
  writeVars: 1,
  iterBlocks: 1,
  doPcaFromVariants: 2,
};

/** The twelve consumers of the package, which are the keys of that table. */
const THE_CONSUMERS = Object.keys(
  THE_PASSES_OF_EACH_CONSUMER,
) as ConsumerName[];

test("the pca of the variants makes two passes when weights are asked for", async () => {
  await init();
  assert.equal(numPassesOf("doPcaFromVariants", { numPrinComps: 10 }), 2);
});

test("the pca of the variants makes one pass when no weights are asked for", async () => {
  await init();
  assert.equal(numPassesOf("doPcaFromVariants", { numPrinComps: 0 }), 1);
});

test("the pca of the variants with no options makes the two passes of its ten components", async () => {
  await init();
  assert.equal(numPassesOf("doPcaFromVariants"), 2);
});

test("the study makes two passes when the grammar gamma approximation is asked for", async () => {
  await init();
  assert.equal(numPassesOf("calcGwas", { useGrammarGammaApprox: true }), 2);
});

test("the study makes one pass when the grammar gamma approximation is not asked for", async () => {
  await init();
  assert.equal(numPassesOf("calcGwas", { useGrammarGammaApprox: false }), 1);
});

test("the study with no options makes the one pass of its exact denominator", async () => {
  await init();
  assert.equal(numPassesOf("calcGwas"), 1);
});

for (const [consumer, numPasses] of Object.entries(
  THE_PASSES_OF_EACH_CONSUMER,
)) {
  test(`numPassesOf says ${consumer} with no options makes ${numPasses}`, async () => {
    await init();
    assert.equal(numPassesOf(consumer as ConsumerName), numPasses);
  });
}

test("a name that is of no consumer is refused, with the names that are", async () => {
  await init();
  assert.throws(() => numPassesOf("calcKinships" as ConsumerName), {
    name: "Error",
    message: /`calcKinships` is not a consumer of popnei.*calcKinship/s,
  });
});

test("every name the refusal gives is a name numPassesOf takes", async () => {
  await init();
  // The twelve names live twice in the binding crate, in the function that
  // takes a name and in the list the message of a refused name is built
  // from, and nothing else holds the two together. A name that the message
  // gives and the function refuses fails the loop below; a name the
  // function takes and the message leaves out fails the comparison with the
  // twelve of this file, which are the twelve of
  // `docs/specs/js_sources.md`.
  const names = theConsumersTheCrateNames();
  assert.deepEqual([...names].sort(), [...THE_CONSUMERS].sort());
  for (const name of names) {
    assert.doesNotThrow(() => numPassesOf(name as ConsumerName));
  }
});

test("a numPrinComps that is not a whole number of 0 or more is refused", async () => {
  await init();
  assert.throws(() => numPassesOf("doPcaFromVariants", { numPrinComps: -1 }), {
    name: "Error",
    message: /`numPrinComps` is a whole number of 0 or more/,
  });
});

test("a useGrammarGammaApprox that is not a boolean is refused", async () => {
  await init();
  assert.throws(
    () => numPassesOf("calcGwas", { useGrammarGammaApprox: "yes" }),
    {
      name: "Error",
      message: /`useGrammarGammaApprox` is /,
    },
  );
});

test("a consumer that is not a name is refused where it was written", async () => {
  await init();
  // A number or an object reaches the core as a pointer into the memory of
  // wasm, which read the bytes of a name where there are none: the wasm
  // trapped with a `RuntimeError` that names neither the argument nor what
  // was given, and the call after it found the module dead.
  assert.throws(() => numPassesOf(42 as unknown as ConsumerName), {
    name: "Error",
    message: /`consumer` is a name, and the number 42 was given/,
  });
  assert.throws(() => numPassesOf({} as unknown as ConsumerName), {
    name: "Error",
    message: /`consumer` is a name, and an object of the type `Object` was given/,
  });
});
