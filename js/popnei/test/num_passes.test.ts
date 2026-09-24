/**
 * That `numPassesOf` says how many times a consumer will read the source
 * before it is started.
 *
 * A consumer is the function that runs a `Variants`, `calcKinship` or the
 * iteration of `iterBlocks`, and a pass is one reading of the source from
 * its start. Every consumer of popnei makes one pass, except the principal
 * components of the variants, which reads the source a second time for the
 * weights of the variants: those need the eigenvectors, which are known
 * when the first pass ends. A page that draws one bar for a whole run asks
 * this function how many passes the bar covers, as
 * `docs/specs/js_sources.md` has it under "How many passes a consumer
 * makes".
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { ConsumerName } from "popnei";
import { init, numPassesOf } from "popnei";

/**
 * The nine consumers that read the source once, which is every one of them
 * but the principal components of the variants.
 */
const THE_CONSUMERS_OF_ONE_PASS: readonly ConsumerName[] = [
  "calcPerVarDistribs",
  "calcPerIndividualStats",
  "calcPairwiseKosmanDists",
  "calcPopDists",
  "calcRogersHuffR2Matrix",
  "calcKinship",
  "calcGwas",
  "writeVars",
  "iterBlocks",
];

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

for (const consumer of THE_CONSUMERS_OF_ONE_PASS) {
  test(`${consumer} makes one pass`, async () => {
    await init();
    assert.equal(numPassesOf(consumer), 1);
  });
}

test("a name that is of no consumer is refused, with the names that are", async () => {
  await init();
  assert.throws(() => numPassesOf("calcKinships" as ConsumerName), {
    name: "Error",
    message: /`calcKinships` is not a consumer of popnei.*calcKinship/s,
  });
});

test("a numPrinComps that is not a whole number of 0 or more is refused", async () => {
  await init();
  assert.throws(() => numPassesOf("doPcaFromVariants", { numPrinComps: -1 }), {
    name: "Error",
    message: /`numPrinComps` is a whole number of 0 or more/,
  });
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
