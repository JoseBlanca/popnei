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

/**
 * The eleven consumers that read the source once when they are asked with no
 * options, which is every one of them but the principal components of the
 * variants.
 *
 * The association study is among them: it reads the source twice only when
 * it is asked for the GRAMMAR-Gamma approximation, which the default does
 * not ask for, and the two tests below are of that option.
 */
const THE_CONSUMERS_OF_ONE_PASS: readonly ConsumerName[] = [
  "calcPerVarDistribs",
  "calcPerIndividualStats",
  "calcPairwiseKosmanDists",
  "calcPopDists",
  "calcPopDiversity",
  "calcRogersHuffR2Matrix",
  "calcLdAndDistPerPop",
  "calcKinship",
  "calcGwas",
  "writeVars",
  "iterBlocks",
];

/**
 * The twelve consumers of the package, which are the eleven above and the
 * principal components of the variants.
 */
const THE_CONSUMERS: readonly ConsumerName[] = [
  ...THE_CONSUMERS_OF_ONE_PASS,
  "doPcaFromVariants",
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

/**
 * The names the refusal of a name that is of no consumer gives, which the
 * crate writes after the colon of that message.
 *
 * @throws {Error} When `calcKinships`, which is the name of no consumer, is
 * taken.
 */
function theNamesOfTheRefusal(): string[] {
  try {
    numPassesOf("calcKinships" as ConsumerName);
  } catch (refused) {
    const message = refused instanceof Error ? refused.message : `${refused}`;
    return message
      .slice(message.lastIndexOf(":") + 1)
      .split(",")
      .map((name) => name.trim());
  }
  throw new Error(
    "popnei: `calcKinships` is the name of no consumer and was taken",
  );
}

test("every name the refusal gives is a name numPassesOf takes", async () => {
  await init();
  // The twelve names live twice in the binding crate, in the function that
  // takes a name and in the list the message of a refused name is built
  // from, and nothing else holds the two together. A name that the message
  // gives and the function refuses fails the loop below; a name the
  // function takes and the message leaves out fails the comparison with the
  // twelve of this file, which are the twelve of
  // `docs/specs/js_sources.md`.
  const names = theNamesOfTheRefusal();
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
