/**
 * How many times a consumer reads the source of a `Variants`.
 *
 * A consumer is a function of this package that runs a `Variants`,
 * `calcKinship` among them, and the iteration of `iterBlocks` is one too. A
 * pass is one reading of the source from its start, through the steps the
 * `Variants` had when the consumer started, and one call of one consumer is
 * a run, of one or two passes.
 *
 * A page that draws a bar over a run asks `numPassesOf` how many passes the
 * bar covers before the run starts, so that it draws one bar for the whole
 * run and not one per pass.
 */

import {
  default_num_prin_comps as defaultNumPrinComps,
  num_passes_of as numPassesOfTheCore,
} from "../wasm/popnei.js";

import { aString, wholeNumberOfZeroOrMore } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";

/**
 * The name of the function of this package that makes the passes, which is
 * what `numPassesOf` is asked about.
 */
export type ConsumerName =
  | "calcPerVarDistribs"
  | "calcPerIndividualStats"
  | "calcPairwiseKosmanDists"
  | "calcPopDists"
  | "calcRogersHuffR2Matrix"
  | "calcKinship"
  | "doPcaFromVariants"
  | "calcGwas"
  | "writeVars"
  | "iterBlocks";

/**
 * How many passes `consumer` will make over the source it is given.
 *
 * Every consumer reads the source once, except `doPcaFromVariants` when it
 * is asked for the weights of the variants, which reads it twice: a weight
 * needs the eigenvectors, and those are known when the first pass ends. Its
 * `numPrinComps` of 0 asks for no weight and reads the source once.
 *
 * `options` is the options object the consumer takes, and only the
 * `numPrinComps` of `doPcaFromVariants` changes the answer; it is checked
 * as that function checks it. A filter of the `Variants` adds no pass: a
 * filter reads the reader of its pass.
 *
 * The progress of every pass of a run carries this same number, from the
 * same function of the binding crate, so what a page draws before the run
 * starts and what it draws while it runs count the same passes.
 *
 * `numPassesOf` has no counterpart in the Python API, which
 * `docs/objectives.md` asks every difference between the two to be written
 * down: what it is for is a page that draws a bar, and Python reads a file
 * by its path in a program that has none.
 *
 * @throws {Error} When `consumer` is not a name, when it is the name of no
 * consumer of the package, when `numPrinComps` is not a whole number of 0 or
 * more and at most 4294967295, and when `init` has not been awaited.
 */
export function numPassesOf(consumer: ConsumerName, options?: object): number {
  theWasmHasToBeLoaded();
  // The name is checked here and not in the core: what the generated code
  // hands the core for a number or an object is a pointer into the memory of
  // wasm and a length, which the module reads the bytes of a name at, and
  // the read traps the module for good where the package owes its user an
  // `Error`.
  const name = aString("consumer", consumer);
  return numPassesOfTheCore(name, numPrinCompsOf(consumer, options));
}

/**
 * The `numPrinComps` the number of passes is asked with: the one of the
 * options of `doPcaFromVariants`, and the core's default for every other
 * consumer, which ignores it.
 *
 * The default is read from the core, as `doPcaFromVariants` reads it, so
 * that the number of passes of a call that says nothing cannot drift from
 * the number of passes that call makes.
 *
 * @throws {Error} When `numPrinComps` is not a whole number of 0 or more.
 */
function numPrinCompsOf(consumer: ConsumerName, options?: object): number {
  const given = (options as { numPrinComps?: unknown } | undefined | null)
    ?.numPrinComps;
  if (consumer !== "doPcaFromVariants" || given === undefined) {
    return defaultNumPrinComps();
  }
  return wholeNumberOfZeroOrMore("numPrinComps", given);
}
