/**
 * The twelve consumers of a `Variants`, each with a call of it over
 * `many.vcf`, and the names the binding crate knows them by.
 *
 * A consumer is a function of the package that reads a `Variants` and gives
 * a result, the iteration of `iterBlocks` and the writer `writeVars` among
 * them, as `docs/glossary.md` has the word. Two tests loop over all twelve:
 * `stop.test.ts`, which stops each of them from inside the function that is
 * told the progress, and `progress.test.ts`, which reads the calls each of
 * them made. Both need the same calls, so they are written here once.
 *
 * The list is held to the crate's twelve by
 * [`theConsumersTheCrateNames`], which reads them out of the message a name
 * of no consumer is refused with: a consumer added to the package and not
 * to this list, or to this list and not to the crate, fails the test of
 * each of those two files that compares the two.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`.
 */

import type { ConsumerName, Variants } from "popnei";
import {
  calcGwas,
  calcKinship,
  calcLdAndDistPerPop,
  calcPairwiseKosmanDists,
  calcPerIndividualStats,
  calcPerVarDistribs,
  calcPopDists,
  calcPopDiversity,
  calcRogersHuffR2Matrix,
  doPcaFromVariants,
  numPassesOf,
  writeVars,
} from "popnei";

/** The names of the 50 individuals of `many.vcf`, `ind00` to `ind49`. */
export const THE_INDIVIDUALS = Array.from(
  { length: 50 },
  (_unused, individual) => `ind${String(individual).padStart(2, "0")}`,
);

/**
 * A number for each individual, for the association study: the values run
 * from 0 to 6 and are not the same for everybody, which a trait a model can
 * be fitted to has to be.
 */
export const THE_TRAIT = Object.fromEntries(
  THE_INDIVIDUALS.map((individual, at) => [individual, at % 7]),
);

/**
 * Two populations of 25 individuals each, for the consumers that are given
 * populations.
 */
export const THE_POPS = {
  one: THE_INDIVIDUALS.slice(0, 25),
  two: THE_INDIVIDUALS.slice(25),
};

/** One consumer of a `Variants`, with a call of it over `many.vcf`. */
export interface TheCallOfAConsumer {
  /** The name `numPassesOf` and the binding crate know it by. */
  name: ConsumerName;
  /** The options `numPassesOf` is asked with, when they change its answer. */
  options?: object;
  /** The call, which reads `variants` to the end of the file. */
  run: (variants: Variants) => void;
}

/**
 * The twelve consumers, each with a call over `many.vcf` that reads the file
 * to its end.
 *
 * `transformToBiallelic` is true where the calculation asks for it, because
 * `many.vcf` holds variants of more than two alleles and the core's default
 * refuses them. The association study is here without a kinship, which is
 * the one pass of its exact denominator; `progress.test.ts` adds the call
 * that asks for the GRAMMAR-Gamma approximation, which is its two passes.
 */
export const THE_CONSUMERS: readonly TheCallOfAConsumer[] = [
  {
    name: "calcPerVarDistribs",
    run: (variants) => {
      calcPerVarDistribs(variants);
    },
  },
  {
    name: "calcPerIndividualStats",
    run: (variants) => {
      calcPerIndividualStats(variants);
    },
  },
  {
    name: "calcPairwiseKosmanDists",
    run: (variants) => {
      calcPairwiseKosmanDists(variants);
    },
  },
  {
    name: "calcPopDists",
    run: (variants) => {
      calcPopDists(variants, THE_POPS, {
        measures: ["fst"],
        jackknifeGroup: null,
        minNumIndividuals: 1,
      });
    },
  },
  {
    name: "calcPopDiversity",
    run: (variants) => {
      calcPopDiversity(variants, { pops: THE_POPS, minNumIndividuals: 1 });
    },
  },
  {
    name: "calcRogersHuffR2Matrix",
    run: (variants) => {
      calcRogersHuffR2Matrix(variants);
    },
  },
  {
    name: "calcLdAndDistPerPop",
    run: (variants) => {
      calcLdAndDistPerPop(variants, { pops: THE_POPS });
    },
  },
  {
    name: "calcKinship",
    run: (variants) => {
      calcKinship(variants, { transformToBiallelic: true });
    },
  },
  {
    name: "doPcaFromVariants",
    options: { numPrinComps: 10 },
    run: (variants) => {
      doPcaFromVariants(variants, {
        numPrinComps: 10,
        transformToBiallelic: true,
      });
    },
  },
  {
    name: "calcGwas",
    run: (variants) => {
      calcGwas(variants, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        transformToBiallelic: true,
      });
    },
  },
  {
    name: "writeVars",
    run: (variants) => {
      writeVars(variants);
    },
  },
  {
    name: "iterBlocks",
    run: (variants) => {
      for (const _block of variants.iterBlocks()) {
        // Every block is read, so the pass reads the file to its end.
      }
    },
  },
];

/**
 * The names of the consumers the binding crate knows, which it writes after
 * the colon of the message a name of no consumer is refused with.
 *
 * The names live twice in the crate, in the function that takes a name and
 * in the list that message is built from, and the test of
 * `num_passes.test.ts` that runs each of them is what holds those two
 * together. What this gives the other files is a list that no file of
 * TypeScript wrote, so a list of this package that names another twelve is
 * caught.
 *
 * @throws {Error} When `calcKinships`, which is the name of no consumer, is
 * taken.
 */
export function theConsumersTheCrateNames(): string[] {
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
