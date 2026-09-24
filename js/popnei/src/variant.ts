/**
 * The handle a user holds: a source of variants, its individuals and the
 * steps that were put on it, and the counts of a pass over it.
 *
 * It lives in the memory of wasm, which the garbage collector of JavaScript
 * does not see, so what a user holds is freed by hand: `free()` on the
 * handle, and, for one pass over the variants, the `finally` of the
 * generator that the blocks of `iterBlocks` come from.
 */

import type {
  Blocks as PassOfTheCore,
  PassCounts,
  Steps as StepsOfTheCore,
  VarsSource,
  VcfSource,
} from "../wasm/popnei.js";
import { Steps } from "../wasm/popnei.js";

import {
  aNumber,
  distanceInBasePairs,
  namesOf,
  whatWasGiven,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import type { Block, Field } from "./block.js";
import { blockOf } from "./block.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { FilteringStats, Step } from "./filters.js";

/**
 * A source of variants in the memory of wasm: a VCF that `openVcf` opened,
 * or a vars file that `openVars` did.
 *
 * The two answer the same calls, so everything a `Variants` does works over
 * either of them.
 */
export type SourceOfVariants = VcfSource | VarsSource;

/**
 * The source of a `Variants` and the steps that were put on it, which is
 * what a pass over it is built from.
 */
export interface SourceAndSteps {
  source: SourceOfVariants;
  steps: StepsOfTheCore;
}

/**
 * The counts of one pass over a source of variants.
 *
 * Every consumer of a `Variants` gives them back with its result, and so do
 * the blocks of `iterBlocks`. A pass is one reading of the source from its
 * start, through the steps the `Variants` had when it started, so these
 * counts are of that reading alone and of what the `Variants` holds
 * afterwards nothing reaches them.
 */
export interface PassStats {
  /** How many variants the consumer took, after the steps. */
  numVars: number;

  /**
   * How many variants each filter of the pass was given and kept, under the
   * kind of the filter, `"missing_data"`, `"maf"`, `"obs_het"` or `"ld"`,
   * in the order of the steps. It is empty for a pass with no filter.
   */
  filtering: Record<string, FilteringStats>;
}

/**
 * How far a pass over the source has got, which a page draws a bar from.
 *
 * A pass is one reading of the source from its start, and a run is one call
 * of one consumer, `calcKinship` or the iteration of `iterBlocks`, with the
 * passes it makes. `Variants.onProgress` is where the function that is told
 * these four numbers is set, and it says when the calls are made.
 */
export interface Progress {
  /** How many bytes of the file this pass has read, `numBytes` at most. */
  bytesRead: number;

  /** How many bytes the file holds. */
  numBytes: number;

  /** Which pass of the run is reading, 1 for the first. */
  pass: number;

  /** How many passes the run makes, `numPassesOf` of its consumer. */
  numPasses: number;
}

/**
 * The blocks of one pass, one after another, and the counts of it.
 *
 * It is what `iterBlocks` gives: the blocks in a `for ... of`, and a
 * `passStats` that says how many variants have come out of it and what each
 * filter of the pass has been given and kept.
 */
export interface Blocks extends IterableIterator<Block> {
  /**
   * The counts of the pass as it stands.
   *
   * Read when the pass is over, they are of everything it gave. Read
   * between two blocks, they are of the blocks that came out so far, which
   * can be fewer variants than the filters of the pass have kept: the
   * reader that cuts the blocks to the size the user asked for keeps the
   * variants of the next block.
   */
  readonly passStats: PassStats;
}

/**
 * The counts that the core gives for one pass, as the numbers a user reads.
 *
 * The core gives the filters of the chain of readers, the outermost first,
 * and a user reads them in the order of the steps, which is the one the
 * filters were put on the `Variants` in and the reverse of the chain's. The
 * counts hold memory of wasm until they are read, so they are freed here.
 */
export function passStatsOf(counts: PassCounts): PassStats {
  try {
    const kinds = counts.kinds();
    const varsProcessed = counts.vars_processed();
    const varsKept = counts.vars_kept();
    const filtering: Record<string, FilteringStats> = {};
    for (let filter = kinds.length - 1; filter >= 0; filter -= 1) {
      const kind = kinds[filter];
      const processed = varsProcessed[filter];
      const kept = varsKept[filter];
      if (kind === undefined || processed === undefined || kept === undefined) {
        throw new Error(
          `popnei: the counts of this pass hold ${kinds.length} filters and ` +
            `not the counts of every one of them`,
        );
      }
      filtering[kind] = { varsProcessed: processed, varsKept: kept };
    }
    return { numVars: counts.num_vars(), filtering };
  } finally {
    counts.free();
  }
}

/**
 * The kind of an argument whose value is the threshold of a filter, which
 * is in `arg_thresholds`, of one whose value is the names of the
 * individuals to keep, which are in `arg_individuals`, and of one whose
 * value is a window of base pairs, which is in `arg_distances`. They are
 * the three numbers `arg_kinds` of the binding crate gives.
 */
const A_THRESHOLD = 0;
const THE_NAMES_OF_INDIVIDUALS = 1;
const A_DISTANCE = 2;

/** The steps of the core as the steps a user reads, in their order. */
function stepsOf(steps: StepsOfTheCore): Step[] {
  const kinds = steps.kinds();
  const names = steps.arg_names();
  const numArgsPerStep = steps.num_args_per_step();
  // The value of an argument crosses in the array of its kind: the
  // threshold of a filter is one number, and the individuals to keep are
  // their names, one argument after another. Which array each argument is
  // read from is the kind that crosses beside it, and an argument of a
  // kind this version of the package does not know is thrown for and not
  // read as a threshold.
  const argKinds = steps.arg_kinds();
  const numNamesPerArg = steps.num_names_per_arg();
  const thresholds = steps.arg_thresholds();
  const individuals = steps.arg_individuals();
  const distances = steps.arg_distances();
  const ofEachStep: Step[] = [];
  let firstArg = 0;
  let nextThreshold = 0;
  let nextDistance = 0;
  let firstName = 0;
  for (const [step, kind] of kinds.entries()) {
    // The arguments of every step cross flat, the ones of the first step
    // first, and how many each step has is what cuts them apart.
    const numArgs = numArgsPerStep[step];
    if (numArgs === undefined) {
      throw new Error(
        `popnei: these variants hold ${kinds.length} steps and not how many ` +
          "arguments every one of them has",
      );
    }
    const args: Record<string, unknown> = {};
    for (let argument = firstArg; argument < firstArg + numArgs; argument += 1) {
      const name = names[argument];
      const argKind = argKinds[argument];
      const numNames = numNamesPerArg[argument];
      if (name === undefined || argKind === undefined || numNames === undefined) {
        throw new Error(
          `popnei: the step \`${kind}\` of these variants holds ${numArgs} ` +
            `arguments and not the name and the kind of every one of them`,
        );
      }
      if (argKind === A_THRESHOLD) {
        const threshold = thresholds[nextThreshold];
        if (threshold === undefined) {
          throw new Error(
            `popnei: the argument \`${name}\` of the step \`${kind}\` of these ` +
              "variants is a threshold and has no number",
          );
        }
        args[name] = threshold;
        nextThreshold += 1;
      } else if (argKind === A_DISTANCE) {
        const distance = distances[nextDistance];
        if (distance === undefined) {
          throw new Error(
            `popnei: the argument \`${name}\` of the step \`${kind}\` of these ` +
              "variants is a window of base pairs and has no number",
          );
        }
        args[name] = distance;
        nextDistance += 1;
      } else if (argKind === THE_NAMES_OF_INDIVIDUALS) {
        const kept = individuals.slice(firstName, firstName + numNames);
        if (kept.length !== numNames) {
          throw new Error(
            `popnei: the argument \`${name}\` of the step \`${kind}\` of these ` +
              `variants holds ${numNames} names and not every one of them`,
          );
        }
        args[name] = kept;
        firstName += numNames;
      } else {
        throw new Error(
          `popnei: the argument \`${name}\` of the step \`${kind}\` of these ` +
            `variants is of the kind ${argKind}, which this version of the ` +
            "package does not know; the package and the wasm it was built " +
            "with are of one version",
        );
      }
    }
    firstArg += numArgs;
    ofEachStep.push({ kind, args });
  }
  return ofEachStep;
}

/**
 * The key of the method that hands the source of a `Variants` and its steps
 * to the functions of the package that write it, `writeVars` of
 * `io_vars.ts`.
 *
 * It is a symbol declared here and exported by nothing, so what a user of
 * the package sees is the same as before it was added, and a method with a
 * name would be one more thing they can call.
 */
const THE_SOURCE: unique symbol = Symbol("popnei: the source of a Variants");

/** What each block of `iterBlocks` carries, and how many variants it holds. */
export interface IterBlocksOptions {
  /**
   * What each block carries besides the genotypes, among `"chrom"`,
   * `"pos"`, `"id"`, `"alleles"` and `"qual"`. The chromosome and the
   * position are one field of the reader, so asking for one fills both, and
   * any other name is an `Error`. `["chrom", "pos"]` when it is not given.
   */
  fields?: readonly Field[];
  /**
   * How many variants a block holds, a whole number of 1 or more and at
   * most 4294967295. When it is not given, the size the core works out from
   * the number of individuals.
   */
  numVarsPerBlock?: number;
}

/** What each block carries when `iterBlocks` is not asked for fields. */
const FIELDS_OF_A_BLOCK: Field[] = ["chrom", "pos"];

/**
 * How many iterations of blocks are running, each holding the memory of
 * wasm of its reader and of the block it is building.
 *
 * It is counted inside the generator, from the first `next` of an iteration
 * until it ends, is left with a `break` or throws, and it is for the tests
 * of the package: no entry point exports it.
 */
let openPasses = 0;

/** How many iterations of blocks are running. For the tests. */
export function numberOfOpenPasses(): number {
  return openPasses;
}

/**
 * A source of variants: a VCF with the options it is read with, or a vars
 * file, and the steps that were put on it.
 *
 * It holds no genotypes. A user gets one from `openVcf` or from `openVars`
 * and gives it to as many calculations as they want: each one reads the
 * source again and runs its loop over the variants inside the Rust core, so
 * the dataset is never in memory as a whole. The genotypes come out of it
 * through `iterBlocks` and through nothing else.
 *
 * What is done with it is of two kinds, and what a call gives back says
 * which. A step, a filter of `docs/specs/filters.md`, is a method that adds
 * itself to the list of steps, reads nothing and returns nothing, and
 * `steps` is that list. A consumer, `iterBlocks`, `writeVars` or the
 * function of a calculation, gives something back, and it runs the steps: it
 * makes as many passes over the source as it needs, each one built from the
 * steps the `Variants` has when that pass starts. So a step added between
 * two consumers holds for the second, and one added while a pass runs holds
 * from the next pass.
 *
 * It is pyNei's `Variants` under the word of `docs/glossary.md`: what pyNei
 * calls a sample is here an individual, one organism that was genotyped.
 */
export class Variants {
  /** The file in the memory of wasm, and `null` once `free` took it. */
  #source: SourceOfVariants | null;
  /**
   * The steps in the memory of wasm, which every pass is built from, and
   * `null` once `free` took them.
   */
  #steps: StepsOfTheCore | null;
  #individuals: readonly string[];
  #ploidy: number;

  /**
   * The handle over `source`, which `openVcf` and `openVars` build, with no
   * step on it.
   *
   * The names of the individuals and the ploidy are read here, from the
   * header of the VCF or the schema of the vars file that was read once, so
   * that they answer without the core. The names go to the steps as well,
   * which resolve the names of a filter of individuals against them.
   */
  constructor(source: SourceOfVariants) {
    const steps = new Steps(source.individuals());
    this.#source = source;
    this.#steps = steps;
    // The names of the individuals the next pass gives, which the steps are
    // what says: those of the source until a filter of individuals is put
    // on them. They are kept in JavaScript so that `individuals` answers
    // after `free`.
    this.#individuals = Object.freeze(steps.individuals());
    this.#ploidy = source.ploidy();
  }

  /**
   * The names of the individuals the next pass gives, in its order.
   *
   * They are those of the source, in the order the source has them, until
   * `filterIndividuals` is put on the `Variants`: from then on they are the
   * ones that filter keeps, in the order they were named, which is the
   * order of the genotypes of every block. A pass changes nothing of them,
   * so they are the same read before one and after one, and they answer
   * after `free` as well: they are in JavaScript.
   */
  get individuals(): readonly string[] {
    return this.#individuals;
  }

  /** How many individuals the next pass gives the genotypes of. */
  get numIndividuals(): number {
    return this.#individuals.length;
  }

  /** How many alleles the genotype of one individual holds. */
  get ploidy(): number {
    return this.#ploidy;
  }

  /**
   * The steps that were put on this `Variants`, in their order.
   *
   * Each one has the kind of the step and the arguments it was given. A
   * pass takes the steps that are there when it starts, so this is what the
   * next consumer will run. The array is the user's own: it is built at
   * every read, and pushing a step into it puts no step on the `Variants`.
   *
   * @throws {Error} When the source was freed, and when `init` has not been
   * awaited.
   */
  get steps(): Step[] {
    theWasmHasToBeLoaded();
    return stepsOf(this.#stepsThatWereNotFreed());
  }

  /**
   * Keeps the variants whose missing rate is at most
   * `maxAllowedMissingRate`.
   *
   * The missing rate of a variant is its missing genotypes divided by all
   * the individuals of the dataset, and not by the ones that were called at
   * it. A genotype is missing when one of its alleles at least was not
   * called, so `0/.` in a VCF is a missing genotype, as it is in pyNei and
   * in bcftools.
   *
   * The call adds a step and gives nothing back. What runs it is the next
   * pass over the source, which every consumer makes: a filter added
   * between two of them holds for the second, and one added while a pass
   * runs holds from the pass after it.
   *
   * `maxAllowedMissingRate` has no default, where pyNei's is 0, which keeps
   * only the variants with every genotype called.
   *
   * @throws {Error} When the threshold is not a number from 0 to 1, both
   * included, which names the argument and the value, and when it is not
   * given at all. A second filter of this kind on the same `Variants` is an
   * `Error` too, with the threshold that is set: two thresholds of one kind
   * keep what the stricter of them keeps alone, so the second says that the
   * steps are not what their user thinks, which running a cell of a
   * notebook twice gives, and `steps` is what they hold. It also throws
   * when the variants were freed and when `init` has not been awaited.
   */
  filterByMissingData(maxAllowedMissingRate: number): void {
    theWasmHasToBeLoaded();
    this.#stepsThatWereNotFreed().filter_by_missing_data(
      aNumber("maxAllowedMissingRate", maxAllowedMissingRate),
    );
  }

  /**
   * Keeps the variants whose major allele frequency is at most
   * `maxAllowedMaf`.
   *
   * The major allele frequency of a variant, "maf" in pyNei and in popnei,
   * is the count of its commonest allele divided by its called alleles,
   * where most of the literature and plink2 give those letters to the minor
   * allele. Every allele of a multiallelic variant has its own count, and
   * an allele is counted wherever it was called, in a half called genotype
   * too. A filter at 0.95 takes out the variants that hardly vary among
   * these individuals. A variant with no called allele has no major allele
   * frequency and is not kept, whatever the threshold.
   *
   * It asks for no minimum of called data, as pyNei does not: a variant
   * with one called genotype has the frequency of the alleles of that
   * genotype. A user who does not want the variants that have little called
   * data puts `filterByMissingData` before this one.
   *
   * The call adds a step and gives nothing back.
   *
   * @throws {Error} What `filterByMissingData` throws: a threshold that is
   * not a number from 0 to 1 or is not given, a second filter of this kind,
   * variants that were freed, and `init` that was not awaited.
   */
  filterByMaf(maxAllowedMaf: number): void {
    theWasmHasToBeLoaded();
    this.#stepsThatWereNotFreed().filter_by_maf(
      aNumber("maxAllowedMaf", maxAllowedMaf),
    );
  }

  /**
   * Keeps the variants whose observed heterozygosity is at most
   * `maxAllowedObsHet`.
   *
   * The observed heterozygosity of a variant is its heterozygous genotypes
   * divided by its called ones, where a genotype is heterozygous when it is
   * called and its alleles are not all the same, at any ploidy. It takes
   * out the variants in which too many individuals are heterozygous, which
   * in most datasets are paralogous regions read as one site. A variant
   * with no called genotype has no observed heterozygosity and is not kept,
   * whatever the threshold.
   *
   * It asks for no minimum of called data, as pyNei does not: a variant
   * with one called genotype, heterozygous, has an observed heterozygosity
   * of 1.
   *
   * The call adds a step and gives nothing back.
   *
   * @throws {Error} What `filterByMissingData` throws: a threshold that is
   * not a number from 0 to 1 or is not given, a second filter of this kind,
   * variants that were freed, and `init` that was not awaited.
   */
  filterByObsHet(maxAllowedObsHet: number): void {
    theWasmHasToBeLoaded();
    this.#stepsThatWereNotFreed().filter_by_obs_het(
      aNumber("maxAllowedObsHet", maxAllowedObsHet),
    );
  }

  /**
   * Keeps the variants whose r² against every variant kept no more than
   * `maxDist` base pairs behind them on their chromosome is at most
   * `maxAllowedR2`.
   *
   * r² is the square of the correlation, across the individuals called at
   * both variants, between the dosages of two variants, where the dosage
   * of a genotype is how many of its alleles are not the major allele of
   * its variant. It is 1 when the dosage of an individual at one variant
   * fixes its dosage at the other and 0 when knowing one says nothing
   * about the other, so two variants with a high r² say the same thing
   * about these individuals, and what this filter leaves is a set of
   * variants that says each thing once. A principal component analysis or
   * a kinship over variants that repeat one another counts that stretch of
   * the genome as many times as it has variants, and this is the filter a
   * user puts before them.
   *
   * The variants a candidate is compared with, its window, are those the
   * filter has already kept that are on its chromosome and no more than
   * `maxDist` base pairs behind it, and of two variants above the
   * threshold the one that comes first is the one kept. A variant whose
   * called genotypes hold one dosage, and one with no called genotype at
   * all, is dropped at every threshold, having nothing to tell another
   * variant apart with; a pair whose r² cannot be worked out, the
   * individuals called at both holding one dosage, drops neither of the
   * two.
   *
   * The dosages are read over every individual of the dataset. A user who
   * wants them read over one population puts the filter of individuals
   * before this one. Neither argument has a default.
   *
   * The largest window is 9007199254740991 base pairs, 2^53 - 1, which is
   * `Number.MAX_SAFE_INTEGER`, the largest whole number a number of
   * JavaScript counts to one by one: the core takes a window of up to
   * 2^64 - 1, which is what a user of popnei in Python can write, and above
   * 2^53 - 1 a number of JavaScript counts in twos, so a larger window would
   * reach the core as another number than the one written. No genome comes
   * near it: the largest one known, over 1e11 base pairs in all of its
   * chromosomes together, is smaller by more than four orders of magnitude.
   *
   * The call adds a step and gives nothing back.
   *
   * @throws {Error} When `maxAllowedR2` is not a number from 0 to 1 or is
   * not given, when `maxDist` is not a whole number of base pairs from 1 to
   * 9007199254740991, and when
   * a filter of this kind is set already. It also throws when the variants
   * were freed and when `init` has not been awaited. The variants of each
   * chromosome have to come together and in order of position, which is
   * what this filter alone of popnei asks of a source: a position below
   * the one before it on the same chromosome, and a chromosome that had
   * already ended, are an `Error` thrown by the block of the pass that
   * would have held that variant, and not by this call.
   */
  filterByLd(maxAllowedR2: number, maxDist: number): void {
    theWasmHasToBeLoaded();
    this.#stepsThatWereNotFreed().filter_by_ld(
      aNumber("maxAllowedR2", maxAllowedR2),
      distanceInBasePairs("maxDist", maxDist),
    );
  }

  /**
   * Keeps the genotypes of `individuals` at every variant and drops those of
   * the rest.
   *
   * Every variant stays: the step takes columns of the genotypes away and no
   * row, so it has no entry in the counts of a pass. The individuals are
   * kept in the order they are named here, which is the order of the
   * genotypes of every block and of the rows of every result over
   * individuals, so it is also the way to put a dataset's individuals in the
   * order a user wants. pyNei's `filter_samples` keeps them in the order of
   * the source instead.
   *
   * A step of it is what every step that comes after it sees:
   * `filterByMissingData` before the call divides by all the individuals of
   * the source, and after it by the kept ones alone. `individuals` and
   * `numIndividuals` are the kept ones from the call on, since they are what
   * the next pass gives.
   *
   * The call adds a step and gives nothing back.
   *
   * @throws {Error} When `individuals` is not an array of names, which one
   * name written as a string is: the call would ask for the individuals
   * `i`, `n`, `d` and so on. A name that is not an individual of the source
   * is an `Error` that names it, where pyNei drops it in silence and gives
   * the individuals it did find; a name that is there twice is one too,
   * since one individual is kept once; and so is a call with no name,
   * because variants of nobody are no dataset. A second filter of
   * individuals on the same `Variants` is an `Error` as well: two lists keep
   * the individuals that are in both, which is one list, so the second says
   * that the steps are not what their user thinks. A user who wants two sets
   * of individuals over one file opens it twice. After any of them the steps
   * are as they were. It also throws when the variants were freed and when
   * `init` has not been awaited.
   */
  filterIndividuals(individuals: readonly string[]): void {
    theWasmHasToBeLoaded();
    const kept = namesOf("individuals", individuals, {
      oneOfThem: "individual",
      anExample: "ind00",
    });
    const steps = this.#stepsThatWereNotFreed();
    steps.filter_individuals(kept);
    // The names the next pass gives, which the core is what says: they are
    // kept here as well so that `individuals` answers after `free`.
    this.#individuals = Object.freeze(steps.individuals());
  }

  /**
   * The variants of the source, block by block, from its start.
   *
   * Every call reads the source from its start, so a `Variants` can be
   * given to one calculation after another. The size of the blocks changes
   * nothing but where the cuts fall: the blocks of a source, joined, are
   * the same for any size, and only the last one can be shorter than the
   * rest. When a variant cannot be read, the error is thrown in the place
   * of the block that would have held it, and the variants of that block
   * that were read are lost with it.
   *
   * The pass holds memory of wasm, which is freed when the iteration ends,
   * when it is left with a `break` or when it throws. An iterator that is
   * never iterated at all keeps that memory until the garbage collector
   * reaches it: wasm-bindgen registers what it generates in a
   * `FinalizationRegistry`, which frees it at a moment nobody chooses.
   *
   * What it gives is the blocks one after another and the counts of the
   * pass in their `passStats`: how many variants have come out of it, and
   * what each filter of the `Variants` was given and kept. They are read
   * after the pass is over too, when the memory of wasm it held is back.
   *
   * @throws {Error} When `fields` is not an array of names, when a name of
   * it is not a field of a block, when `numVarsPerBlock` is not a whole
   * number of 1 or more, when the source was freed, and when `init` has not
   * been awaited. Each of them is thrown by this call and not by the first
   * block.
   */
  iterBlocks(options: IterBlocksOptions = {}): Blocks {
    theWasmHasToBeLoaded();
    const source = this.#sourceThatWasNotFreed();
    const steps = this.#stepsThatWereNotFreed();
    const fields = namesOf(
      "fields",
      options.fields === undefined ? FIELDS_OF_A_BLOCK : options.fields,
      { oneOfThem: "field", anExample: "chrom" },
    );
    const numVarsPerBlock =
      options.numVarsPerBlock === undefined
        ? undefined
        : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
    // The steps of the pass are a copy of the list, which the pass takes
    // over and frees. The call it is made for refuses a field that is not
    // one of the five, inside the core and after the copy was made, and the
    // copy is not left behind: an argument of a type this crate exports
    // crosses by value, so the Rust that refuses the field owns it and drops
    // it. Measured on this build: 20000 calls refused for their field left
    // the memory of wasm at the 1310720 bytes it held before them.
    return new BlocksOfOnePass(
      source.blocks(fields, numVarsPerBlock, steps.of_a_pass()),
    );
  }

  /**
   * Sets `told` as the function that is told how far every pass over the
   * source has got, and takes the one that was set off when it is called
   * with nothing.
   *
   * While a consumer runs, the worker is inside wasm and reads no message,
   * so this is how a page learns how a run is going: the source calls
   * `told` at the first read of each pass, at the first read after every
   * 4 MiB of the file that pass has read, and at the read that finds the end
   * of the file, with the bytes read, the bytes the file holds, which pass
   * of the run is reading and how many passes the run makes. A page that
   * draws a bar from those four numbers sees it fill once per pass, so a
   * principal component analysis that reads the file twice does not look
   * broken when the bar goes back to empty.
   *
   * What `told` throws ends the pass where it was reading, and the consumer
   * throws that same value: an application that cancels a run tells its own
   * cancel from a file that could not be read with `===` and without reading
   * a message, and its worker is not ended. Whichever error the read failed
   * with inside popnei is not the one it gets. The `Variants` is then the
   * one it was, and the next run over it reads the file from its start.
   *
   * The function holds until it is set again, and setting it changes nothing
   * about the variants a pass gives. The reads of `openVcf` and `openVars`,
   * the header of a VCF and the schema of a vars file, are told to nobody:
   * they are made before there is a `Variants` to set a function on.
   *
   * It has no counterpart in the Python API, which `docs/objectives.md` asks
   * every difference between the two to be written down: what it is for is a
   * page that draws a bar and a user who presses a button, and Python reads
   * a file by its path in a program that has neither.
   *
   * @throws {Error} When `told` is given and is not a function, when the
   * variants were freed, and when `init` has not been awaited.
   */
  onProgress(told?: (progress: Progress) => void): void {
    theWasmHasToBeLoaded();
    const source = this.#sourceThatWasNotFreed();
    if (told === undefined) {
      source.on_progress(undefined);
      return;
    }
    if (typeof told !== "function") {
      throw new Error(
        "popnei: `told` is the function that is told how far a pass has got, " +
          `and ${whatWasGiven(told)} was given`,
      );
    }
    // The four numbers cross one by one and the object a user reads is built
    // here, as every other result of the package is built in TypeScript.
    source.on_progress(
      (bytesRead: number, numBytes: number, pass: number, numPasses: number) =>
        told({ bytesRead, numBytes, pass, numPasses }),
    );
  }

  /**
   * Gives back the memory of wasm the source and the steps hold.
   *
   * Every call of `iterBlocks` after it throws, and so does a `writeVars`
   * of these variants and every read of `steps`. The names of the
   * individuals and the ploidy still answer: they are in JavaScript. A
   * second call is not an error: it has nothing left to give back.
   */
  free(): void {
    this.#source?.free();
    this.#source = null;
    this.#steps?.free();
    this.#steps = null;
  }

  /**
   * The same as `free`, for `using variants = openVcf(bytes)`, which frees
   * the source at the end of the block it is written in.
   *
   * Every class wasm-bindgen generates has it, and a `Variants` that had
   * only `free` would be left in the memory of wasm by that line.
   */
  [Symbol.dispose](): void {
    this.free();
  }

  /**
   * The source and the steps, for the functions of the package that read or
   * write them. The symbol that names it is this module's, so no user
   * reaches it.
   */
  [THE_SOURCE](): SourceAndSteps {
    return {
      source: this.#sourceThatWasNotFreed(),
      steps: this.#stepsThatWereNotFreed(),
    };
  }

  /** The source, or the `Error` of a source that was freed. */
  #sourceThatWasNotFreed(): SourceOfVariants {
    if (this.#source === null) {
      throw new Error(
        "popnei: these variants were freed, so the source cannot be read again",
      );
    }
    return this.#source;
  }

  /**
   * The steps, or the `Error` of a `Variants` that was freed.
   *
   * What the message says is that they cannot be changed either, because
   * the three filters are what most often reaches it: a user who reads that
   * their steps cannot be read would look for a read they did not make.
   */
  #stepsThatWereNotFreed(): StepsOfTheCore {
    if (this.#steps === null) {
      throw new Error(
        "popnei: these variants were freed, so their steps cannot be read or " +
          "changed any more",
      );
    }
    return this.#steps;
  }
}

/**
 * The source that `value` holds and its steps, when it is a `Variants` that
 * was not freed, for the functions of the package that read or write a
 * source.
 *
 * @throws {Error} When `value` is not a `Variants`, which names what was
 * given, and when it was freed.
 */
export function sourceOfTheVariants(
  argument: string,
  value: unknown,
): SourceAndSteps {
  if (!(value instanceof Variants)) {
    throw new Error(
      `popnei: \`${argument}\` is what openVcf or openVars gives, and ` +
        `${whatWasGiven(value)} was given`,
    );
  }
  return value[THE_SOURCE]();
}

/**
 * The blocks of one pass and the counts of it, which frees the memory of
 * wasm of every block as soon as its columns are copied out, and the pass
 * itself when the iteration is over.
 *
 * It is an iterator of its own and not the generator inside it, because the
 * counts are read after the pass gave its memory back: the counts of the
 * pass are taken from it just before it is freed, and a read after that
 * answers with those. `next`, `return` and `throw` are the generator's, so
 * a `for ... of` left with a `break` frees the pass as it always did.
 */
class BlocksOfOnePass implements Blocks {
  /** The pass in the memory of wasm, and `null` once it was freed. */
  #pass: PassOfTheCore | null;
  /**
   * The counts the pass had when it was freed, which is what `passStats`
   * answers from then on, and `null` while it still answers itself.
   */
  #countsWhenItEnded: PassStats | null = null;
  #blocks: Generator<Block, void, undefined>;

  constructor(pass: PassOfTheCore) {
    this.#pass = pass;
    this.#blocks = this.#blocksOfThePass();
  }

  next(): IteratorResult<Block, void> {
    return this.#blocks.next();
  }

  /** What a `for ... of` calls when it is left with a `break`. */
  return(): IteratorResult<Block, void> {
    return this.#blocks.return(undefined);
  }

  throw(error: unknown): IteratorResult<Block, void> {
    return this.#blocks.throw(error);
  }

  [Symbol.iterator](): Blocks {
    return this;
  }

  get passStats(): PassStats {
    if (this.#pass !== null) {
      return passStatsOf(this.#pass.pass_stats());
    }
    if (this.#countsWhenItEnded === null) {
      throw new Error(
        "popnei: the counts of this pass cannot be read: it gave its memory " +
          "back before they were taken out of it",
      );
    }
    return this.#countsWhenItEnded;
  }

  /**
   * The blocks themselves.
   *
   * The `finally` runs when the iteration ends, when the caller leaves it
   * with a `break`, which calls `return` on the generator, and when a block
   * throws. It does not run for an iterator that was never started, which
   * is why the pass is counted here and not in `iterBlocks`: a generator
   * that never ran its first line has no `finally` to run either.
   */
  *#blocksOfThePass(): Generator<Block, void, undefined> {
    openPasses += 1;
    let theIterationFailed = false;
    try {
      for (;;) {
        const pass = this.#passThatIsRunning();
        const columns = pass.next_block();
        if (columns === undefined) {
          return;
        }
        let block: Block;
        try {
          block = blockOf(columns);
        } finally {
          columns.free();
        }
        yield block;
      }
    } catch (error) {
      theIterationFailed = true;
      throw error;
    } finally {
      openPasses -= 1;
      try {
        this.#endThePass();
      } catch (freeingFailed) {
        // A free of an object of wasm that is still borrowed throws, which
        // a panic inside the core leaves behind, and so does a read of its
        // counts. The error that is on its way out says what went wrong;
        // this one would hide it.
        if (!theIterationFailed) {
          throw freeingFailed;
        }
      }
    }
  }

  /**
   * The counts of the pass, taken out of it, and the memory of wasm it
   * holds, given back.
   */
  #endThePass(): void {
    const pass = this.#pass;
    if (pass === null) {
      return;
    }
    this.#pass = null;
    try {
      this.#countsWhenItEnded = passStatsOf(pass.pass_stats());
    } finally {
      pass.free();
    }
  }

  /**
   * The pass, for the generator inside it, which runs only while the pass
   * is there: what ends it is its own `finally`.
   */
  #passThatIsRunning(): PassOfTheCore {
    if (this.#pass === null) {
      throw new Error(
        "popnei: this pass gave its memory back and has no more blocks",
      );
    }
    return this.#pass;
  }
}
