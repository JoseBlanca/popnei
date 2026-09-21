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
  namesOfFields,
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
   * kind of the filter, `"missing_data"`, `"maf"` or `"obs_het"`, in the
   * order of the steps. It is empty for a pass with no filter.
   */
  filtering: Record<string, FilteringStats>;
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

/** The steps of the core as the steps a user reads, in their order. */
function stepsOf(steps: StepsOfTheCore): Step[] {
  const kinds = steps.kinds();
  const names = steps.arg_names();
  const values = steps.arg_values();
  const numArgsPerStep = steps.num_args_per_step();
  const ofEachStep: Step[] = [];
  let first = 0;
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
    const args: Record<string, number> = {};
    for (let argument = first; argument < first + numArgs; argument += 1) {
      const name = names[argument];
      const value = values[argument];
      if (name === undefined || value === undefined) {
        throw new Error(
          `popnei: the step \`${kind}\` of these variants holds ${numArgs} ` +
            `arguments and not the name and the value of every one of them`,
        );
      }
      args[name] = value;
    }
    first += numArgs;
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
   * that they answer without the core.
   */
  constructor(source: SourceOfVariants) {
    this.#source = source;
    this.#steps = new Steps();
    this.#individuals = Object.freeze(source.individuals());
    this.#ploidy = source.ploidy();
  }

  /** The names of the individuals, in the order the source has them. */
  get individuals(): readonly string[] {
    return this.#individuals;
  }

  /** How many individuals the source holds. */
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
    const fields = namesOfFields(
      "fields",
      options.fields === undefined ? FIELDS_OF_A_BLOCK : options.fields,
    );
    const numVarsPerBlock =
      options.numVarsPerBlock === undefined
        ? undefined
        : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
    // The steps of the pass are a copy of the list, made after every
    // argument was checked so that nothing refused here leaves one behind:
    // the pass takes it over and frees it.
    return new BlocksOfOnePass(
      source.blocks(fields, numVarsPerBlock, steps.of_a_pass()),
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

  /** The steps, or the `Error` of a `Variants` that was freed. */
  #stepsThatWereNotFreed(): StepsOfTheCore {
    if (this.#steps === null) {
      throw new Error(
        "popnei: these variants were freed, so their steps cannot be read again",
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
