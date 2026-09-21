/**
 * The handle a user holds: a source of variants and its individuals.
 *
 * It lives in the memory of wasm, which the garbage collector of JavaScript
 * does not see, so what a user holds is freed by hand: `free()` on the
 * handle, and, for one pass over the variants, the `finally` of the
 * generator that `iterBlocks` gives.
 */

import type { Blocks, VarsSource, VcfSource } from "../wasm/popnei.js";

import {
  namesOfFields,
  whatWasGiven,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import type { Block, Field } from "./block.js";
import { blockOf } from "./block.js";
import { theWasmHasToBeLoaded } from "./core.js";

/**
 * A source of variants in the memory of wasm: a VCF that `openVcf` opened,
 * or a vars file that `openVars` did.
 *
 * The two answer the same calls, so everything a `Variants` does works over
 * either of them.
 */
export type SourceOfVariants = VcfSource | VarsSource;

/**
 * The key of the method that hands the source of a `Variants` to the
 * functions of the package that write it, `writeVars` of `io_vars.ts`.
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
   * How many variants a block holds, a whole number of 1 or more. When it
   * is not given, the size the core works out from the number of
   * individuals.
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
 * file.
 *
 * It holds no genotypes. A user gets one from `openVcf` or from `openVars`
 * and gives it to as many calculations as they want: each one reads the
 * source again and runs its loop over the variants inside the Rust core, so
 * the dataset is never in memory as a whole. The genotypes come out of it
 * through `iterBlocks` and through nothing else.
 *
 * It is pyNei's `Variants` under the word of `docs/glossary.md`: what pyNei
 * calls a sample is here an individual, one organism that was genotyped.
 */
export class Variants {
  /** The file in the memory of wasm, and `null` once `free` took it. */
  #source: SourceOfVariants | null;
  #individuals: readonly string[];
  #ploidy: number;

  /**
   * The handle over `source`, which `openVcf` and `openVars` build.
   *
   * The names of the individuals and the ploidy are read here, from the
   * header of the VCF or the schema of the vars file that was read once, so
   * that they answer without the core.
   */
  constructor(source: SourceOfVariants) {
    this.#source = source;
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
   * @throws {Error} When `fields` is not an array of names, when a name of
   * it is not a field of a block, when `numVarsPerBlock` is not a whole
   * number of 1 or more, when the source was freed, and when `init` has not
   * been awaited. Each of them is thrown by this call and not by the first
   * block.
   */
  iterBlocks(options: IterBlocksOptions = {}): IterableIterator<Block> {
    theWasmHasToBeLoaded();
    const source = this.#sourceThatWasNotFreed();
    const fields = namesOfFields(
      "fields",
      options.fields === undefined ? FIELDS_OF_A_BLOCK : options.fields,
    );
    const numVarsPerBlock =
      options.numVarsPerBlock === undefined
        ? undefined
        : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
    return blocksOfThePass(source.blocks(fields, numVarsPerBlock));
  }

  /**
   * Gives back the memory of wasm the source holds.
   *
   * Every call of `iterBlocks` after it throws. The names of the
   * individuals and the ploidy still answer: they are in JavaScript. A
   * second call is not an error: it has nothing left to give back.
   */
  free(): void {
    this.#source?.free();
    this.#source = null;
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
   * The source, for the functions of the package that write it. The symbol
   * that names it is this module's, so no user reaches it.
   */
  [THE_SOURCE](): SourceOfVariants {
    return this.#sourceThatWasNotFreed();
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
}

/**
 * The source that `value` holds, when it is a `Variants` that was not
 * freed, for the functions of the package that read or write a source.
 *
 * @throws {Error} When `value` is not a `Variants`, which names what was
 * given, and when it was freed.
 */
export function sourceOfTheVariants(
  argument: string,
  value: unknown,
): SourceOfVariants {
  if (!(value instanceof Variants)) {
    throw new Error(
      `popnei: \`${argument}\` is what openVcf or openVars gives, and ` +
        `${whatWasGiven(value)} was given`,
    );
  }
  return value[THE_SOURCE]();
}

/**
 * The blocks of one pass, which frees the memory of wasm of every block as
 * soon as its columns are copied out, and the pass itself when the
 * iteration is over.
 *
 * The `finally` runs when the iteration ends, when the caller leaves it with
 * a `break`, which calls `return` on the generator, and when a block throws.
 * It does not run for an iterator that was never started, which is why the
 * pass is counted here and not in `iterBlocks`: a generator that never ran
 * its first line has no `finally` to run either.
 */
function* blocksOfThePass(pass: Blocks): Generator<Block, void, undefined> {
  openPasses += 1;
  let theIterationFailed = false;
  try {
    for (;;) {
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
      pass.free();
    } catch (freeingFailed) {
      // A free of an object of wasm that is still borrowed throws, which a
      // panic inside the core leaves behind. The error that is on its way
      // out says what went wrong; this one would hide it.
      if (!theIterationFailed) {
        throw freeingFailed;
      }
    }
  }
}
