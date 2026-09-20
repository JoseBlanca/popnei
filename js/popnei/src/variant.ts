/**
 * The handle a user holds: a source of variants and its individuals.
 *
 * It lives in the memory of wasm, which the garbage collector of JavaScript
 * does not see, so what a user holds is freed by hand: `free()` on the
 * handle, and, for one pass over the variants, the `finally` of the
 * generator that `iterBlocks` gives.
 */

import type { Blocks, VcfSource } from "../wasm/popnei.js";

import type { Block } from "./block.js";
import { blockOf } from "./block.js";
import { theWasmHasToBeLoaded } from "./core.js";

/** What each block of `iterBlocks` carries, and how many variants it holds. */
export interface IterBlocksOptions {
  /**
   * What each block carries besides the genotypes, among `"chrom"`,
   * `"pos"`, `"id"`, `"alleles"` and `"qual"`. The chromosome and the
   * position are one field of the reader, so asking for one fills both, and
   * any other name is an `Error`. `["chrom", "pos"]` when it is not given.
   */
  fields?: readonly string[];
  /**
   * How many variants a block holds. When it is not given, the number that
   * gives a block about five million genotypes, never fewer than 100
   * variants and never more than 10000.
   */
  numVarsPerBlock?: number;
}

/** What each block carries when `iterBlocks` is not asked for fields. */
const FIELDS_OF_A_BLOCK = ["chrom", "pos"];

/**
 * How many passes over a source hold memory of wasm that has not been freed.
 *
 * It is for the tests of the package and is not exported to a user: a pass
 * is counted from the call of `iterBlocks` that opened it until its
 * generator ends, is left with a `break` or throws.
 */
let openPasses = 0;

/** How many passes over a source have not been freed. For the tests. */
export function numberOfOpenPasses(): number {
  return openPasses;
}

/**
 * A source of variants: a VCF with the options it is read with.
 *
 * It holds no genotypes. A user gets one from `openVcf` and gives it to as
 * many calculations as they want: each one reads the source again and runs
 * its loop over the variants inside the Rust core, so the dataset is never
 * in memory as a whole. The genotypes come out of it through `iterBlocks`
 * and through nothing else.
 *
 * It is pyNei's `Variants` under the word of `docs/glossary.md`: what pyNei
 * calls a sample is here an individual, one organism that was genotyped.
 */
export class Variants {
  /** The VCF in the memory of wasm, and `null` once `free` took it. */
  #source: VcfSource | null;
  #individuals: readonly string[];
  #ploidy: number;

  /**
   * The handle over `source`, which `openVcf` builds.
   *
   * The names of the individuals and the ploidy are read here, from the
   * header that was read once, so that they answer without the core.
   */
  constructor(source: VcfSource) {
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
   * never iterated keeps that memory until the garbage collector reaches
   * it.
   *
   * @throws {Error} When a name of `fields` is not a field of a block, when
   * `numVarsPerBlock` is 0, when the source was freed, and when `init` has
   * not been awaited. Each of them is thrown by this call and not by the
   * first block.
   */
  iterBlocks(options: IterBlocksOptions = {}): IterableIterator<Block> {
    theWasmHasToBeLoaded();
    const source = this.#sourceThatWasNotFreed();
    const fields = [...(options.fields ?? FIELDS_OF_A_BLOCK)];
    const pass = source.blocks(fields, options.numVarsPerBlock);
    openPasses += 1;
    return blocksOfThePass(pass);
  }

  /**
   * Gives back the memory of wasm the source holds.
   *
   * Every call of `iterBlocks` after it throws. The names of the
   * individuals and the ploidy still answer: they are in JavaScript.
   */
  free(): void {
    this.#source?.free();
    this.#source = null;
  }

  /** The source, or the `Error` of a source that was freed. */
  #sourceThatWasNotFreed(): VcfSource {
    if (this.#source === null) {
      throw new Error(
        "popnei: these variants were freed, so the source cannot be read again",
      );
    }
    return this.#source;
  }
}

/**
 * The blocks of one pass, which frees the memory of wasm of every block as
 * soon as its columns are copied out, and the pass itself when the
 * iteration is over.
 *
 * The `finally` runs when the iteration ends, when the caller leaves it with
 * a `break`, which calls `return` on the generator, and when a block throws.
 */
function* blocksOfThePass(pass: Blocks): Generator<Block, void, undefined> {
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
  } finally {
    pass.free();
    openPasses -= 1;
  }
}
