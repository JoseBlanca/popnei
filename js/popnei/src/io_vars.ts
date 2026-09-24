/** Reading and writing a vars file, the file popnei keeps its variants in. */

import { open_vars as openVarsOfTheCore } from "../wasm/popnei.js";

import { bytes as bytesOf, wholeNumberOfOneOrMore } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats } from "./variant.js";
import { Variants, passStatsOf, sourceOfTheVariants } from "./variant.js";

/** The vars file `writeVars` wrote, and the counts of the pass it made. */
export interface VarsWritten {
  /**
   * The bytes of the whole file, which a page offers as a download: a tab
   * has no filesystem.
   */
  bytes: Uint8Array;

  /**
   * How many variants were written, and how many each filter of the
   * `Variants` was given and kept.
   */
  passStats: PassStats;
}

/** How a vars file is written: how many variants a batch of it holds. */
export interface WriteVarsOptions {
  /**
   * How many variants one batch of the file holds, the last one aside, a
   * whole number of 1 or more and at most 4294967295. When it is not given,
   * the size popnei chooses for the number of individuals of the source,
   * which is the size of its blocks.
   */
  numVarsPerBlock?: number;
}

/**
 * The variants of the vars file in `source`, the bytes of the file.
 *
 * A vars file is one arrow IPC file, also called feather v2, which
 * `writeVars` writes from any source of variants and which pandas, R and
 * polars open as a table with no popnei installed. It is where a user keeps
 * their variants once the VCF has been read, so that the text is parsed once
 * and every later pass reads a file of arrays.
 *
 * What it gives is the handle `openVcf` gives: the names of the individuals
 * and the ploidy, which it reads from the file, and the variants through
 * `iterBlocks`. Only the columns that a pass asks for are decompressed, and
 * the variants themselves are read again at every pass, from the same bytes.
 *
 * It reads the schema of the file and its footer, so these are an `Error`
 * here and not at the first calculation: bytes that are not a vars file, a
 * format version popnei does not read, a column it knows that is of another
 * type, a file with no `gts` column or whose `gts` holds another number of
 * alleles for each variant than the individuals and the ploidy of the file
 * give, and one whose individuals name nobody. A column popnei does not
 * know is read past and is no error, which is what lets a later version of
 * the format add one. What is in the batches is
 * read block by block and refused there: a batch that popnei cannot read, of
 * a file damaged after it was written, and a file whose buffers are
 * compressed with zstd, which no build of popnei carries the code to read;
 * popnei writes lz4 and reads lz4 and no compression.
 *
 * What it returns holds memory of wasm, the bytes of the file among it,
 * until its `free()` is called.
 *
 * A `File` that a user picked in a page is read inside a web worker, which
 * section 11 of `docs/architecture.md` has and this package does not do yet:
 * the source here is the bytes of the file.
 *
 * @throws {Error} When `source` is not a `Uint8Array`, when the bytes are
 * not a vars file popnei can read, and when `init` has not been awaited.
 */
export function openVars(source: Uint8Array): Variants {
  theWasmHasToBeLoaded();
  return new Variants(openVarsOfTheCore(bytesOf("source", source)));
}

/**
 * Every variant of `variants` as the bytes of a vars file, which a page
 * offers as a download: a tab has no filesystem.
 *
 * The source is a VCF or a vars file, whichever `Variants` holds, so a file
 * read with `openVars` is written again with another size of batch, and the
 * variants that are written are the ones the steps of the `Variants` keep.
 *
 * The call reads the whole source once. The file holds the six columns of a
 * VCF, the chromosome, the position, the id, the alleles, the quality and
 * the genotypes, whether or not the user will read them, so that it can
 * stand in for the VCF in any later analysis; a source that has no alleles
 * to give gives a file without that column.
 *
 * The whole file is built in the memory of wasm, which grows and never
 * shrinks, so what the tab holds while the call runs is the source and the
 * file together, and `numVarsPerBlock` is what an application that runs out
 * of memory lowers: it is the size of the block that is read and written at
 * a time. What comes back is the caller's own array, in the heap of
 * JavaScript, which nothing has to free and which the memory of wasm does
 * not hold a second copy of.
 *
 * What it gives back are those bytes and the counts of the pass it made:
 * how many variants were written and what each filter was given and kept.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `numVarsPerBlock` is not a whole number of 1 or more and at most
 * 4294967295, when the source cannot be read, a wrong line of a VCF among
 * the causes, when the memory of the tab does not take the file, and when
 * `init` has not been awaited.
 */
export function writeVars(
  variants: Variants,
  options: WriteVarsOptions = {},
): VarsWritten {
  theWasmHasToBeLoaded();
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const numVarsPerBlock =
    options.numVarsPerBlock === undefined
      ? undefined
      : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
  // The steps of the pass are a copy of the list, made after the argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  //
  // The file comes out of the memory of wasm in pieces, each of them freed
  // there as it is copied here, and the array they are put into is the
  // user's. A file that crossed in one piece would be held twice while it
  // crossed, and the memory of wasm would keep its half of that for as long
  // as the page lives.
  const file = whileTheRunReads(() =>
    source.write_vars(numVarsPerBlock, steps.of_a_pass()),
  );
  try {
    const bytes = new Uint8Array(file.num_bytes());
    let written = 0;
    for (
      let piece = file.next_piece();
      piece !== undefined;
      piece = file.next_piece()
    ) {
      bytes.set(piece, written);
      written += piece.length;
    }
    if (written !== bytes.length) {
      throw new Error(
        `popnei: the vars file says it holds ${bytes.length} bytes and gave ${written}`,
      );
    }
    return { bytes, passStats: passStatsOf(file.pass_stats()) };
  } finally {
    file.free();
  }
}
