/** Reading and writing a vars file, the file popnei keeps its variants in. */

import type { WrittenFile } from "../wasm/popnei.js";
import {
  open_vars as openVarsOfTheCore,
  open_vars_of_a_file as openVarsOfAFileOfTheCore,
} from "../wasm/popnei.js";

import {
  anObjectOfOptions,
  bytesOrFile as bytesOrFileOf,
  whatWasGiven,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { BytesOrFile } from "./io_vcf.js";
import type { PassStats } from "./variant.js";
import { Variants, passStatsOf, sourceOfTheVariants } from "./variant.js";

/** The vars file `writeVars` wrote, and the counts of the pass it made. */
export interface VarsWritten {
  /**
   * The bytes of the whole file, which a page offers as a download: a tab
   * has no file a path names.
   */
  bytes: Uint8Array;

  /**
   * How many variants were written, and how many each filter of the
   * `Variants` was given and kept.
   */
  passStats: PassStats;
}

/**
 * What `writeVars` and `writeVcf` give when they gave the file to `onBytes`
 * in pieces: the counts of the pass that wrote it.
 */
export interface WrittenInPieces {
  /**
   * How many variants were written, and how many each filter of the
   * `Variants` was given and kept.
   */
  passStats: PassStats;
}

/** The function `onBytes` of the two writers, given each piece of the file. */
export type OnBytes = (piece: Uint8Array) => void;

/**
 * How a vars file is written: how many variants a batch of it holds, and
 * whether its bytes go to a function of the application as they are
 * written.
 */
export interface WriteVarsOptions {
  /**
   * How many variants one batch of the file holds, the last one aside, a
   * whole number of 1 or more and at most 4294967295. When it is not given,
   * the size popnei chooses for the number of individuals of the source,
   * which is the size of its blocks.
   */
  numVarsPerBlock?: number;

  /**
   * The function the file is given to while the pass writes it, one piece
   * at a time, instead of the call giving it back whole. Every piece holds
   * 1048576 bytes, 1 MiB, the last one from 1 byte to that, and the pieces
   * in the order they are given are the bytes of the file. Each piece is a
   * new array of the caller's, which it can keep. The function keeps each
   * piece before it returns: one that returns a promise ends the pass with
   * an `Error`, since the writer cannot wait for it.
   */
  onBytes?: OnBytes | undefined;
}

/**
 * The variants of the vars file in `source`, the bytes of the file or the
 * file the user picked in the page.
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
 * What it returns holds memory of wasm until its `free()` is called, the
 * bytes of the file among it when the source is a `Uint8Array`.
 *
 * With a `File` or a `Blob`, popnei reads the ranges of bytes it needs
 * through `FileReaderSync`, a few MiB at a time, so the file is never in
 * the memory of wasm whole: what a pass holds of it is one range and the
 * batch it is reading, which at the size popnei writes is about 10 MB of
 * genotypes for 1000 individuals. A browser has `FileReaderSync` only
 * inside a web worker, so `openVars` of a `File` or a `Blob` on the main
 * thread of a page, or under node, is an `Error` here, at the call. What is
 * held for such a source until its `free()` is called is the handle of the
 * file, its name and its size, and not its bytes.
 *
 * A range that comes back shorter than the one popnei asked for, inside a
 * file of that size, is an `Error` in the middle of a pass and not the end
 * of the file: a browser gives a short range when the file changed on disk
 * after the page got its handle.
 *
 * @throws {Error} When `source` is neither a `Uint8Array` nor a `File` or
 * `Blob`, when a `File` or a `Blob` is given where there is no
 * `FileReaderSync`, which is everywhere but a web worker, when the bytes
 * are not a vars file popnei can read, and when `init` has not been
 * awaited.
 */
export function openVars(source: BytesOrFile): Variants {
  theWasmHasToBeLoaded();
  const file = bytesOrFileOf("source", source);
  return new Variants(
    file instanceof Uint8Array
      ? openVarsOfTheCore(file)
      : openVarsOfAFileOfTheCore(file),
  );
}

/**
 * Every variant of `variants` as the bytes of a vars file, which a page
 * offers as a download, or given to `onBytes` in pieces as they are written.
 *
 * The source is a VCF or a vars file, whichever `Variants` holds, so a file
 * read with `openVars` is written again with another size of batch, and the
 * variants that are written are the ones the steps of the `Variants` keep.
 *
 * The call reads the whole source once. The file holds the six columns of a
 * VCF, the chromosome, the position, the id, the alleles, the quality and
 * the genotypes, and whether each variant passed its FILTER, whether or not
 * the user will read them, so that it can stand in for the VCF in any later
 * analysis; a source that has no alleles to give gives a file without that
 * column, and one that does not keep whether its variants passed, a vars
 * file written before that column existed, gives one without it.
 *
 * Without `onBytes`, the whole file is built in the memory of wasm, which
 * grows and never shrinks, so what the tab holds while the call runs is the
 * source and the file together, and `numVarsPerBlock` is what an application that runs out
 * of memory lowers: it is the size of the block that is read and written at
 * a time. What comes back is the caller's own array, in the heap of
 * JavaScript, which nothing has to free and which the memory of wasm does
 * not hold a second copy of.
 *
 * What it gives back are those bytes and the counts of the pass it made:
 * how many variants were written and what each filter was given and kept.
 *
 * With `onBytes`, the file goes to that function in pieces of 1 MiB while
 * the pass writes it, and the memory of wasm holds one block of it and not
 * the whole file; the call then gives back the counts alone. Where the
 * pieces go is the application's: a worker can write each one to the
 * private file system of the browser with a `FileSystemSyncAccessHandle`,
 * or keep them and make a `Blob` of them. When `onBytes` throws, the pass
 * ends and the call throws what it threw, and one that returns a promise
 * ends it with an `Error`; when the pass fails after some pieces were given, the call throws
 * its error, and the pieces given are the start of a file that has no end.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `options` is not an object, when `numVarsPerBlock` is not a whole number of 1 or more and at most
 * 4294967295, when `onBytes` is not a function, when the source cannot be read, a wrong line of a VCF among
 * the causes, when the memory of the tab does not take the file, when
 * `onBytes` returns a promise, and when `init` has not been awaited; and
 * it throws what `onBytes` threw.
 */
export function writeVars(
  variants: Variants,
  options: WriteVarsOptions & { onBytes: OnBytes },
): WrittenInPieces;
export function writeVars(
  variants: Variants,
  options?: WriteVarsOptions & { onBytes?: undefined },
): VarsWritten;
export function writeVars(
  variants: Variants,
  options?: WriteVarsOptions,
): VarsWritten | WrittenInPieces;
export function writeVars(
  variants: Variants,
  options: WriteVarsOptions = {},
): VarsWritten | WrittenInPieces {
  theWasmHasToBeLoaded();
  anObjectOfOptions("writeVars", options, ["numVarsPerBlock", "onBytes"]);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const numVarsPerBlock =
    options.numVarsPerBlock === undefined
      ? undefined
      : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
  const onBytes = theFunctionOfThePieces(options.onBytes);
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
    source.write_vars(numVarsPerBlock, steps.of_a_pass(), onBytes),
  );
  return theFileOrTheCountsOf(file, onBytes, "vars file");
}

/**
 * `onBytes` as the binding crate takes it, and nothing when the application
 * gave none.
 *
 * @throws {Error} When `onBytes` is given and is not a function.
 */
export function theFunctionOfThePieces(onBytes: unknown): OnBytes | undefined {
  if (onBytes === undefined) {
    return undefined;
  }
  if (typeof onBytes !== "function") {
    throw new Error(
      "popnei: `onBytes` is a function that is given each piece of the file, " +
        `and ${whatWasGiven(onBytes)} was given`,
    );
  }
  return onBytes as OnBytes;
}

/**
 * What a writer gives back: the counts of the pass alone when the file went
 * to `onBytes`, and otherwise the bytes of the file with them.
 */
export function theFileOrTheCountsOf(
  file: WrittenFile,
  onBytes: OnBytes | undefined,
  what: string,
): VarsWritten | WrittenInPieces {
  if (onBytes === undefined) {
    return theBytesAndTheCountsOf(file, what);
  }
  try {
    return { passStats: passStatsOf(file.pass_stats()) };
  } finally {
    file.free();
  }
}

/**
 * The bytes of a file the core wrote in the memory of wasm, put together
 * into one array of the caller's, and the counts of the pass that wrote it;
 * the file is freed in wasm whatever happens. `what` names the file in the
 * error of one that gave fewer bytes than it said it holds.
 *
 * @throws {Error} When the pieces hold another number of bytes than the
 * file says, which is a defect of popnei.
 */
export function theBytesAndTheCountsOf(
  file: WrittenFile,
  what: string,
): { bytes: Uint8Array; passStats: PassStats } {
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
        `popnei: the ${what} says it holds ${bytes.length} bytes and gave ${written}`,
      );
    }
    return { bytes, passStats: passStatsOf(file.pass_stats()) };
  } finally {
    file.free();
  }
}
