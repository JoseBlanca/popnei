/** Reading and writing a vars file, the file popnei keeps its variants in. */

import { open_vars as openVarsOfTheCore } from "../wasm/popnei.js";

import { bytes as bytesOf, wholeNumberOfOneOrMore } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Variants, sourceOfTheVariants } from "./variant.js";

/** How a vars file is written: how many variants a batch of it holds. */
export interface WriteVarsOptions {
  /**
   * How many variants one batch of the file holds, the last one aside, a
   * whole number of 1 or more. When it is not given, the size popnei
   * chooses for the number of individuals of the source, which is the size
   * of its blocks.
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
 * It reads the schema of the file and its footer, so bytes that are not a
 * vars file, one of a format version popnei does not read, one whose columns
 * are not those of a vars file and one whose individuals name nobody are an
 * `Error` here and not at the first calculation. What is in the batches is
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
 * read with `openVars` is written again with another size of batch.
 *
 * The call reads the whole source once. The file holds the six columns of a
 * VCF, the chromosome, the position, the id, the alleles, the quality and
 * the genotypes, whether or not the user will read them, so that it can
 * stand in for the VCF in any later analysis; a source that has no alleles
 * to give gives a file without that column.
 *
 * The whole file is built in the memory of wasm and the `Uint8Array` is a
 * copy of it, so what the tab holds while the call runs is the source and
 * the file together. That memory grows and never shrinks.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `numVarsPerBlock` is not a whole number of 1 or more, when the source
 * cannot be read, a wrong line of a VCF among the causes, and when `init`
 * has not been awaited.
 */
export function writeVars(
  variants: Variants,
  options: WriteVarsOptions = {},
): Uint8Array {
  theWasmHasToBeLoaded();
  const source = sourceOfTheVariants("variants", variants);
  const numVarsPerBlock =
    options.numVarsPerBlock === undefined
      ? undefined
      : wholeNumberOfOneOrMore("numVarsPerBlock", options.numVarsPerBlock);
  return source.write_vars(numVarsPerBlock);
}
