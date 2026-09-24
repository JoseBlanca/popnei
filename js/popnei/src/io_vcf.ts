/** Reading a VCF. */

import {
  default_only_passed as defaultOnlyPassed,
  default_ploidy as defaultPloidy,
  open_vcf as openVcfOfTheCore,
  open_vcf_of_a_file as openVcfOfAFileOfTheCore,
} from "../wasm/popnei.js";

import {
  aBoolean,
  bytesOrFile as bytesOrFileOf,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Variants } from "./variant.js";

/**
 * What a source of variants is opened over: the bytes of the file, or the
 * file itself as the page holds it.
 *
 * A `File`, the handle a page gets when the user picks a file in a form or
 * drops one on it, is a `Blob`, which is a piece of bytes of a page, so a
 * `Blob` an application built itself is read the same way.
 */
export type BytesOrFile = Uint8Array | Blob;

/** How a VCF is read: the ploidy of its genotypes and which variants. */
export interface OpenVcfOptions {
  /**
   * How many alleles every genotype of the file holds, the same for every
   * individual and every variant. 2 when it is not given.
   */
  ploidy?: number;
  /**
   * Whether the variants that failed a filter are left out. True when it is
   * not given.
   */
  onlyPassed?: boolean;
}

/**
 * The variants of the VCF in `source`, the bytes of the file or the file
 * the user picked in the page, plain or gzipped.
 *
 * It reads the header, so bytes that are not a VCF, or whose header popnei
 * cannot read, are an `Error` here and not at the first calculation. The
 * variants themselves are read again at every pass over what it returns.
 *
 * With a `File` or a `Blob`, popnei reads the ranges of bytes it needs
 * through `FileReaderSync`, a few MiB at a time, so the file is never in
 * the memory of wasm whole and a user opens a file larger than the memory
 * of the tab. A browser has `FileReaderSync` only inside a web worker, so
 * `openVcf` of a `File` or a `Blob` on the main thread of a page, or under
 * node, is an `Error` here, at the call; a `Uint8Array` is read wherever it
 * is given, as it is today. What is held for such a source until its
 * `free()` is called is the handle of the file, its name and its size, and
 * not its bytes, and every pass reads the file again.
 *
 * A range that comes back shorter than the one popnei asked for, inside a
 * file of that size, is an `Error` in the middle of a pass and not the end
 * of the file: a browser gives a short range when the file changed on disk
 * after the page got its handle.
 *
 * `ploidy` is how many alleles every genotype of the file holds, and a
 * genotype of any other number of alleles is an `Error` when it is read:
 * popnei does not read a VCF of mixed ploidies, because its calculations are
 * not defined for one. `onlyPassed` leaves out the variants that failed a
 * filter, those whose FILTER column is neither `PASS` nor a dot; a dot says
 * that no filter was applied. With it false every variant of the file is
 * given, and nothing then says which ones had failed. The two defaults are
 * the ones of the core crate, a ploidy of 2 and only the variants that
 * passed.
 *
 * What it returns holds memory of wasm until its `free()` is called.
 *
 * It is pyNei's `vars_from_vcf` under another name, with the ploidy and the
 * filter as arguments, which pyNei has not; pyNei takes the ploidy from the
 * first genotype of the file and gives every variant.
 *
 * @throws {Error} When `source` is neither a `Uint8Array` nor a `File` or
 * `Blob`, when a `File` or a `Blob` is given where there is no
 * `FileReaderSync`, which is everywhere but a web worker, when `ploidy` is
 * not a whole number of 1 or more, when `onlyPassed` is not a boolean, when
 * the bytes are not a VCF popnei can read, when the ploidy is above 255,
 * and when `init` has not been awaited.
 */
export function openVcf(
  source: BytesOrFile,
  options: OpenVcfOptions = {},
): Variants {
  theWasmHasToBeLoaded();
  const ploidy =
    options.ploidy === undefined
      ? defaultPloidy()
      : wholeNumberOfOneOrMore("ploidy", options.ploidy);
  const onlyPassed =
    options.onlyPassed === undefined
      ? defaultOnlyPassed()
      : aBoolean("onlyPassed", options.onlyPassed);
  const file = bytesOrFileOf("source", source);
  return new Variants(
    file instanceof Uint8Array
      ? openVcfOfTheCore(file, ploidy, onlyPassed)
      : openVcfOfAFileOfTheCore(file, ploidy, onlyPassed),
  );
}
