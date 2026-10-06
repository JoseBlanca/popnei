/** Reading and writing a VCF. */

import {
  default_only_passed as defaultOnlyPassed,
  open_vcf as openVcfOfTheCore,
  open_vcf_of_a_file as openVcfOfAFileOfTheCore,
} from "../wasm/popnei.js";

import {
  aBoolean,
  anObjectOfOptions,
  bytesOrFile as bytesOrFileOf,
  wholeNumberOfOneOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { theBytesAndTheCountsOf } from "./io_vars.js";
import type { PassStats } from "./variant.js";
import { Variants, sourceOfTheVariants } from "./variant.js";

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
   * individual and every variant. When it is not given it is read from the
   * file: the number of alleles of its first genotype that is not a single
   * dot.
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
 * not defined for one. When it is not given, `openVcf` reads it from the
 * file before it returns: it is the number of alleles of the first genotype,
 * in the order of the lines and of the individuals, that is not a single
 * dot, which is a missing genotype of any ploidy. A file whose first 4096
 * data lines hold no such genotype is then an `Error` here, which says to
 * give the ploidy, and so is a file with a header and no data line, which
 * says that the file has no variants and its ploidy cannot be inferred. That
 * search reads the file from its start, so a file damaged before its first
 * genotype with alleles, a corrupted member of a bgzipped file among the
 * causes, is an `Error` here as well, and not at the first pass.
 *
 * `onlyPassed` leaves out the variants that failed a filter, those whose
 * FILTER column is neither `PASS` nor a dot; a dot says that no filter was
 * applied. With it false every variant of the file is given, and nothing
 * then says which ones had failed. Only the variants that passed are given
 * when it is not.
 *
 * What it returns holds memory of wasm until its `free()` is called.
 *
 * It is pyNei's `vars_from_vcf` under another name, with the ploidy and the
 * filter as arguments, which pyNei has not; pyNei takes the ploidy from the
 * genotype of the first individual of the first data line, a single dot
 * among them, and gives every variant.
 *
 * @throws {Error} When `source` is neither a `Uint8Array` nor a `File` or
 * `Blob`, when a `File` or a `Blob` is given where there is no
 * `FileReaderSync`, which is everywhere but a web worker, when `ploidy` is
 * not a whole number of 1 or more, when `onlyPassed` is not a boolean, when
 * the bytes are not a VCF popnei can read, when the ploidy is above 255,
 * when no ploidy is given and none can be read from the file, and when
 * `init` has not been awaited.
 */
export function openVcf(
  source: BytesOrFile,
  options: OpenVcfOptions = {},
): Variants {
  theWasmHasToBeLoaded();
  anObjectOfOptions("openVcf", options, ["ploidy", "onlyPassed"]);
  // Left out, it is read from the file by the core.
  const ploidy =
    options.ploidy === undefined
      ? undefined
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

/** The VCF `writeVcf` wrote, and the counts of the pass it made. */
export interface VcfWritten {
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

/** How a VCF is written: bgzipped or as plain text. */
export interface WriteVcfOptions {
  /**
   * Whether the file is compressed with bgzip, which tabix indexes and
   * bcftools asks a region of. True when it is not given: there is no path
   * to read the compression from.
   */
  bgzip?: boolean;
}

/**
 * Every variant of `variants`, after its steps, as the bytes of a VCF,
 * which a page offers as a download: a tab has no filesystem.
 *
 * It is how the variants popnei kept, filtered by missing data, by
 * individual or by any other step, reach plink2, bcftools or a program of
 * the user's own. When the source is a VCF, each variant is written as its
 * line was, every column of it, INFO, FILTER, the phase and the values of
 * each individual other than GT among them, and the header is the source's
 * with a `#CHROM` line of the individuals that were kept, in the order the
 * filter of individuals named them. When that filter took individuals out,
 * AC and AN, counts over individuals that are no longer in the file, are
 * taken out of INFO and their `##INFO` lines out of the header, as
 * `bcftools annotate -x INFO/AC,INFO/AN` does; a filter that keeps every
 * individual, in any order, leaves them.
 *
 * When the source is a vars file, the lines hold what the file holds, with
 * INFO a dot, FORMAT `GT` and the alleles of each genotype joined by `/`,
 * since a vars file keeps no phase, under a header with one `##contig` line
 * for each chromosome whose length the vars file keeps. FILTER is a dot, or
 * `FAIL` for a variant that failed its FILTER in the VCF the vars file was
 * written from, and then the header has a `##FILTER` line for `FAIL`.
 *
 * The lines are written in the order the source gives them. A VCF opened
 * with `onlyPassed` false and written with no step and `{bgzip: false}` is
 * the same bytes, when its lines end in `\n` and none is empty. The call
 * reads the source once, and the whole file is built in the memory of wasm
 * before it crosses, in pieces, into the array that is returned.
 *
 * pyNei has no VCF writer, so nothing is mirrored; the name follows
 * `writeVars`.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `options` is not an object, when `bgzip` is not a boolean, when the source cannot be read, a wrong line of
 * a VCF among the causes, when the memory of the tab does not take the
 * file, and when `init` has not been awaited.
 */
export function writeVcf(
  variants: Variants,
  options: WriteVcfOptions = {},
): VcfWritten {
  theWasmHasToBeLoaded();
  anObjectOfOptions("writeVcf", options, ["bgzip"]);
  const bgzip =
    options.bgzip === undefined ? true : aBoolean("bgzip", options.bgzip);
  const { source, steps, whileTheRunReads } = sourceOfTheVariants(
    "variants",
    variants,
  );
  const file = whileTheRunReads(() =>
    source.write_vcf(bgzip, steps.of_a_pass()),
  );
  return theBytesAndTheCountsOf(file, "VCF");
}
