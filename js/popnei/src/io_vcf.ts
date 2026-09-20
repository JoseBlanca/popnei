/** Reading a VCF. */

import {
  default_only_passed as defaultOnlyPassed,
  default_ploidy as defaultPloidy,
  open_vcf as openVcfOfTheCore,
} from "../wasm/popnei.js";

import { theWasmHasToBeLoaded } from "./core.js";
import { Variants } from "./variant.js";

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
 * The variants of the VCF in `source`, the bytes of the file, plain or
 * gzipped.
 *
 * It reads the header, so bytes that are not a VCF, or whose header popnei
 * cannot read, are an `Error` here and not at the first calculation. The
 * variants themselves are read again at every pass over what it returns.
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
 * @throws {Error} When the bytes are not a VCF popnei can read, when the
 * ploidy is 0 or above 255, and when `init` has not been awaited.
 */
export function openVcf(
  source: Uint8Array,
  options: OpenVcfOptions = {},
): Variants {
  theWasmHasToBeLoaded();
  const ploidy = options.ploidy ?? defaultPloidy();
  const onlyPassed = options.onlyPassed ?? defaultOnlyPassed();
  return new Variants(openVcfOfTheCore(source, ploidy, onlyPassed));
}
