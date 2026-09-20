/**
 * The block: a run of consecutive variants held as arrays.
 *
 * A block is how the genotypes leave popnei. The calculations consume blocks
 * inside the Rust core, and a user who wants the genotypes for an analysis
 * of their own asks a `Variants` for them with `iterBlocks`, which gives one
 * block after another until the source is at its end.
 */

import type { BlockColumns } from "../wasm/popnei.js";

/**
 * The variants of one block, each field a column of the block.
 *
 * A field other than the genotypes is there only when `iterBlocks` was asked
 * for it, and `null` when it was not.
 */
export interface Block {
  /**
   * The genotypes, `numVars` x the individuals of the source x the ploidy
   * alleles, variant after variant and inside a variant individual after
   * individual: the alleles of the individual `i` of the variant `v` of a
   * source of `n` individuals and ploidy `p` are the `p` numbers that start
   * at `(v * n + i) * p`.
   *
   * 0 is the reference allele and 1 and above the alternative ones, in the
   * order in which the source declares them, and -1 an allele that was not
   * called. It is a copy of what the core filled, which nothing writes into
   * again.
   */
  gts: Int8Array;

  /** How many variants the block holds. */
  numVars: number;

  /** The name of the chromosome of each variant. */
  chrom: string[] | null;

  /**
   * The position of each variant, 1 based as in a VCF.
   *
   * The positions are float64, which holds a position of up to 2^53
   * exactly: an array of unsigned 64 bit numbers would give a `BigInt` for
   * each, which does not mix with the ordinary numbers of JavaScript in
   * arithmetic.
   */
  pos: Float64Array | null;

  /** The id of each variant, `null` for a variant that has none. */
  id: (string | null)[] | null;

  /**
   * The alleles of each variant, the reference one first.
   *
   * Each one is the text the source gave, so a symbolic allele, `<DEL>`, and
   * the allele of an overlapping deletion, `*`, are alleles like any other.
   */
  alleles: string[][] | null;

  /**
   * The quality of each variant, phred scaled as the QUAL of a VCF: 30 is
   * one chance in a thousand that there is no variant at that site. It is
   * NaN for a variant whose source gives no quality.
   */
  qual: Float32Array | null;
}

/**
 * The block of the columns of one block of the core.
 *
 * Every column is read once, which is what takes it out of the memory of
 * wasm, and the caller frees `columns` as soon as this returns.
 */
export function blockOf(columns: BlockColumns): Block {
  const alleles = columns.alleles();
  const numAllelesPerVar = columns.num_alleles_per_var();
  return {
    gts: columns.gts(),
    numVars: columns.num_vars(),
    chrom: columns.chrom() ?? null,
    pos: columns.pos() ?? null,
    // The core gives an empty text for a variant with no id, and a user of
    // the package reads `null`, as a Python user reads `None`.
    id: columns.id()?.map((id) => (id === "" ? null : id)) ?? null,
    alleles:
      alleles === undefined || numAllelesPerVar === undefined
        ? null
        : allelesOfEachVariant(alleles, numAllelesPerVar),
    qual: columns.qual() ?? null,
  };
}

/**
 * The alleles of a block cut into the alleles of each variant.
 *
 * The core holds them in one buffer, the alleles of one variant after those
 * of the variant before it, and that is how they cross: an array of arrays
 * of strings is not one of the types wasm-bindgen carries.
 */
function allelesOfEachVariant(
  alleles: string[],
  numAllelesPerVar: Uint32Array,
): string[][] {
  const ofEachVariant: string[][] = [];
  let first = 0;
  for (const numAlleles of numAllelesPerVar) {
    ofEachVariant.push(alleles.slice(first, first + numAlleles));
    first += numAlleles;
  }
  return ofEachVariant;
}
