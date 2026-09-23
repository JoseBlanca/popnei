/**
 * How strongly two variants go together: the matrix of r² of every pair of
 * the variants of a dataset.
 *
 * Two variants are in linkage disequilibrium when the genotype of one tells
 * something about the genotype of the other, which happens when they sit
 * close enough on a chromosome that few recombinations have separated them.
 * r² measures it: each variant becomes one number per individual, its
 * dosage, how many alleles of the genotype are not the major allele of that
 * variant, and r² is the square of the correlation between the two dosage
 * vectors, 0 when one variant says nothing about the other and 1 when it
 * says everything. `docs/specs/ld.md` has what is computed.
 */

import { default_max_num_vars as defaultMaxNumVars } from "../wasm/popnei.js";

import { wholeNumberOfOneOrMore } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/** How the matrix of r² is calculated: the variants it is taken of. */
export interface CalcRogersHuffR2MatrixOptions {
  /**
   * How many variants the calculation takes before it refuses, a whole
   * number of 1 or more. 5000 when it is not given, which is the default of
   * the core crate and 200 MB of matrix.
   *
   * The matrix holds one r² for each pair of the variants of the pass, 8
   * bytes each, so it grows with the square of them: 5000 variants are 200
   * MB and 100000 are 80 GB. It is the one calculation of popnei whose
   * result grows with the square of its input, and a page holds 4 GB of
   * everything that is open in it at a time, so a pass of more variants
   * than this is an `Error` and not a matrix the tab is asked for the
   * memory of. A user who wants the matrix of more variants, and has the
   * memory, raises it; a user who has more variants than memory puts a
   * filter on the `Variants` first.
   */
  maxNumVars?: number;
}

/**
 * The r² of every pair of the variants of a dataset, with the chromosome
 * and the position of each of them and the counts of the pass.
 *
 * It is the `R2Matrix` of the Python package with the names of TypeScript:
 * a `Float64Array` where Python has a numpy array, and an array of strings
 * where Python has a tuple.
 */
export interface R2Matrix {
  /** How many variants the pass gave, which is the side of the matrix. */
  readonly numVars: number;
  /**
   * The r² of every pair, `numVars` x `numVars` row after row: the r² of
   * the variants `i` and `j` is the value at `i * numVars + j`, and the
   * value at `j * numVars + i` is the same one.
   *
   * A pair that has no r² is NaN, which is a pair holding a variant whose
   * called genotypes have one dosage only, and the diagonal of a variant
   * that has two dosages at least is 1. The individuals whose genotype is
   * missing at either variant of a pair are left out of that pair, so each
   * pair has its own number of individuals.
   */
  readonly r2: Float64Array;
  /**
   * The name of the chromosome of each variant, in the order the pass gave
   * them, which is the order of the rows and of the columns of the matrix.
   */
  readonly chroms: readonly string[];
  /**
   * The position of each variant, 1 based as in a VCF. The distance of a
   * pair is the difference of two of them, and a pair whose variants are on
   * two chromosomes has no distance and an r² like any other pair.
   */
  readonly poss: Float64Array;
  /**
   * How many variants the calculation took, after the steps, and how many
   * variants each filter of the pass was given and kept.
   */
  readonly passStats: PassStats;
}

/**
 * The r² of every pair of the variants of `variants`.
 *
 * Each variant becomes one number per individual, its dosage: how many
 * alleles of the genotype are not the major allele of that variant, which
 * is the most frequent among its called alleles and the lowest numbered of
 * two that are equally frequent. r² is the square of the correlation
 * between the dosages of the two variants of a pair, over the individuals
 * whose genotype is missing at neither of them, so each pair has its own
 * number of individuals; a variant whose called genotypes hold one dosage
 * only has no correlation with anything and is NaN in its whole row, in its
 * whole column and in its diagonal cell.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, so the matrix
 * is of the variants the filters kept, and the `variants` are as they were
 * afterwards.
 *
 * It is `calc_rogers_huff_r2_matrix` of the Python package, which it
 * mirrors, and pyNei's function of that name with two differences a caller
 * sees: pyNei gives r and this gives its square, which loses the sign and
 * nothing else, and pyNei leaves an individual with a missing genotype in
 * the pair with a dosage of -1 where this leaves it out, which is what
 * plink2 does and what makes plink2 the reference program of this
 * calculation.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `maxNumVars` is not a whole number of 1 or more, when the pass gives more
 * variants than `maxNumVars`, whose message has both numbers and the memory
 * the matrix would have needed, when the pass gives no variant, whose
 * message says whether the source held none or the steps kept none and how
 * many variants each filter was given and kept, when the source cannot be
 * read, a wrong line of a VCF among the causes, when a position of the
 * source is above 2^53, which a number of JavaScript rounds, when the
 * memory of the tab does not take the matrix, 8 bytes a pair, and when
 * `init` has not been awaited.
 */
export function calcRogersHuffR2Matrix(
  variants: Variants,
  options: CalcRogersHuffR2MatrixOptions = {},
): R2Matrix {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  // The default is the core's, as the ploidy of `openVcf` is, so that
  // Python and TypeScript cannot drift apart on how many variants the
  // matrix is taken of when the caller says nothing.
  const maxNumVars =
    options.maxNumVars === undefined
      ? defaultMaxNumVars()
      : wholeNumberOfOneOrMore("maxNumVars", options.maxNumVars);
  // The steps of the pass are a copy of the list, made after the arguments
  // were checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = source.calc_rogers_huff_r2_matrix(
    maxNumVars,
    steps.of_a_pass(),
  );
  try {
    // Each array is moved out of the result as it is read, and not cloned.
    // The generated code still copies the values into an array of the
    // JavaScript heap and frees the memory of wasm after it, so the matrix
    // is held twice while that copy is made, 400 MB at 5000 variants, and
    // once afterwards.
    return {
      numVars: calculated.num_vars(),
      r2: thePartOfTheMatrix(calculated.r2(), "r2"),
      chroms: thePartOfTheMatrix(calculated.chroms(), "chroms"),
      poss: thePartOfTheMatrix(calculated.poss(), "poss"),
      passStats: passStatsOf(calculated.pass_stats()),
    };
  } finally {
    calculated.free();
  }
}

/**
 * One array of the matrix on its way out of WebAssembly.
 *
 * Each of them leaves the memory of wasm the first time it is asked for, so
 * the call after that gives nothing; `calcRogersHuffR2Matrix` reads each of
 * them once, and the `Error` here is a defect of this package and not
 * something a caller can do. It is `theValuesOf` of `pca.ts` for the three
 * arrays of this result, one of which is an array of names.
 *
 * @throws {Error} When the array had been read already.
 */
function thePartOfTheMatrix<Part>(part: Part | undefined, name: string): Part {
  if (part === undefined) {
    throw new Error(
      `popnei: \`${name}\` was read twice out of the memory of WebAssembly, ` +
        "which is a defect of popnei; please report it",
    );
  }
  return part;
}
