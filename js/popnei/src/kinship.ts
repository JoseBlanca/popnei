/**
 * The kinship of every pair of the individuals of a dataset: how much more
 * of their genome the two share than two individuals drawn at random from
 * the same panel would.
 *
 * It is the matrix of VanRaden 2008, which plink2's `--make-rel` computes
 * and which GCTA is built on. An entry off the diagonal is twice the
 * coancestry of its pair, about 0.5 for full sibs or for a parent and a
 * child, about 0.25 for half sibs and near 0 for two individuals with no
 * recent ancestor in common; an entry on the diagonal is 1 plus the
 * inbreeding of that individual. Entries below 0 are ordinary and mean a
 * pair less alike than the average pair of the panel, because the whole
 * matrix is measured against that average.
 *
 * `calcKinship` takes it from the variants of a dataset, and `Kinship` is
 * what it gives and also what a user builds from the matrix another program
 * wrote, so that it can be passed on. `docs/specs/kinship.md` has what the
 * calculation does with a missing genotype and with a variant that has no
 * variance.
 */

import { default_transform_to_biallelic as defaultTransformToBiallelic } from "../wasm/popnei.js";

import { aBoolean, namesOf } from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * How far a matrix a user built may be from its own transpose before it is
 * refused, as a part of its largest absolute entry.
 *
 * A kinship is symmetric, and a matrix that is not says that its rows and
 * its columns are not the same individuals in the same order, which every
 * entry of the result would then be wrong about. The tolerance costs a
 * matrix that came from a program nothing: both matrices plink2 wrote for
 * the reference panels of `docs/specs/kinship.md` are symmetric to the bit,
 * a largest `|m - m'|` of 0.
 */
const FROM_ITS_OWN_TRANSPOSE = 1e-9;

/**
 * The kinship of every pair of a set of individuals, with the names of those
 * individuals and how many variants it was built from.
 *
 * `matrix` holds the individuals x individuals values row after row: the
 * entry of the individuals `i` and `j` is the value at `i * N + j`, and the
 * value at `j * N + i` is the same one. `individuals` names the rows and the
 * columns, in that order.
 *
 * It is the `Kinship` of the Python package with the names of TypeScript: a
 * `Float64Array` and an array of names where Python has a pandas frame
 * indexed by those names. A user who has the matrix of plink2 or the one a
 * pedigree gave them builds one with the constructor and passes it on.
 */
export class Kinship {
  /** The individuals x individuals values, row after row. */
  readonly matrix: Float64Array;

  /** The names of the individuals, the rows and the columns of the matrix. */
  readonly individuals: readonly string[];

  /**
   * How many variants the matrix was built from: those that had variance
   * among these individuals. A variant whose called genotypes all have one
   * dosage is in no sum and in no denominator.
   */
  readonly numVars: number;

  /**
   * How many variants the pass gave, used or not, and what each filter of
   * the `Variants` was given and kept, or `undefined` for a kinship a user
   * built by hand, since no pass produced it.
   */
  readonly passStats: PassStats | undefined;

  /**
   * The kinship of the pairs of `individuals`, which is what `calcKinship`
   * builds and what a user builds from a matrix of their own.
   *
   * @throws {Error} When `matrix` is not a `Float64Array` of one value for
   * each pair of `individuals`, which is `individuals.length` squared of
   * them, when `individuals` is not an array of names or names one twice,
   * and when the matrix is further from its own transpose than 1e-9 of its
   * largest absolute entry, which says that its rows and its columns are
   * not the same individuals in the same order.
   */
  constructor(
    matrix: Float64Array,
    individuals: readonly string[],
    numVars: number,
    passStats?: PassStats,
  ) {
    const names = namesOf("individuals", individuals, {
      oneOfThem: "individual",
      anExample: "ind00",
    });
    const numIndividuals = names.length;
    if (!(matrix instanceof Float64Array)) {
      throw new Error(
        "popnei: the matrix of a kinship is a Float64Array of one value for " +
          "each pair of its individuals",
      );
    }
    if (matrix.length !== numIndividuals * numIndividuals) {
      throw new Error(
        `popnei: the kinship of ${numIndividuals} individuals is ` +
          `${numIndividuals * numIndividuals} values, and ${matrix.length} ` +
          "were given",
      );
    }
    const named = new Set(names);
    if (named.size !== numIndividuals) {
      throw new Error(
        "popnei: the individuals of a kinship are named once each, and " +
          `${numIndividuals - named.size} of the names given are there twice`,
      );
    }
    theMatrixIsSymmetric(matrix, numIndividuals);
    this.matrix = matrix;
    this.individuals = Object.freeze([...names]);
    this.numVars = numVars;
    this.passStats = passStats;
  }

  /**
   * The rows and the columns of `individuals`, in the order they are named
   * here, with `numVars` and `passStats` as they are.
   *
   * It takes a part of this matrix out and calculates nothing: the
   * frequencies, the means and the denominators of every entry are those of
   * the individuals the matrix was built from. A kinship of some individuals
   * alone, with their own frequencies, is `calcKinship` with its
   * `individuals`, and the two differ: on the reference panel of
   * `docs/specs/kinship.md` by up to 0.129 over the same 40 individuals.
   *
   * It is pyNei's `Kinship.filter_samples`.
   *
   * @throws {Error} When `individuals` is not an array of names, when a name
   * is of nobody in the matrix, and when a name is there twice.
   */
  filterIndividuals(individuals: readonly string[]): Kinship {
    const names = namesOf("individuals", individuals, {
      oneOfThem: "individual",
      anExample: "ind00",
    });
    const at = new Map(this.individuals.map((name, row) => [name, row]));
    const rows = names.map((name) => {
      const row = at.get(name);
      if (row === undefined) {
        throw new Error(
          `popnei: \`${name}\` is not an individual of this kinship, which ` +
            `is of ${this.individuals.length} of them`,
        );
      }
      return row;
    });
    const numIndividuals = this.individuals.length;
    const kept = new Float64Array(rows.length * rows.length);
    for (const [row, rowOfTheMatrix] of rows.entries()) {
      for (const [column, columnOfTheMatrix] of rows.entries()) {
        kept[row * rows.length + column] = this.matrix[
          rowOfTheMatrix * numIndividuals + columnOfTheMatrix
        ] as number;
      }
    }
    return new Kinship(kept, names, this.numVars, this.passStats);
  }
}

/**
 * That `matrix`, `numIndividuals` x `numIndividuals` row after row, is as
 * far from its own transpose as a kinship is allowed to be.
 *
 * @throws {Error} When two entries that are the same pair differ by more
 * than 1e-9 of the largest absolute entry of the matrix, naming the pair and
 * the two values.
 */
function theMatrixIsSymmetric(
  matrix: Float64Array,
  numIndividuals: number,
): void {
  let largest = 0;
  for (const entry of matrix) {
    const size = Math.abs(entry);
    if (size > largest) {
      largest = size;
    }
  }
  const allowed = largest * FROM_ITS_OWN_TRANSPOSE;
  for (let row = 0; row < numIndividuals; row += 1) {
    for (let column = row + 1; column < numIndividuals; column += 1) {
      const entry = matrix[row * numIndividuals + column] as number;
      const mirrored = matrix[column * numIndividuals + row] as number;
      if (!(Math.abs(entry - mirrored) <= allowed)) {
        throw new Error(
          `popnei: a kinship is symmetric, and the pair of the individuals ` +
            `${row} and ${column} is ${entry} in one half of the matrix and ` +
            `${mirrored} in the other`,
        );
      }
    }
  }
}

/** How the kinship of the variants of a dataset is taken. */
export interface CalcKinshipOptions {
  /**
   * The names of the individuals the matrix is of, in the order it has
   * them, and every individual the pass gives, in its order, when it is not
   * given.
   *
   * Every frequency, mean and denominator is of those individuals: the
   * kinship of some of them is not the rows and the columns of the kinship
   * of the whole panel, which is what `filterIndividuals` on the result
   * gives. On the reference panel of `docs/specs/kinship.md` the two differ
   * by up to 0.129 over the same 40 individuals, and the first uses 1195
   * variants where the whole panel uses 1200, the other 5 having no variance
   * among those 40.
   */
  individuals?: readonly string[];
  /**
   * Whether every allele that is not the major one counts the same, which
   * is what gives a variant of more than two alleles a dosage. False when
   * it is not given, and such a variant is then an `Error`.
   */
  transformToBiallelic?: boolean;
}

/**
 * The kinship of every pair of the individuals of `variants`, after its
 * steps.
 *
 * Each variant becomes one number per individual, its dosage: how many
 * alleles of the genotype are not the major allele of that variant, which is
 * the most frequent among its called alleles. The dosages of a variant are
 * then centered and divided by the standard deviation its allele frequency
 * gives it under Hardy Weinberg, `sqrt(ploidy * p * (1 - p))`, and the entry
 * of a pair is the sum over the variants of the two standardized dosages
 * multiplied, divided by how many of those variants have a called genotype
 * in both individuals. That divisor is not the one `doPcaFromVariants`
 * takes, the standard deviation of the dosages themselves, and the two agree
 * only when the genotypes are in Hardy Weinberg proportions; it is the one
 * that makes an entry twice a coancestry, and it is what plink2 and GCTA
 * use.
 *
 * A genotype with any allele missing takes the mean of the dosages of its
 * variant, so once the variant is centered it pulls its pair nowhere, and it
 * counts in the denominator of no pair it is in. A variant whose called
 * genotypes all have one dosage has no variance and is left out, a variant
 * with one allele and one where every individual is heterozygous among them;
 * `numVars` is how many were used.
 *
 * The call is a consumer of the `variants`: it makes one pass over the
 * source through the steps that are on it when it is called, and nothing of
 * the size of the variants x the individuals is held. Pruning by linkage
 * disequilibrium, `filterByLd`, is the step a user normally puts on the
 * `Variants` first.
 *
 * It is pyNei's `calc_kinship`, whose `samples` is `individuals` here.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed, when
 * `individuals` is not an array of names, when a name of it is of nobody the
 * pass gives, is there twice, or the list is empty, when
 * `transformToBiallelic` is not a boolean, when the source cannot be read, a
 * wrong line of a VCF among the causes, when a variant has more than two
 * alleles among its called genotypes and `transformToBiallelic` is false,
 * when the pass gives no variant or no variant with variance among these
 * individuals, when two individuals have no variant called in both, whose
 * entry would be divided by 0, when the ploidy is above 254 or the
 * individuals are more than 46340, and when `init` has not been awaited.
 */
export function calcKinship(
  variants: Variants,
  options: CalcKinshipOptions = {},
): Kinship {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  const individuals =
    options.individuals === undefined
      ? undefined
      : namesOf("individuals", options.individuals, {
          oneOfThem: "individual",
          anExample: "ind00",
        });
  const transformToBiallelic =
    options.transformToBiallelic === undefined
      ? defaultTransformToBiallelic()
      : aBoolean("transformToBiallelic", options.transformToBiallelic);
  // The steps of the pass are a copy of the list, made after the arguments
  // were checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = source.calc_kinship(
    individuals,
    transformToBiallelic,
    steps.of_a_pass(),
  );
  try {
    // The matrix and the names are moved out of the result as they are
    // read, and not cloned: the generated code copies the values into an
    // array of the JavaScript heap and frees the memory of wasm after it,
    // so the matrix is held twice while that copy is made, 800 MB at 10000
    // individuals, and once afterwards.
    const matrix = calculated.matrix();
    const names = calculated.individuals();
    if (matrix === undefined || names === undefined) {
      throw new Error(
        "popnei: the kinship of this pass gave no matrix or no names, " +
          "which is a defect of popnei; please report it",
      );
    }
    return new Kinship(
      matrix,
      names,
      calculated.num_vars(),
      passStatsOf(calculated.pass_stats()),
    );
  } finally {
    calculated.free();
  }
}
