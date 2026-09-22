/**
 * The files the tests read: the reference files of the repository, and the
 * small VCFs a test writes for a case that no reference file has.
 *
 * The reference VCFs are those of `docs/specs/io_vcf.md`, in
 * `tests/reference/vcf/` at the root of the repository, which
 * `tests/reference/vcf/make_reference.py` writes and bcftools was run on.
 * The reference vars files are those of `docs/specs/io_vars.md`, in
 * `tests/reference/vars/`, which `tests/reference/vars/make_reference.py`
 * writes with pyarrow because popnei cannot write them. The tables of
 * `docs/specs/pca.md` are in `tests/reference/pca/`, which
 * `tests/reference/pca/make_reference.py` writes. The Python tests read the
 * same files.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`.
 */

import { readFile } from "node:fs/promises";

const REFERENCE_VCF_DIR = new URL(
  "../../../tests/reference/vcf/",
  import.meta.url,
);

const REFERENCE_VARS_DIR = new URL(
  "../../../tests/reference/vars/",
  import.meta.url,
);

const REFERENCE_STATS_DIR = new URL(
  "../../../tests/reference/stats/",
  import.meta.url,
);

const REFERENCE_DISTS_DIR = new URL(
  "../../../tests/reference/dists/",
  import.meta.url,
);

/** The bytes of the reference VCF `name`, `cases.vcf` or `many.vcf.gz`. */
export async function referenceVcf(name: string): Promise<Uint8Array> {
  return new Uint8Array(await readFile(new URL(name, REFERENCE_VCF_DIR)));
}

/** The bytes of the reference vars file `name`, `zstd.vars`. */
export async function referenceVars(name: string): Promise<Uint8Array> {
  return new Uint8Array(await readFile(new URL(name, REFERENCE_VARS_DIR)));
}

/**
 * The bytes of the reference file `name` of the stats module, the panel
 * `panel.vcf.gz` and the `panel_pops_bcftools.txt` of its populations.
 *
 * `tests/reference/stats/make_reference.py` writes both, and the same
 * script runs the plink2 and bcftools commands of `docs/specs/stats.md` and
 * keeps their reports beside them.
 */
export async function referenceStats(name: string): Promise<Uint8Array> {
  return new Uint8Array(await readFile(new URL(name, REFERENCE_STATS_DIR)));
}

/**
 * The bytes of the reference file `name` of the Kosman distances,
 * `panel.vcf.gz` or `panel.gdkosman.tsv`.
 *
 * They are the files of "How it is verified" of `docs/specs/dists.md`, in
 * `tests/reference/dists/`, which `tests/reference/dists/make_reference.py`
 * writes: the four datasets as gzipped VCFs and, beside each, the distance
 * and the number of variants that `gd.kosman` of the R package
 * PopGenReport 3.1.3 gives for every pair. The Python tests read the same
 * files.
 */
export async function referenceDists(name: string): Promise<Uint8Array> {
  return new Uint8Array(await readFile(new URL(name, REFERENCE_DISTS_DIR)));
}

/**
 * The bytes of a VCF of three diploid individuals with `dataLines` as its
 * data lines, for a case that no reference file has.
 *
 * Its header is the one of `cases.vcf`, without the `##` lines that no
 * reader of popnei looks at.
 */
/**
 * The bytes of a VCF of three diploid individuals with `numVars` variants
 * of one alternative allele, for a test that needs a file of some megabytes
 * and does not care what is in it.
 */
export function manyVariantsVcf(numVars: number): Uint8Array {
  const lines = Array.from(
    { length: numVars },
    (_unused, variant) =>
      `chr1\t${variant + 1}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1`,
  );
  return vcfOf(lines);
}

export function vcfOf(dataLines: readonly string[]): Uint8Array {
  const header = [
    "##fileformat=VCFv4.4",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3",
  ];
  return new TextEncoder().encode([...header, ...dataLines, ""].join("\n"));
}

const REFERENCE_PCA_DIR = new URL(
  "../../../tests/reference/pca/",
  import.meta.url,
);

/**
 * The bytes of the VCF `name` of `tests/reference/pca/`, `worked.vcf`,
 * which `tests/reference/pca/make_reference.py` writes.
 */
export async function referencePcaVcf(name: string): Promise<Uint8Array> {
  return new Uint8Array(await readFile(new URL(name, REFERENCE_PCA_DIR)));
}

/** A table of numbers: its values row after row, and its two sides. */
export interface Table {
  values: Float64Array;
  numRows: number;
  numCols: number;
}

/**
 * The table `name` of `tests/reference/pca/`, `iris.tsv`, which
 * `tests/reference/pca/make_reference.py` writes with pandas.
 *
 * The first line of such a file names the traits and the first field of
 * every other line names the row, and neither is a value of the table, so
 * both are left out here, as the tests of the core crate leave them out.
 *
 * @throws {Error} When the file has no rows, when its rows do not all hold
 * the same number of fields, or when a field is not a finite number: a
 * table read wrong would be the analysis of something else.
 */
export async function referenceTable(name: string): Promise<Table> {
  const text = await readFile(new URL(name, REFERENCE_PCA_DIR), "utf8");
  const lines = text.split("\n").filter((line) => line.trim() !== "");
  const rows = lines.slice(1).map((line, row) =>
    line
      .split("\t")
      .slice(1)
      .map((field) => {
        const value = Number(field.trim());
        if (!Number.isFinite(value)) {
          throw new Error(
            `${name}: the row ${row} holds the field \`${field}\``,
          );
        }
        return value;
      }),
  );
  const numCols = rows[0]?.length ?? 0;
  if (rows.length === 0 || rows.some((row) => row.length !== numCols)) {
    throw new Error(`${name}: its rows do not all hold ${numCols} values`);
  }
  return {
    values: Float64Array.from(rows.flat()),
    numRows: rows.length,
    numCols,
  };
}
