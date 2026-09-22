/**
 * The files the tests read: the reference files of the repository, and the
 * small VCFs a test writes for a case that no reference file has.
 *
 * The reference VCFs are those of `docs/specs/io_vcf.md`, in
 * `tests/reference/vcf/` at the root of the repository, which
 * `tests/reference/vcf/make_reference.py` writes and bcftools was run on.
 * The reference vars files are those of `docs/specs/io_vars.md`, in
 * `tests/reference/vars/`, which `tests/reference/vars/make_reference.py`
 * writes with pyarrow because popnei cannot write them. The Python tests
 * read the same files.
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
