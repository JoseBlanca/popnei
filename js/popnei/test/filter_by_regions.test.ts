/**
 * The filter by regions from TypeScript: which variants `filterByRegions`
 * keeps, what it counts, the steps of both kinds and what it refuses.
 *
 * The item "The filter by regions" of `docs/specs/filters.md` has the rule
 * and the numbers. The file is `many.vcf` of `docs/specs/io_vcf.md`, 500
 * variants of 50 diploid individuals, read with every variant given, and
 * the BED is `tests/reference/filters/regions.bed`, whose six regions join
 * into five. bcftools 1.24 and plink2 v2.0.0-a.7.7 keep 45 of the variants
 * with it and 455 with the other side, which
 * `tests/reference/filters/make_reference.py` stored by chromosome and
 * position in `regions.txt` and `excluded_regions.txt`. What this file adds
 * to the cargo and the pytest tests is that the filter gives them under
 * WebAssembly, where the rows of a block are read one after another, and
 * that a user reaches it through `filterByRegions`.
 */

import assert from "node:assert/strict";
import { test } from "node:test";
import { gzipSync } from "node:zlib";

import type { PassStats, Variants } from "popnei";
import { init, openVcf } from "popnei";

import { referenceFilters, referenceVcf } from "./reference.ts";

await init();

const MANY_VCF = await referenceVcf("many.vcf");
const REGIONS_BED = await referenceFilters("regions.bed");

/** The variants bcftools and plink2 kept, as `chrom:pos`. */
async function theReference(name: string): Promise<string[]> {
  const text = new TextDecoder().decode(await referenceFilters(`${name}.txt`));
  return text
    .split("\n")
    .filter((line) => line !== "")
    .map((line) => line.replace("\t", ":"));
}

const THE_45 = await theReference("regions");
const THE_455 = await theReference("excluded_regions");

function theDataset(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/** Every variant of one pass as `chrom:pos`, with the counts of the pass. */
function keptBy(
  variants: Variants,
  numVarsPerBlock?: number,
): { kept: string[]; passStats: PassStats } {
  const blocks = variants.iterBlocks({
    fields: ["chrom", "pos"],
    ...(numVarsPerBlock === undefined ? {} : { numVarsPerBlock }),
  });
  const kept: string[] = [];
  for (const block of blocks) {
    if (block.chrom === null || block.pos === null) {
      throw new Error("a block without the chromosome or the position");
    }
    const positions = block.pos;
    kept.push(...block.chrom.map((chrom, row) => `${chrom}:${positions[row]}`));
  }
  return { kept, passStats: blocks.passStats };
}

test("filterByRegions keeps the 45 of bcftools and plink2 and the 455 with exclude", () => {
  assert.equal(THE_45.length, 45);
  assert.equal(THE_455.length, 455);
  for (const numVarsPerBlock of [7, undefined]) {
    const inside = theDataset();
    inside.filterByRegions(REGIONS_BED);
    const ofTheInside = keptBy(inside, numVarsPerBlock);
    assert.deepEqual(ofTheInside.kept, THE_45);
    assert.deepEqual(ofTheInside.passStats.filtering, {
      regions: { varsProcessed: 500, varsKept: 45 },
    });

    const outside = theDataset();
    outside.filterByRegions(REGIONS_BED, { exclude: true });
    const ofTheOutside = keptBy(outside, numVarsPerBlock);
    assert.deepEqual(ofTheOutside.kept, THE_455);
    assert.deepEqual(ofTheOutside.passStats.filtering, {
      excluded_regions: { varsProcessed: 500, varsKept: 455 },
    });
  }
  // The three of chr1 past 2000, and chr2 10250 alone, as the spec names them.
  assert.deepEqual(
    THE_45.filter((variant) => {
      const [chrom, pos] = variant.split(":");
      return (chrom === "chr1" && Number(pos) > 2000) || variant === "chr2:10250";
    }),
    ["chr1:4996", "chr1:5033", "chr1:5070", "chr2:10250"],
  );
});

test("filterByRegions adds a step of each kind with the number of regions", () => {
  const variants = theDataset();
  variants.filterByRegions(REGIONS_BED);
  variants.filterByRegions(new TextEncoder().encode("chr1\t0\t1100\n"), {
    exclude: true,
  });
  assert.deepEqual(variants.steps, [
    { kind: "regions", args: { numRegions: 5 } },
    { kind: "excluded_regions", args: { numRegions: 1 } },
  ]);
  // Of the 45, the three of chr1 up to 1100 are excluded.
  const { kept, passStats } = keptBy(variants);
  assert.equal(kept.length, 42);
  assert.deepEqual(passStats.filtering, {
    regions: { varsProcessed: 500, varsKept: 45 },
    excluded_regions: { varsProcessed: 45, varsKept: 42 },
  });
});

test("a wrong line of the BED is an Error that names the line", () => {
  const variants = theDataset();
  assert.throws(
    () =>
      variants.filterByRegions(
        new TextEncoder().encode("# a comment\nchr1\t0\t10\nchr1\t5\t5\n"),
      ),
    { message: /^line 3 of the BED file: its start, 5, is not below its end, 5/ },
  );
  assert.throws(() => variants.filterByRegions(new TextEncoder().encode("chr1 0 10\n")), {
    message: /line 1 of the BED file: .*separates by tabs/,
  });
  assert.throws(() => variants.filterByRegions(new TextEncoder().encode("track\n")), {
    message: /the BED file holds no region/,
  });
  assert.deepEqual(variants.steps, []);
});

test("a second filterByRegions of a kind is refused and the other kind is not", () => {
  const variants = theDataset();
  variants.filterByRegions(REGIONS_BED);
  assert.throws(() => variants.filterByRegions(REGIONS_BED), {
    message: /filtered by regions already/,
  });
  variants.filterByRegions(REGIONS_BED, { exclude: true });
  assert.equal(variants.steps.length, 2);
});

test("a bed that is not a Uint8Array and an exclude that is not a boolean are refused", () => {
  const variants = theDataset();
  assert.throws(() => variants.filterByRegions("chr1\t0\t10\n" as unknown as Uint8Array), {
    message: /`bed` is the bytes of the file, a Uint8Array/,
  });
  assert.throws(
    () => variants.filterByRegions(REGIONS_BED, { exclude: 1 as unknown as boolean }),
    { message: /`exclude` is true or false/ },
  );
  assert.throws(
    () => variants.filterByRegions(REGIONS_BED, null as unknown as { exclude?: boolean }),
    { message: /the options of `filterByRegions` are an object/ },
  );
  assert.deepEqual(variants.steps, []);
});

test("an option of filterByRegions that is not one is refused and named", () => {
  const variants = theDataset();
  assert.throws(
    () =>
      variants.filterByRegions(REGIONS_BED, {
        excluded: true,
      } as unknown as { exclude?: boolean }),
    {
      message:
        "popnei: `excluded` is not an option of `filterByRegions`, whose options are `exclude`",
    },
  );
  assert.deepEqual(variants.steps, []);
});

test("a gzipped BED cut short is an Error that says the BED could not be read", () => {
  const gzipped = gzipSync(REGIONS_BED);
  const cut = gzipped.subarray(0, gzipped.length - 10);
  const variants = theDataset();
  assert.throws(() => variants.filterByRegions(new Uint8Array(cut)), {
    message: /^the BED could not be read: /,
  });
  assert.deepEqual(variants.steps, []);
});
