/**
 * The density of the variants along the chromosomes from TypeScript.
 *
 * `docs/specs/stats.md` has it under "The density of the variants along the
 * chromosomes", and what it is compared with is tabix 1.24: the count of
 * each window of 1000 base pairs of `tests/reference/vcf/many.vcf.gz`, which
 * `tests/reference/stats/make_reference.py` keeps in `many.density.tsv`.
 * `many.vcf` holds 500 variants of chr1 and chr2, the first of chr1 at 1000
 * and its last at 10213, those of chr2 at 10250 and 19463, and its
 * `##contig` lines have no length.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { VarDensity } from "popnei";
import { calcVarDensity, init, openVars, openVcf, writeVars } from "popnei";

import { referenceStats, referenceVcf } from "./reference.ts";

await init();

const MANY = await referenceVcf("many.vcf");

/**
 * `write.vcf` of the VCF writer: chr1 100, 1000 and 1001 and chr2 1 and
 * 1500 read with the default, with the lengths chr1 2000 and chr2 1500.
 */
const WRITE = await referenceVcf("write.vcf");

/** One window as the tables of the spec give it. */
type Window = [string, number, number, number];

/** Every window of `density`, in its order. */
function windowsOf(density: VarDensity): Window[] {
  return density.chroms.map((chrom, at) => [
    chrom,
    density.start[at] as number,
    density.end[at] as number,
    density.numVars[at] as number,
  ]);
}

/**
 * The windows of `chrom` of `width` base pairs from the position 1, the last
 * one ending at `lastEnd`, with `counts`.
 */
function laidEndToEnd(
  chrom: string,
  width: number,
  lastEnd: number,
  counts: readonly number[],
): Window[] {
  return counts.map((count, at) => [
    chrom,
    at * width + 1,
    at === counts.length - 1 ? lastEnd : (at + 1) * width,
    count,
  ]);
}

/** `count` times `value`. */
function times(count: number, value: number): number[] {
  return Array.from({ length: count }, () => value);
}

/** A VCF of one individual with `contigs` and a variant at each of `lines`. */
function aVcf(contigs: readonly string[], lines: readonly [string, number][]) {
  const text = [
    "##fileformat=VCFv4.3",
    ...contigs,
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta",
    ...lines.map(([chrom, pos]) => `${chrom}\t${pos}\t.\tA\tT\t.\tPASS\t.\tGT\t0/1`),
    "",
  ];
  return new TextEncoder().encode(text.join("\n"));
}

test("the density of many.vcf in windows of 1000 is the one tabix counts", async () => {
  const density = calcVarDensity(openVcf(MANY, { onlyPassed: false }), 1000);
  const expected = [
    ...laidEndToEnd("chr1", 1000, 11000, [1, ...times(9, 27), 6]),
    ...laidEndToEnd("chr2", 1000, 20000, [...times(10, 0), 21, ...times(8, 27), 13]),
  ];
  assert.deepEqual(windowsOf(density), expected);
  const tabix = new TextDecoder()
    .decode(await referenceStats("many.density.tsv"))
    .trim()
    .split("\n")
    .slice(1)
    .map((line) => {
      const [chrom, start, end, count] = line.split("\t");
      return [chrom, Number(start), Number(end), Number(count)];
    });
  assert.deepEqual(windowsOf(density), tabix);
  assert.ok(density.start instanceof Float64Array);
  assert.ok(density.numVars instanceof Uint32Array);
  assert.equal(density.passStats.numVars, 500);
});

test("the density with chromLengths ends each chromosome at its length", () => {
  const density = calcVarDensity(openVcf(MANY, { onlyPassed: false }), 1000, {
    chromLengths: { chr1: 12000, chr2: 19500 },
  });
  assert.deepEqual(windowsOf(density), [
    ...laidEndToEnd("chr1", 1000, 12000, [1, ...times(9, 27), 6, 0]),
    ...laidEndToEnd("chr2", 1000, 19500, [...times(10, 0), 21, ...times(8, 27), 13]),
  ]);
});

test("the density refuses a variant past a length of chromLengths", () => {
  assert.throws(
    () =>
      calcVarDensity(openVcf(MANY, { onlyPassed: false }), 1000, {
        chromLengths: { chr1: 10000, chr2: 20000 },
      }),
    (error: Error) =>
      error.message.includes("10028") &&
      error.message.includes("10000") &&
      error.message.includes("`chromLengths`"),
  );
});

test("the density refuses a width and lengths that are no whole numbers of 1 or more", () => {
  const variants = openVcf(MANY);
  for (const windowSize of [0, -1, 2.5, Number.NaN, 2 ** 53]) {
    assert.throws(() => calcVarDensity(variants, windowSize), /`windowSize`/);
  }
  assert.throws(
    () => calcVarDensity(variants, 1000, { chromLengths: { chr1: 0 } }),
    /`chromLengths.chr1`/,
  );
  assert.throws(
    () =>
      calcVarDensity(variants, 1000, {
        chromLengths: [12000] as unknown as Record<string, number>,
      }),
    /`chromLengths`/,
  );
  assert.throws(
    () =>
      calcVarDensity(variants, 1000, {
        chromLenghts: {},
      } as unknown as { chromLengths: Record<string, number> }),
    /chromLenghts/,
  );
});

test("the density refuses a window that ends past 2^53 and takes one that ends on it", () => {
  const onIt = aVcf([`##contig=<ID=chr1,length=${2 ** 53}>`], [["chr1", 2 ** 53]]);
  const density = calcVarDensity(openVcf(onIt), 2 ** 52);
  assert.deepEqual(windowsOf(density), [
    ["chr1", 1, 2 ** 52, 0],
    ["chr1", 2 ** 52 + 1, 2 ** 53, 1],
  ]);
  const pastIt = aVcf(["##contig=<ID=chr1,length=9007199254740993>"], [["chr1", 5]]);
  assert.throws(
    () => calcVarDensity(openVcf(pastIt), 2 ** 52),
    /9007199254740993/,
  );
});

test("the density gives the chromosomes in the order of chromLengths", () => {
  // chr2 first, which is neither the order of the names nor that of the
  // variants of write.vcf.
  const density = calcVarDensity(openVcf(WRITE), 500, {
    chromLengths: { chr2: 1500, chr1: 2000 },
  });
  assert.deepEqual(windowsOf(density), [
    ...laidEndToEnd("chr2", 500, 1500, [1, 0, 1]),
    ...laidEndToEnd("chr1", 500, 2000, [1, 1, 1, 0]),
  ]);
});

test("the density takes the lengths in the order JavaScript gives the keys, whole numbers first", () => {
  const lines: [string, number][] = [
    ["X", 1],
    ["10", 1],
    ["2", 1],
  ];
  const density = calcVarDensity(openVcf(aVcf([], lines)), 100, {
    chromLengths: { X: 100, "10": 100, "2": 100 },
  });
  assert.deepEqual(density.chroms, ["2", "10", "X"]);
});

test("the density refuses chromLengths that is a Map or any object but a plain one", () => {
  const variants = openVcf(WRITE);
  for (const given of [
    new Map([["chr1", 10]]),
    new Date(0),
    Object.create({ chr1: 10 }),
  ]) {
    assert.throws(
      () =>
        calcVarDensity(variants, 500, {
          chromLengths: given as unknown as Record<string, number>,
        }),
      /`chromLengths`/,
    );
  }
  const noPrototype = Object.assign(Object.create(null), { chr1: 2000 });
  const density = calcVarDensity(variants, 500, { chromLengths: noPrototype });
  assert.deepEqual(density.chroms.slice(0, 1), ["chr1"]);
});

test("the density of a vars file reads the lengths it keeps", () => {
  const vars = writeVars(openVcf(WRITE)).bytes;
  const density = calcVarDensity(openVars(vars), 500);
  assert.deepEqual(windowsOf(density), [
    ...laidEndToEnd("chr1", 500, 2000, [1, 1, 1, 0]),
    ...laidEndToEnd("chr2", 500, 1500, [1, 0, 1]),
  ]);
});

test("the density refuses the last window of a chromosome with no length that ends past 2^53", () => {
  // The window of 2^52 + 1 that holds 2^53 ends at 2^53 + 2.
  const pastIt = aVcf([], [["chr1", 2 ** 53]]);
  assert.throws(
    () => calcVarDensity(openVcf(pastIt), 2 ** 52 + 1),
    /9007199254740994/,
  );
});
