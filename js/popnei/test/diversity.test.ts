/**
 * How much variety each population holds, from TypeScript: the alleles each
 * population called, the private ones among them, the variants that vary in
 * it and F_IS, over the panel and over the cases a user can reach.
 *
 * `docs/specs/diversity.md` has the five statistics and, under "How it is
 * verified" of each, the program its numbers come from. The literals here
 * are the ones `tests/test_diversity.py` asserts, which is what goal 1 of
 * `docs/objectives.md` asks for: the two packages give a user the same
 * numbers. The Python suite reads them from the files of
 * `tests/reference/diversity/`, which `adegenet` 2.1.11, `poppr` 2.9.8 and
 * `scikit-allel` 1.3.13 wrote; node runs none of those programs, so they are
 * written here as literals, as section 11 of `docs/architecture.md` has it
 * for the tests of this package.
 *
 * The dataset is the panel of `docs/specs/stats.md`,
 * `tests/reference/stats/panel.vcf.gz`: 1200 biallelic diploid variants of
 * 200 individuals, 3 in 100 genotypes missing whole, in the three
 * populations `p0`, `p1` and `p2` of 48, 68 and 84 individuals that
 * `panel_pops_bcftools.txt` beside it holds. At `minNumIndividuals` 20 all
 * 1200 variants count for all three.
 *
 * The draw of a common number of called alleles is refused and not
 * asserted: the three standardized columns and the folded spectrum are work
 * package 3 of `docs/plans/diversity.md`, and until they are there a call
 * that gives `numCalledAlleles`, or that asks for the spectrum, is an
 * `Error`, since the NaN and the 0 such a call would read are what a result
 * gives for a draw no population reached. One test holds those two refusals
 * and fails the day the draw arrives, which is what it is for.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PopDiversity, Variants } from "popnei";
import { calcPopDiversity, init, openVcf } from "popnei";

import { referenceStats, vcfOf } from "./reference.ts";

await init();

/** The bytes of the panel, gzipped, read once for the whole file. */
const PANEL = await referenceStats("panel.vcf.gz");

/** The 1200 variants of the panel. */
const PANEL_NUM_VARS = 1200;

/**
 * How many called genotypes a population needs at a variant for the variant
 * to count for it, which is what "How it is verified" of every item of the
 * spec ran the reference programs at.
 */
const PANEL_MIN_NUM_INDIVIDUALS = 20;

/**
 * The individuals of each population of the panel, under the name of the
 * population, the names in the order `p0`, `p1`, `p2`, which is the order
 * every result of a call with these populations is in.
 */
const PANEL_POPS = popsOf(
  new TextDecoder().decode(await referenceStats("panel_pops_bcftools.txt")),
);

/** The three populations in the order of the keys of `PANEL_POPS`. */
const PANEL_POP_NAMES = ["p0", "p1", "p2"];

/**
 * The alleles the three populations called over the 1200 variants and the
 * mean of each, which is the allelic richness, from `adegenet`.
 */
const PANEL_NUM_ALLELES = [2373, 2377, 2384];
const PANEL_MEAN_NUM_ALLELES = [
  1.9775, 1.9808333333333332, 1.9866666666666666,
];

/**
 * Of those, the ones no other population of the call called at the same
 * variant, from `poppr`, and their mean over the 1200 variants that counted
 * for every population.
 */
const PANEL_PRIVATE_ALLELES = [0, 0, 1];
const PANEL_MEAN_PRIVATE_ALLELES = [0, 0, 0.0008333333333333334];

/**
 * The variants each population called more than one allele at, from
 * `adegenet`, and that count over the 1200 variants that counted for the
 * population.
 */
const PANEL_VARIABLE_VARS = [1173, 1177, 1184];
const PANEL_VARIABLE_VARS_RATIO = [
  0.9775, 0.9808333333333333, 0.9866666666666667,
];

/**
 * The unbiased F_IS of the three populations, from "How it is verified" of
 * "The inbreeding coefficient F_IS" of the spec. It is the form
 * `calcPopDiversity` gives: one minus the mean observed heterozygosity over
 * the mean unbiased expected one.
 *
 * The item prints those three to ten decimals for a reader and gives them
 * again at the precision they were computed to, which is what is written
 * here: a value of ten decimals stands for anything within 5e-11 of itself,
 * and 1e-12 of 0.0128 is 1.3e-14, so the shorter form could not be compared
 * within the tolerance this file uses. They come from
 * `docs/reports/diversity-method/panel.py`, which computes the five
 * quantities of the spec in Python, so this assertion says that popnei
 * agrees with that Python and nothing more. What checks F_IS against a
 * program outside popnei is the plain form of it, which
 * `tests/test_diversity.py` builds from `calcPerVarDistribs` and compares
 * with `scikit-allel`.
 */
const PANEL_FIS = [
  -0.012758486763377208, -0.018110713076467277, -0.018458583231322434,
];

/**
 * How far a value of popnei may be from the number it is compared with,
 * which is what every item of the spec asks for its floats: 1e-12 of the
 * value. Every number here carries all its digits, and the two sides add the
 * same per variant values in different orders, so the last bits differ.
 */
const OF_A_REFERENCE_VALUE = 1e-12;

/**
 * The four statistics of this work package, which are every one but the
 * folded spectrum: that one needs a draw, and the draw is work package 3 of
 * `docs/plans/diversity.md`.
 */
const WITH_NO_SPECTRUM = [
  "num_alleles",
  "private_alleles",
  "variable_vars_ratio",
  "fis",
] as const;

/**
 * Four variants of the three diploid individuals `ind1`, `ind2` and `ind3`,
 * for the counts that the panel cannot show: there every population counted
 * every variant, so the variants of a population and the variants of every
 * population are the same 1200 and a divisor read from the wrong one of the
 * two gives the same number.
 *
 * The genotypes are `0/0 1/1 1/1`, `0/0 0/0 0/0`, `0/1 0/0 0/0` and
 * `0/0 ./. ./.`: the last one is called in `ind1` alone, so a population of
 * `ind2` and `ind3` counts three variants where a population of `ind1`
 * counts four.
 */
const SOME_VARS_MISSED = vcfOf([
  "chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1",
  "chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0",
  "chr1\t3\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0",
  "chr1\t4\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t./.",
]);

/**
 * The individuals of each population of a file of two columns, the name of
 * an individual and the name of its population, under the names of the
 * populations in order.
 *
 * The order of the keys is the order of every array of a result, so the
 * names are sorted and not taken as the file has them: what the literals
 * above are written in is `p0`, `p1`, `p2`.
 */
function popsOf(text: string): Record<string, string[]> {
  const ofEachName = new Map<string, string[]>();
  for (const line of text.split("\n")) {
    if (line === "") {
      continue;
    }
    const [individual, pop] = line.split("\t");
    if (individual === undefined || pop === undefined) {
      throw new Error(
        `the line \`${line}\` is not an individual and a population`,
      );
    }
    const individuals = ofEachName.get(pop) ?? [];
    individuals.push(individual);
    ofEachName.set(pop, individuals);
  }
  const pops: Record<string, string[]> = {};
  for (const pop of [...ofEachName.keys()].sort()) {
    pops[pop] = ofEachName.get(pop) as string[];
  }
  return pops;
}

/** The 1200 variants of the panel. */
function panel(): Variants {
  return openVcf(PANEL);
}

/**
 * The diversity of the three populations of the panel, with the threshold
 * the reference programs were run at, and the variants freed when it
 * returns.
 */
function ofThePanel(
  options: Parameters<typeof calcPopDiversity>[1] = {},
): PopDiversity {
  const variants = panel();
  try {
    return calcPopDiversity(variants, {
      pops: PANEL_POPS,
      minNumIndividuals: PANEL_MIN_NUM_INDIVIDUALS,
      ...options,
    });
  } finally {
    variants.free();
  }
}

/** That `found` holds the numbers of `expected`, each within
 * `OF_A_REFERENCE_VALUE` of its value. */
function assertValues(
  found: Float64Array | Uint32Array | null | undefined,
  expected: readonly number[],
  what: string,
): void {
  if (found === null || found === undefined) {
    throw new Error(`the result holds no ${what}`);
  }
  assert.equal(found.length, expected.length, `the populations of ${what}`);
  for (const [pop, value] of expected.entries()) {
    const ours = found[pop] as number;
    assert.ok(
      Math.abs(ours - value) <=
        OF_A_REFERENCE_VALUE * Math.max(Math.abs(ours), Math.abs(value)),
      `the ${what} of ${PANEL_POP_NAMES[pop]} is ${ours} and not ${value}`,
    );
  }
}

/** The counts of a statistic that was asked for, or the `Error` of a result
 * that holds none. */
function counted<T>(counts: T | null, statistic: string): T {
  if (counts === null) {
    throw new Error(`the pass gives no ${statistic}`);
  }
  return counts;
}

test("the alleles of the panel and its F_IS are the numbers of the reference programs", () => {
  // All 1200 variants of the panel count for all three populations at a
  // threshold of 20 called genotypes, so the counts below are over the
  // whole panel and the private alleles are over the same 1200 variants.
  const diversity = ofThePanel({ stats: WITH_NO_SPECTRUM });

  assert.deepEqual(diversity.pops, PANEL_POP_NAMES);
  assert.deepEqual(
    counted(diversity.numAlleles, "count of alleles").total,
    Uint32Array.from(PANEL_NUM_ALLELES),
  );
  assert.deepEqual(
    counted(diversity.privateAlleles, "count of private alleles").total,
    Uint32Array.from(PANEL_PRIVATE_ALLELES),
  );
  assert.deepEqual(
    counted(diversity.variableVarsRatio, "count of variable variants").total,
    Uint32Array.from(PANEL_VARIABLE_VARS),
  );
  assertValues(diversity.fis, PANEL_FIS, "F_IS");
  assert.deepEqual(
    diversity.numVars.withData,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
  assert.equal(diversity.numVarsEveryPop, PANEL_NUM_VARS);
  assert.equal(diversity.passStats.numVars, PANEL_NUM_VARS);
  assert.deepEqual(diversity.passStats.filtering, {});
});

test("the means and the ratios of the panel are the totals over the variants that counted", () => {
  // The mean alleles and the ratio of variable variants are over the
  // variants that counted for the population, and the mean private alleles
  // over the variants that counted for every population.
  const diversity = ofThePanel({ stats: WITH_NO_SPECTRUM });

  assertValues(
    counted(diversity.numAlleles, "count of alleles").mean,
    PANEL_MEAN_NUM_ALLELES,
    "mean alleles",
  );
  assertValues(
    counted(diversity.privateAlleles, "count of private alleles").mean,
    PANEL_MEAN_PRIVATE_ALLELES,
    "mean private alleles",
  );
  assertValues(
    counted(diversity.variableVarsRatio, "count of variable variants").ratio,
    PANEL_VARIABLE_VARS_RATIO,
    "ratio of variable variants",
  );
});

test("a statistic that was not asked for has no value", () => {
  // The counts of the variants are there whatever was asked for: they are
  // what every mean and every ratio is over.
  const diversity = ofThePanel({ stats: ["fis"] });

  assertValues(diversity.fis, PANEL_FIS, "F_IS");
  assert.equal(diversity.numAlleles, null);
  assert.equal(diversity.privateAlleles, null);
  assert.equal(diversity.variableVarsRatio, null);
  assert.equal(diversity.foldedSfs, null);
  assert.deepEqual(
    diversity.numVars.withData,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
});

test("one population holds every allele it called as a private one", () => {
  // With no `pops` there is one population of every individual, and with no
  // other population to hold them every allele it called is private.
  const diversity = ofThePanel({
    pops: undefined,
    stats: ["num_alleles", "private_alleles"],
  });

  assert.deepEqual(diversity.pops, ["pop"]);
  const total = counted(diversity.numAlleles, "count of alleles").total;
  assert.ok((total[0] as number) > 0);
  assert.deepEqual(
    counted(diversity.privateAlleles, "count of private alleles").total,
    total,
  );
  assert.equal(diversity.numVarsEveryPop, PANEL_NUM_VARS);
});

test("a call that names no statistic gives the four that need no draw", () => {
  // `calcPopDiversity(variants)` with no options at all: the four statistics
  // that need no draw and no spectrum, over one population of every
  // individual. The bins of the spectrum are counts of the rarer allele in a
  // draw of `numCalledAlleles`, so while the default was all five this call
  // refused itself and named an option the user had not written.
  const variants = panel();
  let diversity: PopDiversity;
  try {
    diversity = calcPopDiversity(variants);
  } finally {
    variants.free();
  }

  assert.deepEqual(diversity.pops, ["pop"]);
  const alleles = counted(diversity.numAlleles, "count of alleles");
  assert.ok((alleles.total[0] as number) > 0);
  assert.ok(
    (counted(diversity.privateAlleles, "count of private alleles")
      .total[0] as number) > 0,
  );
  assert.ok(
    (counted(diversity.variableVarsRatio, "ratio of variable variants")
      .total[0] as number) > 0,
  );
  assert.ok(!Number.isNaN((diversity.fis as Float64Array)[0] as number));
  // The one statistic that needs a draw is not in the default, and naming it
  // without a draw is refused as it was.
  assert.equal(diversity.foldedSfs, null);
  assert.deepEqual(diversity.numVars.withData, Uint32Array.of(PANEL_NUM_VARS));
});

test("a draw of a common number of called alleles is refused until it is calculated", () => {
  // THE DAY THE DRAW ARRIVES THIS TEST FAILS, AND THAT IS WHAT IT IS FOR.
  // Nothing of the draw is calculated: the three `inDraw` columns are NaN,
  // `numVars.inDraw` is 0 and `foldedSfs` is `null`, which are the values a
  // result gives for a draw no population reached and for a statistic nobody
  // asked for, so a call that asked for a draw would read them as an answer.
  // Both calls are refused instead.
  //
  // Work package 3 of `docs/plans/diversity.md` calculates it, and whoever
  // does it takes the two refusals out of `js/popnei/src/diversity.ts` and
  // replaces this test with the values of "How it is verified" of the spec:
  // at `numCalledAlleles` 20 the mean alleles of the panel in the draw are
  // 1.9283948650041205, 1.9219209943237829 and 1.9197370843937562, its
  // ratios of variable variants in the draw are those three less one, and
  // `numVars.inDraw` is 1200 for each of the three populations.
  assert.throws(
    () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: 20 }),
    { message: /`numCalledAlleles` is 20 and the draw/ },
  );
  assert.throws(() => ofThePanel({ stats: ["folded_sfs"] }), {
    message: /folded site frequency spectrum/,
  });
  assert.throws(() => ofThePanel({ stats: ["fis", "folded_sfs"] }), {
    message: /folded site frequency spectrum/,
  });
});

test("a numCalledAlleles that is no draw at all is refused as the wrong argument it is", () => {
  // The draw is refused above whatever its value, and a value that is no
  // draw is refused before that, with the rule it broke: a number of
  // JavaScript reaches a whole number of the core as 32 bits with no error,
  // so a draw of 20.5 alleles would be a draw of 20 and one of -1 a draw of
  // 4294967295, and a draw of one allele finds one allele whatever the
  // population holds.
  for (const given of [20.5, -1, 0, 1]) {
    assert.throws(
      () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: given }),
      { message: /a whole number of 2 or more/ },
      `numCalledAlleles: ${given}`,
    );
  }
});

test("the mean private alleles are over the variants that counted for every population", () => {
  // The four variants of `SOME_VARS_MISSED`, in two populations, `p1` of
  // `ind1` and `p2` of `ind2` and `ind3`, at a threshold of one called
  // genotype. `p2` called nothing at the fourth variant, so that variant
  // counts for `p1` and not for `p2`: `p1` counts 4 variants, `p2` counts 3,
  // and 3 counted for both.
  //
  // The private alleles are over the variants that counted for every
  // population, the fourth left out, where the allele `p1` called would be
  // private only because `p2` has no data there. `p1` holds a private allele
  // at the first variant, the `0` that `p2` did not call, and at the third,
  // the `1`, which makes 2 over 3 variants, 0.6666666666666666; with the
  // variants of `p1` itself as the divisor it would be 2 over 4, 0.5. `p2`
  // holds the `1` of the first variant, 1 over 3.
  const variants = openVcf(SOME_VARS_MISSED);

  const diversity = calcPopDiversity(variants, {
    pops: { p1: ["ind1"], p2: ["ind2", "ind3"] },
    stats: ["num_alleles", "private_alleles", "variable_vars_ratio"],
    minNumIndividuals: 1,
  });

  assert.deepEqual(diversity.pops, ["p1", "p2"]);
  assert.deepEqual(diversity.numVars.withData, Uint32Array.of(4, 3));
  assert.equal(diversity.numVarsEveryPop, 3);
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  assert.deepEqual(theirPrivate.total, Uint32Array.of(2, 1));
  assert.deepEqual(
    theirPrivate.mean,
    Float64Array.of(0.6666666666666666, 0.3333333333333333),
  );
  // The alleles called and the variable variants are over the variants of
  // each population, the fourth among those of `p1`: `p1` called one allele
  // at the first, the second and the fourth and two at the third, and it is
  // the third that makes it variable.
  const alleles = counted(diversity.numAlleles, "count of alleles");
  assert.deepEqual(alleles.total, Uint32Array.of(5, 3));
  assert.deepEqual(alleles.mean, Float64Array.of(1.25, 1));
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  assert.deepEqual(variable.total, Uint32Array.of(1, 0));
  assert.deepEqual(variable.ratio, Float64Array.of(0.25, 0));
  variants.free();
});

test("a population no variant counted for has no mean, no ratio and no F_IS", () => {
  // At a threshold of 60 called genotypes the 48 individuals of `p0` can
  // never reach it, so no variant counts for that population: its totals are
  // 0, its mean, its ratio and its F_IS are NaN, and the other two
  // populations are given as they are. The variants that counted for every
  // population are 0 too, so no population has a mean private allele, `p1`
  // and `p2` included, whose own counts are above 0: that is the divisor of
  // the private alleles and not the variants of each population.
  const diversity = ofThePanel({
    stats: WITH_NO_SPECTRUM,
    minNumIndividuals: 60,
  });

  const alleles = counted(diversity.numAlleles, "count of alleles");
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  const fis = diversity.fis as Float64Array;
  const p0 = diversity.pops.indexOf("p0");
  assert.equal(diversity.numVars.withData[p0], 0);
  assert.equal(alleles.total[p0], 0);
  assert.equal(variable.total[p0], 0);
  assert.ok(Number.isNaN(alleles.mean[p0] as number), "the mean of p0");
  assert.ok(Number.isNaN(variable.ratio[p0] as number), "the ratio of p0");
  assert.ok(Number.isNaN(fis[p0] as number), "the F_IS of p0");
  assert.equal(diversity.numVarsEveryPop, 0);
  for (const pop of PANEL_POP_NAMES) {
    const which = diversity.pops.indexOf(pop);
    assert.ok(
      Number.isNaN(theirPrivate.mean[which] as number),
      `the mean private alleles of ${pop}`,
    );
  }
  for (const pop of ["p1", "p2"]) {
    const which = diversity.pops.indexOf(pop);
    assert.ok((diversity.numVars.withData[which] as number) > 0, pop);
    assert.ok((alleles.total[which] as number) > 0, pop);
    assert.ok(!Number.isNaN(fis[which] as number), pop);
  }
});

test("a pass that gives no variant is refused and says where the variants went", () => {
  // The whole message is asserted, because it is the one the Python function
  // gives for the same source: the two languages say the same thing, and
  // Python writes the path of the file before it, which the bytes a
  // TypeScript user gives have not.
  const ofNoVariant = openVcf(vcfOf([]));
  try {
    assert.throws(
      () => calcPopDiversity(ofNoVariant, { stats: WITH_NO_SPECTRUM }),
      {
        message:
          "the pass gave no variant and its source holds none: a statistic " +
          "of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    ofNoVariant.free();
  }
  // A major allele frequency is at least one over the alleles of a variant,
  // so a threshold of 0 keeps none of the four variants, and the message
  // names the filter with what it was given and kept.
  const filtered = openVcf(SOME_VARS_MISSED);
  try {
    filtered.filterByMaf(0);
    assert.throws(
      () => calcPopDiversity(filtered, { stats: WITH_NO_SPECTRUM }),
      {
        message:
          "the pass gave no variant: its source gave 4 and the steps kept " +
          "none of them, the `maf` filter was given 4 and kept 0; a " +
          "statistic of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    filtered.free();
  }
});

test("a population that names an individual of no pass is refused", () => {
  // `pops` is read as `calcPerVarDistribs` reads it: a name that is not an
  // individual of the pass names the population and the name.
  assert.throws(
    () =>
      ofThePanel({
        stats: WITH_NO_SPECTRUM,
        pops: { p0: ["s000", "s999"] },
      }),
    { message: /s999/ },
  );
});

test("stats names the statistics there are", () => {
  // Both messages are the core's, so a Python user reads the same sentence:
  // a name of no statistic is refused with the five there are, and a `stats`
  // that names none with what such a pass would do, read every variant of the
  // source and compute nothing of them.
  assert.throws(() => ofThePanel({ stats: ["allelic_richness"] as never }), {
    message: /`allelic_richness` is not one of the statistics/,
  });
  assert.throws(() => ofThePanel({ stats: [] }), {
    message: /`stats` names no statistic and a pass that computes none/,
  });
});

test("variants that is not a Variants is refused", () => {
  assert.throws(() => calcPopDiversity(PANEL as never), {
    message: /openVcf/,
  });
});
