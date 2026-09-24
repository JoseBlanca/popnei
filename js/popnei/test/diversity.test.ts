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
 * The draw of a common number of called alleles is not asserted here: the
 * standardized values and the folded spectrum are work package 3 of
 * `docs/plans/diversity.md`, and what this file asserts of
 * `numCalledAlleles` is that it is taken and that the call is refused when
 * no draw can be made.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PopDiversity, Variants } from "popnei";
import { calcPopDiversity, init, openVcf } from "popnei";

import { referenceStats } from "./reference.ts";

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

test("the totals and the F_IS do not read the draw", () => {
  // A draw changes the standardized values alone: the alleles called, the
  // private ones, the variable variants and F_IS are the same numbers with
  // `numCalledAlleles` 20 as with none.
  const withNoDraw = ofThePanel({ stats: WITH_NO_SPECTRUM });

  const ofADrawOf20 = ofThePanel({
    stats: WITH_NO_SPECTRUM,
    numCalledAlleles: 20,
  });

  assert.deepEqual(
    counted(ofADrawOf20.numAlleles, "count of alleles").total,
    counted(withNoDraw.numAlleles, "count of alleles").total,
  );
  assert.deepEqual(
    counted(ofADrawOf20.privateAlleles, "count of private alleles").total,
    counted(withNoDraw.privateAlleles, "count of private alleles").total,
  );
  assert.deepEqual(
    counted(ofADrawOf20.variableVarsRatio, "count of variable variants").total,
    counted(withNoDraw.variableVarsRatio, "count of variable variants").total,
  );
  assert.deepEqual(ofADrawOf20.fis, withNoDraw.fis);
});

test("the spectrum asked for with no draw is refused", () => {
  // The bins of a folded spectrum are the counts of the rarer allele in a
  // draw, so the spectrum needs `numCalledAlleles`. It is among the
  // statistics of a call that names none, so such a call is refused too.
  assert.throws(() => ofThePanel({ stats: ["folded_sfs"] }), {
    message: /numCalledAlleles/,
  });
  assert.throws(() => ofThePanel(), {
    message: /folded site frequency spectrum/,
  });
});

test("a draw of fewer than two alleles is refused", () => {
  // A draw of one allele finds one allele whatever the population holds, so
  // every standardized value of it would say nothing.
  assert.throws(
    () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: 1 }),
    { message: /numCalledAlleles/ },
  );
});

test("a numCalledAlleles that is no whole number is refused", () => {
  // A number of JavaScript reaches a whole number of the core as 32 bits
  // with no error, so a draw of 20.5 alleles would be a draw of 20 and a
  // draw of -1 one of 4294967295.
  assert.throws(
    () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: 20.5 }),
    { message: /numCalledAlleles/ },
  );
  assert.throws(
    () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: -1 }),
    { message: /numCalledAlleles/ },
  );
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
  // A name of no statistic is refused with the five there are, and a call
  // that names none would make a pass over the whole source for nothing.
  assert.throws(() => ofThePanel({ stats: ["allelic_richness"] as never }), {
    message: /num_alleles/,
  });
  assert.throws(() => ofThePanel({ stats: [] }), {
    message: /names no statistic/,
  });
});

test("variants that is not a Variants is refused", () => {
  assert.throws(() => calcPopDiversity(PANEL as never), {
    message: /openVcf/,
  });
});
