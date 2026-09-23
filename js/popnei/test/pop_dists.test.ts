/**
 * The distances between populations from TypeScript: `calcPopDists` and the
 * `PopDists` it gives.
 *
 * `docs/specs/dists.md` has the seven measures and the numbers this file
 * asserts, and the files the reference programs wrote them into are in
 * `tests/reference/pop_dists/`. The literals are the spec's and are the ones
 * `tests/test_pop_dists.py` asserts, which is what goal 1 of
 * `docs/objectives.md` asks for: the two packages give a user the same
 * numbers.
 *
 * Three programs are compared with. plink2 v2.0.0-a.7.7 gives Hudson's F_ST
 * of both panels, which it prints to six digits, so the comparison is within
 * 1e-6 absolute. ADMIXTOOLS 2.0.10 gives f_2 and its jackknife standard
 * error of the biallelic panel, which it writes to seventeen digits, so that
 * comparison is within 1e-12 relative. mmod 1.3.3 under R 4.6.1 gives Jost's
 * D, Nei's G_ST and the standardized G''_ST of both panels with another
 * estimator of each, so those three are an agreement within 5e-4 and not an
 * equality.
 *
 * pyNei is the fourth, and it is a Python library that cannot be installed
 * here: it is the one program that computes the estimator of Jost's D that
 * popnei computes, `tests/test_pop_dists.py` runs it and matches it within
 * 1e-12 relative, and what this file asserts are the ten digits of the spec
 * that that run gave.
 *
 * ADMIXTOOLS was run on the biallelic panel at three lengths of the
 * resampling groups and `calcPopDists` refuses a pass of fewer than 20 of
 * them, so of the three runs only the one at 55 000 base pairs, which cuts
 * the panel into 22 groups, can be asked for through this package; the other
 * two, at 100 000 and at 250 000, are cargo tests of
 * `crates/popnei/src/pop_dists.rs`.
 *
 * All seven measures are calculated, so the refusal of one that popnei has
 * no value for has nothing left to refuse.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type {
  CalcPopDistsOptions,
  PopDistMeasure,
  PopDists,
  Variants,
} from "popnei";
import {
  calcPairwiseKosmanDists,
  calcPopDists,
  init,
  openVcf,
} from "popnei";

import {
  popsOfTheFile,
  referenceDists,
  referencePopDists,
  referenceStats,
  vcfOf,
} from "./reference.ts";

await init();

/**
 * The biallelic panel of "How it is verified" of the spec: 1200 variants of
 * 200 diploid individuals over two chromosomes of 600 each, 3 in 100
 * genotypes missing whole, in three populations of 48, 68 and 84.
 */
const PANEL_VCF = await referenceDists("panel.vcf.gz");
const PANEL_NUM_VARS = 1200;
const PANEL_POPS = popsOfTheFile(await referenceStats("panel_pops.txt"));

/**
 * The multiallelic panel of the same section, written for this item: 120
 * microsatellite loci of 90 diploid individuals in three populations of 30,
 * six alleles to a locus, 4 in 100 genotypes missing whole. It is the panel
 * that shows that the arithmetic does not assume two alleles.
 */
const MICRO_VCF = await referencePopDists("micro.vcf.gz");
const MICRO_NUM_VARS = 120;
const MICRO_POPS = popsOfTheFile(await referencePopDists("micro_pops.txt"));

/**
 * The three pairs of both panels, in the order of the distance vector, which
 * is the order the populations are named in.
 */
const PAIRS = ["p0-p1", "p0-p2", "p1-p2"];

/**
 * Hudson's F_ST of the three pairs of each panel, which plink2 v2.0.0-a.7.7
 * printed with `--fst popcat method=hudson`, from the F_ST item of the spec.
 */
const PANEL_FST = [0.104962, 0.102736, 0.109621];
const MICRO_FST = [0.0642281, 0.0694106, 0.0699954];

/**
 * plink2 prints six digits, so half a unit of the last of them is the
 * distance the numbers can be apart, which is what the spec asks for.
 */
const FST_TOLERANCE = 1e-6;

/**
 * The length in base pairs that cuts the biallelic panel into 22 resampling
 * groups, ten of 55 variants and one of 50 on each of its two chromosomes.
 * It is the only one of the three lengths ADMIXTOOLS was run at that reaches
 * the 20 groups `calcPopDists` asks for.
 */
const PANEL_JACKKNIFE_GROUP = 55000;
const PANEL_NUM_GROUPS = 22;

/**
 * How many variants each of those groups holds, in the order they were
 * started: the panel has its variants 1000 base pairs apart, so 55 of them
 * fall in the first 55 000 of a chromosome and the eleventh group of each
 * chromosome holds the 50 that are left.
 */
const PANEL_VARS_OF_EACH_GROUP = [
  ...Array.from({ length: 10 }, () => 55),
  50,
  ...Array.from({ length: 10 }, () => 55),
  50,
];

/**
 * f_2 and its jackknife standard error of the three pairs of the biallelic
 * panel at that length, which ADMIXTOOLS 2.0.10 wrote into
 * `tests/reference/pop_dists/panel.f2.min20.tsv`, from the f_2 item of the
 * spec and from "The standard errors".
 */
const PANEL_F2 = [
  0.041181109098151751, 0.039890789655075108, 0.042798563747209556,
];
const PANEL_F2_STANDARD_ERRORS = [
  0.0020502481330704485, 0.0016837006366670754, 0.0019859616713111257,
];

/**
 * Both libraries add the same numbers in a different order, so the last bits
 * are what they can differ by, which is the tolerance the spec gives.
 */
const F2_TOLERANCE = 1e-12;

/**
 * The threshold of called genotypes at which the pairs of the biallelic
 * panel part, from "How it is verified" of the spec: the variants that count
 * are 688 for p0-p1 and p0-p2, whose p0 has 48 individuals, and all 1200 for
 * p1-p2. Every variant counts for every pair at the default of 20.
 */
const PARTING_MIN_NUM_INDIVIDUALS = 47;
const PANEL_NUM_VARS_OF_EACH_PAIR = [688, 688, 1200];

/**
 * Jost's D of the three pairs of each panel at the default threshold of 20
 * called genotypes, which pyNei's `calc_jost_dest_pop_dists` gives and which
 * the Jost's D item of the spec prints to ten digits. pyNei is the one
 * program that computes the estimator popnei computes; it is not installable
 * here, so what this package checks is that it gives the numbers
 * `tests/test_pop_dists.py` got out of a live pyNei.
 */
const PANEL_DEST = [0.063543463, 0.0612981314, 0.0656705213];
const MICRO_DEST = [0.1661307946, 0.1822136338, 0.1819999162];

/**
 * The same three of the biallelic panel at the 47 where its pairs part. The
 * D of p1-p2 is the one above, since that pair keeps all 1200 variants.
 */
const PANEL_DEST_AT_THE_PARTING = [0.0595097904, 0.0612873957, 0.0656705213];

/**
 * Nei's G_ST and the standardized G''_ST of the three pairs of each panel at
 * that same default, which the G_ST item of the spec prints to ten digits as
 * popnei's own numbers: no program outside popnei computes this estimator of
 * either, so what they are checked against is mmod below.
 */
const PANEL_GST = [0.0554614481, 0.0542279036, 0.0580557529];
const PANEL_GST_STANDARDIZED = [0.161959631, 0.1578689665, 0.1682042511];
const MICRO_GST = [0.0331595308, 0.0360169279, 0.0363566716];
const MICRO_GST_STANDARDIZED = [0.2196573038, 0.2390740032, 0.2393928219];

/**
 * The six sets of literals above are the spec's ten digits, so half a unit
 * of the tenth is how far a value may be from the one written there.
 */
const TEN_DIGITS_TOLERANCE = 5e-11;

/**
 * Jost's D, Nei's G_ST and the standardized G''_ST of the three pairs of each
 * panel as `pairwise_D`, `pairwise_Gst_Nei` and `pairwise_Gst_Hedrick` of
 * mmod 1.3.3 under R 4.6.1 give them, which
 * `tests/reference/pop_dists/panel.mmod.tsv` and `micro.mmod.tsv` hold.
 * `pairwise_Gst_Hedrick` computes the standardized G''_ST of Meirmans and
 * Hedrick (2011) and not the G'_ST its name suggests, which "How it is
 * verified" of the G_ST item of the spec shows from its source.
 */
const PANEL_MMOD = {
  dest: [0.0634704859, 0.0612312792, 0.0656071139],
  gst: [0.055389669, 0.0541614926, 0.0579922831],
  gst_standardized: [0.1617736271, 0.1576967932, 0.1680418436],
};
const MICRO_MMOD = {
  dest: [0.1662094246, 0.1819580763, 0.1816462134],
  gst: [0.0331789783, 0.0359529575, 0.0362672006],
  gst_standardized: [0.219761268, 0.238738698, 0.2389275805],
};

/**
 * mmod computes another estimator of the same three quantities: it leaves the
 * observed heterozygosity term out of both corrections and uses 2n/(2n - 1)
 * where popnei, which is pyNei and Nei and Chesser (1983), uses n/(n - 1) and
 * subtracts H_obs/(2n). So the check is an agreement and not an equality, and
 * 5e-4 is the tolerance the two items of the spec give: the furthest of these
 * eighteen numbers, a G''_ST of the multiallelic panel, is 4.7e-4 away.
 */
const MMOD_TOLERANCE = 5e-4;

/** The seven measures, in the order `PopDistMeasure` has them, all of which
 * popnei calculates. */
const EVERY_MEASURE = [
  "fst",
  "f2",
  "chord",
  "da",
  "dest",
  "gst",
  "gst_standardized",
];

/**
 * The distances between the populations of `pops` over the variants of the
 * VCF in `vcf`, with the `Variants` freed whatever happens.
 */
function popDistsOf(
  vcf: Uint8Array,
  pops: Record<string, readonly string[]>,
  askedFor: CalcPopDistsOptions,
  filter?: (variants: Variants) => void,
): PopDists {
  const variants = openVcf(vcf);
  try {
    filter?.(variants);
    return calcPopDists(variants, pops, askedFor);
  } finally {
    variants.free();
  }
}

/**
 * Every value of `found` within `tolerance` of the one beside it in
 * `expected`, relative to the expected value or absolute.
 */
function assertWithin(
  found: Float64Array | null | undefined,
  expected: readonly number[],
  tolerance: number,
  relative: boolean,
  what: string,
): void {
  assert.ok(found !== null && found !== undefined, what);
  assert.equal(found.length, expected.length, what);
  for (const [pair, theirs] of expected.entries()) {
    const ours = found[pair] as number;
    const apart = relative
      ? Math.abs(ours - theirs) / Math.abs(theirs)
      : Math.abs(ours - theirs);
    assert.ok(
      apart <= tolerance,
      `${what} of the pair ${PAIRS[pair] as string}: ${ours} and not ${theirs}`,
    );
  }
}

test("the F_ST of both panels is the one plink2 gives", () => {
  // The multiallelic panel is the one that says that every allele counts as
  // itself: the same formula on the major allele against the rest gives
  // 0.0707950 for p0 and p1 there, and on the reference allele against the
  // rest 0.0598407, neither of which is what plink2 printed.
  //
  // The two panels are cut into resampling groups in different ways, a
  // length in base pairs for the one that has positions along two
  // chromosomes and one group for each locus for the microsatellites, which
  // have no linkage to speak of. Neither changes an F_ST: the groups are
  // what the standard errors are built from.
  for (const { vcf, pops, expected, group, numVars, name } of [
    {
      vcf: PANEL_VCF,
      pops: PANEL_POPS,
      expected: PANEL_FST,
      group: PANEL_JACKKNIFE_GROUP,
      numVars: PANEL_NUM_VARS,
      name: "panel.vcf.gz",
    },
    {
      vcf: MICRO_VCF,
      pops: MICRO_POPS,
      expected: MICRO_FST,
      group: "variant" as const,
      numVars: MICRO_NUM_VARS,
      name: "micro.vcf.gz",
    },
  ]) {
    const dists = popDistsOf(vcf, pops, {
      jackknifeGroup: group,
      measures: ["fst"],
    });

    const ofThePanel = `the F_ST of ${name}`;
    assertWithin(
      dists.fst?.distVector,
      expected,
      FST_TOLERANCE,
      false,
      ofThePanel,
    );
    assert.deepEqual(dists.pops, ["p0", "p1", "p2"], ofThePanel);
    assert.deepEqual(dists.fst?.names, ["p0", "p1", "p2"], ofThePanel);
    assert.equal(dists.passStats.numVars, numVars, ofThePanel);
  }
});

test("the f_2 of the panel and its standard error are the ones ADMIXTOOLS gives", () => {
  // The run at 55 000 base pairs is the one of the three that is above the
  // 20 groups the function asks for, and the standard errors are what checks
  // the delete-m jackknife against a program that computes the same one.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: PANEL_JACKKNIFE_GROUP,
    measures: ["f2", "fst"],
  });

  assertWithin(
    dists.f2?.distVector,
    PANEL_F2,
    F2_TOLERANCE,
    true,
    "the f_2 of the panel",
  );
  assertWithin(
    dists.f2?.standardErrors,
    PANEL_F2_STANDARD_ERRORS,
    F2_TOLERANCE,
    true,
    "the standard error of the f_2 of the panel",
  );
  // The F_ST of the same call, which shares the pass and the groups, has a
  // standard error of the same jackknife, and no program prints one to
  // compare it with: what is asserted is that it is there and is a number of
  // the size of an F_ST of these populations.
  const ofTheFst = dists.fst?.standardErrors;
  assert.ok(ofTheFst !== null && ofTheFst !== undefined);
  assert.ok([...ofTheFst].every((error) => error > 0 && error < 0.01));
});

test("the dest of both panels is the one pyNei gives", () => {
  // pyNei is the one program that computes the estimator popnei computes,
  // the Nei and Chesser correction that GenAlEx prints, and
  // `tests/test_pop_dists.py` runs it against these same numbers within
  // 1e-12 relative. mmod's D, which the test below compares with, is another
  // estimator of the same quantity and is 7.3e-5 away on the biallelic panel
  // and 3.5e-4 on the multiallelic one.
  //
  // Both panels are read at the default of 20 called genotypes, where every
  // variant counts for every pair of both of them. The multiallelic one is
  // the one that says that the arithmetic does not assume two alleles: its
  // 120 loci have six alleles each.
  for (const { vcf, pops, expected, group, name } of [
    {
      vcf: PANEL_VCF,
      pops: PANEL_POPS,
      expected: PANEL_DEST,
      group: PANEL_JACKKNIFE_GROUP,
      name: "panel.vcf.gz",
    },
    {
      vcf: MICRO_VCF,
      pops: MICRO_POPS,
      expected: MICRO_DEST,
      group: "variant" as const,
      name: "micro.vcf.gz",
    },
  ]) {
    const dists = popDistsOf(vcf, pops, {
      jackknifeGroup: group,
      measures: ["dest"],
    });

    assertWithin(
      dists.dest?.distVector,
      expected,
      TEN_DIGITS_TOLERANCE,
      false,
      `Jost's D of ${name}`,
    );
  }
});

test("the dest where the pairs part is over each pair's own variants", () => {
  // At 47 called genotypes the two pairs of p0 lose the variants at which
  // p0, of 48 individuals, has fewer than 47 called. Each pair is a mean over
  // its own variants: 688 for p0-p1 and p0-p2 and all 1200 for p1-p2, whose D
  // is therefore the one at the default of 20 and does not move.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["dest"],
    minNumIndividuals: PARTING_MIN_NUM_INDIVIDUALS,
  });

  assertWithin(
    dists.dest?.distVector,
    PANEL_DEST_AT_THE_PARTING,
    TEN_DIGITS_TOLERANCE,
    false,
    "Jost's D of the panel at 47 called genotypes",
  );
  assert.deepEqual([...dists.numVars], PANEL_NUM_VARS_OF_EACH_PAIR);
});

test("the dest, the gst and the gst_standardized of both panels agree with mmod", () => {
  // mmod computes another estimator of each of the three, so this is an
  // agreement within 5e-4 and not an equality, and it says that popnei
  // computes these three quantities and not other statistics. What pins the
  // estimator is the comparison with pyNei above, which
  // `tests/test_pop_dists.py` makes exact to 1e-12 relative. Tightening this
  // tolerance, or moving the arithmetic towards mmod, breaks that one.
  //
  // The same call also gives the ten digits the spec prints for popnei's own
  // G_ST and G''_ST, which no program outside popnei computes and which are
  // here so that a change to the arithmetic that stays inside 5e-4 of mmod is
  // still caught.
  for (const { vcf, pops, ofMmod, ofPopnei, name } of [
    {
      vcf: PANEL_VCF,
      pops: PANEL_POPS,
      ofMmod: PANEL_MMOD,
      ofPopnei: { gst: PANEL_GST, gst_standardized: PANEL_GST_STANDARDIZED },
      name: "panel.vcf.gz",
    },
    {
      vcf: MICRO_VCF,
      pops: MICRO_POPS,
      ofMmod: MICRO_MMOD,
      ofPopnei: { gst: MICRO_GST, gst_standardized: MICRO_GST_STANDARDIZED },
      name: "micro.vcf.gz",
    },
  ]) {
    const dists = popDistsOf(vcf, pops, {
      jackknifeGroup: null,
      measures: ["dest", "gst", "gst_standardized"],
    });
    const ofTheResult = {
      dest: dists.dest,
      gst: dists.gst,
      gst_standardized: dists.gstStandardized,
    };

    for (const [measure, expected] of Object.entries(ofMmod)) {
      assertWithin(
        ofTheResult[measure as keyof typeof ofTheResult]?.distVector,
        expected,
        MMOD_TOLERANCE,
        false,
        `the ${measure} of ${name} against mmod`,
      );
    }
    for (const [measure, expected] of Object.entries(ofPopnei)) {
      assertWithin(
        ofTheResult[measure as keyof typeof ofTheResult]?.distVector,
        expected,
        TEN_DIGITS_TOLERANCE,
        false,
        `the ${measure} of ${name}`,
      );
    }
  }
});

test("the groups the variants fell into carry their chromosome and positions", () => {
  // The 22 groups of the biallelic panel: 11 on each of its two chromosomes,
  // the first of each holding the positions 1000 to 55000, which are 55
  // variants of the 1000 apart the panel has them at. A group is anchored on
  // its own first variant, so the first group of chr2 starts at that
  // chromosome's first position and not where the last group of chr1 ended.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: PANEL_JACKKNIFE_GROUP,
    measures: ["f2"],
  });

  assert.equal(dists.groupIds.length, PANEL_NUM_GROUPS);
  assert.equal(dists.numGroups, PANEL_NUM_GROUPS);
  assert.deepEqual(dists.groupIds[0], { chrom: "chr1", start: 1000, end: 55000 });
  assert.deepEqual(dists.groupIds[10], {
    chrom: "chr1",
    start: 551000,
    end: 600000,
  });
  assert.deepEqual(dists.groupIds[11], {
    chrom: "chr2",
    start: 1000,
    end: 55000,
  });
  assert.deepEqual(dists.groupIds[21], {
    chrom: "chr2",
    start: 551000,
    end: 600000,
  });
});

test("f2Groups holds the f_2 of every pair within every group", () => {
  // It is 22 groups by 3 pairs. The f_2 of a pair over the whole panel is
  // the mean of its variants, so the f_2 of each group weighted by the
  // variants of that group, 55 in twenty of them and 50 in the other two,
  // comes back to it: that is what ties the table to the numbers ADMIXTOOLS
  // gave.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: PANEL_JACKKNIFE_GROUP,
    measures: ["f2"],
  });

  const f2Groups = dists.f2Groups;
  assert.ok(f2Groups !== null);
  assert.equal(dists.numPairs, 3);
  assert.equal(f2Groups.length, PANEL_NUM_GROUPS * 3);
  for (const [pair, overThePanel] of PANEL_F2.entries()) {
    let weighted = 0;
    for (let group = 0; group < PANEL_NUM_GROUPS; group += 1) {
      const ofTheGroup = f2Groups[group * 3 + pair] as number;
      assert.ok(
        Number.isFinite(ofTheGroup),
        `${PAIRS[pair] as string} in the group ${group}`,
      );
      weighted += ofTheGroup * (PANEL_VARS_OF_EACH_GROUP[group] as number);
    }
    weighted /= PANEL_NUM_VARS;
    assert.ok(
      Math.abs(weighted - overThePanel) / overThePanel < F2_TOLERANCE,
      `${PAIRS[pair] as string}: ${weighted} and not ${overThePanel}`,
    );
  }
});

test("no jackknifeGroup gives no standard error and no group", () => {
  // The measures are the same numbers as with the groups, since the groups
  // are only what the standard errors are resampled over.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["f2", "fst"],
  });

  assertWithin(
    dists.f2?.distVector,
    PANEL_F2,
    F2_TOLERANCE,
    true,
    "the f_2 of the panel",
  );
  assert.equal(dists.f2?.standardErrors, null);
  assert.equal(dists.fst?.standardErrors, null);
  assert.equal(dists.f2?.squareStandardErrors(), null);
  assert.equal(dists.f2Groups, null);
  assert.equal(dists.numGroups, 0);
  assert.deepEqual(dists.groupIds, []);
});

test("a call with no jackknifeGroup is refused and names the argument", () => {
  // The length of the resampling groups has no default and the call fails
  // without it: the right length depends on the linkage disequilibrium of
  // the populations being compared, which popnei cannot know. TypeScript
  // refuses the call at build time, and a web application may call the same
  // function from JavaScript, which has no declarations, so it is refused at
  // run time too.
  const fromUntypedJavaScript = calcPopDists as unknown as (
    variants: Variants,
    pops: Record<string, readonly string[]>,
    options?: unknown,
  ) => PopDists;
  const variants = openVcf(PANEL_VCF);
  try {
    assert.throws(
      () => fromUntypedJavaScript(variants, PANEL_POPS, { measures: ["fst"] }),
      { message: /jackknifeGroup/ },
    );
    assert.throws(() => fromUntypedJavaScript(variants, PANEL_POPS), {
      message: /jackknifeGroup/,
    });
  } finally {
    variants.free();
  }
});

test("a jackknifeGroup that is no length and no word is refused", () => {
  // A length of 0 base pairs is no stretch of a chromosome: a user who wants
  // each variant in a group of its own writes `"variant"` and one who wants
  // no standard error writes `null`, and neither is a length. A word that is
  // not `"variant"` is refused the same way.
  for (const given of [0, -1, 2.5, "variants"]) {
    assert.throws(
      () =>
        popDistsOf(PANEL_VCF, PANEL_POPS, {
          jackknifeGroup: given as number,
          measures: ["fst"],
        }),
      { message: /jackknifeGroup/ },
      `the jackknifeGroup ${String(given)}`,
    );
  }
});

test("fewer than twenty groups are refused with how many there were", () => {
  // A length of 100 000 base pairs cuts the biallelic panel into 12 groups.
  // Each group is left out in turn, so a standard error built from a handful
  // of them says more about where the cuts fell than about the populations,
  // and a user who chose a length too long for their data is told how many
  // groups it gave, which is what says how much shorter to make them.
  assert.throws(
    () =>
      popDistsOf(PANEL_VCF, PANEL_POPS, {
        jackknifeGroup: 100000,
        measures: ["fst"],
      }),
    { message: /12.*20|20.*12/ },
  );
});

test("fewer than two populations are refused", () => {
  // Every one of the seven measures is of a pair of populations, so one
  // population makes no pair and there is nothing to give.
  assert.throws(
    () =>
      popDistsOf(
        PANEL_VCF,
        { p0: PANEL_POPS["p0"] as string[] },
        { jackknifeGroup: null, measures: ["fst"] },
      ),
    { message: /pops/ },
  );
});

test("pops with no population at all is refused with the rule of the call", () => {
  // `pops` is a required argument of `calcPopDists`, so the advice the
  // statistics give for an empty one, to leave it out for one population of
  // every individual, cannot be followed here: what the user has to do is
  // name two populations.
  assert.throws(
    () => popDistsOf(PANEL_VCF, {}, { jackknifeGroup: null, measures: ["fst"] }),
    { message: /two/ },
  );
});

test("a name that is not an individual of the pass is refused", () => {
  // The names are looked up among the individuals the pass gives, which are
  // those of the source after a filter of individuals when the `Variants`
  // has one, so only the pass knows them.
  assert.throws(
    () =>
      popDistsOf(
        PANEL_VCF,
        { p0: ["s000", "nobody"], p1: PANEL_POPS["p1"] as string[] },
        { jackknifeGroup: null, measures: ["fst"] },
      ),
    { message: /nobody/ },
  );
});

test("a source with no variant is refused and says that the source has none", () => {
  // A VCF whose header names three individuals and that has no data line.
  // The whole message is asserted, because it is the one the Python function
  // gives for the same source: the two languages say the same thing, and
  // Python writes the path of the file before it, which the bytes a
  // TypeScript user gives have not.
  const variants = openVcf(vcfOf([]));
  try {
    assert.throws(
      () =>
        calcPopDists(
          variants,
          { one: ["ind1"], two: ["ind2", "ind3"] },
          { jackknifeGroup: null, measures: ["fst"], minNumIndividuals: 1 },
        ),
      {
        message:
          "the pass gave no variant and its source holds none: a statistic " +
          "of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    variants.free();
  }
});

test("steps that keep no variant are refused with the counts of each filter", () => {
  // Four variants, each with one genotype of the three missing, and a filter
  // of missing data that keeps the variants with no missing genotype. The
  // counts of a pass that could not finish are otherwise lost, so the
  // message carries them.
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t./.",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
      "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t0/1\t1/1",
      "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
    ]),
  );
  try {
    variants.filterByMissingData(0);
    assert.throws(
      () =>
        calcPopDists(
          variants,
          { one: ["ind1"], two: ["ind2", "ind3"] },
          { jackknifeGroup: null, measures: ["fst"], minNumIndividuals: 1 },
        ),
      {
        message:
          "the pass gave no variant: its source gave 4 and the steps kept " +
          "none of them, the `missing_data` filter was given 4 and kept 0; " +
          "a statistic of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    variants.free();
  }
});

test("the variants that count are counted for each pair on its own", () => {
  // The biallelic panel at a `minNumIndividuals` of 47, which is above the
  // 42 called genotypes that p0, of 48 individuals, has at its emptiest
  // variant and below the 60 and 76 of p1 and p2: a variant p0 is short of
  // is lost by the two pairs p0 is in and kept by p1-p2. At the default of
  // 20 every variant counts for every pair, which is the other half of the
  // test: a fixture in which the counts are always the same number cannot
  // tell a count per pair from one count for all of them.
  const parted = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["fst"],
    minNumIndividuals: PARTING_MIN_NUM_INDIVIDUALS,
  });
  assert.deepEqual([...parted.numVars], PANEL_NUM_VARS_OF_EACH_PAIR);
  assert.ok(parted.numVars instanceof Int32Array);

  const atTheDefault = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["fst"],
  });
  assert.deepEqual(
    [...atTheDefault.numVars],
    [PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS],
  );
  // The F_ST of p1-p2 is over the same 1200 variants at both thresholds and
  // is the same number; the two pairs that lost variants moved.
  assert.equal(
    atTheDefault.fst?.distVector[2],
    parted.fst?.distVector[2],
  );
  assert.notEqual(
    atTheDefault.fst?.distVector[0],
    parted.fst?.distVector[0],
  );
});

test("the pairs are in the order the populations were named in", () => {
  // The populations stay in the order of the keys of `pops`, which is where
  // popnei parts from pyNei: pyNei sorts their names. No value changes with
  // the order, so named p2, p1, p0 the pairs are p2-p1, p2-p0 and p1-p0, and
  // their F_ST is that of p1-p2, of p0-p2 and of p0-p1.
  const turnedAround = {
    p2: PANEL_POPS["p2"] as string[],
    p1: PANEL_POPS["p1"] as string[],
    p0: PANEL_POPS["p0"] as string[],
  };

  const dists = popDistsOf(PANEL_VCF, turnedAround, {
    jackknifeGroup: null,
    measures: ["fst"],
  });

  assert.deepEqual(dists.pops, ["p2", "p1", "p0"]);
  assertWithin(
    dists.fst?.distVector,
    [PANEL_FST[2] as number, PANEL_FST[1] as number, PANEL_FST[0] as number],
    FST_TOLERANCE,
    false,
    "the F_ST of the panel named the other way round",
  );
});

test("the counts of the pass hold the variants and the filters", () => {
  // The filter keeps the variants whose major allele frequency is below
  // 0.95, and the count of the variants of the pass is what it kept.
  const dists = popDistsOf(
    PANEL_VCF,
    PANEL_POPS,
    { jackknifeGroup: null, measures: ["fst"] },
    (variants) => {
      variants.filterByMaf(0.95);
    },
  );

  const kept = dists.passStats.filtering["maf"];
  assert.ok(kept !== undefined);
  assert.equal(kept.varsProcessed, PANEL_NUM_VARS);
  assert.ok(kept.varsKept > 0 && kept.varsKept < PANEL_NUM_VARS);
  assert.equal(dists.passStats.numVars, kept.varsKept);
});

test("a measure that was not asked for is null in the result", () => {
  // Asking for one measure gives the same number as asking for both: the
  // pass is what costs and each measure is a division at the end of it.
  const ofOne = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["f2"],
  });
  const ofBoth = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: null,
    measures: ["fst", "f2"],
  });

  assert.equal(ofOne.fst, null);
  assert.ok(ofOne.f2 !== null);
  assert.ok(ofBoth.f2 !== null);
  assert.ok(ofBoth.fst !== null);
  assert.deepEqual([...ofOne.f2.distVector], [...ofBoth.f2.distVector]);
  assert.deepEqual(
    [
      ofBoth.chord,
      ofBoth.da,
      ofBoth.dest,
      ofBoth.gst,
      ofBoth.gstStandardized,
    ],
    [null, null, null, null, null],
  );
});

test("no measure of the seven is refused", () => {
  // Every measure has a value, so no `measures`, which asks for all seven,
  // is taken and the refusal of one popnei has no value for has nothing
  // left to refuse. The list of the ones that have a value is the core's,
  // and a measure added to it without a formula beside it would be refused
  // here rather than handed to a user as a vector of NaN.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, { jackknifeGroup: null });
  const ofEveryMeasure = [
    dists.fst,
    dists.f2,
    dists.chord,
    dists.da,
    dists.dest,
    dists.gst,
    dists.gstStandardized,
  ];

  assert.equal(ofEveryMeasure.length, EVERY_MEASURE.length);
  ofEveryMeasure.forEach((values, measure) => {
    assert.ok(values !== null, `the ${EVERY_MEASURE[measure]} of the panel`);
  });
});

test("a measure that is of none of the seven is refused with the seven", () => {
  assert.throws(
    () =>
      popDistsOf(PANEL_VCF, PANEL_POPS, {
        jackknifeGroup: null,
        measures: ["fsts" as PopDistMeasure],
      }),
    { message: /fst/ },
  );
  assert.throws(
    () =>
      popDistsOf(PANEL_VCF, PANEL_POPS, {
        jackknifeGroup: null,
        measures: [] as PopDistMeasure[],
      }),
    { message: /measures/ },
  );
  assert.throws(
    () =>
      popDistsOf(PANEL_VCF, PANEL_POPS, {
        jackknifeGroup: null,
        measures: 2 as unknown as PopDistMeasure[],
      }),
    { message: /measures/ },
  );
});

test("squareStandardErrors is the square matrix of the vector", () => {
  // The standard error of a pair is in both of its cells and the diagonal is
  // NaN, where the square matrix of the distances has 0: a distance of a
  // population with itself is 0 and known, and how far that 0 would move is
  // nothing the calculation gives.
  const dists = popDistsOf(PANEL_VCF, PANEL_POPS, {
    jackknifeGroup: PANEL_JACKKNIFE_GROUP,
    measures: ["f2"],
  });

  const square = dists.f2?.squareStandardErrors();
  assert.ok(square !== null && square !== undefined);
  assert.equal(square.length, 9);
  assert.equal(square[0 * 3 + 1], square[1 * 3 + 0]);
  assert.ok(
    Math.abs(
      ((square[0 * 3 + 1] as number) - (PANEL_F2_STANDARD_ERRORS[0] as number)) /
        (PANEL_F2_STANDARD_ERRORS[0] as number),
    ) <= F2_TOLERANCE,
  );
  for (const place of [0, 1, 2]) {
    assert.ok(Number.isNaN(square[place * 3 + place]));
  }
});

test("the Kosman distances between individuals have no standard errors", () => {
  // `standardErrors` is the field the distances between populations added,
  // and the Kosman distances leave it null: nothing of theirs changes.
  const variants = openVcf(MICRO_VCF);
  try {
    const dists = calcPairwiseKosmanDists(variants);
    assert.equal(dists.standardErrors, null);
    assert.equal(dists.squareStandardErrors(), null);
  } finally {
    variants.free();
  }
});

test("the variants are as they were after the call", () => {
  // The call is a consumer: it makes one pass over the source through the
  // steps the `Variants` has and puts nothing on it, so a second call reads
  // the source again and gives the same numbers.
  const variants = openVcf(PANEL_VCF);
  try {
    const first = calcPopDists(variants, PANEL_POPS, {
      jackknifeGroup: null,
      measures: ["fst"],
    });
    assert.deepEqual(variants.steps, []);
    const second = calcPopDists(variants, PANEL_POPS, {
      jackknifeGroup: null,
      measures: ["fst"],
    });

    assert.deepEqual(
      [...(second.fst as { distVector: Float64Array }).distVector],
      [...(first.fst as { distVector: Float64Array }).distVector],
    );
  } finally {
    variants.free();
  }
});
