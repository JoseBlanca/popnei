/**
 * The association study from TypeScript: `calcGwas` and the result it gives.
 *
 * "The linear model" and "The worked example" of `docs/specs/gwas.md` have
 * the numbers. The worked example is 3 variants of 6 diploid individuals
 * with one covariate, written as a VCF here, and pyNei at commit ef0ca6e
 * gave its null model and its three rows; it reads no reference file and
 * nothing of it is rounded away. The panel is
 * `tests/reference/kinship/panel_called.vcf.gz`, 200 individuals and 1200
 * biallelic variants with every genotype called, with the trait `cont` and
 * the covariates `cov1` and `cov2` of `tests/reference/gwas/phenotypes.csv`,
 * and the six variants asserted here are what plink2 v2.0.0-a.7.7 wrote for
 * it.
 *
 * The calculation is tested in the core crate, over all 1200 variants of the
 * panel. What these tests say is that the study reaches TypeScript with the
 * numbers it has in the core, that the individuals that are tested are the
 * ones that have a phenotype, and that what popnei refuses is thrown as an
 * `Error` with the message it has in Rust.
 *
 * Nothing here computes an expected value with popnei: every number comes
 * from the spec, which has them from plink2 and from pyNei.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { calcGwas, init, openVcf } from "popnei";
import type { GwasResult, Variants } from "popnei";

import { referenceGwas, referenceKinship } from "./reference.ts";

await init();

/**
 * How far a number of the worked example may be from pyNei's, as a share of
 * that number: 7e-15.
 *
 * "The worked example" of `docs/specs/gwas.md` asks for 1e-12, and "How it
 * is verified" of "What every model shares" asks that such a number be
 * lowered until it fails and then set two or three times above where it
 * broke. It was, under node on 23 September 2026: it breaks at 2e-15, where
 * the `beta` of `v1` is 2.31e-15 away, and that is the worst of the twelve
 * numbers it holds there, the three of the null model and the nine of the
 * three rows. So this is 3 times the worst measured. It holds one more
 * number in another test, the frequency of `v0` over five individuals,
 * which is 0.4 exactly in both libraries.
 *
 * The cargo test of the same twelve numbers measures 1.24e-15 on faer and
 * 7.4e-16 on Accelerate and is set at 3e-15. This build is faer as well,
 * since that is the linear algebra of wasm, and it sits 1.9 times further
 * from pyNei than faer natively does: a wasm module has no instruction that
 * adds a product to a sum in one rounding, which the native build uses where
 * the processor has it, so the two round the same sums differently.
 *
 * No number of the example is near 0: the smallest is the `beta` of `v1`,
 * 0.3125, against an `se` of 1.3, so a share of the number and a share of
 * the scale of what is estimated are the same bound here, and the rule of
 * the spec that a `beta` is measured against its `se` changes nothing.
 */
const OF_THE_WORKED_EXAMPLE = 7e-15;

/**
 * How far a `beta` or an `se` of the panel may be from plink2's, as a share
 * of the `se` plink2 printed for that variant: 1e-5.
 *
 * It is what "How it is verified" of "The linear model" of the spec asks,
 * and it is a share of the `se` and not of the value, because a study is
 * mostly null and a `beta` that cancelled to near 0 is no guide to its own
 * error. plink2 prints six significant digits, so the number it is compared
 * against is itself rounded by up to 5e-7 absolute, which of the 1.2e-6 this
 * bound comes to at the smallest `se` of the six is 41 per cent. Measured
 * under node on 23 September 2026, the worst of the six is 3.74e-6 of the
 * `se` of its variant, the `se` of `var0052`, and the worst p-value is
 * 2.82e-6 of itself, at `var0629`; the cargo test of the core measured the
 * same two to two digits on faer and on Accelerate, and in all three every
 * `allele_freq` is plink2's `A1_FREQ` exactly.
 *
 * None of the six has a `beta` above 1, where six significant digits round a
 * value by up to 5e-6 absolute instead of 5e-7 and the printing alone can
 * pass this bound.
 */
const OF_PLINK2 = 1e-5;

/**
 * How far a p-value of the panel may be from plink2's, as a share of it.
 * plink2 prints it to six significant digits like the rest, so the printing
 * alone can take half of this.
 */
const OF_PLINK2_P_VALUE = 1e-5;

/**
 * How far an `allele_freq` of the panel may be from plink2's `A1_FREQ`: 1e-6
 * absolute, since it is a frequency and lies between 0 and 1.
 */
const OF_PLINK2_FREQUENCY = 1e-6;

/** How many variants the panel holds, and how many individuals. */
const PANEL_NUM_VARS = 1200;
const PANEL_NUM_INDIVIDUALS = 200;

/**
 * The worked example of the spec as a VCF: six diploid individuals, `i0` to
 * `i5`, and three variants. The `./.` of `i3` at `v1` is the genotype that
 * takes the mean dosage of its variant, 0.8, and `v2`, where every
 * individual is heterozygous, has no variance and no answer.
 */
const WORKED_EXAMPLE = vcfOf(
  ["i0", "i1", "i2", "i3", "i4", "i5"],
  [
    { id: "v0", genotypes: "0/0\t0/1\t1/1\t0/0\t0/1\t1/1" },
    { id: "v1", genotypes: "0/0\t0/1\t1/1\t./.\t0/1\t0/0" },
    { id: "v2", genotypes: "0/1\t0/1\t0/1\t0/1\t0/1\t0/1" },
  ],
);

/** The trait of the worked example, one measurement for each individual. */
const THE_TRAIT: Record<string, number> = {
  i0: 2,
  i1: 3,
  i2: 5,
  i3: 4,
  i4: 4,
  i5: 7,
};

/** Its one covariate beside the intercept. */
const THE_COVARIATE: Record<string, Record<string, number>> = {
  cov: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
};

/**
 * The null model of the worked example, which pyNei gave at commit ef0ca6e:
 * an intercept of 3.6666666666666683, an effect of the covariate of 1.0 and
 * a residual variance of 3.333333333333334, which is its residual sum of
 * squares of 13.333333333333336 over the 4 degrees of freedom of 6
 * individuals and 2 coefficients.
 */
const THE_NULL = {
  intercept: 3.666_666_666_666_668_3,
  cov: 1.0,
  residualVariance: 3.333_333_333_333_334,
};

/**
 * The three rows of the worked example, which pyNei gave: `v2` keeps its
 * frequency and has no test, since every individual is heterozygous there.
 */
const THE_ROWS: {
  id: string;
  alleleFreq: number;
  beta: number;
  se: number;
  pValue: number;
}[] = [
  {
    id: "v0",
    alleleFreq: 0.5,
    beta: 1.5,
    se: 0.600_925_212_577_332,
    pValue: 0.088_004_892_382_756,
  },
  {
    id: "v1",
    alleleFreq: 0.4,
    beta: 0.3125,
    se: 1.305_204_592_306_424,
    pValue: 0.826_200_867_452_417,
  },
  {
    id: "v2",
    alleleFreq: 0.5,
    beta: Number.NaN,
    se: Number.NaN,
    pValue: Number.NaN,
  },
];

/**
 * What plink2 v2.0.0-a.7.7 wrote for six variants of the panel in
 * `tests/reference/gwas/plink2.panel_called.glm.linear.tsv`, with `cov1` and
 * `cov2` as covariates: `A1_FREQ`, `BETA`, `SE` and `P`, the four columns the
 * spec compares. Five of them are the causal variants of `causal_vars.csv`
 * and `var0000` is not causal.
 */
const OF_PLINK2_SIX: {
  id: string;
  alleleFreq: number;
  beta: number;
  se: number;
  pValue: number;
}[] = [
  {
    id: "var0000",
    alleleFreq: 0.21,
    beta: -0.424136,
    se: 0.139354,
    pValue: 0.00265846,
  },
  {
    id: "var0052",
    alleleFreq: 0.3225,
    beta: -0.697724,
    se: 0.122348,
    pValue: 4.28981e-8,
  },
  {
    id: "var0629",
    alleleFreq: 0.175,
    beta: -0.813852,
    se: 0.161809,
    pValue: 1.10646e-6,
  },
  {
    id: "var0751",
    alleleFreq: 0.41,
    beta: -0.0963977,
    se: 0.130636,
    pValue: 0.461451,
  },
  {
    id: "var1137",
    alleleFreq: 0.2625,
    beta: -0.171137,
    se: 0.145987,
    pValue: 0.242511,
  },
  {
    id: "var1188",
    alleleFreq: 0.4675,
    beta: -0.655393,
    se: 0.138912,
    pValue: 4.51958e-6,
  },
];

/** The bytes of the panel and its phenotypes, read once for every test. */
const PANEL_VCF = await referenceKinship("panel_called.vcf.gz");
const PHENOTYPES = theColumnsOfTheFile(await referenceGwas("phenotypes.csv"));

/**
 * A VCF of the individuals `names` with one data line for each of
 * `variants`, every variant declaring the two alleles `A` and `C` and
 * carrying its id and a position of its own.
 */
function vcfOf(
  names: readonly string[],
  variants: readonly { id: string; genotypes: string }[],
): Uint8Array {
  const lines = variants.map(
    ({ id, genotypes }, at) =>
      `chr1\t${(at + 1) * 1000}\t${id}\tA\tC\t.\tPASS\t.\tGT\t${genotypes}`,
  );
  const header = [
    "##fileformat=VCFv4.4",
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
  ];
  return new TextEncoder().encode([...header, ...lines, ""].join("\n"));
}

/**
 * The columns of `phenotypes.csv` as one object of the name of an individual
 * to its value for each of them: `IID` names the rows and `cont`, `binom`,
 * `cov1`, `cov2` and `pop` are the columns.
 *
 * @throws {Error} When a line does not hold every column, or when a value is
 * not a number: the study would then be of something else.
 */
function theColumnsOfTheFile(
  text: string,
): Record<string, Record<string, number>> {
  const lines = text.split("\n").filter((line) => line.trim() !== "");
  const header = (lines[0] as string).split(",");
  const columns: Record<string, Record<string, number>> = {};
  for (const name of header.slice(1)) {
    columns[name] = {};
  }
  for (const line of lines.slice(1)) {
    const fields = line.split(",");
    if (fields.length !== header.length) {
      throw new Error(`phenotypes.csv: the line \`${line}\` has no value`);
    }
    const individual = fields[0] as string;
    for (const [at, name] of header.slice(1).entries()) {
      const value = Number(fields[at + 1]);
      if (!Number.isFinite(value)) {
        throw new Error(
          `phenotypes.csv: the ${name} of ${individual} is \`${fields[at + 1]}\``,
        );
      }
      (columns[name] as Record<string, number>)[individual] = value;
    }
  }
  return columns;
}

/** The study of `bytes` with `options` as a user would give them. */
function gwasOf(
  bytes: Uint8Array,
  options: Parameters<typeof calcGwas>[1],
): GwasResult {
  const variants: Variants = openVcf(bytes, { onlyPassed: false });
  try {
    return calcGwas(variants, options);
  } finally {
    variants.free();
  }
}

/** The study of the worked example, with `phenotype` as a test wants it. */
function theWorkedExample(
  phenotype: Record<string, number> = THE_TRAIT,
): GwasResult {
  return gwasOf(WORKED_EXAMPLE, {
    phenotype,
    trait: "continuous",
    covariates: THE_COVARIATE,
  });
}

/** The study of the panel with `cov1` and `cov2`, which plink2 was given. */
function theStudyOfThePanel(): GwasResult {
  return gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.cont as Record<string, number>,
    trait: "continuous",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
  });
}

/** The row of the variant `id` in the result of a study. */
function rowOf(result: GwasResult, id: string): number {
  const at = result.stats.id?.indexOf(id) ?? -1;
  assert.ok(at >= 0, `${id} is not a variant of the study`);
  return at;
}

/** `found` within `tolerance` of `expected` as a share of `expected`. */
function assertWithin(
  found: number,
  expected: number,
  tolerance: number,
  what: string,
): void {
  const difference = Math.abs(found - expected);
  assert.ok(
    difference <= tolerance * Math.abs(expected),
    `${what} is ${found} and the literal is ${expected}, ${difference} away, ` +
      `which is ${difference / Math.abs(expected)} of it against the ` +
      `${tolerance} allowed`,
  );
}

/**
 * `found` within `tolerance` times `scale`, which for a `beta` and for an
 * `se` is the `se` of that variant.
 */
function assertWithinTheScale(
  found: number,
  expected: number,
  scale: number,
  tolerance: number,
  what: string,
): void {
  const difference = Math.abs(found - expected);
  assert.ok(
    difference <= tolerance * scale,
    `${what} is ${found} and plink2 gives ${expected}, ${difference} away, ` +
      `which is ${difference / scale} of the ${scale} it is uncertain by ` +
      `against the ${tolerance} allowed`,
  );
}

test("the worked example gives pyNei's null model", () => {
  const result = theWorkedExample();

  assert.equal(result.nullModel.model, "lm");
  assert.equal(result.test, "wald");
  assert.equal(result.trait, "continuous");
  assert.equal(result.nullModel.numIndividuals, 6);
  assert.equal(result.usedGrammarGammaApprox, false);
  assert.deepEqual(result.individuals, ["i0", "i1", "i2", "i3", "i4", "i5"]);
  assert.deepEqual(Object.keys(result.nullModel.covariateEffects), [
    "intercept",
    "cov",
  ]);
  assertWithin(
    result.nullModel.covariateEffects["intercept"] as number,
    THE_NULL.intercept,
    OF_THE_WORKED_EXAMPLE,
    "the intercept",
  );
  assertWithin(
    result.nullModel.covariateEffects["cov"] as number,
    THE_NULL.cov,
    OF_THE_WORKED_EXAMPLE,
    "the effect of the covariate",
  );
  assertWithin(
    result.nullModel.residualVariance as number,
    THE_NULL.residualVariance,
    OF_THE_WORKED_EXAMPLE,
    "the residual variance",
  );
  assert.equal(result.nullModel.geneticVariance, undefined);
  assert.equal(result.nullModel.heritability, undefined);
});

test("the worked example gives pyNei's three rows and no answer for the third", () => {
  const result = theWorkedExample();

  assert.deepEqual(result.stats.id, ["v0", "v1", "v2"]);
  assert.deepEqual(result.stats.chrom, ["chr1", "chr1", "chr1"]);
  assert.deepEqual([...(result.stats.pos as Float64Array)], [1000, 2000, 3000]);
  assert.equal(result.passStats.numVars, 3);
  assert.deepEqual(result.passStats.filtering, {});
  for (const { id, alleleFreq, beta, se, pValue } of THE_ROWS) {
    const at = rowOf(result, id);
    assertWithin(
      result.stats.alleleFreq[at] as number,
      alleleFreq,
      OF_THE_WORKED_EXAMPLE,
      `the frequency of ${id}`,
    );
    if (Number.isNaN(beta)) {
      assert.ok(Number.isNaN(result.stats.beta[at] as number), `beta of ${id}`);
      assert.ok(Number.isNaN(result.stats.se[at] as number), `se of ${id}`);
      assert.ok(
        Number.isNaN(result.stats.pValue[at] as number),
        `the p-value of ${id}`,
      );
      continue;
    }
    assertWithin(
      result.stats.beta[at] as number,
      beta,
      OF_THE_WORKED_EXAMPLE,
      `the effect of ${id}`,
    );
    assertWithin(
      result.stats.se[at] as number,
      se,
      OF_THE_WORKED_EXAMPLE,
      `the standard error of ${id}`,
    );
    assertWithin(
      result.stats.pValue[at] as number,
      pValue,
      OF_THE_WORKED_EXAMPLE,
      `the p-value of ${id}`,
    );
  }
});

test("the six variants of the panel are plink2's effect, error and p-value", () => {
  const result = theStudyOfThePanel();

  assert.equal(result.nullModel.numIndividuals, PANEL_NUM_INDIVIDUALS);
  assert.equal(result.passStats.numVars, PANEL_NUM_VARS);
  assert.equal(result.stats.beta.length, PANEL_NUM_VARS);
  assert.deepEqual(Object.keys(result.nullModel.covariateEffects), [
    "intercept",
    "cov1",
    "cov2",
  ]);
  for (const { id, alleleFreq, beta, se, pValue } of OF_PLINK2_SIX) {
    const at = rowOf(result, id);
    const frequency = result.stats.alleleFreq[at] as number;
    assert.ok(
      Math.abs(frequency - alleleFreq) <= OF_PLINK2_FREQUENCY,
      `the frequency of ${id} is ${frequency} and plink2 gives ${alleleFreq}`,
    );
    assertWithinTheScale(
      result.stats.beta[at] as number,
      beta,
      se,
      OF_PLINK2,
      `the effect of ${id}`,
    );
    assertWithinTheScale(
      result.stats.se[at] as number,
      se,
      se,
      OF_PLINK2,
      `the standard error of ${id}`,
    );
    assertWithin(
      result.stats.pValue[at] as number,
      pValue,
      OF_PLINK2_P_VALUE,
      `the p-value of ${id}`,
    );
  }
});

test("an individual with no phenotype is not tested and the frequencies are of the rest", () => {
  const { i5: _left_out, ...withoutI5 } = THE_TRAIT;

  const result = theWorkedExample(withoutI5);

  assert.deepEqual(result.individuals, ["i0", "i1", "i2", "i3", "i4"]);
  assert.equal(result.nullModel.numIndividuals, 5);
  // The dosages of `v0` over the five are 0, 1, 2, 0, 1, whose mean is 0.8,
  // and `allele_freq` is that mean over the ploidy: 0.4 where the six give
  // 0.5. The frequency of a variant is of the individuals that are tested,
  // which is what "What it gives" of the spec says of the mean, the
  // frequency and whether a variant varies at all.
  assertWithin(
    result.stats.alleleFreq[rowOf(result, "v0")] as number,
    0.4,
    OF_THE_WORKED_EXAMPLE,
    "the frequency of v0 over the five",
  );
});

test("a binomial trait is refused with the model that is being written", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
        trait: "binomial",
        covariates: THE_COVARIATE,
      }),
    { message: /logistic regression, which is being written/ },
  );
});

test("a trait of another name is refused with the two names", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "quantitative" as "continuous",
      }),
    { message: /`continuous`.*`binomial`.*`quantitative`/s },
  );
});

test("an individual of the phenotype that the dataset has not is refused by its name", () => {
  assert.throws(
    () => theWorkedExample({ ...THE_TRAIT, i9: 5 }),
    { message: /`i9` has a phenotype/ },
  );
});

test("a covariate that does not cover a tested individual is refused", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: { cov: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0 } },
      }),
    { message: /the covariate `cov` has no value for `i5`/ },
  );
});

test("a covariate that is a copy of another is refused as collinear", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: {
          cov: THE_COVARIATE["cov"] as Record<string, number>,
          twice: THE_COVARIATE["cov"] as Record<string, number>,
        },
      }),
    { message: /3 columns, the intercept among them, and only 2 of them are/ },
  );
});

test("the options of the linear mixed model are refused by name", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        kinship: undefined as unknown as Parameters<
          typeof calcGwas
        >[1]["kinship"],
        useGrammarGammaApprox: true,
      }),
    { message: /`useGrammarGammaApprox` belongs to the linear mixed model/ },
  );
});

test("an option of the mixed model written as undefined is not given", () => {
  // Spreading an object of options over a call leaves `undefined` for every
  // one that was not filled in, and a user who writes the documented
  // default explicitly is asking for nothing.
  const result = gwasOf(WORKED_EXAMPLE, {
    phenotype: THE_TRAIT,
    trait: "continuous",
    covariates: THE_COVARIATE,
    kinship: undefined,
    useGrammarGammaApprox: undefined,
  });

  assert.equal(result.nullModel.model, "lm");
  assert.equal(result.usedGrammarGammaApprox, false);
});

test("the wald test is the linear model's own and the score test is refused", () => {
  const asked = gwasOf(WORKED_EXAMPLE, {
    phenotype: THE_TRAIT,
    trait: "continuous",
    covariates: THE_COVARIATE,
    test: "wald",
  });
  const byDefault = theWorkedExample();

  // What a user reads in the result is what they may write back into the
  // call, which is the one thing a refusal of `test` broke.
  assert.equal(byDefault.test, "wald");
  assert.equal(asked.test, "wald");
  assert.deepEqual([...asked.stats.beta], [...byDefault.stats.beta]);
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: THE_COVARIATE,
        test: "score",
      }),
    { message: /only test is the t test of the effect it fitted/ },
  );
});

test("a test of another name is refused with the two names", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: THE_COVARIATE,
        test: "rao" as "wald",
      }),
    { message: /`wald`.*`score`.*`rao`/s },
  );
});

test("a study asked for with no options at all says what to write", () => {
  const variants: Variants = openVcf(WORKED_EXAMPLE, { onlyPassed: false });
  try {
    assert.throws(
      () => (calcGwas as (variants: Variants) => GwasResult)(variants),
      { message: /^popnei: a study is asked for with the trait/ },
    );
  } finally {
    variants.free();
  }
});

test("a covariate that is not finite is refused by its name", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: {
          cov: { ...THE_COVARIATE["cov"], i2: Number.POSITIVE_INFINITY },
        },
      }),
    { message: /the covariate `cov` at `i2` is Infinity/ },
  );
});

/**
 * The call of each case of `refusals_of_both_layers.json`, over the worked
 * example.
 *
 * The Python suite holds the same cases under the same names, and each
 * suite writes the call in its own language: the names of the options
 * differ between the two.
 */
function theCallsThatAreRefused(): Record<string, () => GwasResult> {
  const cov = THE_COVARIATE["cov"] as Record<string, number>;
  const { i5: _withoutI5, ...ofFive } = cov;
  const { i5: _alsoWithoutI5, ...ofThree } = THE_TRAIT;
  const study = (options: Record<string, unknown>) => () =>
    gwasOf(WORKED_EXAMPLE, {
      phenotype: THE_TRAIT,
      trait: "continuous",
      covariates: THE_COVARIATE,
      ...options,
    } as Parameters<typeof calcGwas>[1]);
  return {
    "a kinship": study({ kinship: "a matrix" }),
    "the grammar gamma approximation": study({ useGrammarGammaApprox: true }),
    "the score test": study({ test: "score" }),
    "a test of another name": study({ test: "rao" }),
    "a trait of another name": study({ trait: "quantitative" }),
    "a binomial trait": study({
      phenotype: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
      trait: "binomial",
    }),
    "a covariate named intercept": study({ covariates: { intercept: cov } }),
    "a covariate that has no value for a tested individual": study({
      covariates: { cov: ofFive },
    }),
    "a covariate that is not finite": study({
      covariates: { cov: { ...cov, i2: Number.POSITIVE_INFINITY } },
    }),
    "a covariate that is a copy of another": study({
      covariates: { cov, twice: cov },
    }),
    "an individual of the phenotype that the variants have not": study({
      phenotype: { ...THE_TRAIT, i9: 5 },
    }),
    "fewer individuals than the design has columns plus two": study({
      phenotype: { i0: 2, i1: 3, i2: 5 },
    }),
  };
}

test("both layers refuse the same calls", async () => {
  const listed = JSON.parse(
    await referenceGwas("refusals_of_both_layers.json"),
  ) as { refusals: { case: string; match: string }[] };
  const calls = theCallsThatAreRefused();

  assert.ok(listed.refusals.length > 0, "the file lists no refusal");
  for (const { case: name, match } of listed.refusals) {
    const call = calls[name];
    assert.ok(
      call !== undefined,
      `\`${name}\` is in refusals_of_both_layers.json and this suite has no ` +
        "call for it: a refusal both layers make is written in both",
    );
    assert.throws(call, new RegExp(match), name);
  }
  assert.deepEqual(
    Object.keys(calls).sort(),
    listed.refusals.map(({ case: name }) => name).sort(),
    "this suite makes a call that the file does not list",
  );
});

test("a covariate named intercept is refused with the name it collides with", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: { intercept: THE_COVARIATE["cov"] as Record<string, number> },
      }),
    { message: /a covariate is named `intercept`/ },
  );
});
