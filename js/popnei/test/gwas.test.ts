/**
 * The association study from TypeScript: `calcGwas` and the result it gives.
 *
 * "The linear model", "The linear mixed model" and "The worked example" of
 * `docs/specs/gwas.md` have the numbers. The worked example is 3 variants of
 * 6 diploid individuals with one covariate, written as a VCF here, and pyNei
 * at commit ef0ca6e gave its null model and its three rows; it reads no
 * reference file and nothing of it is rounded away. The panel is
 * `tests/reference/kinship/panel_called.vcf.gz`, 200 individuals and 1200
 * biallelic variants with every genotype called, with the trait `cont` and
 * the covariates `cov1` and `cov2` of `tests/reference/gwas/phenotypes.csv`.
 * The six variants asserted for the linear model are what plink2
 * v2.0.0-a.7.7 wrote for it, and the six of the mixed model are what GMMAT
 * 1.5.0 wrote, over the kinship that `plink2 --make-rel` wrote and that
 * neither popnei nor pyNei calculated.
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
import { gunzipSync } from "node:zlib";

import { calcGwas, init, Kinship, openVcf } from "popnei";
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

/**
 * How far `1 / se**2` of the score test may be from GMMAT's `VAR`, as a
 * share of it, and how far a p-value may be from GMMAT's in `log10`: 1e-5
 * and 1e-4, which is what "How it is verified" of "The linear mixed model"
 * of the spec asks of every variant and what the cargo tests hold the same
 * six literals to.
 *
 * `gmmat.panel_called.lmm.score.tsv` is printed to six significant digits,
 * which rounds a value by up to 5e-6 of itself, so half of the first bound
 * can go on GMMAT's printing alone. Measured under node on 24 September
 * 2026, the worst of the six is 1.837e-6 of `VAR`, at `var0052`, and the
 * worst p-value is 4.165e-5 in `log10`, at `var0629`; the cargo test of the
 * core measures 1.843e-6 and 4.164e-5 on faer natively and 1.842e-6 and
 * 4.164e-5 on Accelerate, so nothing of the distance from GMMAT is
 * WebAssembly's own rounding. That p-value is not the printing either:
 * 4.16e-5 in `log10` is 9.6e-5 of the p-value, where six digits round it by
 * 5e-6, and it is the two fits landing 1.2e-6 apart in the genetic variance
 * at a p-value of 4.8e-5, where the tail of the chi square turns a small
 * move of the statistic into a larger one of the p-value.
 */
const OF_GMMAT_VARIANCE = 1e-5;
const OF_GMMAT_P_VALUE = 1e-4;

/**
 * How far each of the two variances of the null model may be from GMMAT's
 * `glmmkin`: 1e-5 absolute, which is the spec's and is how far two
 * restricted maximum likelihood searches land apart.
 *
 * Measured under node on 24 September 2026, the worst of the three numbers
 * this test holds to it is the genetic variance, 1.273e-6 from GMMAT's,
 * which is 13 per cent of what is allowed; the cargo test measures 1.218e-6
 * on faer natively. Where that 1.2e-6 comes from is the kinship and not the
 * search: GMMAT was given the six printed digits of plink2's matrix and
 * this suite reads the float64 beside them.
 */
const OF_GMMAT_NULL_MODEL = 1e-5;

/**
 * What GMMAT 1.5.0's `glmm.score` gave for six variants of the panel with
 * every genotype called, from `tests/reference/gwas/gmmat.panel_called.lmm.
 * score.tsv`: the variance of the score, which is `x' p x` and which popnei
 * gives as `1 / se**2`, and the p-value.
 *
 * GMMAT was given both covariates and the kinship plink2 wrote for this
 * panel, which is what this suite gives popnei. Five of the six are the
 * causal variants of `causal_vars.csv` and `var0000` is not causal.
 */
const OF_GMMAT_SIX: { id: string; variance: number; pValue: number }[] = [
  { id: "var0000", variance: 29.8774, pValue: 0.360_526 },
  { id: "var0052", variance: 43.8076, pValue: 0.001_189_85 },
  { id: "var0629", variance: 31.7241, pValue: 4.810_05e-5 },
  { id: "var0751", variance: 44.5825, pValue: 0.004_392_26 },
  { id: "var1137", variance: 43.3724, pValue: 0.013_926 },
  { id: "var1188", variance: 47.3766, pValue: 0.001_073_4 },
];

/**
 * The two variances GMMAT's `glmmkin` fitted for the panel, from
 * `tests/reference/gwas/gmmat.null_models.tsv`, which the reference script
 * writes at full precision: the variance of the random effect of the
 * kinship, which GMMAT calls `tau`, and what is left over, its `sigma2`.
 */
const OF_GMMAT_NULL = { geneticVariance: 1.221_616_675_296_99, residualVariance: 0.342_359_482_266_917 };

/** The bytes of the panel and its phenotypes, read once for every test. */
const PANEL_VCF = await referenceKinship("panel_called.vcf.gz");
const PHENOTYPES = theColumnsOfTheFile(await referenceGwas("phenotypes.csv"));

/**
 * The kinship that `plink2 --make-rel square bin` wrote for the panel, at
 * full precision, with the individuals plink2 wrote beside it.
 *
 * It is the one the mixed model is given, as "How it is verified" of "What
 * every model shares" of the spec asks: a kinship that came from neither
 * popnei nor pyNei, and the one `tests/reference/gwas/make_reference.py`
 * gave GMMAT.
 */
const PANEL_KINSHIP = theKinshipOfThePanel(
  await referenceKinship("panel_called.plink2.rel.bin.gz"),
  new TextDecoder().decode(
    await referenceKinship("panel_called.plink2.rel.id"),
  ),
);

/**
 * The four lists of `tests/reference/gwas/refusals_of_both_layers.json`,
 * which the Python suite walks as well: the calls both layers refuse, the
 * values both read as a number, the values that mean an individual with no
 * phenotype in both, and the ones TypeScript alone refuses. That file says
 * what each list is and why the last one is there.
 */
const OF_BOTH_LAYERS = JSON.parse(
  await referenceGwas("refusals_of_both_layers.json"),
) as {
  refusals: { case: string; match: string }[];
  coercions: { case: string }[];
  no_phenotype_in_both_layers: { case: string }[];
  refused_in_typescript_alone: { case: string; match: string }[];
};

/**
 * The kinship of the panel out of the two files plink2 wrote: the little
 * endian float64 of `--make-rel square bin`, gzipped, and the names of the
 * individuals in the order of its rows, under one header line of `#IID`.
 */
function theKinshipOfThePanel(gzipped: Uint8Array, ids: string): Kinship {
  const bytes = gunzipSync(gzipped);
  const values = new Float64Array(
    bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength),
  );
  const individuals = ids
    .split("\n")
    .slice(1)
    .map((line) => line.trim())
    .filter((line) => line !== "");
  assert.equal(individuals.length, PANEL_NUM_INDIVIDUALS);
  assert.equal(values.length, PANEL_NUM_INDIVIDUALS * PANEL_NUM_INDIVIDUALS);
  return new Kinship(values, individuals, PANEL_NUM_VARS);
}

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

/**
 * The study of the panel with the kinship and the score test, which is what
 * GMMAT was given: both covariates and the matrix plink2 wrote.
 */
function theMixedStudyOfThePanel(): GwasResult {
  return gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.cont as Record<string, number>,
    trait: "continuous",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
    kinship: PANEL_KINSHIP,
    test: "score",
  });
}

/**
 * A kinship of the six individuals of the worked example, or of the five
 * that are not `without`: the identity, which is the relatedness of
 * individuals with no recent ancestor in common.
 *
 * The calls it is written for are refused before any model is fitted, so
 * what the matrix holds only has to be a kinship.
 */
function theKinshipOfTheWorkedExample(without?: string): Kinship {
  const names = ["i0", "i1", "i2", "i3", "i4", "i5"].filter(
    (name) => name !== without,
  );
  const matrix = new Float64Array(names.length * names.length);
  for (let row = 0; row < names.length; row += 1) {
    matrix[row * names.length + row] = 1;
  }
  return new Kinship(matrix, names, 3);
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

test("the six variants of the panel are gmmat's score test under a kinship", () => {
  const result = theMixedStudyOfThePanel();

  assert.equal(result.nullModel.model, "lmm");
  assert.equal(result.test, "score");
  assert.equal(result.nullModel.numIndividuals, PANEL_NUM_INDIVIDUALS);
  assert.equal(result.stats.beta.length, PANEL_NUM_VARS);
  assert.ok(
    Math.abs(
      (result.nullModel.geneticVariance as number) -
        OF_GMMAT_NULL.geneticVariance,
    ) <= OF_GMMAT_NULL_MODEL,
    `the genetic variance is ${result.nullModel.geneticVariance} and GMMAT ` +
      `gives ${OF_GMMAT_NULL.geneticVariance}`,
  );
  assert.ok(
    Math.abs(
      (result.nullModel.residualVariance as number) -
        OF_GMMAT_NULL.residualVariance,
    ) <= OF_GMMAT_NULL_MODEL,
    `the residual variance is ${result.nullModel.residualVariance} and ` +
      `GMMAT gives ${OF_GMMAT_NULL.residualVariance}`,
  );
  // The heritability is the one number of the null model built from the two
  // variances rather than read off the fit, and it is `undefined` for every
  // model but this one.
  const heritability =
    OF_GMMAT_NULL.geneticVariance /
    (OF_GMMAT_NULL.geneticVariance + OF_GMMAT_NULL.residualVariance);
  assert.ok(
    Math.abs((result.nullModel.heritability as number) - heritability) <=
      OF_GMMAT_NULL_MODEL,
    `the heritability is ${result.nullModel.heritability} and the two ` +
      `variances of GMMAT give ${heritability}`,
  );
  for (const { id, variance, pValue } of OF_GMMAT_SIX) {
    const at = rowOf(result, id);
    const se = result.stats.se[at] as number;
    assertWithin(1 / (se * se), variance, OF_GMMAT_VARIANCE, `1 / se² of ${id}`);
    const found = result.stats.pValue[at] as number;
    const apart = Math.abs(Math.log10(found / pValue));
    assert.ok(
      apart <= OF_GMMAT_P_VALUE,
      `the p-value of ${id} is ${found} and GMMAT gives ${pValue}, ${apart} ` +
        `apart in log10 against the ${OF_GMMAT_P_VALUE} allowed`,
    );
  }
});

test("the kinship is read in the order the source has the individuals", () => {
  // The phenotype, the rows of the design and the dosages of a block are
  // read together row by row, and the kinship is the relatedness of those
  // rows, so a matrix left in the order the user built it in would put one
  // individual's relatedness against another's genotypes. No message can
  // catch this one, both matrices being kinships of the same 200
  // individuals: what says it is that the two studies are equal. The Python
  // suite makes the same pair of calls.
  const ofThePanel = PANEL_KINSHIP.individuals;
  const names = [...ofThePanel].reverse();
  const matrix = new Float64Array(names.length * names.length);
  for (const [row, ofTheRow] of names.entries()) {
    for (const [column, ofTheColumn] of names.entries()) {
      matrix[row * names.length + column] = PANEL_KINSHIP.matrix[
        ofThePanel.indexOf(ofTheRow) * ofThePanel.length +
          ofThePanel.indexOf(ofTheColumn)
      ] as number;
    }
  }

  const backwards = gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.cont as Record<string, number>,
    trait: "continuous",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
    kinship: new Kinship(matrix, names, PANEL_NUM_VARS),
    test: "score",
  });
  const inTheSourcesOrder = theMixedStudyOfThePanel();

  assert.deepEqual([...backwards.stats.beta], [...inTheSourcesOrder.stats.beta]);
  assert.deepEqual([...backwards.stats.se], [...inTheSourcesOrder.stats.se]);
});

test("a tested individual the kinship has not is refused by name", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        covariates: THE_COVARIATE,
        kinship: theKinshipOfTheWorkedExample("i2"),
      }),
    { message: /`i2` is tested and is not one of the 5 individuals/ },
  );
});

test("the grammar gamma approximation of a mixed model is being written", () => {
  assert.throws(
    () =>
      gwasOf(PANEL_VCF, {
        phenotype: PHENOTYPES.cont as Record<string, number>,
        trait: "continuous",
        kinship: PANEL_KINSHIP,
        useGrammarGammaApprox: true,
      }),
    { message: /GRAMMAR-Gamma approximation is being written/ },
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
 * The study of the worked example with `options` written over its trait and
 * its covariate, as a function that makes it.
 *
 * It is what the cases of `refusals_of_both_layers.json` are written with,
 * and the option it is given is typed as it reaches a user's editor and not
 * as `calcGwas` declares it: a case of that file is a call popnei has to
 * answer for, a kinship of `"a matrix"` among them.
 */
function theStudyWith(options: Record<string, unknown>): () => GwasResult {
  return () =>
    gwasOf(WORKED_EXAMPLE, {
      phenotype: THE_TRAIT,
      trait: "continuous",
      covariates: THE_COVARIATE,
      ...options,
    } as Parameters<typeof calcGwas>[1]);
}

/**
 * The call of each case of `refusals` of `refusals_of_both_layers.json`,
 * over the worked example.
 *
 * The Python suite holds the same cases under the same names, and each
 * suite writes the call in its own language: the names of the options
 * differ between the two.
 */
function theCallsThatAreRefused(): Record<string, () => GwasResult> {
  const cov = THE_COVARIATE["cov"] as Record<string, number>;
  const { i5: _withoutI5, ...ofFive } = cov;
  const { i5: _alsoWithoutI5, ...ofThree } = THE_TRAIT;
  return {
    "a kinship that is not a kinship": theStudyWith({ kinship: "a matrix" }),
    "the grammar gamma approximation with no kinship": theStudyWith({
      useGrammarGammaApprox: true,
    }),
    "the grammar gamma approximation with a kinship": theStudyWith({
      kinship: theKinshipOfTheWorkedExample(),
      useGrammarGammaApprox: true,
    }),
    "a tested individual the kinship has not": theStudyWith({
      kinship: theKinshipOfTheWorkedExample("i2"),
    }),
    "the score test": theStudyWith({ test: "score" }),
    "a test of another name": theStudyWith({ test: "rao" }),
    "a trait of another name": theStudyWith({ trait: "quantitative" }),
    "a binomial trait": theStudyWith({
      phenotype: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
      trait: "binomial",
    }),
    "a covariate named intercept": theStudyWith({
      covariates: { intercept: cov },
    }),
    "a covariate that has no value for a tested individual": theStudyWith({
      covariates: { cov: ofFive },
    }),
    "a covariate that is not finite": theStudyWith({
      covariates: { cov: { ...cov, i2: Number.POSITIVE_INFINITY } },
    }),
    "a covariate that is a copy of another": theStudyWith({
      covariates: { cov, twice: cov },
    }),
    "an individual of the phenotype that the variants have not": theStudyWith({
      phenotype: { ...THE_TRAIT, i9: 5 },
    }),
    "fewer individuals than the design has columns plus two": theStudyWith({
      phenotype: { i0: 2, i1: 3, i2: 5 },
    }),
    "a phenotype that is a name": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: "tall" },
    }),
    "a phenotype that is the empty string": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: "" },
    }),
    "a covariate that is a name": theStudyWith({
      covariates: { cov: { ...cov, i2: "north" } },
    }),
  };
}

/**
 * The call of each case of `coercions` of that file: the worked example with
 * one of its values written as something that is not a number.
 *
 * Each of them is the worked example and nothing else, so each gives the
 * null model that example gives, which is what both suites assert of them:
 * `Number` here and `float` in Python read the same number out of the
 * string and out of the boolean.
 */
function theCallsThatAreCoerced(): Record<string, () => GwasResult> {
  const cov = THE_COVARIATE["cov"] as Record<string, number>;
  const written = (
    values: Record<string, number>,
    how: (of: number) => unknown,
  ) =>
    Object.fromEntries(
      Object.entries(values).map(([individual, value]) => [
        individual,
        how(value),
      ]),
    );
  return {
    "a phenotype written as strings": theStudyWith({
      phenotype: written(THE_TRAIT, String),
    }),
    "a covariate written as strings": theStudyWith({
      covariates: { cov: written(cov, String) },
    }),
    // The covariate of the worked example is 0 and 1, which is what a
    // boolean is read as: false, true, false, true, false, true.
    "a covariate written as booleans": theStudyWith({
      covariates: { cov: written(cov, (of) => of === 1) },
    }),
  };
}

/**
 * The call of each case of `no_phenotype_in_both_layers` of that file: a
 * value that both layers read as an individual with no phenotype, which is
 * left untested.
 *
 * A value means no phenotype exactly where `float` of it gives NaN, which is
 * NaN itself and the string `nan`. The second is the one a user does not
 * expect, and it is the one a table of traits written by a program that
 * prints NaN as text arrives with.
 */
function theCallsWithNoPhenotypeForOne(): Record<string, () => GwasResult> {
  return {
    "a phenotype that is NaN": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: Number.NaN },
    }),
    "a phenotype that is the string nan": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: "nan" },
    }),
  };
}

/**
 * The call of each case of `refused_in_typescript_alone` of that file: a
 * phenotype written as what Python has for an individual with no phenotype
 * and TypeScript has not.
 *
 * The Python suite holds the same cases and asserts the other half of each,
 * that the individual is left out of the study.
 */
function theCallsThatTypescriptAloneRefuses(): Record<string, () => GwasResult> {
  return {
    "a phenotype that is null": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: null },
    }),
    "a phenotype that is undefined": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: undefined },
    }),
  };
}

/**
 * The call of each case of `listed`, which is one list of
 * `refusals_of_both_layers.json`, with the case this suite has no call for
 * and the call the list does not hold failing here: that is what binds the
 * file to the two suites.
 */
function theCallsOf<OfTheList extends { case: string }>(
  listed: OfTheList[],
  calls: Record<string, () => GwasResult>,
  list: string,
): [OfTheList, () => GwasResult][] {
  assert.ok(listed.length > 0, `\`${list}\` of the file lists nothing`);
  assert.deepEqual(
    Object.keys(calls).sort(),
    listed.map(({ case: name }) => name).sort(),
    `the cases of \`${list}\` of refusals_of_both_layers.json and the calls ` +
      "this suite writes are not the same: a case of the file is answered " +
      "for in both suites, which is what makes it bind",
  );
  return listed.map((ofTheList) => [
    ofTheList,
    calls[ofTheList.case] as () => GwasResult,
  ]);
}

test("both layers refuse the same calls", () => {
  for (const [{ case: name, match }, call] of theCallsOf(
    OF_BOTH_LAYERS.refusals,
    theCallsThatAreRefused(),
    "refusals",
  )) {
    assert.throws(call, new RegExp(match), name);
  }
});

test("both layers read a value that is not a number as the number it holds", () => {
  for (const [{ case: name }, call] of theCallsOf(
    OF_BOTH_LAYERS.coercions,
    theCallsThatAreCoerced(),
    "coercions",
  )) {
    const result = call();

    assert.equal(result.nullModel.numIndividuals, 6, name);
    assertWithin(
      result.nullModel.covariateEffects["intercept"] as number,
      THE_NULL.intercept,
      OF_THE_WORKED_EXAMPLE,
      `the intercept of the study with ${name}`,
    );
    assertWithin(
      result.nullModel.covariateEffects["cov"] as number,
      THE_NULL.cov,
      OF_THE_WORKED_EXAMPLE,
      `the effect of the covariate of the study with ${name}`,
    );
    assertWithin(
      result.nullModel.residualVariance as number,
      THE_NULL.residualVariance,
      OF_THE_WORKED_EXAMPLE,
      `the residual variance of the study with ${name}`,
    );
  }
});

test("both layers leave an individual with no phenotype untested", () => {
  // A value means no phenotype exactly where Python's `float` of it gives
  // NaN, which is the rule the spec settled on 23 September 2026 by the
  // oracle, and the string `nan` is the case nobody guesses: `float('nan')`
  // is NaN, so it is an individual that is not tested and not a refusal.
  // The Python suite asserts the same five individuals of the same calls.
  for (const [{ case: name }, call] of theCallsOf(
    OF_BOTH_LAYERS.no_phenotype_in_both_layers,
    theCallsWithNoPhenotypeForOne(),
    "no_phenotype_in_both_layers",
  )) {
    const result = call();

    assert.deepEqual(result.individuals, ["i0", "i1", "i3", "i4", "i5"], name);
    assert.equal(result.nullModel.numIndividuals, 5, name);
  }
});

test("what says no phenotype in python is refused here, by the individual", () => {
  // `Number` turns `null` into 0 and `undefined` into NaN, and `float` of
  // either raises, so both are refused here where pandas reads them as its
  // missing value and drops the individual. A key the object has not is
  // what says that an individual has no phenotype here, and the Python
  // suite asserts of these two that its layer tests five individuals.
  for (const [{ case: name, match }, call] of theCallsOf(
    OF_BOTH_LAYERS.refused_in_typescript_alone,
    theCallsThatTypescriptAloneRefuses(),
    "refused_in_typescript_alone",
  )) {
    assert.throws(call, new RegExp(match), name);
  }
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
