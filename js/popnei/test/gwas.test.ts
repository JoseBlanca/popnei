/**
 * The association study from TypeScript: `calcGwas` and the result it gives.
 *
 * "The linear model", "The linear mixed model", "The logistic model", "The
 * logistic mixed model" and
 * "The worked example" of `docs/specs/gwas.md` have the numbers. The worked
 * example is 3 variants of
 * 6 diploid individuals with one covariate, written as a VCF here, and pyNei
 * at commit ef0ca6e gave its null model and its three rows; it reads no
 * reference file and nothing of it is rounded away. The panel is
 * `tests/reference/kinship/panel_called.vcf.gz`, 200 individuals and 1200
 * biallelic variants with every genotype called, with the trait `cont` and
 * the covariates `cov1` and `cov2` of `tests/reference/gwas/phenotypes.csv`.
 * The six variants asserted for the linear model are what plink2
 * v2.0.0-a.7.7 wrote for it, the six of the logistic model are what plink2
 * wrote for the binomial trait `binom` of the same file, and the six of each
 * mixed model are what GMMAT 1.5.0 wrote, for `cont` and for `binom`, over
 * the kinship that
 * `plink2 --make-rel` wrote and that
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
 * What plink2 v2.0.0-a.7.7 wrote for six variants of the panel with the
 * binomial trait `binom` and the covariates `cov1` and `cov2`, in
 * `tests/reference/gwas/plink2.panel_called.glm.logistic.hybrid.tsv`: the
 * effect as a log odds ratio, which is the logarithm of the `OR` plink2
 * prints, the `LOG(OR)_SE` beside it and the p-value.
 *
 * The `beta` carries more digits than plink2 prints because it is the
 * logarithm of the six digits of the odds ratio. Five of the six are the
 * causal variants of `causal_vars.csv` and `var0000` is not causal; none of
 * the six is the variant plink2 fell back to a penalized regression for,
 * which is `var0006` and which the test below is about.
 */
const OF_PLINK2_LOGISTIC_SIX: {
  id: string;
  beta: number;
  se: number;
  pValue: number;
}[] = [
  {
    id: "var0000",
    beta: -0.572_578_694_541_525_8,
    se: 0.261_917,
    pValue: 0.028_808_1,
  },
  {
    id: "var0052",
    beta: -0.852_823_096_430_041_7,
    se: 0.248_979,
    pValue: 0.000_614_166,
  },
  {
    id: "var0629",
    beta: -0.949_570_924_908_968,
    se: 0.323_553,
    pValue: 0.003_337_39,
  },
  {
    id: "var0751",
    beta: -0.265_916_666_783_602_6,
    se: 0.219_207,
    pValue: 0.225_098,
  },
  {
    id: "var1137",
    beta: -0.427_448_481_503_581_95,
    se: 0.252_402,
    pValue: 0.090_356_1,
  },
  {
    id: "var1188",
    beta: -0.830_184_139_078_324,
    se: 0.264_071,
    pValue: 0.001_667_71,
  },
];

/**
 * How far a `beta` and an `se` of the six above may be from plink2's, as a
 * share of the `se` of that variant, and how far a p-value may be, as a
 * share of itself: 1e-5, 1e-4 and 5e-3, which is what "How it is verified"
 * of "The logistic model" holds the same six literals to and what pyNei
 * holds them to.
 *
 * They are wider than the 5e-6 plink2's six printed digits round a value by
 * because plink2 stops its logistic fit earlier than popnei does: what they
 * measure is the distance between two fits. Measured under node on 24
 * September 2026, the worst effect is 2.688e-6 of the `se` of its variant,
 * `var0052`, 27 per cent of what is allowed; the worst standard error is
 * 5.301e-5 of itself, `var1137`, 53 per cent of its bound; and the worst
 * p-value is 1.878e-4 of itself, `var1137` again, 4 per cent of its bound.
 * The same six measure 2.688e-6, 5.301e-5 and 1.878e-4 natively on faer and
 * on Accelerate, so none of the distance from plink2 is WebAssembly's own
 * rounding.
 */
const OF_PLINK2_LOGISTIC_BETA = 1e-5;
const OF_PLINK2_LOGISTIC_SE = 1e-4;
const OF_PLINK2_LOGISTIC_P_VALUE = 5e-3;

/**
 * What R 4.6.1's `anova(glm, test = "Rao")` wrote for six variants of the
 * panel with every genotype called, one logistic regression per variant with
 * `cov1` and `cov2` as covariates, in
 * `tests/reference/gwas/r.panel_called.glm.score.tsv`: the score statistic,
 * which popnei gives as `(beta / se)**2`, and its p-value.
 *
 * They are the same six variants the Wald literals above are of, five of
 * them the causal variants of `causal_vars.csv` and `var0000` not causal.
 * The score test fits nothing per variant, so `var0006`, whose Wald fit runs
 * away, has an answer here and is not left out of anything.
 *
 * R reports the statistic and popnei reports `beta` and `se`, and the two
 * libraries count the dosages of a variant from different alleles, R from
 * the alternative one and popnei from the one that is not the major one
 * among the tested individuals. That turns the sign of `beta` over for some
 * variants and leaves the statistic and the p-value as they are.
 */
const OF_R_LOGISTIC_SCORE_SIX: {
  id: string;
  statistic: number;
  pValue: number;
}[] = [
  { id: "var0000", statistic: 4.938_245, pValue: 0.026_268_700 },
  { id: "var0052", statistic: 12.484_427, pValue: 0.000_410_359 },
  { id: "var0629", statistic: 9.165_576, pValue: 0.002_466_100 },
  { id: "var0751", statistic: 1.480_401, pValue: 0.223_711_736 },
  { id: "var1137", statistic: 2.911_424, pValue: 0.087_954_199 },
  { id: "var1188", statistic: 10.382_961, pValue: 0.001_271_835 },
];

/**
 * How far the score statistic of one of the six may be from R's, absolute,
 * and how far its p-value may be, in `log10`: 1e-3 and 1e-3, which is what
 * "How it is verified" of "The logistic model" holds the same six literals
 * to and what the cargo test of them holds.
 *
 * R writes both at full precision, so what these measure is where R's fit
 * stopped against where popnei's did, and not the width of a printed digit:
 * R's glm converges to 1e-8 in the deviance. The first is absolute where the
 * statistics of the six run from 1.48 to 12.48, which on a panel whose
 * statistics were far larger would fail a right answer rather than pass a
 * wrong one.
 *
 * Measured under node on 24 September 2026, the worst statistic is 6.024e-4,
 * `var0629`, 60 per cent of what is allowed, and the worst p-value is
 * 1.430e-4 in `log10`, `var0629` again, 14 per cent of its bound. The cargo
 * test of the same six measures 6.024e-4 and 1.430e-4 on faer natively and
 * on Accelerate, so none of the distance from R is WebAssembly's own
 * rounding.
 */
const OF_R_LOGISTIC_STATISTIC = 1e-3;
const OF_R_LOGISTIC_P_VALUE_IN_LOG10 = 1e-3;

/**
 * How far `1 / se**2` of a mixed model's score test may be from GMMAT's
 * `VAR`, as a share of it, and how far a p-value may be from GMMAT's in
 * `log10`: 1e-5 and 1e-4, which is what "How it is verified" of "The linear
 * mixed model" and of "The logistic mixed model" of the spec ask of every
 * variant and what the cargo tests hold the same six literals of each model
 * to. Both models take them, since both are scored against `glmm.score`.
 *
 * `gmmat.panel_called.lmm.score.tsv` and its `glmm` counterpart are printed
 * to six significant digits, which rounds a value by up to 5e-6 of itself,
 * so half of the first bound can go on GMMAT's printing alone.
 *
 * Measured under node on 24 September 2026 for the linear mixed model, the
 * worst of the six is 1.837e-6 of `VAR`, at `var0052`, and the worst p-value
 * is 4.165e-5 in `log10`, at `var0629`; the cargo test of the core measures
 * 1.843e-6 and 4.164e-5 on faer natively and 1.842e-6 and 4.164e-5 on
 * Accelerate, so nothing of the distance from GMMAT is WebAssembly's own
 * rounding. That p-value is not the printing either: 4.16e-5 in `log10` is
 * 9.6e-5 of the p-value, where six digits round it by 5e-6, and it is the
 * two fits landing 1.2e-6 apart in the genetic variance at a p-value of
 * 4.8e-5, where the tail of the chi square turns a small move of the
 * statistic into a larger one of the p-value.
 *
 * For the logistic mixed model, measured the same day: the worst of its six
 * is 2.6832488e-6 of `VAR`, at `var0751`, and the worst p-value is 5.2446e-6
 * in `log10`, at `var0052`, where the cargo test measures 2.6832488e-6 and
 * 5.2446e-6 on both backends natively. So WebAssembly agrees with the native
 * builds to eight digits there as well.
 */
const OF_GMMAT_VARIANCE = 1e-5;
const OF_GMMAT_P_VALUE = 1e-4;

/**
 * How far a number of a null model fitted over a kinship may be from
 * GMMAT's `glmmkin`: 1e-5 absolute, which is the spec's and is how far two
 * searches for the variance land apart, the restricted maximum likelihood
 * one of the linear mixed model and the penalized quasi-likelihood one of
 * the logistic mixed model.
 *
 * Seven numbers of the two models are held to it, all of them fitted over
 * the panel with every genotype called, which is the only panel this suite
 * studies. Measured under node on 24 September 2026, the worst is the
 * logistic mixed model's variance of the kinship effect, 6.033e-6 from
 * GMMAT's and 60 per cent of what is allowed, which pytest measures at the
 * same 6.033e-6 natively on both backends; its three covariate effects are
 * far nearer, the intercept 1.045e-6 away, `cov2` 1.119e-6 and `cov1`
 * 3.390e-7. The linear mixed model's three are nearer as well: its genetic
 * variance is 1.273e-6 away, 13 per cent of what is allowed, where the
 * cargo test measures 1.218e-6 on faer natively, its residual variance
 * 1.073e-6 and its heritability 7.138e-7.
 *
 * The two models are that far away for different reasons. The logistic
 * model's 6.0e-6 is the two programs and not the route popnei takes:
 * `docs/reports/glmm-method/README.md` measured pyNei's fit and the cheaper
 * one alike at 6.3e-6 of GMMAT's variance. The linear model's 1.2e-6 is the
 * kinship and not the search: GMMAT was given the six printed digits of
 * plink2's matrix and this suite reads the float64 beside them.
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
 * What GMMAT 1.5.0's `glmm.score` gave for six variants of the panel with
 * every genotype called under the logistic mixed model, from
 * `tests/reference/gwas/gmmat.panel_called.glmm.score.tsv`: the variance of
 * the score, which is `x' p x` and which popnei gives as `1 / se**2`, and
 * the p-value.
 *
 * They are the literals of "How it is verified" of "The logistic mixed
 * model" of the spec, and they are not the linear mixed model's above: the
 * trait is `binom` and not `cont`, so every one of the twelve numbers
 * differs. GMMAT was given both covariates and the kinship plink2 wrote for
 * this panel, which is what this suite gives popnei.
 */
const OF_GMMAT_GLMM_SIX: { id: string; variance: number; pValue: number }[] = [
  { id: "var0000", variance: 6.486_64, pValue: 0.702_659 },
  { id: "var0052", variance: 8.988_34, pValue: 0.030_670_3 },
  { id: "var0629", variance: 6.499_56, pValue: 0.089_510_4 },
  { id: "var0751", variance: 10.563, pValue: 0.014_233_1 },
  { id: "var1137", variance: 9.098, pValue: 0.093_808 },
  { id: "var1188", variance: 9.050_95, pValue: 0.026_273_7 },
];

/**
 * The four numbers GMMAT's `glmmkin` fitted for the logistic mixed model of
 * the panel, from the `glmm` row of
 * `tests/reference/gwas/gmmat.null_models.tsv`, at full precision: the
 * variance of the random effect of the kinship, which GMMAT calls `tau`,
 * and the effects of the intercept, of `cov1` and of `cov2`.
 *
 * Its `sigma2` is 1 and is not read: a logistic model has no free residual
 * variance, so popnei gives `undefined` for it and for the heritability
 * built from the two, which this suite asserts instead.
 */
const OF_GMMAT_GLMM_NULL = {
  geneticVariance: 1.508_056_730_382_11,
  intercept: -1.416_463_894_533_21,
  cov1: 0.753_476_451_040_228,
  cov2: 1.583_209_937_850_49,
};

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
  refusals: { case: string; match?: string; match_in_typescript?: string }[];
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
 * The study of the binomial trait of the panel with `cov1` and `cov2`, which
 * is what plink2 was given for its logistic regression and, with `test` set
 * to `score`, what R fitted one regression per variant for.
 */
function theLogisticStudyOfThePanel(test?: "wald" | "score"): GwasResult {
  return gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.binom as Record<string, number>,
    trait: "binomial",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
    test,
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
 * The study of the panel with the kinship and the default test of the linear
 * mixed model, which is the Wald one, `useGrammarGammaApprox` as it is
 * given.
 *
 * It is the study the Python suite compares the approximation against, so
 * the three numbers the approximation gives here are the same three
 * measured there, on the other backend.
 */
function theMixedStudyOfThePanelWith(
  useGrammarGammaApprox: boolean,
): GwasResult {
  return gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.cont as Record<string, number>,
    trait: "continuous",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
    kinship: PANEL_KINSHIP,
    useGrammarGammaApprox,
  });
}

/**
 * The study of the panel with the binomial trait, both covariates and the
 * kinship, which is what GMMAT was given for its logistic mixed model.
 *
 * No test is asked for: the logistic mixed model has the score test alone,
 * and that is what the default takes, where the other three default to the
 * Wald test.
 *
 * `useGrammarGammaApprox` asks for the approximate denominator in place of
 * the exact one, which only the test of the approximation does: the
 * comparison with GMMAT is against the exact answer.
 */
function theLogisticMixedStudyOfThePanel(
  useGrammarGammaApprox = false,
): GwasResult {
  return gwasOf(PANEL_VCF, {
    phenotype: PHENOTYPES.binom as Record<string, number>,
    trait: "binomial",
    covariates: {
      cov1: PHENOTYPES.cov1 as Record<string, number>,
      cov2: PHENOTYPES.cov2 as Record<string, number>,
    },
    kinship: PANEL_KINSHIP,
    useGrammarGammaApprox,
  });
}

/**
 * A binomial trait of the six individuals of the worked example that its
 * covariate, 0, 1, 0, 1, 0, 1, does not separate.
 *
 * The two calls that have to reach a logistic mixed model take this one: the
 * trait 0, 1, 0, 1, 0, 1 is the covariate itself, so every fit of it walks
 * towards an infinite effect and is refused before the mixed model is
 * reached at all.
 */
const A_BINOMIAL_TRAIT_THE_COVARIATE_DOES_NOT_SEPARATE = {
  i0: 0,
  i1: 0,
  i2: 1,
  i3: 1,
  i4: 0,
  i5: 1,
};

/**
 * A kinship of the six individuals of the worked example with `i0` and `i1`
 * given a relatedness of 100, which no covariance has.
 *
 * The 2 by 2 block of that pair has an eigenvalue of -99, and a weight of
 * the logistic mixed model is at most 0.25, so the reciprocals put 4 at
 * least on every diagonal entry of the covariance of the working trait: a
 * variance of the kinship effect above about 0.04 takes that covariance
 * below 0 and the Cholesky factorization refuses it. The fit starts at half
 * the variance of the first working trait, which is far above that, so the
 * first linearization is where it stops.
 *
 * What a user reaches it with is a kinship built from variants with many
 * genotypes missing, where every pair is counted over its own variants; this
 * is that matrix pushed far enough to fail on six individuals.
 */
function aKinshipThatIsNotACovariance(): Kinship {
  const names = ["i0", "i1", "i2", "i3", "i4", "i5"];
  const matrix = new Float64Array(names.length * names.length);
  for (let row = 0; row < names.length; row += 1) {
    matrix[row * names.length + row] = 1;
  }
  matrix[1] = 100;
  matrix[names.length] = 100;
  return new Kinship(matrix, names, 3);
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

test("the six variants of the panel are plink2's logistic effect, error and p-value", () => {
  const result = theLogisticStudyOfThePanel();

  assert.equal(result.nullModel.model, "glm");
  assert.equal(result.test, "wald");
  assert.equal(result.trait, "binomial");
  assert.equal(result.nullModel.numIndividuals, PANEL_NUM_INDIVIDUALS);
  assert.equal(result.stats.beta.length, PANEL_NUM_VARS);
  // A binomial trait has no residual variance, its variance being decided
  // by its mean, and a study with no kinship has no genetic variance.
  assert.equal(result.nullModel.residualVariance, undefined);
  assert.equal(result.nullModel.geneticVariance, undefined);
  assert.equal(result.nullModel.heritability, undefined);
  for (const { id, beta, se, pValue } of OF_PLINK2_LOGISTIC_SIX) {
    const at = rowOf(result, id);
    assertWithinTheScale(
      result.stats.beta[at] as number,
      beta,
      se,
      OF_PLINK2_LOGISTIC_BETA,
      `the effect of ${id}`,
    );
    assertWithinTheScale(
      result.stats.se[at] as number,
      se,
      se,
      OF_PLINK2_LOGISTIC_SE,
      `the standard error of ${id}`,
    );
    assertWithin(
      result.stats.pValue[at] as number,
      pValue,
      OF_PLINK2_LOGISTIC_P_VALUE,
      `the p-value of ${id}`,
    );
  }
});

test("the one variant whose logistic fit runs away has no answer here either", () => {
  // `var0006` separates the individuals that have the condition from those
  // that have not, so its effect has no finite value to reach: plink2 falls
  // back to a penalized regression there and popnei gives NaN. What catches
  // it is the mark for an effect past 30 in absolute value, at round 29 on
  // both linear algebra backends, which is what "What it gives" of "The
  // logistic model" of the spec measured; the two marks for a value that is
  // not finite fire at no fixture of any of the three suites.
  const result = theLogisticStudyOfThePanel();

  const withoutAnAnswer = (result.stats.id as readonly string[]).filter(
    (_id, at) => Number.isNaN(result.stats.pValue[at]),
  );
  assert.deepEqual(withoutAnAnswer, ["var0006"]);
  const at = rowOf(result, "var0006");
  assert.ok(Number.isNaN(result.stats.beta[at] as number));
  assert.ok(Number.isNaN(result.stats.se[at] as number));
  // The frequency of such a variant is still there, as it is for a variant
  // with no variance: what it has not is a test.
  assert.ok(!Number.isNaN(result.stats.alleleFreq[at] as number));
});

test("the logistic score test of the panel is fitted and the wald test is the default", () => {
  // `test` takes a value that is not the default here, which is the one
  // thing this suite has no other call for on this model, and the two tests
  // share the null fit and nothing else: the coefficients are the same and
  // the rows are not.
  const byDefault = theLogisticStudyOfThePanel();
  const score = theLogisticStudyOfThePanel("score");

  assert.equal(score.nullModel.model, "glm");
  assert.equal(score.test, "score");
  assert.equal(byDefault.test, "wald");
  assert.deepEqual(Object.keys(score.nullModel.covariateEffects), [
    "intercept",
    "cov1",
    "cov2",
  ]);
  for (const name of ["intercept", "cov1", "cov2"]) {
    assert.equal(
      score.nullModel.covariateEffects[name],
      byDefault.nullModel.covariateEffects[name],
      `the ${name} of the two tests is the same null fit`,
    );
  }
  // The score test fits nothing per variant, so `var0006` has an answer
  // there where the Wald test has none.
  const separating = score.stats.pValue[rowOf(score, "var0006")] as number;
  assert.ok(!Number.isNaN(separating));
});

test("the six variants of the panel are r's logistic score statistic and p-value", () => {
  // The per variant arithmetic of the score test is its own, and until this
  // test was written nothing under node asserted a number of it: the test
  // above compares its null fit with the Wald study's, which is popnei
  // against popnei, and then only that `var0006` is not NaN. A reviewer
  // multiplied this test's `beta` and its statistic by 1.5 in the core,
  // rebuilt the WebAssembly and saw all 328 node tests pass while two cargo
  // tests failed.
  const result = theLogisticStudyOfThePanel("score");

  assert.equal(result.nullModel.model, "glm");
  assert.equal(result.test, "score");
  assert.equal(result.trait, "binomial");
  assert.equal(result.nullModel.numIndividuals, PANEL_NUM_INDIVIDUALS);
  assert.equal(result.stats.beta.length, PANEL_NUM_VARS);
  for (const { id, statistic, pValue } of OF_R_LOGISTIC_SCORE_SIX) {
    const at = rowOf(result, id);
    const ours =
      (result.stats.beta[at] as number) / (result.stats.se[at] as number);
    const ofPopnei = ours * ours;
    const difference = Math.abs(ofPopnei - statistic);
    assert.ok(
      difference <= OF_R_LOGISTIC_STATISTIC,
      `the score statistic of ${id} is ${ofPopnei} and R gives ${statistic}, ` +
        `${difference} away against the ${OF_R_LOGISTIC_STATISTIC} allowed`,
    );
    const found = result.stats.pValue[at] as number;
    const apart = Math.abs(Math.log10(found / pValue));
    assert.ok(
      apart <= OF_R_LOGISTIC_P_VALUE_IN_LOG10,
      `the p-value of ${id} is ${found} and R gives ${pValue}, ${apart} ` +
        `apart in log10 against the ${OF_R_LOGISTIC_P_VALUE_IN_LOG10} allowed`,
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

test("the wald test of a binomial trait with a kinship is refused", () => {
  assert.throws(
    () =>
      gwasOf(WORKED_EXAMPLE, {
        phenotype: A_BINOMIAL_TRAIT_THE_COVARIATE_DOES_NOT_SEPARATE,
        trait: "binomial",
        covariates: THE_COVARIATE,
        kinship: theKinshipOfTheWorkedExample(),
        test: "wald",
      }),
    { message: /one mixed model for every variant/ },
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

test("the six variants of the panel are gmmat's logistic mixed score test", () => {
  const result = theLogisticMixedStudyOfThePanel();

  assert.equal(result.nullModel.model, "glmm");
  assert.equal(result.trait, "binomial");
  // The score test is the default of this model and the only test it has,
  // where the other three default to the Wald test.
  assert.equal(result.test, "score");
  assert.equal(result.nullModel.numIndividuals, PANEL_NUM_INDIVIDUALS);
  assert.equal(result.stats.beta.length, PANEL_NUM_VARS);
  const fitted: typeof OF_GMMAT_GLMM_NULL = {
    geneticVariance: result.nullModel.geneticVariance as number,
    intercept: result.nullModel.covariateEffects["intercept"] as number,
    cov1: result.nullModel.covariateEffects["cov1"] as number,
    cov2: result.nullModel.covariateEffects["cov2"] as number,
  };
  for (const [name, expected] of Object.entries(OF_GMMAT_GLMM_NULL)) {
    const found = fitted[name as keyof typeof OF_GMMAT_GLMM_NULL];
    assert.ok(
      Math.abs(found - expected) <= OF_GMMAT_NULL_MODEL,
      `the ${name} of the null model is ${found} and GMMAT gives ` +
        `${expected}, against the ${OF_GMMAT_NULL_MODEL} allowed`,
    );
  }
  // A logistic model has no free residual variance, its trait's variance
  // being decided by its mean, so neither it nor the heritability built from
  // the two is there.
  assert.equal(result.nullModel.residualVariance, undefined);
  assert.equal(result.nullModel.heritability, undefined);
  for (const { id, variance, pValue } of OF_GMMAT_GLMM_SIX) {
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

test("a kinship that does not tell the two variances apart gives none of them", () => {
  // The identity is what a user passes to mean no relatedness, and with it
  // the model is the ordinary linear one whatever the split between the
  // genetic variance and the residual one: the criterion of the search is
  // flat over the whole grid and which point wins is rounding. The three
  // fields are `undefined`, which is what the linear model gives for the
  // two it has not, and not 0 and not NaN. The rows are still the worked
  // example's, because the test is scale free. It is the meanwhile of
  // Open 3 of `docs/specs/gwas.md`, and the test below is its other end.
  const result = gwasOf(WORKED_EXAMPLE, {
    phenotype: THE_TRAIT,
    trait: "continuous",
    covariates: THE_COVARIATE,
    kinship: theKinshipOfTheWorkedExample(),
  });

  assert.equal(result.nullModel.model, "lmm");
  assert.equal(result.nullModel.geneticVariance, undefined);
  assert.equal(result.nullModel.residualVariance, undefined);
  assert.equal(result.nullModel.heritability, undefined);
  assert.equal(result.nullModel.numIndividuals, 6);
  for (const { id, beta, pValue } of THE_ROWS) {
    const at = rowOf(result, id);
    if (Number.isNaN(beta)) {
      assert.ok(Number.isNaN(result.stats.beta[at] as number), `beta of ${id}`);
      continue;
    }
    assertWithin(
      result.stats.beta[at] as number,
      beta,
      OF_THE_WORKED_EXAMPLE,
      `the effect of ${id}`,
    );
    assertWithin(
      result.stats.pValue[at] as number,
      pValue,
      OF_THE_WORKED_EXAMPLE,
      `the p-value of ${id}`,
    );
  }
});

test("the panel tells the two variances apart and gives all three", () => {
  // The other end of the test above: a rule that gave nothing wherever a
  // mixed model was fitted would pass that one and take the heritability
  // away from every user. The numbers are GMMAT's and are asserted where
  // the six literals are.
  const result = theMixedStudyOfThePanel();

  assert.notEqual(result.nullModel.geneticVariance, undefined);
  assert.notEqual(result.nullModel.residualVariance, undefined);
  assert.notEqual(result.nullModel.heritability, undefined);
});

/**
 * How far the GRAMMAR-Gamma approximation may be from the exact answer on the
 * panel, under the linear mixed model and under the logistic one: the median
 * of `log10(p_approx/p_exact)` within 0.1 of 0, the largest of those within
 * 1.5, and `beta` within 0.5 of itself.
 *
 * The three are pyNei's own numbers in `test_grammar_gamma_approx`, which
 * "How it is verified" of the approximation in `docs/specs/gwas.md` carries,
 * and the Python suite asserts the same three on the same study. They are
 * not bounds measured here and lowered until they broke: the spec gives
 * them, and a run either fills them or does not.
 *
 * The third is nearly full, and that is the method and not a defect. The
 * factor is one number standing in for a quantity that differs from variant
 * to variant, and on this panel the 100 ratios it is the mean of run from
 * 0.312 to 0.670 around a mean of 0.517. Measured on 24 September 2026: in
 * WebAssembly, which is faer, the median is -5.1885e-4, the largest 0.51185
 * and the worst `beta` 0.48966; natively the same three are -5.1885e-4,
 * 0.51185 and 0.48966 on Accelerate and on faer. So the `beta` bound is 98
 * per cent spent and the other two are nowhere near theirs.
 *
 * The same three under the logistic mixed model, on the same panel and the
 * same day: in WebAssembly the median is -1.0467e-3, the largest 0.52373 and
 * the worst `beta` 0.39335, and the Python suite measures the same three to
 * six digits on Accelerate. That model spends less of the `beta` bound than
 * the linear one and a little more of the other two.
 */
const OF_THE_APPROXIMATION_MEDIAN = 0.1;
const OF_THE_APPROXIMATION_LARGEST = 1.5;
const OF_THE_APPROXIMATION_BETA = 0.5;

/**
 * How far the worst effect of the panel has to move under the approximation
 * for the approximation to have been made at all: 0.05 of it.
 *
 * The three bounds above are ceilings, and a run that made no approximation
 * would pass every one of them, an exact answer being at no distance from
 * itself. It is not a hypothesis: a reviewer replaced the approximation with
 * the exact denominator in the core on 24 September 2026, and this suite and
 * the Python one stayed green under both mixed models. So the worst effect
 * is held above a floor as well, which is the assertion a study that had
 * quietly stopped approximating would fail.
 *
 * The floor is 0.05 where the worst effect measured that day is 0.48966
 * under the linear mixed model and 0.39335 under the logistic one, a factor
 * of eight of room. The p-values carry no floor: they are where the
 * approximation moves least, the median being 5.2e-4 and 1.0e-3 in log10,
 * and a floor near those would go red on a panel the approximation happens
 * to suit.
 */
const THE_APPROXIMATION_MOVES_THE_EFFECT = 0.05;

/**
 * How far the study `approximated` sits from the study `exact` over the
 * panel: how many variants both answered, the median and the largest of
 * `log10(p_approx / p_exact)` over those, and the worst `beta` as a share of
 * itself.
 *
 * A variant the exact test leaves with no answer is answered under the
 * approximation, which "Open 2's threshold under the approximation" of the
 * spec says: the approximate denominator is a positive factor times a sum of
 * squares and holds no cancellation. So the two are compared where both have
 * an answer, and the count of those is returned for the caller to assert.
 */
function theDistanceFromTheExactAnswer(
  approximated: GwasResult,
  exact: GwasResult,
): { answered: number; median: number; largest: number; worst: number } {
  const moved: number[] = [];
  const apart: number[] = [];
  for (let at = 0; at < PANEL_NUM_VARS; at += 1) {
    const ofTheApproximation = approximated.stats.pValue[at] as number;
    const ofTheExact = exact.stats.pValue[at] as number;
    if (!Number.isFinite(ofTheApproximation) || !Number.isFinite(ofTheExact)) {
      continue;
    }
    moved.push(Math.log10(ofTheApproximation / ofTheExact));
    const betaOfTheApproximation = approximated.stats.beta[at] as number;
    const betaOfTheExact = exact.stats.beta[at] as number;
    apart.push(
      Math.abs(betaOfTheApproximation - betaOfTheExact) /
        Math.abs(betaOfTheExact),
    );
  }
  const sorted = [...moved].sort((one, other) => one - other);
  const half = sorted.length / 2;
  return {
    answered: moved.length,
    median: ((sorted[half - 1] as number) + (sorted[half] as number)) / 2,
    largest: Math.max(...moved.map((value) => Math.abs(value))),
    worst: Math.max(...apart),
  };
}

test("the approximation of the panel is near the exact answer", () => {
  const approximated = theMixedStudyOfThePanelWith(true);
  const exact = theMixedStudyOfThePanelWith(false);

  assert.equal(approximated.usedGrammarGammaApprox, true);
  assert.equal(exact.usedGrammarGammaApprox, false);
  assert.equal(approximated.stats.pValue.length, PANEL_NUM_VARS);
  const { answered, median, largest, worst } = theDistanceFromTheExactAnswer(
    approximated,
    exact,
  );
  // This panel has no variant that the exact test leaves without an answer.
  assert.equal(answered, PANEL_NUM_VARS);
  assert.ok(
    Math.abs(median) <= OF_THE_APPROXIMATION_MEDIAN,
    `the median of log10(p_approx / p_exact) over the panel is ${median}`,
  );
  assert.ok(
    largest <= OF_THE_APPROXIMATION_LARGEST,
    `the largest |log10(p_approx / p_exact)| over the panel is ${largest}`,
  );
  assert.ok(
    worst <= OF_THE_APPROXIMATION_BETA,
    `the worst \`beta\` of the panel is ${worst} of itself away from the exact one`,
  );
  // Every bound above is a ceiling on how far the two runs lie apart, and a
  // study that had silently stopped approximating would pass all of them,
  // being the exact answer compared with itself. This is what would fail
  // there.
  assert.ok(
    worst >= THE_APPROXIMATION_MOVES_THE_EFFECT,
    `the worst \`beta\` of the panel is ${worst} of itself away from the ` +
      `exact one, and an approximation that was made moves it by ` +
      `${THE_APPROXIMATION_MOVES_THE_EFFECT} at least`,
  );
});

test("the logistic approximation of the panel is near the exact answer", () => {
  // Both mixed models take the approximation, and until this test the
  // layers above the core ran it on the continuous trait alone. The study
  // is the one GMMAT was given, the binomial trait with both covariates and
  // the kinship plink2 wrote, and what it is checked against is popnei's own
  // exact answer: GMMAT computes the exact denominator, so no program
  // outside popnei answers the approximated question.
  const approximated = theLogisticMixedStudyOfThePanel(true);
  const exact = theLogisticMixedStudyOfThePanel(false);

  assert.equal(approximated.nullModel.model, "glmm");
  assert.equal(approximated.usedGrammarGammaApprox, true);
  assert.equal(exact.usedGrammarGammaApprox, false);
  assert.equal(approximated.stats.pValue.length, PANEL_NUM_VARS);
  const { answered, median, largest, worst } = theDistanceFromTheExactAnswer(
    approximated,
    exact,
  );
  assert.equal(answered, PANEL_NUM_VARS);
  assert.ok(
    Math.abs(median) <= OF_THE_APPROXIMATION_MEDIAN,
    `the median of log10(p_approx / p_exact) over the panel is ${median}`,
  );
  assert.ok(
    largest <= OF_THE_APPROXIMATION_LARGEST,
    `the largest |log10(p_approx / p_exact)| over the panel is ${largest}`,
  );
  assert.ok(
    worst <= OF_THE_APPROXIMATION_BETA,
    `the worst \`beta\` of the panel is ${worst} of itself away from the exact one`,
  );
  assert.ok(
    worst >= THE_APPROXIMATION_MOVES_THE_EFFECT,
    `the worst \`beta\` of the panel is ${worst} of itself away from the ` +
      `exact one, and an approximation that was made moves it by ` +
      `${THE_APPROXIMATION_MOVES_THE_EFFECT} at least`,
  );
});

test("the approximation is asked for and the result says it was used", () => {
  // It is the one place a user can tell the two apart: the columns of
  // `stats` have the same names and the same shape either way, and a study
  // that had quietly made the exact test would look the same.
  const approximated = theMixedStudyOfThePanelWith(true);

  assert.equal(approximated.usedGrammarGammaApprox, true);
  assert.equal(approximated.nullModel.model, "lmm");
  // The pass that estimates the factor reads the same variants as the one
  // that tests them, and it is the counts of the first that come back, so a
  // study that read the source twice reports the variants once.
  assert.equal(approximated.passStats.numVars, PANEL_NUM_VARS);
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
    "a tested individual the kinship has not": theStudyWith({
      kinship: theKinshipOfTheWorkedExample("i2"),
    }),
    "the score test": theStudyWith({ test: "score" }),
    "a test of another name": theStudyWith({ test: "rao" }),
    "a trait of another name": theStudyWith({ trait: "quantitative" }),
    "the wald test of a binomial trait with a kinship": theStudyWith({
      phenotype: A_BINOMIAL_TRAIT_THE_COVARIATE_DOES_NOT_SEPARATE,
      trait: "binomial",
      kinship: theKinshipOfTheWorkedExample(),
      test: "wald",
    }),
    "a kinship the logistic mixed model cannot factor": theStudyWith({
      phenotype: A_BINOMIAL_TRAIT_THE_COVARIATE_DOES_NOT_SEPARATE,
      trait: "binomial",
      kinship: aKinshipThatIsNotACovariance(),
    }),
    // The trait of the worked example is 2, 3, 5, 4, 4, 7, so the first
    // tested individual is the one the message names.
    "a binomial trait that is neither 0 nor 1": theStudyWith({
      trait: "binomial",
    }),
    // The covariate of the worked example is 0, 1, 0, 1, 0, 1, which is this
    // phenotype, so it separates the individuals that have the condition
    // from the ones that have not and the null fit walks towards an infinite
    // effect for it instead of settling.
    "a binomial null model that walks towards an infinite coefficient":
      theStudyWith({
        phenotype: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
        trait: "binomial",
      }),
    // The same trait and the same covariate with a kinship, which makes it a
    // logistic mixed model. That model starts from the plain logistic null
    // above, so it is that null that runs away, and the message names the
    // model the user asked for and not the one the starting fit is.
    "a binomial null model with a kinship that runs away": theStudyWith({
      phenotype: { i0: 0, i1: 1, i2: 0, i3: 1, i4: 0, i5: 1 },
      trait: "binomial",
      kinship: theKinshipOfTheWorkedExample(),
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
    "a trait that is the same in every individual": theStudyWith({
      phenotype: Object.fromEntries(
        Object.keys(THE_TRAIT).map((individual) => [individual, 4]),
      ),
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
    "a phenotype that is an infinity": theStudyWith({
      phenotype: { ...THE_TRAIT, i2: Number.POSITIVE_INFINITY },
    }),
    "a covariate that has no value at an individual": theStudyWith({
      covariates: { cov: { ...cov, i2: Number.NaN } },
    }),
    "the covariates explaining the whole of the trait": theStudyWith({
      covariates: { itself: THE_TRAIT },
      kinship: theKinshipOfTheWorkedExample(),
    }),
    "a kinship that is not symmetric": theStudyWith({
      kinship: aKinshipWrittenInto(),
    }),
    "an approximation that is not a boolean": theStudyWith({
      useGrammarGammaApprox: "no",
    }),
  };
}

/**
 * A kinship of the worked example with one cell of a pair written into
 * after it was built, which is what the core is there to catch.
 *
 * The constructor of `Kinship` refuses a matrix that is not symmetric, and
 * the `Float64Array` it keeps is the caller's to write into, so the check it
 * made says nothing about what a study is given later. The eigendecomposition
 * reads the lower triangle alone, so such a matrix was being read as that
 * half mirrored.
 */
function aKinshipWrittenInto(): Kinship {
  const kinship = theKinshipOfTheWorkedExample();
  kinship.matrix[1] = 0.5;
  return kinship;
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
  // The digits of any script are the digits `float` reads, so a trait
  // written with the full width ones of a spreadsheet is the same trait.
  const inFullWidth = (of: number) =>
    String(of).replaceAll(/\d/gu, (digit) =>
      String.fromCodePoint(0xff10 + Number(digit)),
    );
  return {
    "a phenotype written in full width digits": theStudyWith({
      phenotype: written(THE_TRAIT, inFullWidth),
    }),
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
  for (const [ofTheList, call] of theCallsOf(
    OF_BOTH_LAYERS.refusals,
    theCallsThatAreRefused(),
    "refusals",
  )) {
    // A case whose refusal comes from the machinery of its own language
    // carries one expression for each, since neither message can hold the
    // other's words; what binds there is that both layers refuse the call.
    const match = ofTheList.match ?? ofTheList.match_in_typescript;
    assert.ok(
      match !== undefined,
      `\`${ofTheList.case}\` has neither \`match\` nor \`match_in_typescript\``,
    );
    assert.throws(call, new RegExp(match), ofTheList.case);
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
