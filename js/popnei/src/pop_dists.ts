/**
 * The distances between populations, and the result they come in.
 *
 * Two populations are far apart when the alleles of their individuals are not
 * the same alleles in the same proportions, and there are seven measures of
 * how far, each answering a different question. `calcPopDists` calculates
 * the ones a user asks for in one pass over the variants, for every pair of
 * the populations they name, and gives them in a `PopDists`: one `Distances`
 * for each measure, with the standard error of each pair beside its value.
 *
 * The seven are the names of `PopDistMeasure`. Two of them are calculated
 * today, Hudson's F_ST and f_2; the other five are refused, and the work
 * packages 2 and 3 of `docs/plans/dists-pops.md` add them.
 *
 * `docs/specs/dists.md` has each measure, what it answers, the program it is
 * verified against and the numbers the tests assert.
 */

import {
  default_min_num_individuals as defaultMinNumIndividuals,
  pop_dist_measures_that_have_a_value as measuresThatHaveAValue,
} from "../wasm/popnei.js";

import {
  namesOf,
  popsOfTheObject,
  whatWasGiven,
  wholeNumberOfZeroOrMore,
} from "./arguments.js";
import { theWasmHasToBeLoaded } from "./core.js";
import { Distances } from "./dists.js";
import type { PassStats, Variants } from "./variant.js";
import { passStatsOf, sourceOfTheVariants } from "./variant.js";

/**
 * The seven measures of how far apart two populations are, each named as the
 * field of the result that holds it is named in Python.
 */
const THE_MEASURES = [
  "fst",
  "f2",
  "chord",
  "da",
  "dest",
  "gst",
  "gst_standardized",
] as const;

/**
 * The name of one of the seven measures of how far apart two populations
 * are: Hudson's F_ST, how much of the diversity of the two taken together
 * lies between them rather than within them; f_2, how much allele frequency
 * they have drifted apart by, in the units it was measured in, which is what
 * makes it add up along a tree; the chord distance of Cavalli-Sforza and
 * Edwards and Nei's D_A, its square, the two that are Euclidean and that a
 * tree or a principal coordinate analysis is built from; Jost's D, how
 * different the alleles the two hold are; and Nei's G_ST with the
 * standardized G''_ST, which is G_ST divided by the largest value it could
 * reach with the diversity the two hold.
 */
export type PopDistMeasure = (typeof THE_MEASURES)[number];

/** One resampling group: its chromosome and the positions it holds. */
export interface PopDistGroup {
  /** The name of its chromosome. */
  readonly chrom: string;
  /** The position of its first variant, 1 based as in a VCF. */
  readonly start: number;
  /** The position of its last variant, that one included. */
  readonly end: number;
}

/** What `calcPopDists` calculates, and how it cuts the variants. */
export interface CalcPopDistsOptions {
  /**
   * How the variants are cut into the groups the standard errors are
   * resampled over. It has no default and the call fails without it.
   *
   * A number is a length in base pairs of a chromosome, 1 or more;
   * `"variant"` makes each variant a group of its own; and `null` asks for
   * no standard errors, which is the only one of the three that does not
   * need the chromosome and the position of the variants.
   *
   * A length needs the variants of each chromosome to come together and in
   * order of position, and a source whose variants go back is an `Error`
   * that names the chromosome and the two positions: the cut compares the
   * position of a variant with the first position of the group being
   * filled, so a variant that goes back joins that group instead of
   * starting one and the groups are not the stretches that were asked for.
   * `"variant"` and `null` take a source in any order.
   *
   * A group has to be longer than the distance over which two variants still
   * carry the same history, because two groups that share it are not the
   * independent draws the standard error takes them for, and there have to
   * be at least 20 groups, which is an `Error` below. The number of the
   * f-statistics literature is 5 centimorgans, about 5 million base pairs in
   * humans, and it does not carry over by itself: linkage disequilibrium
   * runs 6.1 to 12.5 centimorgans in cultivated tomato and falls off within
   * 18 thousand base pairs in its wild relative *S. pimpinellifolium*, and
   * one length in base pairs is several different lengths in centimorgans
   * along one chromosome. A user who does not know the decay distance of
   * their own panel measures it from the curve of r^2 against distance. A
   * few hundred microsatellite loci scattered over a genome have no linkage
   * to speak of and take `"variant"`; a panel of linked variants must not
   * use it.
   *
   * `"variant"` is also what the memory of the pass grows with. popnei keeps
   * six numbers, 48 bytes, for each pair of populations and each group, so a
   * group of each variant makes that 48 bytes for each pair and each
   * variant: 1200 variants of 20 populations, which are 190 pairs, are
   * 10.9 MB, and a million variants are 144 MB for 3 populations and 9.1 GB
   * for 20, which no tab gives. A length in base pairs, whose groups are as
   * many as the stretches of the chromosomes, does not grow with the
   * variants. A tab that has not the memory is an `Error` and not a wrong
   * number.
   */
  jackknifeGroup: number | "variant" | null;

  /**
   * Which of the seven measures to calculate, all of them when it is not
   * given, since the pass is what costs and each measure is a division at
   * the end of it. A result holds `null` for one nobody asked for.
   *
   * Five of the seven are not calculated yet, so asking for one of them, and
   * asking for all of them by leaving this out, is an `Error` that names the
   * two that are.
   */
  measures?: readonly PopDistMeasure[];

  /**
   * How many called genotypes a population needs at a variant for that
   * variant to count for a pair, 20 when it is not given. The test is made
   * for each pair on its own, so a population with fewer loses that variant
   * in every pair it is in and the pairs it is not in keep it.
   */
  minNumIndividuals?: number;
}

/**
 * The measures of how far apart the populations are, one `Distances` for
 * each that was asked for and `null` for the ones that were not.
 *
 * Every measure is over the same pairs, in the order (0, 1), (0, 2), ...,
 * (1, 2), ..., over the populations in the order `pops` has them, which is
 * the order the keys of the `pops` that was given iterate in.
 */
export interface PopDists {
  /**
   * The names of the populations, in the order the keys of `pops` iterate
   * in, which is the order of the pairs of every measure.
   */
  readonly pops: readonly string[];

  /** Hudson's F_ST of every pair, and `null` when it was not asked for. */
  readonly fst: Distances | null;

  /** f_2 of every pair, and `null` when it was not asked for. */
  readonly f2: Distances | null;

  /** The chord distance of every pair, and `null` when it was not asked for. */
  readonly chord: Distances | null;

  /** Nei's D_A of every pair, and `null` when it was not asked for. */
  readonly da: Distances | null;

  /** Jost's D of every pair, and `null` when it was not asked for. */
  readonly dest: Distances | null;

  /** Nei's G_ST of every pair, and `null` when it was not asked for. */
  readonly gst: Distances | null;

  /**
   * The standardized G''_ST of every pair, and `null` when it was not asked
   * for.
   */
  readonly gstStandardized: Distances | null;

  /**
   * How many variants counted for each pair, in the order of the pairs.
   *
   * A variant counts for a pair when both of its populations have at least
   * `minNumIndividuals` called genotypes there, so a variant one population
   * is short of is lost by the pairs that population is in and kept by the
   * others: two pairs are means over different variants, which is why each
   * pair carries its own count.
   */
  readonly numVars: Int32Array;

  /**
   * f_2 within each resampling group, `numGroups` x `numPairs` values row by
   * row, and `null` when no standard errors were asked for.
   *
   * f_3 and f_4, the statistics of three and of four populations that
   * admixture graphs are fitted with, are sums and differences of these, so
   * they can be built from this array without reading the genotypes again.
   */
  readonly f2Groups: Float64Array | null;

  /**
   * How many resampling groups the variants fell into, which is the rows of
   * `f2Groups` and the length of `groupIds`, and 0 when no standard errors
   * were asked for.
   */
  readonly numGroups: number;

  /**
   * How many pairs the populations make, which is the columns of `f2Groups`
   * and the length of the `distVector` of every measure.
   */
  readonly numPairs: number;

  /**
   * One group for each of `numGroups`, in the order the groups were started:
   * the chromosome of its variants and the position of its first and of its
   * last one, both included.
   */
  readonly groupIds: readonly PopDistGroup[];

  /**
   * How many variants the pass gave, after the steps of the `Variants`, and
   * what each filter of it was given and kept.
   */
  readonly passStats: PassStats;
}

/**
 * How far apart every pair of the populations of `pops` is, in one pass over
 * `variants`.
 *
 * The measures are functions of the same three counts of a population at a
 * variant, how often each allele was called there, how many genotypes were
 * called whole and how many of those are heterozygous, so a user who wants
 * to compare two of them pays for one reading of the variants and not two.
 * Every allele counts as itself: a variant of three alleles is not collapsed
 * to the commonest one against the rest, which is what lets the same call
 * serve microsatellites and single nucleotide polymorphisms.
 *
 * It is a consumer of the `variants`: it makes one pass over the source
 * through the steps the `Variants` has when it is called, and the `Variants`
 * is as it was afterwards.
 *
 * `pops` is an object of population name to the names of its individuals,
 * which are looked up among the individuals the pass gives. Fewer than two
 * populations is an `Error`, since every measure is of a pair. An individual
 * in two populations is taken and counted in each of them, and one in none
 * takes no part.
 *
 * A pair that counted no variant has NaN for every measure, which is not an
 * error: the other pairs may have values. f_2 and F_ST can come out
 * negative, for one variant and for a whole dataset, and popnei does not
 * clamp them: it is the correction doing its work, and a user who sees a
 * small negative F_ST has two populations this dataset cannot tell apart.
 *
 * It is Python's `calc_pop_dists` under the names of this package, and the
 * measures of pyNei's `calc_jost_dest_pop_dists`, which calculates Jost's D
 * and nothing else, with these differences: popnei has one function for all
 * the measures, because the pass is the cost and the counts are shared; the
 * populations are in the order of the keys of `pops`, where pyNei sorts
 * their names, which changes no value; there is no `numThreads`, since wasm
 * has one thread; there is no `alleles`, since popnei counts the alleles
 * each variant has and lines nothing up across blocks; and
 * `minNumIndividuals` is pyNei's `min_num_samples` under the word of the
 * glossary.
 *
 * @throws {Error} When `variants` is not a `Variants` or was freed; when
 * `options` holds no `jackknifeGroup`, which TypeScript refuses at build
 * time and which a page that calls from JavaScript reaches, and when that
 * `jackknifeGroup` is neither `null`, nor `"variant"`, nor a whole length in
 * base pairs of 1 or more; when `measures` is not an array of names, when a
 * name of it is of no measure, when it names none at all and when it names
 * one of the two that are not calculated yet; when `pops` is not an object
 * of names to arrays of names, when a population names an individual the
 * pass does not give, names one twice or names none, and when it holds fewer
 * than two populations; when `minNumIndividuals` is not a whole number of 0
 * or more; when the source cannot be read, a wrong line of a VCF among the
 * causes; when the pass gives no variant, whether the source holds none or
 * the steps kept none; when the variants fall into fewer than 20 resampling
 * groups, with how many they fell into; when the memory of the tab does not
 * take the six sums popnei keeps for each pair and group; and when `init`
 * has not been awaited.
 */
export function calcPopDists(
  variants: Variants,
  pops: Record<string, readonly string[]>,
  options: CalcPopDistsOptions,
): PopDists {
  theWasmHasToBeLoaded();
  const { source, steps } = sourceOfTheVariants("variants", variants);
  // The resampling groups are checked first, because the argument has no
  // default and a call that left it out has nothing else worth telling its
  // writer about.
  const group = theJackknifeGroup(options);
  const measures = theMeasures(options.measures);
  const thePops = popsOfTheObject(pops);
  const minNumIndividuals =
    options.minNumIndividuals === undefined
      ? defaultMinNumIndividuals()
      : wholeNumberOfZeroOrMore("minNumIndividuals", options.minNumIndividuals);
  // The steps of the pass are a copy of the list, made after every argument
  // was checked so that nothing refused here leaves one behind: the call
  // takes it over and frees it.
  const calculated = source.calc_pop_dists(
    steps.of_a_pass(),
    thePops.names,
    thePops.individuals,
    thePops.numIndividualsPerPop,
    measures,
    group.perVariant,
    group.basePairs,
    minNumIndividuals,
  );
  // Every array is copied out of the memory of wasm as it is read, and the
  // result holds that memory until it is freed, which is here: what the user
  // gets are the copies.
  try {
    const popNames = calculated.pop_names();
    const numPairs = calculated.num_pairs();
    const values = calculated.values();
    const standardErrors = calculated.standard_errors();
    const numVars = calculated.num_vars();
    if (
      popNames === undefined ||
      values === undefined ||
      numVars === undefined
    ) {
      throw new Error(
        "popnei: the pass gave no populations, no measures or no counts of " +
          "the variants of each pair",
      );
    }
    const passStats = passStatsOf(calculated.pass_stats());
    const ofEachMeasure = new Map<PopDistMeasure, Distances>();
    for (const [which, measure] of measures.entries()) {
      const first = which * numPairs;
      const pastTheLast = first + numPairs;
      ofEachMeasure.set(
        measure,
        new Distances(
          values.subarray(first, pastTheLast),
          popNames,
          passStats,
          standardErrors === undefined
            ? null
            : standardErrors.subarray(first, pastTheLast),
        ),
      );
    }
    return {
      pops: Object.freeze(popNames),
      fst: ofEachMeasure.get("fst") ?? null,
      f2: ofEachMeasure.get("f2") ?? null,
      chord: ofEachMeasure.get("chord") ?? null,
      da: ofEachMeasure.get("da") ?? null,
      dest: ofEachMeasure.get("dest") ?? null,
      gst: ofEachMeasure.get("gst") ?? null,
      gstStandardized: ofEachMeasure.get("gst_standardized") ?? null,
      numVars,
      f2Groups: calculated.f2_groups() ?? null,
      numGroups: calculated.num_groups(),
      numPairs,
      groupIds: Object.freeze(
        theGroups(
          calculated.group_chroms(),
          calculated.group_starts(),
          calculated.group_ends(),
        ),
      ),
      passStats,
    };
  } finally {
    calculated.free();
  }
}

/**
 * The resampling groups of a pass, in the order they were started, each with
 * the name of its chromosome and the two positions it runs between.
 *
 * They cross as the three arrays this takes, one for each of those, because
 * an array of objects is not one of the types wasm-bindgen carries.
 *
 * @throws {Error} When the three do not hold one value for each group, which
 * is a defect of popnei.
 */
function theGroups(
  chroms: string[] | undefined,
  starts: Float64Array | undefined,
  ends: Float64Array | undefined,
): PopDistGroup[] {
  if (chroms === undefined || starts === undefined || ends === undefined) {
    throw new Error(
      "popnei: the pass gave the resampling groups without their " +
        "chromosomes or without their positions",
    );
  }
  return chroms.map((chrom, group) => {
    const start = starts[group];
    const end = ends[group];
    if (start === undefined || end === undefined) {
      throw new Error(
        `popnei: the resampling group ${group} of the chromosome ${chrom} ` +
          "has no first or no last position",
      );
    }
    return { chrom, start, end };
  });
}

/**
 * The names of the measures a user asked for, each once and in the order
 * they named them, and the seven of them when they named none.
 *
 * Which names there are is the binding crate's rule, as the names of the
 * statistics are, and which of them have a value today is the core's: what
 * is refused here is what is no array of names, an array of none, which
 * would make a pass over the whole source for nothing, and a measure that
 * has no value yet.
 *
 * @throws {Error} When `measures` is not an array of names, when it names
 * none, and when it names one of the two that are not calculated yet.
 */
function theMeasures(
  measures: readonly PopDistMeasure[] | undefined,
): PopDistMeasure[] {
  const asked =
    measures === undefined
      ? [...THE_MEASURES]
      : (namesOf("measures", measures, {
          oneOfThem: "measure",
          anExample: "fst",
        }) as PopDistMeasure[]);
  if (asked.length === 0) {
    throw new Error(
      "popnei: `measures` names no measure, and a result holds the ones " +
        "that were asked for: leave `measures` out for every measure there is",
    );
  }
  const askedFor = [...new Set(asked)];
  // Which of the seven have a value today is the core's, so that the
  // measures work package 3 of `docs/plans/dists-pops.md` adds are added
  // there and not here as well.
  const haveAValue = measuresThatHaveAValue() as PopDistMeasure[];
  // A name that is of none of the seven goes on to the binding crate, which
  // refuses it with the seven: which names there are is the core's rule too.
  const notWrittenYet = askedFor.filter(
    (measure) =>
      THE_MEASURES.includes(measure) && !haveAValue.includes(measure),
  );
  if (notWrittenYet.length > 0) {
    throw new Error(
      `popnei: ${named(notWrittenYet)} ` +
        `${notWrittenYet.length === 1 ? "is" : "are"} not calculated yet, ` +
        `and what popnei calculates today is ` +
        `${named(haveAValue)}: ask for those`,
    );
  }
  return askedFor;
}

/**
 * How the variants are cut into the resampling groups, out of the
 * `jackknifeGroup` a user wrote: whether each variant is a group of its own,
 * and the length in base pairs otherwise, which is nothing when they asked
 * for no standard errors.
 *
 * The argument has no default and the call fails without it, which
 * TypeScript refuses at build time; a page may call the same function from
 * JavaScript, which has no declarations, so it is refused here too. Whether
 * a length is one a group can have is the binding crate's rule, which writes
 * the three kinds in its message.
 *
 * @throws {Error} When `options` is no object, when it holds no
 * `jackknifeGroup`, and when that `jackknifeGroup` is neither `null`, nor a
 * number, nor the word `"variant"`.
 */
function theJackknifeGroup(options: CalcPopDistsOptions): {
  perVariant: boolean;
  basePairs: number | undefined;
} {
  const noGroup = () =>
    new Error(
      "popnei: `jackknifeGroup` says how the variants are cut into the " +
        "resampling groups the standard errors are built from, and it has " +
        "no default: write a length in base pairs of 1 or more, " +
        '`"variant"` for a group of each variant, or `null` for no standard ' +
        "error",
    );
  if (
    typeof options !== "object" ||
    options === null ||
    !("jackknifeGroup" in options)
  ) {
    throw noGroup();
  }
  const given: unknown = options.jackknifeGroup;
  if (given === null) {
    return { perVariant: false, basePairs: undefined };
  }
  if (given === "variant") {
    return { perVariant: true, basePairs: undefined };
  }
  if (typeof given !== "number") {
    throw new Error(
      "popnei: `jackknifeGroup` is a length in base pairs of 1 or more, " +
        '`"variant"` for a group of each variant, or `null` for no standard ' +
        `error, and ${whatWasGiven(given)} was given`,
    );
  }
  return { perVariant: false, basePairs: given };
}

/**
 * `measures` in one sentence, each in backticks, the last one after an
 * "and".
 */
function named(measures: readonly string[]): string {
  const inBackticks = measures.map((measure) => `\`${measure}\``);
  if (inBackticks.length === 1) {
    return inBackticks[0] as string;
  }
  return `${inBackticks.slice(0, -1).join(", ")} and ${
    inBackticks[inBackticks.length - 1] as string
  }`;
}
