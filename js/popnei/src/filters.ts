/**
 * The filters that keep the variants of a dataset that pass a threshold.
 *
 * A filter is a step of a `Variants`: a method of it that adds itself to the
 * list of steps and returns nothing, and that every pass over the source
 * runs inside the Rust core. The three threshold filters are
 * `filterByMissingData`, over the missing rate of a variant, `filterByMaf`,
 * over its major allele frequency, and `filterByObsHet`, over its observed
 * heterozygosity, and each of them keeps the variants whose number is at
 * most the threshold it was given. The fourth, `filterByLd`, compares a
 * variant with the variants kept behind it on its chromosome instead of
 * with a number of its own, and keeps the ones that do not repeat what a
 * variant near them already said. `filterIndividuals` keeps individuals and
 * not variants: it takes the genotypes of the individuals a user names, at
 * every variant, and drops those of the rest.
 *
 * What is here is what a user reads of them: the step that a filter is in
 * the steps of a `Variants`, and the counts that a pass holds for each
 * filter of it.
 */

/**
 * How many variants one filter of one pass was given and kept.
 *
 * They are the counts of that pass alone, one reading of the source from its
 * start: a `Variants` can be given to any number of consumers, and each of
 * them counts its own.
 */
export interface FilteringStats {
  /**
   * The variants the filter was given, which are those that the filter
   * before it in the steps kept, and all of them for the first filter.
   */
  varsProcessed: number;

  /**
   * Those of them that passed its threshold, which are the ones the filter
   * after it was given.
   */
  varsKept: number;
}

/**
 * One step of a `Variants`: what every pass over its source runs.
 *
 * A filter is the only kind of step there is: one of the three over a number
 * of a variant, the one by linkage disequilibrium, or the one that keeps the
 * individuals a user names.
 */
export interface Step {
  /**
   * What the step does: `"missing_data"`, `"maf"`, `"obs_het"`, `"ld"` or
   * `"individuals"`. The kind of a filter of the variants is the name its
   * counts have in the counts of a pass, where the filter of individuals
   * has no entry, since it takes no variant away.
   */
  kind: string;

  /**
   * What the step was given, under the names of the arguments of the method
   * that added it, `{maxAllowedMaf: 0.95}` for a filter of one threshold,
   * `{maxAllowedR2: 0.3, maxDist: 10000}` for the filter by linkage
   * disequilibrium, which takes two, and `{individuals: ["ind05", "ind00"]}`
   * for the filter of individuals, whose names are in the order they were
   * given.
   *
   * The values are what the arguments of that method take, a number for a
   * threshold and for a window and an array of names for the individuals,
   * so they arrive here as an `unknown` and a user who does arithmetic with
   * a threshold narrows it first: `const maf = step.args["maxAllowedMaf"];
   * if (typeof maf === "number") ...`.
   */
  args: Record<string, unknown>;
}
