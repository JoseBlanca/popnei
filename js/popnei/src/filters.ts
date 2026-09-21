/**
 * The filters that keep the variants of a dataset that pass a threshold.
 *
 * A filter is a step of a `Variants`: a method of it that adds itself to the
 * list of steps and returns nothing, and that every pass over the source
 * runs inside the Rust core. The three threshold filters are
 * `filterByMissingData`, over the missing rate of a variant, `filterByMaf`,
 * over its major allele frequency, and `filterByObsHet`, over its observed
 * heterozygosity, and each of them keeps the variants whose number is at
 * most the threshold it was given.
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
 * A filter is the only kind of step there is.
 */
export interface Step {
  /**
   * What the step does, which is the name its counts have in the counts of
   * a pass: `"missing_data"`, `"maf"` or `"obs_het"`.
   */
  kind: string;

  /**
   * What the step was given, under the names of the arguments of the method
   * that added it, `{maxAllowedMaf: 0.95}`.
   *
   * The values are what the argument of that method takes, so that the
   * steps of the later filters, which take other arguments than a
   * threshold, fit in it. The three that are there take a number.
   */
  args: Record<string, unknown>;
}
