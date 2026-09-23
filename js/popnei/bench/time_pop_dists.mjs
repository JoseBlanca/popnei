/**
 * How long the distances between populations take in WebAssembly under node,
 * and how long reading the same file alone takes.
 *
 * Nothing had timed `calcPopDists` in wasm when this was written. "Speed" of
 * docs/specs/dists.md names the dataset it is to be timed on, 100000
 * variants x 1000 individuals, biallelic, 3 in 100 genotypes missing, read
 * from a vars file, with the individuals cut into 3 populations and into 20,
 * so that how the cost grows with the pairs is measured and not guessed: 3
 * populations make 3 pairs and 20 make 190. That section has no number to
 * reach yet, and this is the clock in wasm that gives the first one.
 * `crates/popnei/benches/time_pop_dists.py` is the same clock on the native
 * build.
 *
 * It times two passes over one source. The calculation is one call of
 * `calcPopDists`, which reads every block and gives the seven measures of
 * every pair of populations with the standard error of each beside it. The
 * reading alone is `iterBlocks({fields: ["chrom", "pos"]})`, the same pass
 * carrying the two fields that calculation asks its reader for, with nothing
 * done to the genotypes but adding up how many alleles came out. The
 * difference of the two is what the calculation costs beyond the reader,
 * which is the number "Speed" asks for.
 *
 * The standard errors are resampled over groups of variants, each group a
 * stretch of one chromosome, and the fourth argument is how long one stretch
 * is in base pairs, 1000000 when it is not given. The dataset above is two
 * chromosomes of 50000 variants 1000 base pairs apart, which 1000000 cuts
 * into 100 groups. popnei refuses a pass that falls into fewer than 20
 * groups, so a shorter file needs a shorter stretch: the 5000 variants of
 * `big5000.vars`, 1000 base pairs apart, fall into 5 groups with 1000000,
 * which is refused, and into 25 with 200000.
 *
 * A populations file holds one line for each individual, the name of the
 * individual and the name of its population with a tab between them, which
 * is what `crates/popnei/benches/make_pops.py` writes. The populations are
 * given to `calcPopDists` in the order in which they first appear in the
 * file, and the pairs are in that order too, so a population whose name is a
 * whole number must not be used: JavaScript keeps such a name at the front
 * of an object whatever order it was put in.
 *
 * Reading the vars file into a `Uint8Array` and `openVars` are outside both
 * timings, as they are outside the native ones: `openVars` copies those
 * bytes into the memory of wasm and reads the schema and the footer, and
 * every pass reads them again from their start. Reading the populations file
 * is outside both as well, and is node's own work. What the two reads and
 * `openVars` took is printed under the runs. wasm has no threads, so there
 * is one.
 *
 *     node js/popnei/bench/time_pop_dists.mjs <path to a vars file>
 *         <path to a populations file> <runs> [group length in base pairs]
 *
 * It needs `npm run build` in `js/popnei` before it: that is what compiles
 * the core to WebAssembly in release and writes `wasm/` and `dist/`, which
 * the name `popnei` resolves to here. One pass of each kind runs before the
 * timed ones and is not timed, and the two kinds are run one after the other
 * inside each round, so that a machine that grows busier over the runs falls
 * on both.
 *
 * Beside each timed calculation it prints the variants the pass took, the
 * pairs of the result, how many resampling groups the variants fell into and
 * Hudson's F_ST of the first pair with its standard error, which together
 * show that the call read the genotypes and resampled them. A call whose
 * pairs are not the pairs the populations of the file make is refused with a
 * message, because a run over the wrong pairs is a time nobody can use.
 */

import { readFile } from "node:fs/promises";

import { calcPopDists, init, openVars } from "popnei";

/** How long one resampling group is when the command line does not say it,
 * which cuts the dataset of "Speed" into 100 groups. */
const DEFAULT_GROUP_BASE_PAIRS = 1000000;

const USAGE =
  "node js/popnei/bench/time_pop_dists.mjs <path to a vars file> " +
  "<path to a populations file> <runs> [group length in base pairs]";

/**
 * The populations of a populations file, the name of each population to the
 * names of its individuals, in the order the populations first appear.
 *
 * @throws {Error} When a line is not the name of an individual, a tab and
 * the name of its population, with the file and the line in the message.
 */
function popsOfTheFile(text, path) {
  const ofEachPop = new Map();
  for (const [which, line] of text.split(/\r?\n/).entries()) {
    if (line === "") {
      continue;
    }
    const fields = line.split("\t");
    if (fields.length !== 2 || fields[0] === "" || fields[1] === "") {
      throw new Error(
        `${path}, line ${which + 1}: a line of a populations file is the ` +
          "name of an individual, a tab and the name of its population, " +
          `and this one is \`${line}\``,
      );
    }
    const [individual, pop] = fields;
    const ofThePop = ofEachPop.get(pop);
    if (ofThePop === undefined) {
      ofEachPop.set(pop, [individual]);
    } else {
      ofThePop.push(individual);
    }
  }
  if (ofEachPop.size === 0) {
    throw new Error(`${path}: the populations file holds no individual`);
  }
  return Object.fromEntries(ofEachPop);
}

/**
 * One call of `calcPopDists`, as what it gave and the variants it took.
 *
 * What it gave is the line printed beside the time: the pairs, the
 * resampling groups, and Hudson's F_ST of the first pair with its standard
 * error.
 *
 * @throws {Error} When the pairs of the result are not `numPairsOfThePops`,
 * the pairs the populations make, and when the result carries no F_ST, no
 * standard error or no value for the first pair.
 */
function theCalculation(variants, pops, groupBasePairs, numPairsOfThePops) {
  const popDists = calcPopDists(variants, pops, {
    jackknifeGroup: groupBasePairs,
  });
  if (popDists.numPairs !== numPairsOfThePops) {
    throw new Error(
      `the ${Object.keys(pops).length} populations of the file make ` +
        `${numPairsOfThePops} pairs, and the call gave ${popDists.numPairs}: ` +
        "the run is not over the pairs it was asked for",
    );
  }
  const fst = popDists.fst;
  if (fst === null || fst.standardErrors === null) {
    throw new Error(
      "the call gave no F_ST, or gave it with no standard error, so it did " +
        "not do the work this times",
    );
  }
  const value = fst.distVector[0];
  const standardError = fst.standardErrors[0];
  if (value === undefined || standardError === undefined) {
    throw new Error("the call gave no value for the first pair");
  }
  const [firstPop, secondPop] = popDists.pops;
  return [
    `${popDists.numPairs} pairs, ${popDists.numGroups} resampling groups, ` +
      `fst of ${firstPop} and ${secondPop} ${value.toFixed(6)} ` +
      `+- ${standardError.toFixed(6)}`,
    popDists.passStats.numVars,
  ];
}

/**
 * One pass carrying the chromosome and the position, which are the fields
 * `calcPopDists` asks its reader for when the groups are stretches of base
 * pairs, as the alleles that came out of it and the variants it took.
 *
 * The alleles are added up so that a reader which stopped filling the
 * genotypes, or filled them only when somebody read them, would show.
 */
function theReadingAlone(variants) {
  const blocks = variants.iterBlocks({ fields: ["chrom", "pos"] });
  let alleles = 0;
  for (const block of blocks) {
    alleles += block.gts.length;
  }
  return [`${alleles} alleles`, blocks.passStats.numVars];
}

/**
 * The same pass with nothing read out of the block at all, which is the
 * baseline to subtract.
 *
 * `theReadingAlone` above reads `block.gts`, and a column of a block in
 * wasm is copied when it crosses: over the dataset of "Speed" that is
 * 200 million bytes into fresh typed arrays, and reaching the block at all
 * builds a string for the chromosome of every variant, 100000 of them.
 * `calcPopDists` copies none of that: no block of it leaves the core, and
 * it makes one string for each resampling group, 100 of them. So
 * subtracting `theReadingAlone` takes out more than the reader and leaves
 * the calculation looking faster than it is. This pass reads only what the
 * reader counted, so the bytes stay in wasm, and it is the honest
 * subtrahend. The performance review of 23 September 2026 found this.
 */
function theReadingWithNothingRead(variants) {
  const blocks = variants.iterBlocks({ fields: ["chrom", "pos"] });
  let numBlocks = 0;
  for (const _block of blocks) {
    numBlocks += 1;
  }
  return [`${numBlocks} blocks`, blocks.passStats.numVars];
}

/** The best, the median and the worst of `times`, in seconds. */
function saidAbout(what, times) {
  const sorted = [...times].sort((one, other) => one - other);
  const median = sorted[(sorted.length - 1) >> 1];
  return (
    `${what}: best ${sorted[0].toFixed(3)} s, ` +
    `median ${median.toFixed(3)} s, ` +
    `worst ${sorted[sorted.length - 1].toFixed(3)} s`
  );
}

const [path, popsPath, runsAsText, groupAsText] = process.argv.slice(2);
if (path === undefined || popsPath === undefined || runsAsText === undefined) {
  console.log(USAGE);
  process.exit(1);
}
const runs = Number.parseInt(runsAsText, 10);
const groupBasePairs =
  groupAsText === undefined
    ? DEFAULT_GROUP_BASE_PAIRS
    : Number.parseInt(groupAsText, 10);
if (!Number.isInteger(runs) || runs < 1) {
  console.log(`the runs are a whole number of 1 or more\n${USAGE}`);
  process.exit(1);
}
if (!Number.isInteger(groupBasePairs) || groupBasePairs < 1) {
  console.log(
    `the length of a resampling group is a whole number of base pairs of 1 ` +
      `or more\n${USAGE}`,
  );
  process.exit(1);
}

let started = performance.now();
const pops = popsOfTheFile(await readFile(popsPath, "utf8"), popsPath);
const popsReadIn = (performance.now() - started) / 1000;
const numPops = Object.keys(pops).length;
const numPairsOfThePops = (numPops * (numPops - 1)) / 2;
console.log(
  `${path}, the ${numPops} populations of ${popsPath}, ` +
    `${numPairsOfThePops} pairs, resampling groups of ${groupBasePairs} ` +
    `base pairs, popnei in wasm under node, ${runs} runs, one thread`,
);

await init();
started = performance.now();
const bytes = new Uint8Array(await readFile(path));
const readIn = (performance.now() - started) / 1000;
started = performance.now();
const variants = openVars(bytes);
const openedIn = (performance.now() - started) / 1000;

console.log(
  "the first run of each kind, which is not timed: " +
    `${theCalculation(variants, pops, groupBasePairs, numPairsOfThePops)}`,
);
theReadingAlone(variants);
theReadingWithNothingRead(variants);

const calculations = [];
const readings = [];
const bareReadings = [];
for (let run = 1; run <= runs; run += 1) {
  for (const [what, passOfTheRun, times] of [
    [
      "the calculation",
      () => theCalculation(variants, pops, groupBasePairs, numPairsOfThePops),
      calculations,
    ],
    ["the reading alone", () => theReadingAlone(variants), readings],
    [
      "the reading with nothing read",
      () => theReadingWithNothingRead(variants),
      bareReadings,
    ],
  ]) {
    started = performance.now();
    const [counted, numVars] = passOfTheRun();
    const took = (performance.now() - started) / 1000;
    times.push(took);
    console.log(
      `run ${run}, ${what}: ${took.toFixed(3)} s, ${numVars} variants, ${counted}`,
    );
  }
}
console.log(saidAbout("the calculation", calculations));
console.log(saidAbout("the reading alone", readings));
console.log(saidAbout("the reading with nothing read", bareReadings));
const difference = Math.min(...calculations) - Math.min(...bareReadings);
console.log(
  `the calculation with the reading taken out, on the bests: ` +
    `${difference.toFixed(3)} s, which subtracts the reading with nothing ` +
    `read: subtracting the reading alone instead would give ` +
    `${(Math.min(...calculations) - Math.min(...readings)).toFixed(3)} s, ` +
    `and that one takes out the copies out of wasm as well as the reader`,
);
console.log(
  `reading the vars file into a Uint8Array, which is in neither: ` +
    `${readIn.toFixed(3)} s; openVars, which is in neither: ` +
    `${openedIn.toFixed(3)} s; reading the populations file, which is in ` +
    `neither: ${popsReadIn.toFixed(3)} s`,
);
variants.free();
