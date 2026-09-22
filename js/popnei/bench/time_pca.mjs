// How long the principal component analysis of the variants of a vars file
// takes under node on the WebAssembly of this package.
//
// It times one thing: `openVars` over the bytes of a vars file and
// `doPcaFromVariants` over what that gives, which is the analysis and the
// reading of the file inside WebAssembly. The bytes are read off the disc
// before the clock starts, because that read is node's and not popnei's: a
// page gets them from a `File` or from the network and the package takes
// them as a `Uint8Array` either way. The time to beat is in "Speed" of
// docs/specs/pca.md: 7 s for 100000 variants of 1000 individuals without
// the vector instructions of WebAssembly and 5 s with them.
//
// It is what task 4.2 of docs/plans/pca.md timed the package with, and
// docs/reports/pca-measurement.md has the numbers with the commands that
// produced them. `crates/popnei/benches/pca_vars.rs` times the same
// analysis natively.
//
//     node bench/time_pca.mjs <path to a vars file> [--runs n]
//         [--num-prin-comps n]
//
// from `js/popnei`, after `npm run build`, which builds the WebAssembly
// this reads. `--runs` is 5 and `--num-prin-comps` is 0 when they are not
// given; with 0 there is no second reading of the file and no weight.
//
// One run comes before the timed ones whose time is not taken, so that the
// timed runs pay neither the growing of the memory of WebAssembly nor the
// code that node compiles on the first call. Every run opens the file
// again and frees what it opened, so no run reads what another left
// behind.
//
// It prints the wall time of each run with what the analysis gave and the
// memory the node process holds after it, which is the memory of
// WebAssembly and everything else node has, and then the best, the
// median and the worst of the times. The best is what the report states,
// since every other process on the machine can only make a run longer.
//
// The vars file of 100000 variants of 1000 individuals is written as
// `crates/popnei/benches/pca_vars.rs` says, and the one of 20000 by the
// same commands with `20000` after the path of the VCF.

import { readFile } from "node:fs/promises";

import { doPcaFromVariants, init, openVars } from "popnei";

// How many times the analysis is timed when the command line does not say.
const DEFAULT_RUNS = 5;

// How many components the weights are given for when the command line does
// not say. With 0 there is no second reading of the file.
const DEFAULT_NUM_PRIN_COMPS = 0;

const USAGE = `\
node bench/time_pca.mjs <path to a vars file> [--runs n] [--num-prin-comps n]

It times \`openVars\` and \`doPcaFromVariants\` over that file under node, on
the WebAssembly that \`npm run build\` leaves in wasm/.

  --runs n              how many times it analyses the file, 5 by default
  --num-prin-comps n    how many components the weights are given for, 0 by
                        default, and with 0 the file is read once
`;

/** What the command line asked for, or the message that says what it should
 * have been. */
function theArguments(argv) {
  let path;
  let runs = DEFAULT_RUNS;
  let numPrinComps = DEFAULT_NUM_PRIN_COMPS;
  for (let at = 0; at < argv.length; at += 1) {
    const argument = argv[at];
    if (argument === "--runs" || argument === "--num-prin-comps") {
      at += 1;
      const number = Number(argv[at]);
      if (!Number.isInteger(number) || number < 0) {
        throw new Error(`${argument} takes a whole number of 0 or more`);
      }
      if (argument === "--runs") {
        runs = number;
      } else {
        numPrinComps = number;
      }
    } else if (argument === "--help" || argument === "-h") {
      throw new Error(USAGE);
    } else if (argument.startsWith("-")) {
      throw new Error(`\`${argument}\` is not an argument\n\n${USAGE}`);
    } else {
      path = argument;
    }
  }
  if (path === undefined) {
    throw new Error(`no file was given\n\n${USAGE}`);
  }
  if (runs < 1) {
    throw new Error("--runs is 1 or more");
  }
  return { path, runs, numPrinComps };
}

/**
 * One whole analysis of the bytes of a vars file: how long it took and the
 * line that says what it gave.
 *
 * The percentage of the variance of the first component is in that line,
 * which is what says the projections were computed.
 */
function oneAnalysis(bytes, numPrinComps) {
  const started = performance.now();
  const variants = openVars(bytes);
  let result;
  try {
    result = doPcaFromVariants(variants, { numPrinComps });
  } finally {
    variants.free();
  }
  const took = (performance.now() - started) / 1000;
  const first = result.explainedVariancePercent[0] ?? Number.NaN;
  const did =
    `${result.passStats.numVars} variants, ${result.usedVars.length} of them ` +
    `with variance, ${result.individuals.length} individuals, ` +
    `${result.numComps} components, the first with ${first.toFixed(3)} per 100 ` +
    `of the variance, weights for ${result.numPrinComps} of them`;
  return { took, did };
}

/** The memory the node process holds, in GB, which is the memory of
 * WebAssembly and everything else node has. */
function memoryOfTheProcess() {
  return `${(process.memoryUsage().rss / 1e9).toFixed(2)} GB`;
}

/** The time at `part` of the times sorted from the shortest to the longest. */
function sortedTime(times, part) {
  return [...times].sort((one, other) => one - other)[part];
}

/** The median of the times: the middle one of an odd number of them, and the
 * middle of the two middle ones of an even number. */
function medianTime(times) {
  const half = Math.floor(times.length / 2);
  const upper = sortedTime(times, half);
  if (times.length % 2 === 1) {
    return upper;
  }
  return (sortedTime(times, half - 1) + upper) / 2;
}

async function main() {
  const { path, runs, numPrinComps } = theArguments(process.argv.slice(2));
  const bytes = new Uint8Array(await readFile(path));
  await init();
  console.log(
    `${path}, ${bytes.length} bytes, node ${process.version}, ${runs} runs, ` +
      `weights for ${numPrinComps} components`,
  );
  const times = [];
  for (let run = 0; run <= runs; run += 1) {
    const { took, did } = oneAnalysis(bytes, numPrinComps);
    if (run === 0) {
      console.log(
        `the first run, which is not timed: ${did}, in ${took.toFixed(3)} s`,
      );
      continue;
    }
    times.push(took);
    console.log(
      `run ${run}: ${took.toFixed(3)} s, ${did}, the process holds ` +
        memoryOfTheProcess(),
    );
  }
  console.log(
    `best ${sortedTime(times, 0).toFixed(3)} s, ` +
      `median ${medianTime(times).toFixed(3)} s, ` +
      `worst ${sortedTime(times, times.length - 1).toFixed(3)} s`,
  );
}

await main();
