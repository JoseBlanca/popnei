/**
 * How long the Kosman distances of every pair take in WebAssembly under
 * node, and how long reading the same file alone takes.
 *
 * It is what task 3.1 of docs/plans/dists-kosman.md timed the wasm build
 * with, against the 1.43 s of "Speed" of docs/specs/dists.md, which is of
 * the calculation with the reading of the file taken out.
 * `docs/reports/dists-kosman-measurement.md` has the numbers and the load
 * averages they were taken at, and
 * `crates/popnei/benches/time_kosman_dists.py` is the same clock on the
 * native build.
 *
 * It times two passes over one source. The calculation is one call of
 * `calcPairwiseKosmanDists`, which reads every block and gives the distance
 * of every pair. The reading alone is `iterBlocks({fields: []})`, the same
 * pass with the genotypes as the only field the blocks carry, which is what
 * the calculation asks its reader for, with nothing done to the genotypes
 * but adding up how many alleles came out. The difference of the two is
 * what the calculation costs beyond the reader.
 *
 * Reading the file into a `Uint8Array` and `openVars` are outside both
 * timings, as `open_vars` is outside the native ones: `openVars` copies
 * those bytes into the memory of wasm and reads the schema and the footer,
 * and every pass reads them again from their start. What each of the two
 * took is printed under the runs. wasm has no threads, so there is one.
 *
 *     node bench/time_kosman_dists.mjs <path to a vars file> <runs>
 *
 * It is run from `js/popnei`, and it needs `npm run build` before it: that
 * is what compiles the core to WebAssembly in release and writes `wasm/`
 * and `dist/`, which the name `popnei` resolves to here. One pass of each
 * kind runs before the timed ones and is not timed, and the two kinds are
 * run one after the other inside each round, so that a machine that grows
 * busier over the runs falls on both.
 */

import { readFile } from "node:fs/promises";

import { calcPairwiseKosmanDists, init, openVars } from "popnei";

/** One call of `calcPairwiseKosmanDists`, as its pairs and its variants. */
function theCalculation(variants) {
  const distances = calcPairwiseKosmanDists(variants);
  return [distances.distVector.length, distances.passStats.numVars];
}

/**
 * One pass with the genotypes as the only field, as the alleles that came
 * out of it and the variants it took.
 *
 * The alleles are added up so that a reader which stopped filling the
 * genotypes, or filled them only when somebody read them, would show.
 */
function theReadingAlone(variants) {
  const blocks = variants.iterBlocks({ fields: [] });
  let alleles = 0;
  for (const block of blocks) {
    alleles += block.gts.length;
  }
  return [alleles, blocks.passStats.numVars];
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

const [path, runsAsText] = process.argv.slice(2);
if (path === undefined || runsAsText === undefined) {
  console.log("node bench/time_kosman_dists.mjs <path to a vars file> <runs>");
  process.exit(1);
}
const runs = Number.parseInt(runsAsText, 10);
console.log(`${path}, popnei in wasm under node, ${runs} runs, one thread`);

await init();
let started = performance.now();
const bytes = new Uint8Array(await readFile(path));
const readIn = (performance.now() - started) / 1000;
started = performance.now();
const variants = openVars(bytes);
const openedIn = (performance.now() - started) / 1000;

console.log(
  `the first run of each kind, which is not timed: ${theCalculation(variants)}`,
);
theReadingAlone(variants);

const calculations = [];
const readings = [];
for (let run = 1; run <= runs; run += 1) {
  for (const [what, passOfTheRun, times] of [
    ["the calculation", theCalculation, calculations],
    ["the reading alone", theReadingAlone, readings],
  ]) {
    started = performance.now();
    const [counted, numVars] = passOfTheRun(variants);
    const took = (performance.now() - started) / 1000;
    times.push(took);
    console.log(
      `run ${run}, ${what}: ${took.toFixed(3)} s, ${numVars} variants, ${counted}`,
    );
  }
}
console.log(saidAbout("the calculation", calculations));
console.log(saidAbout("the reading alone", readings));
const difference = Math.min(...calculations) - Math.min(...readings);
console.log(
  `the calculation with the reading taken out, on the bests: ` +
    `${difference.toFixed(3)} s`,
);
console.log(
  `reading the file into a Uint8Array, which is in neither: ` +
    `${readIn.toFixed(3)} s; openVars, which is in neither: ` +
    `${openedIn.toFixed(3)} s`,
);
variants.free();
