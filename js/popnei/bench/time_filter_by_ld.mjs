// What the filter by linkage disequilibrium costs a whole pass over a vars
// file under node, on the WebAssembly of this package.
//
// The filter keeps the variants that do not repeat what a variant kept near
// them already said: r² is the square of the correlation between the dosages
// of two variants, where the dosage of an individual is how many alleles of
// its genotype are not the major allele of that variant, and a variant is
// dropped when its r² is above `--max-allowed-r2` against any variant the
// filter has kept within `--max-dist` base pairs behind it on its
// chromosome. Those kept variants are its window. In a browser there is no
// BLAS and no thread: the products run on faer with the 128 bit vector
// instructions of WebAssembly that `.cargo/config.toml` passes, on the one
// thread a page has.
//
//     node bench/time_filter_by_ld.mjs <path to a vars file> [--runs n]
//         [--max-allowed-r2 r] [--max-dist d]
//
// from `js/popnei`, after `npm run build`, which builds the WebAssembly this
// reads. `--runs` is 5, `--max-allowed-r2` is 0.3 and `--max-dist` is 250000
// when they are not given.
//
// It times two whole passes in each run, one after the other, as
// `crates/popnei/benches/filter_vars.rs` does natively. The first has the
// filter on it and the second has no filter, and both read every block to
// the end with the genotypes as the only field, which is what a calculation
// over them asks for and what the filter reads. The difference of the two is
// what the filter costs. They run back to back inside each run, so that a
// machine that grows busier over the runs falls on both.
//
// A time for this filter says nothing on its own, because what it costs is
// set by how many variants its window holds, and that is of the dataset and
// not of the filter. So the run that is not timed prints what the window
// held, the mean and the largest of the kept variants that lie within
// `--max-dist` base pairs behind a kept variant. It is counted outside every
// clock, from the chromosome and the position of the variants that untimed
// pass gave, and no timed run asks for either column or counts any of it.
//
// Reading the bytes off the disc and `openVars` are outside both timings, as
// they are in `time_kosman_dists.mjs`: node reads the file where a page
// would get the bytes from a `File` or from the network, and `openVars`
// copies them into the memory of WebAssembly and reads the schema and the
// footer once for every run that follows. What it took is printed above the
// runs. The file is opened twice, because a filter is a step of a
// `Variants` that nothing takes off again: one `Variants` carries the filter
// and the other has no step, and each run reads both of them from their
// start. So the memory of WebAssembly holds the bytes of the file twice.
//
// The alleles of every block are added up inside both clocks, so that the
// pass with no filter, which reads no genotype of its own, cannot come out
// short because the genotypes were never there.
//
// WebAssembly never gives memory back to the host, so what the module holds
// is the high water mark of everything that has run in it. The byte length
// of that memory is printed after each `openVars`, after the runs that are
// not timed and after every timed run.

import { readFile } from "node:fs/promises";

import { init, openVars } from "popnei";

// The WebAssembly of the core, which is what the byte length of the memory
// of the module is read from. `init` of the package loads it, and this call
// gives back what that one left.
import theWasmOfTheCore from "../wasm/popnei.js";

// How many times both passes are timed when the command line does not say.
const DEFAULT_RUNS = 5;

// The largest r² against a variant of its window that keeps a variant, when
// the command line does not say.
const DEFAULT_MAX_ALLOWED_R2 = 0.3;

// How many base pairs behind a variant its window reaches, when the command
// line does not say.
const DEFAULT_MAX_DIST = 250000;

const USAGE = `\
node bench/time_filter_by_ld.mjs <path to a vars file> [--runs n]
     [--max-allowed-r2 r] [--max-dist d]

It times a whole pass over that file with the filter by linkage
disequilibrium on it, and beside it the same pass with no filter, under node
on the WebAssembly that \`npm run build\` leaves in wasm/. Both read every
block to the end with the genotypes as the only field. The difference of the
two is what the filter costs.

  --runs n             how many times both passes are timed, 5 by default
  --max-allowed-r2 r   the largest r² against a variant of its window that
                       keeps a variant, a number from 0 to 1, 0.3 by default
  --max-dist d         how many base pairs behind a variant its window
                       reaches, 1 or more, 250000 by default

One pass of each kind comes before the timed ones and is not timed; the one
with the filter also counts what the window held, which no timed pass does.
Each run prints its two times, how many variants each pass gave, how many
the filter was given and kept, and the byte length of the memory of
WebAssembly, which never shrinks.
`;

/** What the command line asked for, or the message that says what it should
 * have been. */
function theArguments(argv) {
  let path;
  let runs = DEFAULT_RUNS;
  let maxAllowedR2 = DEFAULT_MAX_ALLOWED_R2;
  let maxDist = DEFAULT_MAX_DIST;
  for (let at = 0; at < argv.length; at += 1) {
    const argument = argv[at];
    if (argument === "--runs" || argument === "--max-dist") {
      at += 1;
      const number = Number(argv[at]);
      if (!Number.isInteger(number) || number < 1) {
        throw new Error(`${argument} takes a whole number of 1 or more`);
      }
      if (argument === "--runs") {
        runs = number;
      } else {
        maxDist = number;
      }
    } else if (argument === "--max-allowed-r2") {
      at += 1;
      const number = Number(argv[at]);
      if (!Number.isFinite(number) || number < 0 || number > 1) {
        throw new Error(`${argument} takes a number from 0 to 1`);
      }
      maxAllowedR2 = number;
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
  return { path, runs, maxAllowedR2, maxDist };
}

/**
 * One whole pass over `variants` with the genotypes as the only field, read
 * to its last block and timed.
 *
 * The alleles are added up inside the clock, so that a pass which reads no
 * genotype of its own cannot come out short because the genotypes were never
 * filled. What each filter of the pass was given and kept is read after the
 * clock stops, off the counts the pass has already built.
 */
function onePass(variants) {
  const started = performance.now();
  const blocks = variants.iterBlocks({ fields: [] });
  let alleles = 0;
  for (const block of blocks) {
    alleles += block.gts.length;
  }
  const took = (performance.now() - started) / 1000;
  const stats = blocks.passStats;
  let did = `${stats.numVars} variants, ${alleles} alleles`;
  for (const [kind, counts] of Object.entries(stats.filtering)) {
    did +=
      `, the ${kind} filter was given ${counts.varsProcessed} and kept ` +
      `${counts.varsKept}`;
  }
  return { took, did };
}

/**
 * Where every variant of one pass over `variants` lies, as the chromosome
 * and the position of each of them in the order the pass gave them.
 *
 * It is what the window is counted from, and it is asked for in the pass
 * that is not timed and in no other: the two columns are fields a timed pass
 * does not ask for and holding them is work no timed pass does.
 */
function whereTheVariantsLie(variants) {
  const lie = [];
  for (const block of variants.iterBlocks({ fields: ["chrom", "pos"] })) {
    for (let at = 0; at < block.numVars; at += 1) {
      lie.push({ chrom: block.chrom[at], pos: block.pos[at] });
    }
  }
  return lie;
}

/**
 * What the window held over a pass whose variants are `lie`: for each of
 * them, how many of them lie within `maxDist` base pairs behind it on its
 * chromosome, as the mean of those counts and the largest of them.
 *
 * This is not what the filter computes. The filter drops a variant as soon
 * as one variant of its window is too close a match to it, so it works out
 * fewer values of r² than these counts; what the counts say is how many
 * variants the window held, which is what a time for the filter has to be
 * read against.
 *
 * The filter refuses a variant that does not come after the one before it,
 * so the variants of a chromosome come together and in the order of their
 * positions: the ones within `maxDist` behind the variant at a place are
 * those from the first that is near enough up to that place, and one index
 * that never goes back walks the whole pass.
 */
function theWindowOver(lie, maxDist) {
  if (lie.length === 0) {
    return undefined;
  }
  let firstOfTheWindow = 0;
  let heldTogether = 0;
  let largest = 0;
  for (let place = 0; place < lie.length; place += 1) {
    const variant = lie[place];
    while (firstOfTheWindow < place) {
      const behind = lie[firstOfTheWindow];
      if (
        behind.chrom === variant.chrom &&
        variant.pos - behind.pos <= maxDist
      ) {
        break;
      }
      firstOfTheWindow += 1;
    }
    const held = place - firstOfTheWindow;
    heldTogether += held;
    largest = Math.max(largest, held);
  }
  return {
    mean: heldTogether / lie.length,
    largest,
    ofTheVars: lie.length,
  };
}

/** The byte length of the memory of WebAssembly, and what it is in MB. It
 * never shrinks: a module gives no page back to the host. */
function theMemoryOfTheWasm(wasm) {
  const bytes = wasm.memory.buffer.byteLength;
  return `${bytes} bytes of WebAssembly, ${(bytes / 1e6).toFixed(0)} MB`;
}

/** The memory the node process holds, in GB, which is the memory of
 * WebAssembly and everything else node has. */
function theMemoryOfTheProcess() {
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

/** The line that says what `what` came to over the runs: its best, its
 * median and its worst. */
function theSpreadOf(what, times) {
  return (
    `${what}: best ${sortedTime(times, 0).toFixed(3)} s, ` +
    `median ${medianTime(times).toFixed(3)} s, ` +
    `worst ${sortedTime(times, times.length - 1).toFixed(3)} s`
  );
}

async function main() {
  const { path, runs, maxAllowedR2, maxDist } = theArguments(
    process.argv.slice(2),
  );
  const bytes = new Uint8Array(await readFile(path));
  await init();
  const wasm = await theWasmOfTheCore();
  console.log(
    `${path}, ${bytes.length} bytes, popnei in wasm under node ` +
      `${process.version}, one thread, ${runs} runs, the filter by linkage ` +
      `disequilibrium at ${maxAllowedR2} over ${maxDist} base pairs`,
  );
  const startedOpening = performance.now();
  const filtered = openVars(bytes);
  filtered.filterByLd(maxAllowedR2, maxDist);
  const unfiltered = openVars(bytes);
  const openedIn = (performance.now() - startedOpening) / 1000;
  console.log(
    `the two openVars, which are in neither timing: ${openedIn.toFixed(3)} ` +
      `s together, ${theMemoryOfTheWasm(wasm)}`,
  );

  const firstFiltered = onePass(filtered);
  const firstUnfiltered = onePass(unfiltered);
  console.log(
    `the first pass of each kind, which is not timed: with the filter ` +
      `${firstFiltered.took.toFixed(3)} s, ${firstFiltered.did}; with no ` +
      `filter ${firstUnfiltered.took.toFixed(3)} s, ${firstUnfiltered.did}`,
  );
  const window = theWindowOver(whereTheVariantsLie(filtered), maxDist);
  console.log(
    window === undefined
      ? "the window of that pass: it gave no variant"
      : `the window of a pass with that filter, counted outside every ` +
          `clock: ${window.mean.toFixed(1)} variants on average and ` +
          `${window.largest} at most, of the ${window.ofTheVars} variants it ` +
          `gave, within ${maxDist} base pairs behind each of them`,
  );
  console.log(
    `after the passes that are not timed: ${theMemoryOfTheWasm(wasm)}, and ` +
      `the process holds ${theMemoryOfTheProcess()}`,
  );

  const withTheFilter = [];
  const withNoFilter = [];
  for (let run = 1; run <= runs; run += 1) {
    const one = onePass(filtered);
    const other = onePass(unfiltered);
    withTheFilter.push(one.took);
    withNoFilter.push(other.took);
    console.log(
      `run ${run}: with the filter ${one.took.toFixed(3)} s, ${one.did}; ` +
        `with no filter ${other.took.toFixed(3)} s, ${other.did}; ` +
        `${theMemoryOfTheWasm(wasm)}, the process holds ` +
        theMemoryOfTheProcess(),
    );
  }
  console.log(theSpreadOf("with the filter", withTheFilter));
  console.log(theSpreadOf("with no filter", withNoFilter));
  // The two passes read the same file, so the one with no filter is the
  // shorter of the two and the difference is a time; a machine that was busy
  // during one of them and not the other could turn it around, and 0 is what
  // that gives.
  const difference = Math.max(
    sortedTime(withTheFilter, 0) - sortedTime(withNoFilter, 0),
    0,
  );
  console.log(
    `the filter, on the bests: ${difference.toFixed(3)} s over the pass with ` +
      `no filter`,
  );
  console.log(`after the last run: ${theMemoryOfTheWasm(wasm)}`);
  filtered.free();
  unfiltered.free();
}

await main();
