// How long the matrix of r² of every pair of the variants of a vars file
// takes under node on the WebAssembly of this package.
//
// r² is the square of the correlation between the dosages of two variants,
// where the dosage of an individual is how many alleles of its genotype are
// not the major allele of that variant, and `calcRogersHuffR2Matrix` gives
// it for every pair of the variants of one pass, as a square matrix of
// `Float64Array`. In a browser there is no BLAS and no thread: the matrix
// products run on faer with the 128 bit vector instructions of WebAssembly
// that `.cargo/config.toml` passes, on the one thread a page has.
//
//     node bench/time_r2_matrix.mjs <path to a vars file> [--runs n]
//         [--max-num-vars n]
//
// from `js/popnei`, after `npm run build`, which builds the WebAssembly this
// reads. `--runs` is 5 and `--max-num-vars` is 5000 when they are not given,
// which is the cap a user gets when they name none and 200 MB of matrix.
//
// It times two passes over the same source in each run, as
// `crates/popnei/benches/r2_matrix.rs` does natively. The first is
// `calcRogersHuffR2Matrix`, which reads the source and computes the matrix.
// The second is `iterBlocks` with the chromosome and the position asked for
// beside the genotypes, which are the fields the matrix asks its reader for,
// with nothing done to the genotypes but adding up how many alleles came
// out. The difference of the two is the products and the arithmetic over
// them, free of the reading.
//
// The pass gives every variant of the file: the call refuses a pass of more
// variants than `--max-num-vars` instead of stopping at that number, and
// nothing a TypeScript user has cuts a pass short, so the file itself is
// what holds the variants the matrix is of. The native benchmark puts
// `Reblock` in front of the file and gives one block on, which TypeScript
// has no way to do.
//
// Reading the bytes off the disc and `openVars` are outside both timings, as
// they are in `time_kosman_dists.mjs`: node reads the file where a page
// would get the bytes from a `File` or from the network, and `openVars`
// copies them into the memory of WebAssembly and reads the schema and the
// footer once for every run that follows. What it took is printed above the
// runs. One `Variants` is opened and every run reads it again from its
// start, so no run pays another's copy of the bytes.
//
// The memory of the matrix is inside the time of the matrix and not outside
// it, in WebAssembly and in the heap of JavaScript both: the call asks for 8
// bytes a pair, 200 MB at 5000 variants, and wasm-bindgen writes that into
// the `Float64Array` of the JavaScript heap a user gets. The binding copies
// nothing inside WebAssembly, the core handing its own allocation over, so
// the module holds the matrix once.
//
// WebAssembly never gives memory back to the host, so what the module holds
// is the high water mark of everything that has run in it. The byte length
// of that memory is printed after `openVars`, after the run that is not
// timed and after every timed run: it is reached by the run that is not
// timed and no timed run adds to it.
//
// It prints the checksum the native benchmark prints, what the values of the
// matrix that are numbers add to and how many of them are not a number, so
// that the browser and the native build can be shown to give the same
// matrix. It is taken after the clock stops. A pair has no r² when either of
// its variants has one dosage among its called genotypes, and that pair is
// NaN, so the two numbers together read every value of the matrix.
//
// The vars file of the first 5000 variants of the 100000 the native
// benchmark reads is written from the first 5000 variants of the VCF those
// were written from:
//
//     head -n 5005 big.vcf > big5000.vcf
//     python -c "import popnei; popnei.write_vars( \
//         popnei.open_vcf('big5000.vcf'), 'big5000.vars')"
//
// where 5005 is the 5 lines of the header of that VCF and 5000 variants
// after them. A VCF that `make_big_vcf.py` writes for 5000 variants holds
// other genotypes: its generator draws every variant of the file at once, so
// no such file is the beginning of another.

import { readFile } from "node:fs/promises";

import { calcRogersHuffR2Matrix, init, openVars } from "popnei";

// The WebAssembly of the core, which is what the byte length of the memory
// of the module is read from. `init` of the package loads it, and this call
// gives back what that one left.
import theWasmOfTheCore from "../wasm/popnei.js";

// How many times the matrix is timed when the command line does not say.
const DEFAULT_RUNS = 5;

// How many variants the matrix takes before it refuses when the command line
// does not say, which is the cap the core gives a caller who names none.
const DEFAULT_MAX_NUM_VARS = 5000;

// The fields the matrix asks its reader for beside the genotypes, which the
// pass that computes nothing asks for too. The genotypes are in every block.
const THE_FIELDS_OF_THE_MATRIX = ["chrom", "pos"];

const USAGE = `\
node bench/time_r2_matrix.mjs <path to a vars file> [--runs n]
     [--max-num-vars n]

It times \`calcRogersHuffR2Matrix\` over that file under node, on the
WebAssembly that \`npm run build\` leaves in wasm/, and beside it a pass that
reads the same variants with the same fields and computes nothing, so that
what the products and the arithmetic over them take is the difference of the
two.

  --runs n            how many times both passes are timed, 5 by default
  --max-num-vars n    how many variants the matrix takes before it refuses,
                      5000 by default, which is the cap a user gets when
                      they name none and 200 MB of matrix

Every variant of the file goes into the matrix: the call refuses a pass of
more variants than the cap instead of stopping at it. One run of each kind
comes before the timed ones and is not timed. Each run prints its two times
and their difference, how many variants the matrix is of, what the values of
it that are numbers add to and how many are not a number, which is the
checksum the native benchmark prints, and the byte length of the memory of
WebAssembly, which never shrinks.
`;

/** What the command line asked for, or the message that says what it should
 * have been. */
function theArguments(argv) {
  let path;
  let runs = DEFAULT_RUNS;
  let maxNumVars = DEFAULT_MAX_NUM_VARS;
  for (let at = 0; at < argv.length; at += 1) {
    const argument = argv[at];
    if (argument === "--runs" || argument === "--max-num-vars") {
      at += 1;
      const number = Number(argv[at]);
      if (!Number.isInteger(number) || number < 1) {
        throw new Error(`${argument} takes a whole number of 1 or more`);
      }
      if (argument === "--runs") {
        runs = number;
      } else {
        maxNumVars = number;
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
  return { path, runs, maxNumVars };
}

/**
 * What a matrix of r² holds, as the line the run prints: how many variants
 * it is of, what the values of it that are numbers add to and how many of
 * them are not a number.
 *
 * It is taken after the caller has stopped its clock, and the two counts
 * together read every value of the matrix.
 */
function theChecksumOf(matrix) {
  let notANumber = 0;
  let addsTo = 0;
  for (const value of matrix.r2) {
    if (Number.isNaN(value)) {
      notANumber += 1;
      continue;
    }
    addsTo += value;
  }
  return (
    `${matrix.numVars} variants, the ${matrix.r2.length} values of the ` +
    `matrix add to ${addsTo.toFixed(6)} over the ${notANumber} of them that ` +
    `are not a number`
  );
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
  const { path, runs, maxNumVars } = theArguments(process.argv.slice(2));
  const bytes = new Uint8Array(await readFile(path));
  await init();
  const wasm = await theWasmOfTheCore();
  console.log(
    `${path}, ${bytes.length} bytes, popnei in wasm under node ` +
      `${process.version}, one thread, ${runs} runs, the matrix of at most ` +
      `${maxNumVars} variants`,
  );
  console.log(
    "the memory of the matrix, 8 bytes for each pair of variants, is asked " +
      "for inside the time of the matrix and not before it",
  );
  const startedOpening = performance.now();
  const variants = openVars(bytes);
  const openedIn = (performance.now() - startedOpening) / 1000;
  console.log(
    `openVars, which is in neither timing: ${openedIn.toFixed(3)} s, ` +
      theMemoryOfTheWasm(wasm),
  );

  /** One matrix of r² over the variants, timed. */
  const theMatrix = () => {
    const started = performance.now();
    const matrix = calcRogersHuffR2Matrix(variants, { maxNumVars });
    const took = (performance.now() - started) / 1000;
    return { took, did: theChecksumOf(matrix) };
  };

  /** The same variants read with the same fields and nothing computed,
   * timed. The alleles are added up so that a reader which filled the
   * genotypes only when somebody read them would show. */
  const theReadingAlone = () => {
    const started = performance.now();
    const blocks = variants.iterBlocks({ fields: THE_FIELDS_OF_THE_MATRIX });
    let alleles = 0;
    for (const block of blocks) {
      alleles += block.gts.length;
    }
    const took = (performance.now() - started) / 1000;
    return {
      took,
      did: `${blocks.passStats.numVars} variants, ${alleles} alleles`,
    };
  };

  const firstMatrix = theMatrix();
  const firstReading = theReadingAlone();
  console.log(
    `the first run of each kind, which is not timed: the matrix in ` +
      `${firstMatrix.took.toFixed(3)} s, reading in ` +
      `${firstReading.took.toFixed(3)} s, ${firstMatrix.did}, and the pass ` +
      `that computed nothing read ${firstReading.did}`,
  );
  console.log(
    `after the run that is not timed: ${theMemoryOfTheWasm(wasm)}, and the ` +
      `process holds ${theMemoryOfTheProcess()}`,
  );

  const ofTheMatrix = [];
  const ofTheReading = [];
  const ofTheDifference = [];
  for (let run = 1; run <= runs; run += 1) {
    const matrix = theMatrix();
    const reading = theReadingAlone();
    // The two passes read the same variants through the same source, so the
    // reading is the shorter of the two and the difference is a time; a
    // machine that was doing something else during one of them and not the
    // other could turn it around, and 0 is what that gives.
    const difference = Math.max(matrix.took - reading.took, 0);
    ofTheMatrix.push(matrix.took);
    ofTheReading.push(reading.took);
    ofTheDifference.push(difference);
    console.log(
      `run ${run}: the matrix ${matrix.took.toFixed(3)} s, reading ` +
        `${reading.took.toFixed(3)} s, the matrix less the reading ` +
        `${difference.toFixed(3)} s, ${matrix.did}, and the pass that ` +
        `computed nothing read ${reading.did}, ${theMemoryOfTheWasm(wasm)}, ` +
        `the process holds ${theMemoryOfTheProcess()}`,
    );
  }
  console.log(theSpreadOf("the matrix", ofTheMatrix));
  console.log(theSpreadOf("reading", ofTheReading));
  console.log(theSpreadOf("the matrix less the reading", ofTheDifference));
  console.log(`after the last run: ${theMemoryOfTheWasm(wasm)}`);
  variants.free();
}

await main();
