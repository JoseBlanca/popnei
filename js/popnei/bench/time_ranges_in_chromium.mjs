/**
 * What the size of a range costs popnei in a browser: one pass over a VCF of
 * a few hundred MB in Chromium, over a `File` of the page read one range at
 * a time, at each of the sizes of range it is asked for, and over the same
 * file held whole as a `Uint8Array`, which is what an application does
 * today.
 *
 *     node bench/time_ranges_in_chromium.mjs [--copies n] [--runs n]
 *         [--sizes 262144,1048576,4194304,16777216] [--calls n]
 *
 * from `js/popnei`, after `npm run build`. `--copies` is 2570, `--runs` 5,
 * `--calls` 2000 and `--sizes` the four of "Speed" of
 * `docs/specs/js_sources.md` when they are not given.
 *
 * The file is the header of `tests/reference/vcf/many.vcf` with the body of
 * that file after it `--copies` times, built in the worker:
 * 617 + 116729 * `--copies` bytes and 500 * `--copies` variants of 50
 * diploid individuals. 2570 copies make 299994147 bytes and 1285000
 * variants.
 *
 * The size of a range is `NUM_BYTES_PER_RANGE` of
 * `crates/popnei-js/src/source.rs`, a constant of the build and no argument
 * of the API, so there is one build of the WebAssembly for each size: this
 * script writes the size into that constant, runs `npm run build:wasm`,
 * takes the times, and then, whatever happened and on an interrupt as well,
 * writes the file back as it found it and builds the WebAssembly of `wasm/`
 * from it again. Both are needed: the source alone put back leaves `wasm/`
 * built at the last size measured, `git status` says nothing of it, and the
 * next `npm run test:browser` fails at the case of a file of more than one
 * range and `npm test` at the counts of a pass, neither of them naming this
 * script. It is why it is a script of `bench/` and not a test: nothing of a
 * test run edits the source of popnei.
 *
 * Each measurement gets a page and a worker of its own, so the memory of
 * wasm it reports is of that measurement alone: the memory of a module never
 * shrinks, so what it holds when a pass has ended is the most it held while
 * the pass ran, and a second pass in the same worker would read the high
 * water mark of the first.
 *
 * What it prints is a table of the runs and, under it, the best run of each
 * size, which is what `docs/reports/js-sources-measurement.md` carries,
 * and the load averages of the machine before and after, since a machine
 * that is busy makes every row of the table slower.
 */

import { spawn, spawnSync } from "node:child_process";
import { writeFileSync } from "node:fs";
import { readFile, writeFile } from "node:fs/promises";
import { loadavg } from "node:os";
import { fileURLToPath } from "node:url";

import { chromium } from "@playwright/test";

/** Where the constant that sets the size of a range is. */
const SOURCE_RS = new URL(
  "../../../crates/popnei-js/src/source.rs",
  import.meta.url,
);

/** The line of that constant, whose value this script writes. */
const THE_CONSTANT = /^const NUM_BYTES_PER_RANGE: u64 = .+;$/m;

/** The page that starts the worker every measurement runs in. */
const PAGE = "/js/popnei/bench/ranges_page.html";

/**
 * The port the server of the measurement listens on. It is one more than the
 * port of the browser tests, so that the two can run at the same time.
 */
const PORT = 8974;

/** What the arguments of the command are when they are not given. */
const BY_DEFAULT = {
  copies: 2570,
  runs: 5,
  sizes: [256 * 1024, 1024 * 1024, 4 * 1024 * 1024, 16 * 1024 * 1024],
  calls: 2000,
};

/** The sizes of range the cost of one call into JavaScript is measured at. */
const BYTES_OF_A_CALL = [1024, 256 * 1024, 1024 * 1024, 4 * 1024 * 1024, 16 * 1024 * 1024];

/** The arguments of the command, over the ones above. */
function theArguments(argv) {
  const given = { ...BY_DEFAULT };
  for (let at = 0; at < argv.length; at += 2) {
    const name = argv[at]?.replace(/^--/, "");
    const value = argv[at + 1];
    if (!(name in given) || value === undefined) {
      throw new Error(`\`${argv[at]} ${value}\` is no argument of this script`);
    }
    given[name] =
      name === "sizes" ? value.split(",").map(Number) : Number(value);
  }
  return given;
}

/** Where `npm` is run from, the directory of the package. */
const JS_POPNEI = fileURLToPath(new URL("..", import.meta.url));

/** Runs `command` in `js/popnei` and waits for it. */
function runs(command, args) {
  return new Promise((ran, failed) => {
    const child = spawn(command, args, { cwd: JS_POPNEI, stdio: "inherit" });
    child.on("error", failed);
    child.on("exit", (code) => {
      if (code === 0) {
        ran();
      } else {
        failed(new Error(`\`${command} ${args.join(" ")}\` exited ${code}`));
      }
    });
  });
}

/** Builds the WebAssembly with `numBytesPerRange` as the size of a range. */
async function buildsTheWasmAt(numBytesPerRange, theSource) {
  if (!THE_CONSTANT.test(theSource)) {
    throw new Error(
      `\`${fileURLToPath(SOURCE_RS)}\` has no line \`const ` +
        "NUM_BYTES_PER_RANGE: u64 = ...;\`, so this script cannot set the " +
        "size of a range: the constant was renamed or written in another way",
    );
  }
  const patched = theSource.replace(
    THE_CONSTANT,
    `const NUM_BYTES_PER_RANGE: u64 = ${numBytesPerRange};`,
  );
  await writeFile(SOURCE_RS, patched);
  await runs("npm", ["run", "build:wasm"]);
}

/**
 * Writes `theSource` back into `source.rs` and builds the WebAssembly of
 * `wasm/` from it, so that the source of popnei and the artifact the tests
 * of the package load are of the same size of range.
 */
async function putsTheSourceBackAndBuilds(theSource) {
  await writeFile(SOURCE_RS, theSource);
  await runs("npm", ["run", "build:wasm"]);
}

/** Asks the worker of a page of its own for one measurement. */
async function measures(browser, ask) {
  const page = await browser.newPage();
  page.on("pageerror", (error) => {
    process.stderr.write(`the page threw: ${error.message}\n`);
  });
  try {
    await page.goto(`http://127.0.0.1:${PORT}${PAGE}`);
    return await page.evaluate(
      (theAsk) => window.measureInTheWorker(theAsk),
      ask,
    );
  } finally {
    await page.close();
  }
}

/** The best of the runs of one measurement, in milliseconds. */
function bestOf(measured) {
  return Math.min(...measured.runs.map((run) => run.ms));
}

/** What the run that was fastest held in the memory of wasm. */
function memoryOfTheBest(measured) {
  const best = measured.runs.reduce((one, other) =>
    other.ms < one.ms ? other : one,
  );
  return best;
}

/** Starts the server the page and the WebAssembly are read from. */
function startsTheServer() {
  const server = spawn("node", ["test/browser/server.ts", String(PORT)], {
    cwd: fileURLToPath(new URL("..", import.meta.url)),
    stdio: ["ignore", "inherit", "inherit"],
  });
  return server;
}

/** Waits until the server answers the page, or gives up. */
async function theServerAnswers() {
  for (let attempt = 0; attempt < 100; attempt += 1) {
    try {
      const answer = await fetch(`http://127.0.0.1:${PORT}${PAGE}`);
      if (answer.ok) {
        return;
      }
    } catch {
      // The server has not opened its port yet.
    }
    await new Promise((later) => setTimeout(later, 100));
  }
  throw new Error(`the server did not answer \`${PAGE}\` on port ${PORT}`);
}

const given = theArguments(process.argv.slice(2));
const theSource = await readFile(SOURCE_RS, "utf8");

// Ctrl-C ends the process where the `finally` below does not run, and what
// that would leave behind is `source.rs` at the size being measured and the
// WebAssembly of `wasm/` built from it. Both are put back here, with the
// synchronous calls a handler can make.
for (const signal of ["SIGINT", "SIGTERM"]) {
  process.on(signal, () => {
    process.stderr.write(
      `\n${signal}: putting \`source.rs\` back and building the ` +
        "WebAssembly of `wasm/` from it again\n",
    );
    writeFileSync(SOURCE_RS, theSource);
    spawnSync("npm", ["run", "build:wasm"], {
      cwd: JS_POPNEI,
      stdio: "inherit",
    });
    process.exit(1);
  });
}

const server = startsTheServer();
const browser = await chromium.launch();
const table = [];
let calls;
let ofTheWholeFile;

process.stdout.write(`the load averages before: ${loadavg().join(", ")}\n`);
try {
  await theServerAnswers();
  for (const numBytesPerRange of given.sizes) {
    await buildsTheWasmAt(numBytesPerRange, theSource);
    const measured = await measures(browser, {
      what: "passes",
      over: "file",
      numCopies: given.copies,
      numRuns: given.runs,
    });
    table.push({ numBytesPerRange, ...measured });
    process.stdout.write(
      `${numBytesPerRange} bytes per range: ` +
        `${measured.runs.map((run) => run.ms.toFixed(0)).join(", ")} ms, ` +
        `${measured.runs[0].numVars} variants\n`,
    );
    if (calls === undefined) {
      calls = await measures(browser, {
        what: "calls",
        numCalls: given.calls,
        numBytesPerCall: BYTES_OF_A_CALL,
        numCopies: given.copies,
      });
      for (const run of calls.runs) {
        process.stdout.write(
          `one call of ${run.numBytesPerCall} bytes: ` +
            `${(run.msPerCall * 1000).toFixed(1)} us over the file of many ` +
            `pieces, ${(run.msPerCallOfOneBlock * 1000).toFixed(1)} us over ` +
            `the same bytes as one, ${run.numCalls} calls\n`,
        );
      }
      process.stdout.write(
        `the whole file in one call: ` +
          `${calls.wholeFile.map((run) => run.ms.toFixed(0)).join(", ")} ms\n`,
      );
    }
  }
  ofTheWholeFile = await measures(browser, {
    what: "passes",
    over: "bytes",
    numCopies: given.copies,
    numRuns: given.runs,
  });
  process.stdout.write(
    `the whole file as a Uint8Array: ` +
      `${ofTheWholeFile.runs.map((run) => run.ms.toFixed(0)).join(", ")} ms, ` +
      `${ofTheWholeFile.runs[0].numVars} variants\n`,
  );
} finally {
  await browser.close();
  server.kill();
  await putsTheSourceBackAndBuilds(theSource);
}

process.stdout.write(`the load averages after: ${loadavg().join(", ")}\n\n`);
const numBytes = ofTheWholeFile.numBytes;
process.stdout.write(
  `one pass over ${numBytes} bytes, ${ofTheWholeFile.numVars} variants, ` +
    `the best of ${given.runs} runs\n`,
);
for (const measured of table) {
  const best = memoryOfTheBest(measured);
  process.stdout.write(
    `${String(measured.numBytesPerRange).padStart(9)} bytes per range: ` +
      `${bestOf(measured).toFixed(0).padStart(6)} ms, ` +
      `${(numBytes / bestOf(measured) / 1000).toFixed(1)} MB/s, ` +
      `the memory of wasm ${best.memoryAfterThePass} bytes, ` +
      `${best.memoryBefore} before it\n`,
  );
}
const bestWhole = memoryOfTheBest(ofTheWholeFile);
process.stdout.write(
  `the whole file as a Uint8Array: ` +
    `${bestOf(ofTheWholeFile).toFixed(0)} ms, ` +
    `${(numBytes / bestOf(ofTheWholeFile) / 1000).toFixed(1)} MB/s, ` +
    `the memory of wasm ${bestWhole.memoryAfterThePass} bytes, ` +
    `${bestWhole.memoryBefore} before it\n`,
);
process.stdout.write(`\n${JSON.stringify({ table, calls, ofTheWholeFile })}\n`);
