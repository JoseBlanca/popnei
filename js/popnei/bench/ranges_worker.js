/**
 * The web worker that times one pass of popnei over a VCF of a few hundred
 * MB, once over a `File` the page holds, which popnei reads one range of
 * bytes at a time, and once over the same file as a `Uint8Array`, which is
 * what an application does today; and that times the call into JavaScript a
 * range costs, on its own.
 *
 * `FileReaderSync`, which reads a range of a `File`, exists only inside a
 * web worker, so the whole measurement runs in here.
 * `bench/time_ranges_in_chromium.mjs` starts Chromium, opens
 * `ranges_page.html`, which starts this worker, and prints what it answers.
 *
 * The file is the header of `tests/reference/vcf/many.vcf` with the body of
 * that file after it as many times as the driver asks for, which is how
 * `test/browser/cases/many_ranges.ts` builds its file: the body is one block
 * of memory that the `File` reads again for each copy, so the worker holds
 * 117346 bytes of it and the file is as large as it is asked to be. What it
 * costs to read is the same as any other VCF of that size, except that the
 * bytes come from the memory of the browser and not from a disc, which no
 * measurement of a browser can help.
 *
 * The consumer is `calcPerIndividualStats`: it makes one pass, reads the
 * genotypes and adds two counts per individual, so what the clock holds is
 * the reading of the file and not a calculation over it.
 *
 * It is written in JavaScript and not in TypeScript because it is not built:
 * the server of the browser tests serves it as it is.
 */

import { calcPerIndividualStats, init, openVcf } from "../dist/web.js";
import loadTheWasm from "../wasm/popnei.js";

/** Where the header of `many.vcf` ends, which is the line of its columns. */
const NUM_BYTES_OF_THE_HEADER = 617;

/** The variants of one copy of the body of `many.vcf`, its 25 that failed
 * their FILTER among them, which `onlyPassed` false asks a VCF for. */
const NUM_VARS_PER_COPY = 500;

/**
 * How many bytes the calls of one size of range read altogether at most, and
 * how few calls a size is measured over however large it is. The browser
 * gives back an `ArrayBuffer` for each call, and a worker that asked for tens
 * of GB of them is refused the next read.
 */
const MOST_BYTES_OF_A_MEASUREMENT = 256 * 1024 * 1024;
const LEAST_CALLS_OF_A_MEASUREMENT = 16;

/** How many times the whole file is read into an array of bytes in one
 * call, which is what an application does before its first pass today. */
const NUM_RUNS_OF_THE_WHOLE_FILE = 3;

/** The WebAssembly of popnei, whose memory this worker reads. */
const loading = (async () => {
  await init();
  return loadTheWasm();
})();

/** The bytes of `many.vcf`, fetched once for every measurement. */
let manyVcf;

/** How many bytes the memory of wasm holds. */
function memoryOfWasm(wasm) {
  return wasm.memory.buffer.byteLength;
}

/**
 * `bytes` copied into a block of memory of their own, which is what a `File`
 * takes.
 */
function aBlockOfItsOwn(bytes) {
  const block = new ArrayBuffer(bytes.byteLength);
  new Uint8Array(block).set(bytes);
  return block;
}

/** The bytes of the file `path` of the repository, from the server. */
async function bytesOf(path) {
  const answer = await fetch(path);
  if (!answer.ok) {
    throw new Error(`${path}: the server answered ${answer.status}`);
  }
  return new Uint8Array(await answer.arrayBuffer());
}

/**
 * The header of `many.vcf` with its body after it `numCopies` times, as a
 * `File` of the page and as the same bytes in one array.
 *
 * The array is built only when it is asked for: it is a few hundred MB that
 * the measurement over the `File` is there not to need.
 */
async function theFileOf(numCopies, alsoAsBytes) {
  manyVcf ??= await bytesOf("/tests/reference/vcf/many.vcf");
  const header = manyVcf.subarray(0, NUM_BYTES_OF_THE_HEADER);
  const body = manyVcf.subarray(NUM_BYTES_OF_THE_HEADER);
  const numBytes = header.byteLength + body.byteLength * numCopies;
  const pieces = [aBlockOfItsOwn(header)];
  const oneBody = aBlockOfItsOwn(body);
  for (let copy = 0; copy < numCopies; copy += 1) {
    pieces.push(oneBody);
  }
  const file = new File(pieces, "many_repeated.vcf");
  if (file.size !== numBytes) {
    throw new Error(`the file is ${file.size} bytes and not ${numBytes}`);
  }
  if (!alsoAsBytes) {
    return { file, bytes: null, numBytes };
  }
  const bytes = new Uint8Array(numBytes);
  bytes.set(header, 0);
  for (let copy = 0; copy < numCopies; copy += 1) {
    bytes.set(body, header.byteLength + body.byteLength * copy);
  }
  return { file, bytes, numBytes };
}

/**
 * One pass with `calcPerIndividualStats` over `source`, timed, with what the
 * memory of wasm held before it, after the source was opened and after the
 * pass.
 */
function onePassOver(source, wasm) {
  const before = memoryOfWasm(wasm);
  const atTheStart = performance.now();
  const variants = openVcf(source, { onlyPassed: false });
  const opened = performance.now();
  const afterOpen = memoryOfWasm(wasm);
  let numVars;
  try {
    numVars = calcPerIndividualStats(variants).passStats.numVars;
  } finally {
    variants.free();
  }
  // The memory of wasm never shrinks, so what it holds when the pass has
  // ended is the most it held while the pass ran.
  return {
    ms: performance.now() - atTheStart,
    msToOpen: opened - atTheStart,
    numVars,
    memoryBefore: before,
    memoryAfterOpen: afterOpen,
    memoryAfterThePass: memoryOfWasm(wasm),
  };
}

/**
 * Times `numRuns` passes over the file of `numCopies` copies, over the
 * `File` when `over` is `"file"` and over the bytes of the file in one array
 * when it is `"bytes"`.
 */
async function timesThePasses({ over, numCopies, numRuns }) {
  const wasm = await loading;
  const { file, bytes, numBytes } = await theFileOf(numCopies, over === "bytes");
  const source = over === "bytes" ? bytes : file;
  const runs = [];
  for (let run = 0; run < numRuns; run += 1) {
    runs.push(onePassOver(source, wasm));
  }
  return {
    over,
    numBytes,
    numVars: numCopies * NUM_VARS_PER_COPY,
    memoryOfWasmAtTheStart: runs[0].memoryBefore,
    runs,
  };
}

/**
 * Times `numCalls` reads of a range of `numBytes` bytes of a file, each one
 * `Blob.prototype.slice` and one `FileReaderSync.readAsArrayBuffer`, which
 * is what popnei calls into JavaScript for each range it reads.
 *
 * The ranges walk the file from its start and begin again at its start when
 * they reach its end, so no two calls in a row read the same bytes. What it
 * leaves out of what popnei pays for a range is the copy of the
 * `ArrayBuffer` into the memory of wasm, which the timings of a pass hold.
 *
 * Every size is timed twice, once over the file of many pieces the worker
 * builds and once over the same bytes as one block, which is what a file the
 * user picked is: over the file of pieces alone, what the browser spends
 * walking the pieces a range falls in would count as the cost of the call.
 *
 * A size of range whose `numCalls` calls would read more than
 * `MOST_BYTES_OF_A_MEASUREMENT` gets as many calls as fit in that: the
 * browser hands back an `ArrayBuffer` for each call and refuses to read at
 * all once the worker has asked for tens of GB of them.
 */
async function timesTheCalls({ numCalls, numBytesPerCall, numCopies }) {
  await loading;
  const { file, bytes } = await theFileOf(numCopies, true);
  const ofOneBlock = new File([bytes.buffer], "many_repeated.vcf");
  const reader = new FileReaderSync();
  const runs = [];
  for (const numBytes of numBytesPerCall) {
    const calls = Math.max(
      LEAST_CALLS_OF_A_MEASUREMENT,
      Math.min(numCalls, Math.floor(MOST_BYTES_OF_A_MEASUREMENT / numBytes)),
    );
    if (numBytes > file.size) {
      throw new Error(
        `a call of ${numBytes} bytes is more than the ${file.size} bytes of ` +
          "the file, so it would read the file to its end and no further: " +
          "ask for more copies of the body of `many.vcf`",
      );
    }
    const msOf = (over) => {
      let at = 0;
      let readBytes = 0;
      const atTheStart = performance.now();
      for (let call = 0; call < calls; call += 1) {
        if (at + numBytes > over.size) {
          at = 0;
        }
        readBytes += reader.readAsArrayBuffer(
          over.slice(at, at + numBytes),
        ).byteLength;
        at += numBytes;
      }
      const ms = performance.now() - atTheStart;
      if (readBytes !== calls * numBytes) {
        throw new Error(
          `the ${calls} calls of ${numBytes} bytes read ${readBytes} bytes`,
        );
      }
      return ms;
    };
    const ms = msOf(file);
    const msOfOneBlock = msOf(ofOneBlock);
    runs.push({
      numBytesPerCall: numBytes,
      numCalls: calls,
      ms,
      msPerCall: ms / calls,
      msPerCallOfOneBlock: msOfOneBlock / calls,
    });
  }
  // What an application pays today before its first pass: the whole file
  // from the `File` into an array of bytes, in one call, which is what it
  // then gives `openVcf`. No clock of a pass holds it.
  const wholeFile = [];
  for (let run = 0; run < NUM_RUNS_OF_THE_WHOLE_FILE; run += 1) {
    const atTheStart = performance.now();
    const bytes = new Uint8Array(reader.readAsArrayBuffer(file));
    const ms = performance.now() - atTheStart;
    if (bytes.byteLength !== file.size) {
      throw new Error(
        `the whole file came back as ${bytes.byteLength} bytes and not ` +
          `${file.size}`,
      );
    }
    wholeFile.push({ ms });
  }
  return {
    what: "the calls into JavaScript",
    numBytes: file.size,
    runs,
    wholeFile,
  };
}

/** Runs the measurement the driver asked for. */
async function measures(ask) {
  if (ask.what === "passes") {
    return timesThePasses(ask);
  }
  if (ask.what === "calls") {
    return timesTheCalls(ask);
  }
  throw new Error(`\`${ask.what}\` is no measurement of this worker`);
}

self.addEventListener("message", (event) => {
  const { id, ask } = event.data;
  measures(ask).then(
    (measured) => {
      self.postMessage({ id, measured, failure: null });
    },
    (error) => {
      const failure = error instanceof Error ? error.message : String(error);
      self.postMessage({ id, measured: null, failure });
    },
  );
});
