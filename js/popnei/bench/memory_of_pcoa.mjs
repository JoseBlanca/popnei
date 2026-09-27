/**
 * How far the principal coordinates of the variants of n individuals grow
 * the memory of wasm under node, and whether they run or end the module.
 *
 * It is what task 3.1 of docs/plans/pcoa.md measured the limit of the
 * browser with: "How it runs" of "The principal coordinates of distances"
 * of docs/specs/pca.md counts at most 56.8 bytes for each cell of the
 * individuals x individuals matrix, 8695 individuals in the 4 GiB of a
 * page, and says the number is to be measured through
 * `doPcoaFromVariants` with `correctByLingoes`, as the PCA's was. The
 * numbers, with the command, are in docs/reports/pcoa.md.
 *
 *     node bench/memory_of_pcoa.mjs <num individuals> [--num-vars n]
 *
 * from `js/popnei`, after `npm run build`. The binding refuses more than
 * the limit it holds before the pass starts, so a number of individuals
 * above it is measured only on a build whose limit was lifted, by a local
 * edit of `TENTHS_OF_THE_MATRIX_THE_PRINCIPAL_COORDINATES_HOLD` in
 * `crates/popnei-js/src/pca.rs` that is not committed. Each number of
 * individuals is run in a node process of its own: the memory of wasm
 * never shrinks, and a trap leaves the module unusable.
 *
 * It writes, in memory, a VCF of that many diploid individuals and
 * `--num-vars` biallelic variants, 300 when it is not given. The allele
 * frequency of each variant is drawn between 0.05 and 0.95, each genotype
 * from it as two independent alleles, and 2 in 100 of the genotypes are
 * missing, from a fixed seed so that a number of individuals always gives
 * the same file. `openVcf` opens its bytes as a page opens the bytes of a
 * file, and `doPcoaFromVariants` runs over them with `correctByLingoes`.
 *
 * It prints the number of individuals, whether the analysis ran, trapped
 * (a `WebAssembly.RuntimeError`, which ends the module with no message) or
 * was refused with an `Error` of popnei, and the bytes of the memory of
 * wasm before the VCF was opened, after it was opened and after the
 * analysis. That memory never shrinks, so its size after the analysis is
 * its peak. The growth after the VCF was opened is divided by the 8 n²
 * bytes of an n x n matrix of `f64`, which gives the bytes the analysis
 * held for each cell of that matrix, besides what the open VCF holds.
 */

import { doPcoaFromVariants, init, openVcf } from "popnei";

import loadTheWasm from "../wasm/popnei.js";

// How many variants the VCF holds when the command line does not say.
const DEFAULT_NUM_VARS = 300;

// The share of the genotypes that are missing.
const MISSING_RATE = 0.02;

// The seed of the draws, fixed so that every run of a number of individuals
// analyses the same file.
const SEED = 20260927;

const USAGE = `\
node bench/memory_of_pcoa.mjs <num individuals> [--num-vars n]

It runs doPcoaFromVariants with correctByLingoes over a VCF of that many
diploid individuals and --num-vars biallelic variants, 300 by default, and
prints how far the memory of wasm grew.
`;

/** What the command line asked for, or the message that says what it should
 * have been. */
function theArguments(argv) {
  let numIndividuals;
  let numVars = DEFAULT_NUM_VARS;
  for (let at = 0; at < argv.length; at += 1) {
    const argument = argv[at];
    if (argument === "--num-vars") {
      at += 1;
      numVars = Number(argv[at]);
      if (!Number.isInteger(numVars) || numVars < 1) {
        throw new Error("--num-vars takes a whole number of 1 or more");
      }
    } else if (argument === "--help" || argument === "-h") {
      throw new Error(USAGE);
    } else if (argument.startsWith("-")) {
      throw new Error(`\`${argument}\` is not an argument\n\n${USAGE}`);
    } else {
      numIndividuals = Number(argument);
      if (!Number.isInteger(numIndividuals) || numIndividuals < 2) {
        throw new Error(`the number of individuals is a whole number of 2 or more\n\n${USAGE}`);
      }
    }
  }
  if (numIndividuals === undefined) {
    throw new Error(`no number of individuals was given\n\n${USAGE}`);
  }
  return { numIndividuals, numVars };
}

/** A generator of numbers drawn uniformly in [0, 1), mulberry32, which gives
 * the same numbers for the same seed on every machine. */
function drawsFrom(seed) {
  let state = seed >>> 0;
  return () => {
    state = (state + 0x6d2b79f5) >>> 0;
    let mixed = state;
    mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
    mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
    return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
  };
}

/** The bytes of a VCF of `numIndividuals` diploid individuals and
 * `numVars` biallelic variants, with some genotypes missing. */
function vcfOf(numIndividuals, numVars) {
  const draw = drawsFrom(SEED);
  const names = Array.from({ length: numIndividuals }, (_unused, at) => `i${at}`);
  const lines = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=1>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join("\t")}`,
  ];
  for (let variant = 0; variant < numVars; variant += 1) {
    const frequency = 0.05 + 0.9 * draw();
    const genotypes = new Array(numIndividuals);
    for (let individual = 0; individual < numIndividuals; individual += 1) {
      if (draw() < MISSING_RATE) {
        genotypes[individual] = "./.";
      } else {
        const first = draw() < frequency ? 1 : 0;
        const second = draw() < frequency ? 1 : 0;
        genotypes[individual] = `${first}/${second}`;
      }
    }
    lines.push(
      `1\t${(variant + 1) * 100}\tv${variant}\tA\tC\t.\tPASS\t.\tGT\t${genotypes.join("\t")}`,
    );
  }
  lines.push("");
  return new TextEncoder().encode(lines.join("\n"));
}

/** What the analysis did: it ran, it trapped, or popnei refused it, with
 * the line that says so. */
function theAnalysisOf(variants) {
  const started = performance.now();
  try {
    const result = doPcoaFromVariants(variants, { correctByLingoes: true });
    const took = (performance.now() - started) / 1000;
    return {
      outcome: "ran",
      said:
        `${result.numComps} components, Lingoes' constant ` +
        `${result.lingoesConstant.toPrecision(6)}, ` +
        `${result.passStats.numVars} variants, in ${took.toFixed(1)} s`,
    };
  } catch (error) {
    const took = (performance.now() - started) / 1000;
    const outcome =
      error instanceof WebAssembly.RuntimeError ? "trapped" : "refused";
    return { outcome, said: `${error}, after ${took.toFixed(1)} s` };
  }
}

async function main() {
  const { numIndividuals, numVars } = theArguments(process.argv.slice(2));
  const bytes = vcfOf(numIndividuals, numVars);
  await init();
  const wasm = await loadTheWasm();
  const memoryOfWasm = () => wasm.memory.buffer.byteLength;

  const before = memoryOfWasm();
  const variants = openVcf(bytes);
  const opened = memoryOfWasm();
  const { outcome, said } = theAnalysisOf(variants);
  const after = memoryOfWasm();

  const perCell = (after - opened) / (numIndividuals * numIndividuals);
  console.log(
    `node ${process.version}, ${numIndividuals} individuals, ${numVars} ` +
      `variants, a VCF of ${bytes.length} bytes`,
  );
  console.log(`${outcome}: ${said}`);
  console.log(
    `memory of wasm: ${before} bytes before the VCF was opened, ${opened} ` +
      `after, ${after} after the analysis; it grew ${after - opened} bytes, ` +
      `${perCell.toFixed(2)} bytes for each of the ${numIndividuals}² cells`,
  );
  console.log(
    `row\t${numIndividuals}\t${outcome}\t${before}\t${opened}\t${after}\t${perCell.toFixed(2)}`,
  );
  if (outcome === "ran") {
    variants.free();
  }
}

await main();
