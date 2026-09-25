/**
 * What a source tells the page while a consumer reads it: the calls of
 * `onProgress`, with the bytes the pass has read, the bytes the file holds,
 * which pass of the run is reading and how many passes the run makes.
 *
 * A pass is one reading of a source from its start, a run is one call of one
 * consumer with the passes it makes, and a consumer is a function of the
 * package that reads a `Variants`, the iteration of `iterBlocks` among them,
 * as `docs/glossary.md` has the three words. "What the source tells the
 * page" of `docs/specs/js_sources.md` is what these tests are written from.
 *
 * Three things make a call: the first read of a pass, which says it has read
 * nothing; a read that brings the bytes read since the last call to the
 * 4 MiB of a range; and the end of the run, which makes one call for each of
 * its passes, in the order of their numbers, with the bytes that pass read.
 * The last of the three is what says a pass is over, because no read does: a
 * pass over a vars file stops after its last batch, and a run that fails
 * stops wherever it failed. The files here are `many.vcf` of
 * `docs/specs/io_vcf.md`, 500 variants of 50 individuals in 117346 bytes,
 * and its bgzipped `many.vcf.gz` of 21904, both smaller than one range, so a
 * pass over either is told twice.
 *
 * Nothing here asserts a variant or a statistic: what the counting must not
 * change is the numbers the other tests of the package assert.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { ConsumerName, Progress, Variants } from "popnei";
import {
  calcGwas,
  calcKinship,
  calcLdAndDistPerPop,
  calcPairwiseKosmanDists,
  calcPerIndividualStats,
  calcPerVarDistribs,
  calcPopDists,
  calcPopDiversity,
  calcRogersHuffR2Matrix,
  doPcaFromVariants,
  init,
  numPassesOf,
  openVars,
  openVcf,
  writeVars,
} from "popnei";

import {
  manyVariantsVcf,
  referenceVcf,
  vcfOfDrawnGenotypes,
} from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The bytes of `many.vcf.gz`, which bgzip wrote as four gzip members. */
const MANY_VCF_GZ = await referenceVcf("many.vcf.gz");

/** How many bytes `many.vcf` holds, which `docs/specs/io_vcf.md` gives. */
const BYTES_OF_MANY_VCF = 117346;

/** How many bytes `many.vcf.gz` holds. */
const BYTES_OF_MANY_VCF_GZ = 21904;

/**
 * How many bytes a pass over the vars file written from `many.vcf` reads,
 * which is its footer and its batches and not the whole file.
 *
 * A pass over a vars file reads the last ten bytes, which say how long the
 * footer is, then the footer, which says where the batches are, and then
 * each batch whole. The schema message at the head of the file it never
 * reads, because the footer carries the schema too, so the count ends below
 * the size of the file. It was measured on the file this test writes, with
 * the blocks of the size popnei chooses for 50 individuals.
 */
const BYTES_READ_OF_THE_VARS_FILE = 42552;

/**
 * How many bytes a pass over the vars file of more than one range reads,
 * measured on the file `vcfOfDrawnGenotypes(14000, 600)` writes, 12231602
 * bytes, of which the pass reads all but its schema message.
 *
 * The call that says it is the one the end of the run makes. The reads of
 * that pass after its last range never reach another one, so without that
 * call the page was last told at 8476400 bytes, two thirds of the way, and
 * a bar drawn from the calls stopped there.
 */
const BYTES_READ_OF_THE_LARGE_VARS_FILE = 12225584;

/** The names of the 50 individuals of `many.vcf`, `ind00` to `ind49`. */
const THE_INDIVIDUALS = Array.from(
  { length: 50 },
  (_unused, individual) => `ind${String(individual).padStart(2, "0")}`,
);

/**
 * A number for each individual, for the association study: the values run
 * from 0 to 6 and are not the same for everybody, which a trait a model can
 * be fitted to has to be.
 */
const THE_TRAIT = Object.fromEntries(
  THE_INDIVIDUALS.map((individual, at) => [individual, at % 7]),
);

/** Two populations of 25 individuals each, for the distances between them. */
const THE_POPS = {
  one: THE_INDIVIDUALS.slice(0, 25),
  two: THE_INDIVIDUALS.slice(25),
};

/**
 * The kinship of the 50 individuals of `many.vcf`, for the study that asks
 * for the GRAMMAR-Gamma approximation, which only a mixed model has.
 *
 * It is taken over a `Variants` of its own, so that the pass that builds it
 * is not among the calls the test of a run reads.
 */
const THE_KINSHIP = (() => {
  const variants = openVcf(MANY_VCF);
  try {
    return calcKinship(variants, { transformToBiallelic: true });
  } finally {
    variants.free();
  }
})();

/**
 * The calls the source of `variants` will make, which fill as the consumers
 * of the tests run.
 *
 * The function is set on the `Variants` and holds until it is set again, so
 * one source tells every run over it.
 */
function theCallsOf(variants: Variants): Progress[] {
  const calls: Progress[] = [];
  variants.onProgress((progress) => {
    calls.push(progress);
  });
  return calls;
}

/**
 * That `calls` never says fewer bytes read than the call before it, and
 * never more than the file holds.
 */
function assertTheyRise(calls: readonly Progress[], what: string): void {
  let readSoFar = -1;
  for (const call of calls) {
    assert.ok(
      call.bytesRead >= readSoFar,
      `${what}: a call says ${call.bytesRead} bytes read after one that ` +
        `said ${readSoFar}`,
    );
    assert.ok(
      call.bytesRead <= call.numBytes,
      `${what}: a call says ${call.bytesRead} bytes read of a file of ` +
        `${call.numBytes} bytes`,
    );
    readSoFar = call.bytesRead;
  }
}

test("a pass over a VCF is told from no bytes read to the whole file", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    calcPerVarDistribs(variants);
  } finally {
    variants.free();
  }
  assertTheyRise(calls, "the pass over many.vcf");
  // The first read of the pass and the end of the run: the file is smaller
  // than one range, so no read of it brings the bytes read since the last
  // call to 4 MiB.
  assert.deepEqual(calls, [
    { bytesRead: 0, numBytes: BYTES_OF_MANY_VCF, pass: 1, numPasses: 1 },
    {
      bytesRead: BYTES_OF_MANY_VCF,
      numBytes: BYTES_OF_MANY_VCF,
      pass: 1,
      numPasses: 1,
    },
  ]);
});

test("a pass over a gzipped VCF counts the bytes the file holds on disk", () => {
  const variants = openVcf(MANY_VCF_GZ);
  const calls = theCallsOf(variants);
  try {
    calcPerVarDistribs(variants);
  } finally {
    variants.free();
  }
  assertTheyRise(calls, "the pass over many.vcf.gz");
  // The compressed bytes and not the 117346 of the text they hold: what a
  // page draws a bar against is the file the user picked.
  assert.deepEqual(calls, [
    { bytesRead: 0, numBytes: BYTES_OF_MANY_VCF_GZ, pass: 1, numPasses: 1 },
    {
      bytesRead: BYTES_OF_MANY_VCF_GZ,
      numBytes: BYTES_OF_MANY_VCF_GZ,
      pass: 1,
      numPasses: 1,
    },
  ]);
});

test("a pass over a file of several ranges is told once per range", () => {
  // 150000 variants of 3 individuals, 6188977 bytes, 6.19 MB, which is more
  // than one range of 4 MiB and less than two.
  const vcf = manyVariantsVcf(150000);
  const variants = openVcf(vcf);
  const calls = theCallsOf(variants);
  try {
    calcPerIndividualStats(variants);
  } finally {
    variants.free();
  }
  assertTheyRise(calls, "the pass over a VCF of several ranges");
  assert.ok(
    vcf.length > 4 * 1024 * 1024 && vcf.length < 8 * 1024 * 1024,
    `the VCF holds ${vcf.length} bytes`,
  );
  // The first read, the read that found a range read since it, and the end
  // of the run.
  assert.equal(calls.length, 3, `the calls are ${JSON.stringify(calls)}`);
  assert.equal(calls.at(0)?.bytesRead, 0);
  assert.ok(
    (calls.at(1)?.bytesRead ?? 0) >= 4 * 1024 * 1024,
    `the second call says ${calls.at(1)?.bytesRead} bytes read`,
  );
  assert.equal(calls.at(-1)?.bytesRead, vcf.length);
});

test("a pass over a vars file smaller than a range is told twice", () => {
  const vcf = openVcf(MANY_VCF);
  const file = writeVars(vcf).bytes;
  vcf.free();
  const variants = openVars(file);
  const calls = theCallsOf(variants);
  try {
    calcPerIndividualStats(variants);
  } finally {
    variants.free();
  }
  assertTheyRise(calls, "the pass over a vars file");
  // A pass over a vars file reads the ten last bytes, the footer and each
  // batch whole, every one of them by a read of the length it asks for, so
  // no read of it finds the end of the file. The file holds fewer bytes than
  // one range, so its two calls are the one of its first read and the one
  // the end of the run makes, which says the bytes of its footer and its
  // batches and not the size of the file.
  assert.ok(
    BYTES_READ_OF_THE_VARS_FILE < file.length,
    `the pass read ${BYTES_READ_OF_THE_VARS_FILE} bytes of a file of ${file.length}`,
  );
  assert.deepEqual(calls, [
    { bytesRead: 0, numBytes: file.length, pass: 1, numPasses: 1 },
    {
      bytesRead: BYTES_READ_OF_THE_VARS_FILE,
      numBytes: file.length,
      pass: 1,
      numPasses: 1,
    },
  ]);
});

test("a pass over a vars file of several ranges is told at the end of the run", () => {
  // The vars file of `vars_memory.test.ts`, 14000 variants of 600
  // individuals written in batches of 100, 12231602 bytes, which is more
  // than one range of 4 MiB. A pass over it reads its batches and stops
  // after the last one, with up to a range read since the call before it.
  const vcf = openVcf(vcfOfDrawnGenotypes(14000, 600), { onlyPassed: false });
  const file = writeVars(vcf, { numVarsPerBlock: 100 }).bytes;
  vcf.free();
  const variants = openVars(file);
  const calls = theCallsOf(variants);
  try {
    calcPerIndividualStats(variants);
  } finally {
    variants.free();
  }
  assertTheyRise(calls, "the pass over a vars file of several ranges");
  assert.equal(file.length, 12231602, "the vars file changed size");
  assert.equal(calls.at(0)?.bytesRead, 0);
  assert.equal(calls.at(-1)?.bytesRead, BYTES_READ_OF_THE_LARGE_VARS_FILE);
  assert.ok(
    calls.every((call) => call.numBytes === file.length),
    `a call says another size of the file: ${JSON.stringify(calls)}`,
  );
});

test("the two passes of the pca take their numbers in the order they start", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    doPcaFromVariants(variants, {
      numPrinComps: 10,
      transformToBiallelic: true,
    });
  } finally {
    variants.free();
  }
  assertTheyRise(
    calls.filter((call) => call.pass === 1),
    "the first pass of the pca",
  );
  assertTheyRise(
    calls.filter((call) => call.pass === 2),
    "the second pass of the pca",
  );
  assert.ok(
    calls.every((call) => call.numPasses === 2),
    `a call of the pca says ${JSON.stringify(calls)}`,
  );
  // The four calls come in this order, and not the two of the first pass
  // and then the two of the second: the analysis builds both of its readers
  // before it asks either for a block, and the header of a VCF is read when
  // its reader is built, which is the first read of that pass and the call
  // that says it has read nothing. So the two passes take their numbers, 1
  // and then 2, at the two calls of 0 bytes, and each of them reads the file
  // after that. What a page draws from this rises all the same: the share of
  // the run that is done, `(pass - 1 + bytesRead / numBytes) / numPasses`,
  // goes 0, 0.5, 0.5 and 1.
  assert.deepEqual(
    calls.map((call) => call.pass),
    [1, 2, 1, 2],
  );
  assert.deepEqual(
    calls.map((call) => call.bytesRead),
    [0, 0, BYTES_OF_MANY_VCF, BYTES_OF_MANY_VCF],
  );
});

test("the two passes of the study asked for the approximation take their numbers in the order they start", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    calcGwas(variants, {
      phenotype: THE_TRAIT,
      trait: "continuous",
      kinship: THE_KINSHIP,
      useGrammarGammaApprox: true,
      transformToBiallelic: true,
    });
  } finally {
    variants.free();
  }
  assertTheyRise(
    calls.filter((call) => call.pass === 1),
    "the first pass of the study",
  );
  assertTheyRise(
    calls.filter((call) => call.pass === 2),
    "the second pass of the study",
  );
  assert.ok(
    calls.every((call) => call.numPasses === 2),
    `a call of the study says ${JSON.stringify(calls)}`,
  );
  // Both readers are built before either is asked for a block, and the
  // header of a VCF is read when its reader is built, so the two passes take
  // their numbers, 1 and then 2, at the two calls of 0 bytes. The second
  // pass gives the first block the factor of the approximation is estimated
  // from and is then dropped, so it ends where the end of the run tells the
  // page of it, below the whole file, while the first pass reads the file to
  // its end.
  assert.deepEqual(
    calls.map((call) => call.pass),
    [1, 2, 1, 2],
  );
  assert.deepEqual(
    calls.map((call) => call.bytesRead).slice(0, 2),
    [0, 0],
  );
  assert.equal(
    calls.filter((call) => call.pass === 1).at(-1)?.bytesRead,
    BYTES_OF_MANY_VCF,
  );
});

test("the study that is not asked for the approximation reads the source once", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    calcGwas(variants, {
      phenotype: THE_TRAIT,
      trait: "continuous",
      kinship: THE_KINSHIP,
      transformToBiallelic: true,
    });
  } finally {
    variants.free();
  }
  assert.deepEqual(
    calls.map((call) => ({ pass: call.pass, numPasses: call.numPasses })),
    [
      { pass: 1, numPasses: 1 },
      { pass: 1, numPasses: 1 },
    ],
  );
});

test("the pca that is asked for no weight reads the source once", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    doPcaFromVariants(variants, {
      numPrinComps: 0,
      transformToBiallelic: true,
    });
  } finally {
    variants.free();
  }
  assert.deepEqual(
    calls.map((call) => ({ pass: call.pass, numPasses: call.numPasses })),
    [
      { pass: 1, numPasses: 1 },
      { pass: 1, numPasses: 1 },
    ],
  );
});

test("twelve iterations of blocks over one source are twelve runs of one pass", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    // The twelve are opened together and read one after another: they are
    // twelve runs over one source, and each of them is the pass 1 of 1 of
    // its own run, where twelve passes counted in the source would be the
    // passes 1 to 12.
    const passes = Array.from({ length: 12 }, () =>
      variants.iterBlocks({ numVarsPerBlock: 100 }),
    );
    for (const pass of passes) {
      assert.equal([...pass].length, 5);
    }
  } finally {
    variants.free();
  }
  assert.ok(
    calls.every((call) => call.pass === 1 && call.numPasses === 1),
    `a call of the twelve iterations says ${JSON.stringify(calls)}`,
  );
  assert.equal(calls.length, 24);
});

/**
 * The twelve consumers of the package, each with the options its run is made
 * with and the ones `numPassesOf` is asked with, and the association study
 * twice, once for each number of passes its GRAMMAR-Gamma approximation
 * gives.
 *
 * `transformToBiallelic` is true where the calculation asks for it, because
 * `many.vcf` holds variants of more than two alleles and the core's default
 * refuses them.
 */
const THE_CONSUMERS: readonly {
  name: ConsumerName;
  options?: object;
  run: (variants: Variants) => void;
}[] = [
  {
    name: "calcPerVarDistribs",
    run: (variants) => {
      calcPerVarDistribs(variants);
    },
  },
  {
    name: "calcPerIndividualStats",
    run: (variants) => {
      calcPerIndividualStats(variants);
    },
  },
  {
    name: "calcPairwiseKosmanDists",
    run: (variants) => {
      calcPairwiseKosmanDists(variants);
    },
  },
  {
    name: "calcPopDists",
    run: (variants) => {
      calcPopDists(variants, THE_POPS, {
        measures: ["fst"],
        jackknifeGroup: null,
        minNumIndividuals: 1,
      });
    },
  },
  {
    name: "calcPopDiversity",
    run: (variants) => {
      calcPopDiversity(variants, { pops: THE_POPS, minNumIndividuals: 1 });
    },
  },
  {
    name: "calcRogersHuffR2Matrix",
    run: (variants) => {
      calcRogersHuffR2Matrix(variants);
    },
  },
  {
    name: "calcLdAndDistPerPop",
    run: (variants) => {
      calcLdAndDistPerPop(variants, { pops: THE_POPS });
    },
  },
  {
    name: "calcKinship",
    run: (variants) => {
      calcKinship(variants, { transformToBiallelic: true });
    },
  },
  {
    name: "doPcaFromVariants",
    options: { numPrinComps: 10 },
    run: (variants) => {
      doPcaFromVariants(variants, {
        numPrinComps: 10,
        transformToBiallelic: true,
      });
    },
  },
  {
    name: "calcGwas",
    run: (variants) => {
      calcGwas(variants, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        transformToBiallelic: true,
      });
    },
  },
  {
    name: "calcGwas",
    options: { useGrammarGammaApprox: true },
    run: (variants) => {
      calcGwas(variants, {
        phenotype: THE_TRAIT,
        trait: "continuous",
        kinship: THE_KINSHIP,
        useGrammarGammaApprox: true,
        transformToBiallelic: true,
      });
    },
  },
  {
    name: "writeVars",
    run: (variants) => {
      writeVars(variants);
    },
  },
  {
    name: "iterBlocks",
    run: (variants) => {
      for (const _block of variants.iterBlocks()) {
        // Every block is read, so the pass reads the file to its end.
      }
    },
  },
];

test("every consumer makes the passes numPassesOf says it makes", () => {
  for (const consumer of THE_CONSUMERS) {
    const variants = openVcf(MANY_VCF);
    const calls = theCallsOf(variants);
    try {
      consumer.run(variants);
    } finally {
      variants.free();
    }
    const numPasses = numPassesOf(consumer.name, consumer.options);
    assert.ok(calls.length > 0, `${consumer.name} told the page nothing`);
    assert.equal(
      Math.max(...calls.map((call) => call.pass)),
      numPasses,
      `the passes of ${consumer.name}`,
    );
    assert.ok(
      calls.every((call) => call.numPasses === numPasses),
      `a call of ${consumer.name} says another number of passes`,
    );
  }
});

test("onProgress with nothing takes the function off", () => {
  const variants = openVcf(MANY_VCF);
  const calls = theCallsOf(variants);
  try {
    calcPerIndividualStats(variants);
    const toldOfTheFirstRun = calls.length;
    assert.ok(toldOfTheFirstRun > 0, "the first run told the page nothing");
    variants.onProgress();
    calcPerIndividualStats(variants);
    assert.equal(calls.length, toldOfTheFirstRun);
  } finally {
    variants.free();
  }
});

test("a function that is not a function is refused where it is set", () => {
  const variants = openVcf(MANY_VCF);
  try {
    assert.throws(
      () => variants.onProgress(42 as unknown as (progress: Progress) => void),
      {
        name: "Error",
        message: /`told` is the function/,
      },
    );
  } finally {
    variants.free();
  }
});
