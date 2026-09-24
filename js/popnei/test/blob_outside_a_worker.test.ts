/**
 * That `openVcf` and `openVars` of a file of the page, given where there is
 * no `FileReaderSync`, are refused with a message that names it and the web
 * worker.
 *
 * popnei reads a `File`, the handle a page gets when the user picks a file,
 * and a `Blob`, which a `File` is one of, a range of bytes at a time through
 * `FileReaderSync`, the reader that returns when it has the bytes. A browser
 * gives that reader only inside a web worker, so the same call on the main
 * thread of a page has nothing to read the file with, and
 * `docs/specs/js_sources.md` says under "The cases a reader of the rules
 * would not guess" that it is an `Error` at the call, where the header of the
 * VCF or the footer of the vars file is read.
 *
 * node is that case: node 26.8.2 has `Blob` and `File` and no
 * `FileReaderSync`, which is why this refusal is the one thing of the reading
 * of a file that the tests under node can assert. What a range of a real file
 * gives is asserted in Chromium, by the tests of `test/browser/`.
 *
 * The file is `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
 * individuals in 117346 bytes, and the vars file written from it. Both are
 * opened here as bytes as well, which is what says that the `Error` is about
 * where the call was made and not about what the file holds.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVars, openVcf, writeVars } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`. */
const MANY_VCF = await referenceVcf("many.vcf");

/** How many individuals `many.vcf` holds, which its header names. */
const INDIVIDUALS_OF_MANY_VCF = 50;

/** The bytes of the vars file of the variants of `many.vcf`. */
const MANY_VARS = varsFileOfManyVcf();

/** The vars file written from every variant of `many.vcf`. */
function varsFileOfManyVcf(): Uint8Array {
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  try {
    return writeVars(variants).bytes;
  } finally {
    variants.free();
  }
}

/**
 * A `Blob` of `bytes`, the piece of bytes of a page that a `File` is one of.
 *
 * The bytes are copied into an array whose buffer is its own: the buffer of a
 * `Uint8Array` that a test read a file into is a `SharedArrayBuffer` as far
 * as TypeScript knows, and a `Blob` is not built from one.
 */
function blobOf(bytes: Uint8Array): Blob {
  return new Blob([new Uint8Array(bytes)]);
}

/** A `File` of `bytes` under `name`, which is a `Blob` that has a name. */
function fileOf(bytes: Uint8Array, name: string): File {
  return new File([new Uint8Array(bytes)], name);
}

/**
 * That `error` is popnei's refusal of a file opened where no web worker
 * gives a `FileReaderSync`, which names both.
 *
 * A user who reads the message has to learn what popnei could not find and
 * where it is, so the two names are asserted one by one: a message that named
 * the reader alone would leave them with a name of the web API and nothing to
 * do about it.
 */
function namesTheReaderAndTheWorker(error: unknown): true {
  assert.ok(error instanceof Error, `${String(error)} is not an Error`);
  assert.match(error.message, /FileReaderSync/);
  assert.match(error.message, /web worker/);
  return true;
}

test("node has `Blob` and `File` and no `FileReaderSync`", () => {
  // The premise of every test below: this process is the case of the main
  // thread of a page, where a file of the page can be built and not read.
  assert.equal(typeof globalThis.Blob, "function");
  assert.equal(typeof globalThis.File, "function");
  assert.equal("FileReaderSync" in globalThis, false);
});

test("openVcf of a `Blob` names `FileReaderSync` and the web worker", () => {
  assert.throws(
    () => openVcf(blobOf(MANY_VCF)),
    namesTheReaderAndTheWorker,
  );
});

test("openVcf of a `File` names `FileReaderSync` and the web worker", () => {
  assert.throws(
    () => openVcf(fileOf(MANY_VCF, "many.vcf"), { onlyPassed: false }),
    namesTheReaderAndTheWorker,
  );
});

test("openVars of a `Blob` names `FileReaderSync` and the web worker", () => {
  assert.throws(
    () => openVars(blobOf(MANY_VARS)),
    namesTheReaderAndTheWorker,
  );
});

test("openVars of a `File` names `FileReaderSync` and the web worker", () => {
  assert.throws(
    () => openVars(fileOf(MANY_VARS, "many.vars")),
    namesTheReaderAndTheWorker,
  );
});

test("the bytes of those two files open on the main thread", () => {
  const vcf = openVcf(MANY_VCF);
  try {
    assert.equal(vcf.numIndividuals, INDIVIDUALS_OF_MANY_VCF);
  } finally {
    vcf.free();
  }
  const vars = openVars(MANY_VARS);
  try {
    assert.equal(vars.numIndividuals, INDIVIDUALS_OF_MANY_VCF);
  } finally {
    vars.free();
  }
});
