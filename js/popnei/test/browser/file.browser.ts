/**
 * That popnei reads a file the user picked in the page one range of bytes
 * at a time: a VCF, the same VCF gzipped and a vars file give the variants
 * the tests under node give, a file of more than one range gives the
 * variants of every range of it, a range that comes back short or is
 * refused ends the pass with popnei's error, and a file of 299994147 bytes
 * is passed over with the memory of wasm staying under 24 MiB.
 *
 * A range of a `File` is read with `FileReaderSync`, which a browser has
 * only inside a web worker, so none of these seven can run under node.
 * Playwright starts Chromium and the server of `server.ts`, opens
 * `harness.html` and asks it for a case of `cases/`, which asserts the
 * numbers inside the worker and whose message is what a failure here shows.
 */

import { test } from "@playwright/test";

/** The page that starts the worker every case of this file runs in. */
const HARNESS = "/js/popnei/test/browser/harness.html";

test("a File of many.vcf gives the variants the file has", async ({ page }) => {
  await page.goto(HARNESS);
  await page.evaluate((name) => window.runInTheWorker(name), "vcf_file");
});

test("a File of many.vcf.gz gives those same variants", async ({ page }) => {
  await page.goto(HARNESS);
  await page.evaluate((name) => window.runInTheWorker(name), "vcf_gz_file");
});

test("a File of a vars file written from many.vcf gives those same variants", async ({
  page,
}) => {
  await page.goto(HARNESS);
  await page.evaluate((name) => window.runInTheWorker(name), "vars_file");
});

test("a File of more than one range gives the variants of every range", async ({
  page,
}) => {
  await page.goto(HARNESS);
  await page.evaluate((name) => window.runInTheWorker(name), "many_ranges");
});

test("a range that comes back one byte short ends the pass and names the range", async ({
  page,
}) => {
  await page.goto(HARNESS);
  await page.evaluate((name) => window.runInTheWorker(name), "short_range");
});

test("a range the browser refused ends the pass with what the browser said", async ({
  page,
}) => {
  await page.goto(HARNESS);
  await page.evaluate(
    (name) => window.runInTheWorker(name),
    "browser_refused",
  );
});

test("a File of 299994147 bytes is passed over without the memory of wasm passing 24 MiB", async ({
  page,
}) => {
  await page.goto(HARNESS);
  await page.evaluate(
    (name) => window.runInTheWorker(name),
    "wasm_memory_over_300_mb",
  );
});
