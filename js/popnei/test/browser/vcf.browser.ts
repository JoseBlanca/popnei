/**
 * That popnei runs in a browser: the WebAssembly loads inside a module web
 * worker of Chromium and reads there a VCF held as an array of bytes.
 *
 * Playwright starts the browser and the server of `server.ts`, opens
 * `harness.html` and asks it for the case of `cases/`, which asserts the
 * numbers inside the worker and whose message is what a failure here
 * shows.
 */

import { test } from "@playwright/test";

test("a web worker of Chromium reads many.vcf from an array of bytes and gives its first block", async ({
  page,
}) => {
  await page.goto("/js/popnei/test/browser/harness.html");
  await page.evaluate(() => window.runInTheWorker("vcf_bytes"));
});
