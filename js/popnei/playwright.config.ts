/**
 * How the browser tests of popnei are run: Chromium, headless, over the
 * files of the repository served by `test/browser/server.ts`.
 *
 * `npm run test:browser` runs them. A test of `test/browser/` opens
 * `harness.html`, which starts a module web worker, and the worker loads
 * the WebAssembly of `wasm/` through the entry point of `dist/web.js` and
 * asserts there the numbers the tests under node assert. The browsers
 * Playwright drives are downloaded once, with `npx playwright install
 * chromium`, and Playwright says so when they are missing.
 *
 * Firefox and WebKit are not run: the owner decided on 24 September 2026
 * to start with Chromium and to add them when an application needs them,
 * which is another `projects` entry and their download.
 */

import { defineConfig, devices } from "@playwright/test";

/**
 * The port the server of the tests listens on. It is of the range no
 * program is registered at, and the tests fail with the answer of whatever
 * else holds it.
 */
const PORT = 8973;

export default defineConfig({
  testDir: "./test/browser",
  // The tests under node are `test/*.test.ts` and are run by `npm test`,
  // so a browser test is named for the browser and neither runner takes
  // the files of the other.
  testMatch: "**/*.browser.ts",
  reporter: "list",
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"] } }],
  webServer: {
    command: `node test/browser/server.ts ${PORT}`,
    url: `http://127.0.0.1:${PORT}/js/popnei/test/browser/harness.html`,
    reuseExistingServer: false,
  },
});
