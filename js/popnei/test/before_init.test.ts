/**
 * That a function of the package says what is wrong when it is called
 * before the WebAssembly of the core is loaded.
 *
 * node's test runner gives each test file its own process, and this is the
 * file in which `init` is never called: beside a test that awaits it, this
 * one would pass or fail with the order the tests run in.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { openVcf, version } from "popnei";

import { vcfOf } from "./reference.ts";

test("version throws before init was awaited", () => {
  assert.throws(() => version(), {
    name: "Error",
    message: /await init\(\)/,
  });
});

test("openVcf throws before init was awaited", () => {
  assert.throws(() => openVcf(vcfOf([])), {
    name: "Error",
    message: /await init\(\)/,
  });
});
