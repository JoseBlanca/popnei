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

import type { Variants } from "popnei";
import { openVars, openVcf, version, writeVars } from "popnei";

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

test("openVars throws before init was awaited", () => {
  assert.throws(() => openVars(new Uint8Array([65, 82, 82, 79, 87, 49])), {
    name: "Error",
    message: /await init\(\)/,
  });
});

test("writeVars throws before init was awaited", () => {
  // The handle cannot be made before `init` either, so what is given here
  // is the argument the check of the arguments would refuse: the error of
  // the WebAssembly comes first, which is what this asserts.
  assert.throws(() => writeVars(null as unknown as Variants), {
    name: "Error",
    message: /await init\(\)/,
  });
});
