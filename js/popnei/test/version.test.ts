/**
 * That the version a TypeScript user sees is the one of the core crate.
 *
 * It is the first test of the chain the TypeScript side is made of: the
 * core crate, the binding crate `crates/popnei-js` compiled for
 * `wasm32-unknown-unknown`, the JavaScript that wasm-bindgen generates
 * from it, and this package. When it passes, the four were built together
 * and the package loaded the WebAssembly that was just compiled.
 *
 * It imports `popnei`, the name of this package, which node resolves to
 * the built `dist/node.js` through the field `exports` of `package.json`,
 * so the entry point a user of node gets is the one that is tested.
 */

import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

import { init, version } from "popnei";

import { versionOfTheCore, versionOfTheCoreCrate } from "./manifest.ts";

/** Read the version this package would be published to npm under. */
async function versionOfThePackage(): Promise<string> {
  const manifest = await readFile(
    new URL("../package.json", import.meta.url),
    "utf8",
  );
  return (JSON.parse(manifest) as { version: string }).version;
}

test("the version the package gives is the version of the core crate", async () => {
  await init();
  assert.equal(version(), await versionOfTheCoreCrate());
});

test("the version of the package is the version of the core crate", async () => {
  assert.equal(await versionOfThePackage(), await versionOfTheCoreCrate());
});

test("the version is read from a manifest that writes it with no spaces", () => {
  const manifest = [
    "[workspace]",
    'members = ["crates/popnei"]',
    "",
    "[workspace.dependencies]",
    'popnei = { path = "crates/popnei", version = "0.0.1" }',
    "",
    "[workspace.package]",
    'version="9.9.9"',
    'edition = "2024"',
    "",
  ].join("\n");
  assert.equal(versionOfTheCore(manifest), "9.9.9");
});

test("init loads the WebAssembly once and gives back the same promise", () => {
  const loading = init();
  assert.equal(init(), loading);
});
