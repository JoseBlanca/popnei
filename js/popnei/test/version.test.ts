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

const REPO_DIR = new URL("../../../", import.meta.url);

/**
 * Read the version of `crates/popnei` from the workspace manifest.
 *
 * The crate takes `version.workspace = true`, so the version is written
 * once, in `[workspace.package]` of the `Cargo.toml` of the repository.
 * node has no reader of TOML, the format of that file, so the section is
 * taken as the lines between its heading and the next one, and the line
 * that gives the version is matched in it.
 */
async function versionOfTheCoreCrate(): Promise<string> {
  const manifest = await readFile(new URL("Cargo.toml", REPO_DIR), "utf8");
  const afterTheHeading = manifest.split("\n[workspace.package]\n")[1];
  assert.ok(
    afterTheHeading !== undefined,
    "the workspace manifest has no [workspace.package] section",
  );
  const section = afterTheHeading.split("\n[")[0] ?? "";
  const version = /^version = "([^"]+)"$/m.exec(section)?.[1];
  assert.ok(
    version !== undefined,
    "the [workspace.package] section of the workspace manifest gives no version",
  );
  return version;
}

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
