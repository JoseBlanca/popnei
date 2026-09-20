/**
 * The version of the core crate, read from the cargo manifest of the
 * repository, for the tests that check that every build of popnei
 * publishes that version as its own.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`.
 */

import { readFile } from "node:fs/promises";

const REPO_DIR = new URL("../../../", import.meta.url);

/**
 * The version of the core crate: the `version` of the `[workspace.package]`
 * table of the given cargo manifest.
 *
 * The crate takes `version.workspace = true`, so the version is written
 * once, in that table. node has no reader of TOML, the format of the
 * manifest, so the tables are walked line by line, as
 * `tests/pyodide/smoke.mjs` walks them.
 */
export function versionOfTheCore(manifest: string): string {
  let inWorkspacePackage = false;
  for (const line of manifest.split("\n")) {
    const text = line.trim();
    if (text.startsWith("[")) {
      inWorkspacePackage = text === "[workspace.package]";
    } else if (inWorkspacePackage) {
      const version = /^version\s*=\s*"([^"]+)"/.exec(text)?.[1];
      if (version !== undefined) {
        return version;
      }
    }
  }
  throw new Error("no version in the [workspace.package] of Cargo.toml");
}

/** The version of the core crate, from the manifest of the repository. */
export async function versionOfTheCoreCrate(): Promise<string> {
  return versionOfTheCore(await readFile(new URL("Cargo.toml", REPO_DIR), "utf8"));
}
