// Loads pyodide under node, installs into it the wheel that
// scripts/build_pyodide_wheel.sh left in dist/, and prints the version of
// popnei that answers inside pyodide. It exits with an error when that
// version is not the one of the core crate, which is in
// [workspace.package] of the Cargo.toml of the repository.
//
// README.md, beside this file, says how to run it.

import { readdir, readFile } from "node:fs/promises";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

import { loadPyodide } from "pyodide";

const repoRoot = join(dirname(fileURLToPath(import.meta.url)), "..", "..");

/**
 * The version of the core crate, which every build of popnei publishes as
 * its own: the `version` of the `[workspace.package]` table of the given
 * cargo manifest.
 */
function versionOfTheCore(manifest) {
  let inWorkspacePackage = false;
  for (const line of manifest.split("\n")) {
    const text = line.trim();
    if (text.startsWith("[")) {
      inWorkspacePackage = text === "[workspace.package]";
    } else if (inWorkspacePackage) {
      const version = /^version\s*=\s*"([^"]+)"/.exec(text);
      if (version !== null) {
        return version[1];
      }
    }
  }
  throw new Error("no version in the [workspace.package] of Cargo.toml");
}

/** The one wheel of popnei for pyodide in dist/, as its path and its name. */
async function theWheel() {
  const dist = join(repoRoot, "dist");
  let names;
  try {
    names = await readdir(dist);
  } catch {
    throw new Error(`no ${dist}: run scripts/build_pyodide_wheel.sh first`);
  }
  const wheels = names.filter(
    (name) =>
      name.startsWith("popnei-") &&
      name.includes("pyemscripten") &&
      name.endsWith(".whl"),
  );
  if (wheels.length !== 1) {
    throw new Error(
      `${wheels.length} wheels of popnei with a pyemscripten tag in ${dist},` +
        " and there has to be one: run scripts/build_pyodide_wheel.sh",
    );
  }
  return { path: join(dist, wheels[0]), name: wheels[0] };
}

const expected = versionOfTheCore(
  await readFile(join(repoRoot, "Cargo.toml"), "utf8"),
);
const wheel = await theWheel();

const pyodide = await loadPyodide();
console.log(`pyodide ${pyodide.version}, installing ${wheel.name}`);

// micropip reads a wheel of the file system of emscripten, not of the one
// of node, so the bytes are written there first and the path is given with
// the `emfs:` prefix that asks micropip for a local file.
const inPyodide = `/tmp/${wheel.name}`;
pyodide.FS.writeFile(inPyodide, await readFile(wheel.path));
await pyodide.loadPackage("micropip");
const micropip = pyodide.pyimport("micropip");
await micropip.install(`emfs:${inPyodide}`);

const found = pyodide.runPython("import popnei\npopnei.__version__");
console.log(`popnei.__version__ is ${found}`);
if (found !== expected) {
  console.error(
    `the core crate is ${expected}: the wheel in dist/ is of another build`,
  );
  process.exit(1);
}
