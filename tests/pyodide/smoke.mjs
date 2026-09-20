// Loads pyodide under node, installs into it the wheel that
// scripts/build_pyodide_wheel.sh left in dist/, and checks two things: that
// the version popnei answers with is the one of the core crate, which is in
// [workspace.package] of the Cargo.toml of the repository, and that
// `open_vcf` reads tests/reference/vcf/cases.vcf and cases.vcf.gz there as
// the table of "How it is verified" of docs/specs/io_vcf.md says. It exits
// with an error when anything differs.
//
// README.md, beside this file, says how to run it.

import { readdir, readFile } from "node:fs/promises";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

import { loadPyodide } from "pyodide";

const repoRoot = join(dirname(fileURLToPath(import.meta.url)), "..", "..");

// The four variants of cases.vcf, of three diploid individuals, as the
// table of "How it is verified" of docs/specs/io_vcf.md gives them: the
// position, and then the two alleles of each individual, 0 for the
// reference allele, 1 and 2 for the alternative ones and -1 for an allele
// that was not called.
const EVERY_VARIANT = [
  [100, 0, 0, 0, 1, 1, 1],
  [200, -1, -1, 0, 1, -1, 0],
  [300, 1, 2, 2, 1, 2, 2],
  [400, 0, 0, 0, 0, 0, 0],
];
// The FILTER of the variant at 200 is q10, and of the other three PASS, so
// the default, which gives only the variants that failed no filter, leaves
// that one out.
const PASSED_VARIANTS = [EVERY_VARIANT[0], EVERY_VARIANT[2], EVERY_VARIANT[3]];

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

// The variants that popnei reads from one VCF inside pyodide, in the shape
// of the two tables above: for each variant its position and then its
// genotypes, individual after individual. A block holds several variants,
// and the blocks of a source, one after another, are all of them.
const READ_THE_VARIANTS = `
import json

import popnei


def variants_as_rows(vcf_path, only_passed):
    variants = popnei.open_vcf(vcf_path, only_passed=only_passed)
    rows = []
    for block in variants.iter_blocks():
        for index in range(block.num_vars):
            alleles = [int(allele) for allele in block.gts[index].ravel()]
            rows.append([int(block.pos[index])] + alleles)
    return json.dumps(rows)
`;

const failures = [];

const expectedVersion = versionOfTheCore(
  await readFile(join(repoRoot, "Cargo.toml"), "utf8"),
);
const wheel = await theWheel();

const pyodide = await loadPyodide();
console.log(`pyodide ${pyodide.version}, installing ${wheel.name}`);

// micropip reads a wheel of the file system of emscripten, not of the one
// of node, so the bytes are written there first and the path is given with
// the `emfs:` prefix that asks micropip for a local file.
const wheelInPyodide = `/tmp/${wheel.name}`;
pyodide.FS.writeFile(wheelInPyodide, await readFile(wheel.path));
// The genotypes of a block are a numpy array, so popnei cannot be imported
// before numpy is in pyodide. numpy is a package of pyodide itself and is
// loaded from there; micropip would fetch it for the dependency of the
// wheel anyway, and asking for it here says which numpy answers.
await pyodide.loadPackage(["micropip", "numpy"]);
const micropip = pyodide.pyimport("micropip");
await micropip.install(`emfs:${wheelInPyodide}`);

const foundVersion = pyodide.runPython("import popnei\npopnei.__version__");
const numpyVersion = pyodide.runPython("import numpy\nnumpy.__version__");
console.log(
  `popnei.__version__ is ${foundVersion}, on numpy ${numpyVersion}`,
);
if (foundVersion !== expectedVersion) {
  failures.push(
    `popnei.__version__ is ${foundVersion} and the core crate is` +
      ` ${expectedVersion}: the wheel in dist/ is of another build`,
  );
}

// The reader opens a path of the file system of emscripten, so the two
// files are copied into it, where nothing else of node reaches.
const reference = join(repoRoot, "tests", "reference", "vcf");
pyodide.FS.mkdir("/vcf");
for (const name of ["cases.vcf", "cases.vcf.gz"]) {
  pyodide.FS.writeFile(`/vcf/${name}`, await readFile(join(reference, name)));
}

pyodide.runPython(READ_THE_VARIANTS);

for (const name of ["cases.vcf", "cases.vcf.gz"]) {
  for (const [onlyPassed, expected] of [
    [true, PASSED_VARIANTS],
    [false, EVERY_VARIANT],
  ]) {
    const asked = onlyPassed ? "True" : "False";
    const call = `variants_as_rows("/vcf/${name}", ${asked})`;
    const found = JSON.parse(pyodide.runPython(call));
    const what = `${name} with only_passed=${onlyPassed}`;
    if (JSON.stringify(found) !== JSON.stringify(expected)) {
      failures.push(
        `${what} gives ${JSON.stringify(found)} and the spec says` +
          ` ${JSON.stringify(expected)}`,
      );
    } else {
      console.log(`${what}: ${expected.length} variants, as the spec says`);
    }
  }
}

for (const failure of failures) {
  console.error(failure);
}
if (failures.length > 0) {
  process.exit(1);
}
