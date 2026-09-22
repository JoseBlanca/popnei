// How long the principal component analysis of the variants of a vars file
// takes under pyodide, on the wheel that scripts/build_pyodide_wheel.sh
// leaves in dist/.
//
// It times one thing: `popnei.do_pca_from_variants` over a vars file inside
// the file system of emscripten, which is the analysis and the reading of
// the file, and it is the Python API of popnei, the same call a user makes
// in a notebook in a page. The bytes of the file are written into that file
// system before the clock starts, because that write is node's and not
// popnei's.
//
// It is what task 4.2 of docs/plans/pca.md timed the wheel with, beside the
// same analysis on the wasm package, which `js/popnei/bench/time_pca.mjs`
// times; docs/reports/pca-measurement.md has the numbers with the commands
// that produced them. The two builds are of different targets,
// `wasm32-unknown-emscripten` here and `wasm32-unknown-unknown` there.
//
//     node tests/pyodide/time_pca.mjs <path to a vars file> [runs]
//
// from the root of the repository, after scripts/build_pyodide_wheel.sh.
// `runs` is 3 when it is not given, and one run comes before the timed ones
// whose time is not taken, for the memory that the first analysis grows and
// the code that is compiled in it.
//
// README.md, beside this file, says what has to be installed to run
// anything of this directory.

import { readFile } from "node:fs/promises";
import { readdir } from "node:fs/promises";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

import { loadPyodide } from "pyodide";

const repoRoot = join(dirname(fileURLToPath(import.meta.url)), "..", "..");

// How many times the analysis is timed when the command line does not say.
const DEFAULT_RUNS = 3;

/** The one wheel of popnei for pyodide in dist/, as its path and its name. */
async function theWheel() {
  const dist = join(repoRoot, "dist");
  const names = await readdir(dist);
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

// The analysis, as a user of popnei writes it in Python. It gives the time
// of the call and what it gave, so that the clock is inside pyodide and the
// crossing back into node is not timed.
const THE_ANALYSIS = `
import time

import popnei


def one_analysis(path, num_prin_comps):
    started = time.perf_counter()
    variants = popnei.open_vars(path)
    result = popnei.do_pca_from_variants(variants, num_prin_comps=num_prin_comps)
    took = time.perf_counter() - started
    first = result.explained_variance_percent.iloc[0]
    return (
        f"{took:.3f}|{result.projections.shape[0]} individuals, "
        f"{result.princomps.shape[1]} variants with variance, "
        f"{result.projections.shape[1]} components, the first with "
        f"{first:.3f} per 100 of the variance"
    )
`;

const [path, runs = DEFAULT_RUNS] = process.argv.slice(2);
if (path === undefined) {
  console.error("node tests/pyodide/time_pca.mjs <path to a vars file> [runs]");
  process.exit(1);
}

const wheel = await theWheel();
const pyodide = await loadPyodide();
console.log(`pyodide ${pyodide.version}, installing ${wheel.name}`);
const wheelInPyodide = `/tmp/${wheel.name}`;
pyodide.FS.writeFile(wheelInPyodide, await readFile(wheel.path));
await pyodide.loadPackage(["micropip", "numpy"]);
const micropip = pyodide.pyimport("micropip");
await micropip.install(`emfs:${wheelInPyodide}`);

const varsInPyodide = "/tmp/timed.vars";
const bytes = await readFile(path);
pyodide.FS.writeFile(varsInPyodide, bytes);
console.log(
  `${path}, ${bytes.length} bytes, node ${process.version}, ${runs} runs`,
);

pyodide.runPython(THE_ANALYSIS);
const oneAnalysis = pyodide.globals.get("one_analysis");
const times = [];
for (let run = 0; run <= Number(runs); run += 1) {
  const [took, did] = oneAnalysis(varsInPyodide, 0).split("|");
  if (run === 0) {
    console.log(`the first run, which is not timed: ${did}, in ${took} s`);
    continue;
  }
  times.push(Number(took));
  console.log(`run ${run}: ${took} s, ${did}`);
}
const sorted = [...times].sort((one, other) => one - other);
console.log(
  `best ${sorted[0].toFixed(3)} s, worst ${sorted[sorted.length - 1].toFixed(3)} s`,
);
