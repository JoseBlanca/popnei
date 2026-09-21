// Loads pyodide under node, installs into it the wheel that
// scripts/build_pyodide_wheel.sh left in dist/, and checks four things:
// that the version popnei answers with is the one of the core crate, which
// is in [workspace.package] of the Cargo.toml of the repository; that
// `open_vcf` reads tests/reference/vcf/cases.vcf and cases.vcf.gz there as
// the table of "How it is verified" of docs/specs/io_vcf.md says; that
// `write_vars` writes those variants into a vars file and `open_vars`
// reads the four of them back out of it, which is what says that arrow-rs
// was linked into this wheel and that it writes and decompresses there; and
// that a VCF whose blocks of the size popnei chooses would not fit in what
// a wasm build counts is opened all the same. It exits with an error when
// anything differs.
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

// The variants that popnei reads from one source inside pyodide, in the
// shape of the two tables above: for each variant its position and then its
// genotypes, individual after individual. A block holds several variants,
// and the blocks of a source, one after another, are all of them. The
// snippets below run in this same interpreter and call `rows_of` for the
// variants of a vars file.
const READ_THE_VARIANTS = `
import json

import popnei


def rows_of(variants):
    rows = []
    for block in variants.iter_blocks():
        for index in range(block.num_vars):
            alleles = [int(allele) for allele in block.gts[index].ravel()]
            rows.append([int(block.pos[index])] + alleles)
    return json.dumps(rows)


def variants_as_rows(vcf_path, only_passed):
    return rows_of(popnei.open_vcf(vcf_path, only_passed=only_passed))
`;

// The six bytes that an arrow IPC file begins and ends with. A vars file is
// one, written by arrow-rs, so finding them says that the crates of arrow
// were linked into this wheel and ran inside pyodide.
const ARROW_MARK = "ARROW1";

// The variants of a VCF written into a vars file inside pyodide, at a path
// of the file system of emscripten, which `write_vars` opens as it does a
// path natively, and read back from it with `open_vars`. The file is
// written with the buffers of its batches compressed with lz4, which popnei
// writes in every build, so reading it back is that compression undone in
// wasm.
const A_VARS_FILE_WRITTEN_AND_READ = `
import popnei


def write_a_vars_file(vcf_path, vars_path):
    variants = popnei.open_vcf(vcf_path, only_passed=False)
    popnei.write_vars(variants, vars_path)


def vars_file_as_rows(vars_path):
    return rows_of(popnei.open_vars(vars_path))
`;

// How many individuals the header of the fourth check names, and the ploidy
// it is read with, the largest a reader of popnei takes. A block of the
// size popnei chooses for so many individuals is 100 variants, and its
// genotypes are 100 x 170000 x 255, 4335 million, where a count of things
// in wasm holds 4295 million.
const MANY_INDIVIDUALS = 170000;
const LARGEST_PLOIDY = 255;

// What popnei answers for such a file inside pyodide: opening it reads the
// header and asks for no block, so the individuals and the ploidy come out;
// blocks of 10 variants fit and the file has no variant to put in one; and
// a size that does not fit, 100 and the one popnei chooses, is refused when
// the blocks are asked for and not before. js/popnei/test/open.test.ts has
// the case under node, where a count of things is 64 bits and nothing is
// refused for its size.
const OPEN_A_HEADER_OF_MANY_INDIVIDUALS = `
import json

import popnei


def header_of_many_individuals(vcf_path, num_individuals):
    names = "\\t".join(f"ind{individual}" for individual in range(num_individuals))
    columns = "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t" + names
    with open(vcf_path, "w") as vcf:
        vcf.write("##fileformat=VCFv4.4\\n" + columns + "\\n")


def what_a_header_of_many_individuals_gives(vcf_path, num_individuals, ploidy):
    header_of_many_individuals(vcf_path, num_individuals)
    variants = popnei.open_vcf(vcf_path, ploidy=ploidy)
    what = {
        "num_individuals": variants.num_individuals,
        "ploidy": variants.ploidy,
        "first_individual": variants.individuals[0],
    }
    what["blocks_of_ten"] = [
        block.num_vars for block in variants.iter_blocks(num_vars_per_block=10)
    ]
    for asked_for, size in (("of_a_hundred", 100), ("of_the_default", None)):
        try:
            for block in variants.iter_blocks(num_vars_per_block=size):
                pass
            what[asked_for] = "no error"
        except ValueError as error:
            what[asked_for] = str(error)
    return json.dumps(what)
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

pyodide.runPython(A_VARS_FILE_WRITTEN_AND_READ);
const varsPath = "/vcf/cases.vars";
pyodide.runPython(`write_a_vars_file("/vcf/cases.vcf", "${varsPath}")`);
const varsFile = pyodide.FS.readFile(varsPath);
const decoder = new TextDecoder();
const startsWith = decoder.decode(varsFile.slice(0, ARROW_MARK.length));
const endsWith = decoder.decode(varsFile.slice(-ARROW_MARK.length));
if (startsWith !== ARROW_MARK || endsWith !== ARROW_MARK) {
  failures.push(
    `the vars file written in pyodide starts with "${startsWith}" and ends` +
      ` with "${endsWith}", and an arrow IPC file has "${ARROW_MARK}" at` +
      " both ends",
  );
} else {
  console.log(
    `cases.vcf written as a vars file of ${varsFile.length} bytes, with` +
      ` "${ARROW_MARK}" at both ends`,
  );
}

const fromTheVarsFile = JSON.parse(
  pyodide.runPython(`vars_file_as_rows("${varsPath}")`),
);
if (JSON.stringify(fromTheVarsFile) !== JSON.stringify(EVERY_VARIANT)) {
  failures.push(
    `the vars file read back gives ${JSON.stringify(fromTheVarsFile)} and` +
      ` the spec says ${JSON.stringify(EVERY_VARIANT)}`,
  );
} else {
  console.log(
    `cases.vars read back: ${EVERY_VARIANT.length} variants, as the spec says`,
  );
}

pyodide.runPython(OPEN_A_HEADER_OF_MANY_INDIVIDUALS);
const openCall =
  `what_a_header_of_many_individuals_gives("/vcf/many_individuals.vcf",` +
  ` ${MANY_INDIVIDUALS}, ${LARGEST_PLOIDY})`;
const opened = JSON.parse(pyodide.runPython(openCall));
const whatWasAsked =
  `a header of ${MANY_INDIVIDUALS} individuals of the ploidy ${LARGEST_PLOIDY}`;
if (
  opened.num_individuals !== MANY_INDIVIDUALS ||
  opened.ploidy !== LARGEST_PLOIDY ||
  opened.first_individual !== "ind0"
) {
  failures.push(
    `${whatWasAsked} was opened as ${JSON.stringify(opened)}: opening a file` +
      " reads its header and asks for no block",
  );
} else if (JSON.stringify(opened.blocks_of_ten) !== "[]") {
  failures.push(
    `${whatWasAsked} gives the blocks ${JSON.stringify(opened.blocks_of_ten)}` +
      " for a size of 10, and the file has no variant",
  );
} else if (
  !opened.of_a_hundred.includes("memory") ||
  !opened.of_the_default.includes("memory")
) {
  failures.push(
    `${whatWasAsked} answers "${opened.of_a_hundred}" for blocks of 100 and` +
      ` "${opened.of_the_default}" for the size popnei chooses, and the` +
      " genotypes of both are more than this build counts",
  );
} else {
  console.log(
    `${whatWasAsked}: opened, no block of 10, and 100 needs more memory than` +
      " this build gives",
  );
}

for (const failure of failures) {
  console.error(failure);
}
if (failures.length > 0) {
  process.exit(1);
}
