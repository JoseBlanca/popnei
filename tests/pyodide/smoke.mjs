// Loads pyodide under node, installs into it the wheel that
// scripts/build_pyodide_wheel.sh left in dist/, and checks six things:
// that the version popnei answers with is the one of the core crate, which
// is in [workspace.package] of the Cargo.toml of the repository; that
// `open_vcf` reads tests/reference/vcf/cases.vcf and cases.vcf.gz there as
// the table of "How it is verified" of docs/specs/io_vcf.md says; that
// `write_vars` writes those variants into a vars file and `open_vars`
// reads the four of them back out of it, which is what says that arrow-rs
// was linked into this wheel and that it writes and decompresses there;
// that a VCF whose blocks of the size popnei chooses would not fit in what
// a wasm build counts is opened all the same; and that
// `calc_per_var_distribs` and `calc_per_individual_stats` give the numbers
// of the worked example of docs/specs/stats.md, which is what says that the
// two calculations reach the same values where there is one thread and the
// build is not the native one. It exits with an error when anything
// differs.
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

// The worked example of "How it is verified" of the per variant
// distributions and of the per individual statistics of
// docs/specs/stats.md: six variants of five diploid individuals named `i1`
// to `i5`. popnei has no `Variants.from_gt_array`, so the genotypes are
// written as a VCF, as tests/test_stats.py writes the worked examples it
// runs natively. Every variant is given three alternative alleles, because
// the third of the six holds the alleles 2 and 3.
const THE_WORKED_EXAMPLE = `
import json

import popnei

INDIVIDUALS = ("i1", "i2", "i3", "i4", "i5")

THE_SIX_VARIANTS = (
    "0/0 0/1 0/0 0/0 0/.",
    "0/0 0/1 0/0 ./. 0/.",
    "0/1 2/3 0/1 2/3 ./.",
    "./. ./. ./. ./. ./.",
    "0/0 0/0 0/0 0/0 1/1",
    "0/. ./. ./. ./. ./.",
)


def write_the_worked_example(vcf_path):
    lines = [
        "##fileformat=VCFv4.4",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t"
        + "\\t".join(INDIVIDUALS),
    ]
    for position, genotypes in enumerate(THE_SIX_VARIANTS, start=1):
        columns = ["chr1", str(position), ".", "A", "C,G,T", ".", "PASS", ".", "GT"]
        lines.append("\\t".join(columns + genotypes.split(" ")))
    with open(vcf_path, "w") as vcf:
        vcf.write("\\n".join(lines) + "\\n")


def per_var_distribs_of_the_worked_example(vcf_path, pops_json):
    distribs = popnei.calc_per_var_distribs(
        popnei.open_vcf(vcf_path),
        pops=json.loads(pops_json),
        min_num_individuals=1,
        hist_kwargs={"num_bins": 4},
    )
    what = {
        "num_vars": distribs.pass_stats.num_vars,
        "pops": [str(pop) for pop in distribs.obs_het.mean.index],
        "hist_bin_edges": [float(edge) for edge in distribs.obs_het.hist_bin_edges],
    }
    for stat in ("obs_het", "maf", "exp_het", "unbiased_exp_het"):
        distrib = getattr(distribs, stat)
        what[stat] = {
            "mean": [float(mean) for mean in distrib.mean],
            "hist": [
                [int(count) for count in distrib.hist_counts[pop]]
                for pop in distrib.hist_counts.columns
            ],
        }
    poly = distribs.poly_vars_ratio
    what["poly_vars_ratio"] = {
        "num_poly": [int(count) for count in poly.num_poly],
        "num_variable": [int(count) for count in poly.num_variable],
        "tot_num_variants_with_data": [
            int(count) for count in poly.tot_num_variants_with_data
        ],
        "poly_ratio": [float(ratio) for ratio in poly.poly_ratio],
        "poly_ratio_over_variables": [
            float(ratio) for ratio in poly.poly_ratio_over_variables
        ],
    }
    return json.dumps(what)


def per_individual_stats_of_the_worked_example(vcf_path):
    stats = popnei.calc_per_individual_stats(popnei.open_vcf(vcf_path))
    return json.dumps(
        {
            "num_vars": stats.pass_stats.num_vars,
            "individuals": [str(name) for name in stats.missing_gt_rate.index],
            "missing_gt_rate": [float(rate) for rate in stats.missing_gt_rate],
            "obs_het_rate": [float(rate) for rate in stats.obs_het_rate],
        }
    )
`;

// How many variants the worked example has, which both passes give.
const THE_WORKED_EXAMPLE_NUM_VARS = 6;

// The populations of the worked example, pop1 of i1 and i2 and pop2 of i3,
// i4 and i5, and the one population of every individual, named `pop`, that
// a call with no `pops` gives.
const THE_TWO_POPS = { pop1: ["i1", "i2"], pop2: ["i3", "i4", "i5"] };
const EVERY_INDIVIDUAL = null;

// The edges of the four bins from 0 to 1 that the example asks for.
const HIST_BIN_EDGES = [0, 0.25, 0.5, 0.75, 1];

// The mean and the four histogram counts of each statistic over the two
// populations, pop1 first, from the second table of "How it is verified" of
// the per variant distributions, and the three counts and two ratios of the
// polymorphism ratio from the worked example of its own item.
const OVER_THE_TWO_POPS = {
  pops: ["pop1", "pop2"],
  obs_het: {
    mean: [0.5, 0.25],
    hist: [
      [1, 0, 2, 1],
      [3, 0, 0, 1],
    ],
  },
  maf: {
    mean: [0.6875, 0.729167],
    hist: [
      [0, 1, 0, 3],
      [0, 1, 1, 2],
    ],
  },
  exp_het: {
    mean: [0.375, 0.298611],
    hist: [
      [1, 2, 0, 1],
      [2, 1, 0, 1],
    ],
  },
  unbiased_exp_het: {
    mean: [0.5, 0.383333],
    hist: [
      [1, 0, 2, 1],
      [2, 0, 1, 1],
    ],
  },
  poly_vars_ratio: {
    num_poly: [3, 2],
    num_variable: [3, 2],
    tot_num_variants_with_data: [4, 4],
    poly_ratio: [0.75, 0.5],
    poly_ratio_over_variables: [1, 1],
  },
};

// The same numbers over the five individuals as one population, which the
// same two paragraphs of the spec give for a call with no `pops`.
const OVER_EVERY_INDIVIDUAL = {
  pops: ["pop"],
  obs_het: { mean: [0.395833], hist: [[1, 2, 0, 1]] },
  maf: { mean: [0.699008], hist: [[0, 1, 0, 3]] },
  exp_het: { mean: [0.378107], hist: [[2, 1, 0, 1]] },
  unbiased_exp_het: { mean: [0.430159], hist: [[1, 2, 0, 1]] },
  poly_vars_ratio: {
    num_poly: [4],
    num_variable: [4],
    tot_num_variants_with_data: [4],
    poly_ratio: [1],
    poly_ratio_over_variables: [1],
  },
};

// The four statistics that have a mean and a histogram, and the counts of
// the polymorphism ratio that are whole numbers, which are compared exactly
// while its two ratios are compared as the means are.
const DISTRIBS = ["obs_het", "maf", "exp_het", "unbiased_exp_het"];
const POLY_COUNTS = ["num_poly", "num_variable", "tot_num_variants_with_data"];

// The two rates of each of the five individuals, from "How it is verified"
// of the per individual statistics: i5 has two half called genotypes, which
// are missing, and one called genotype, which is not heterozygous.
const PER_INDIVIDUAL = {
  individuals: ["i1", "i2", "i3", "i4", "i5"],
  missing_gt_rate: [2 / 6, 2 / 6, 2 / 6, 3 / 6, 5 / 6],
  obs_het_rate: [1 / 4, 3 / 4, 1 / 4, 1 / 3, 0],
};

// How far a number popnei gives may be from the one the spec prints. The
// spec prints the means that are not exact to six digits after the point,
// 0.729167 for the major allele frequency of pop2, whose value is
// 0.7291666...
const DIGITS_OF_THE_SPEC = 1e-6;

/**
 * The values of `found` that are further than that from the ones the spec
 * gives, each as a sentence that names where it is, and one sentence when
 * there are not as many values as the spec gives.
 */
function numbersThatDiffer(what, found, expected) {
  if (found.length !== expected.length) {
    return [
      `${what} has ${found.length} values and the spec gives ${expected.length}`,
    ];
  }
  const differences = [];
  for (const [index, value] of found.entries()) {
    if (!(Math.abs(value - expected[index]) <= DIGITS_OF_THE_SPEC)) {
      differences.push(
        `${what} is ${value} at ${index} and the spec says ${expected[index]}`,
      );
    }
  }
  return differences;
}

/**
 * The same for whole numbers and for names, which are compared exactly: one
 * sentence with both lists when they are not equal, and nothing when they
 * are.
 */
function countsThatDiffer(what, found, expected) {
  if (JSON.stringify(found) === JSON.stringify(expected)) {
    return [];
  }
  return [
    `${what} is ${JSON.stringify(found)} and the spec says` +
      ` ${JSON.stringify(expected)}`,
  ];
}

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
// The genotypes of a block are a numpy array and the values of a statistic
// per population are a pandas series, so popnei cannot be imported before
// the two are in pyodide. Both are packages of pyodide itself and are
// loaded from there; micropip would fetch them for the dependencies of the
// wheel anyway, and asking for them here says which ones answer.
await pyodide.loadPackage(["micropip", "numpy", "pandas"]);
const micropip = pyodide.pyimport("micropip");
await micropip.install(`emfs:${wheelInPyodide}`);

const foundVersion = pyodide.runPython("import popnei\npopnei.__version__");
const numpyVersion = pyodide.runPython("import numpy\nnumpy.__version__");
const pandasVersion = pyodide.runPython("import pandas\npandas.__version__");
console.log(
  `popnei.__version__ is ${foundVersion}, on numpy ${numpyVersion} and` +
    ` pandas ${pandasVersion}`,
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

pyodide.runPython(THE_WORKED_EXAMPLE);
const workedExampleVcf = "/vcf/worked_example.vcf";
pyodide.runPython(`write_the_worked_example("${workedExampleVcf}")`);

for (const [asked, pops, expected] of [
  ["the two populations", THE_TWO_POPS, OVER_THE_TWO_POPS],
  ["no pops", EVERY_INDIVIDUAL, OVER_EVERY_INDIVIDUAL],
]) {
  // The populations cross as JSON, whose double quotes the call puts
  // between single ones, and a `pops` of null is Python's `None`.
  const call =
    `per_var_distribs_of_the_worked_example("${workedExampleVcf}",` +
    ` '${JSON.stringify(pops)}')`;
  const found = JSON.parse(pyodide.runPython(call));
  const what = `the worked example of calc_per_var_distribs with ${asked}`;
  const differences = [
    ...countsThatDiffer(
      `${what}: the variants of the pass`,
      found.num_vars,
      THE_WORKED_EXAMPLE_NUM_VARS,
    ),
    ...countsThatDiffer(`${what}: the populations`, found.pops, expected.pops),
    ...numbersThatDiffer(
      `${what}: the edges of the bins`,
      found.hist_bin_edges,
      HIST_BIN_EDGES,
    ),
  ];
  for (const stat of DISTRIBS) {
    differences.push(
      ...numbersThatDiffer(
        `${what}: the mean of ${stat}`,
        found[stat].mean,
        expected[stat].mean,
      ),
      ...countsThatDiffer(
        `${what}: the histogram of ${stat}`,
        found[stat].hist,
        expected[stat].hist,
      ),
    );
  }
  for (const [name, values] of Object.entries(expected.poly_vars_ratio)) {
    const ofThePolymorphism = `${what}: ${name} of the polymorphism ratio`;
    const differ = POLY_COUNTS.includes(name)
      ? countsThatDiffer
      : numbersThatDiffer;
    differences.push(
      ...differ(ofThePolymorphism, found.poly_vars_ratio[name], values),
    );
  }
  if (differences.length > 0) {
    failures.push(...differences);
  } else {
    console.log(
      `${what}: the five statistics over ${expected.pops.length} of them, as` +
        " the spec says",
    );
  }
}

const perIndividual = JSON.parse(
  pyodide.runPython(
    `per_individual_stats_of_the_worked_example("${workedExampleVcf}")`,
  ),
);
const ofTheIndividuals = "the worked example of calc_per_individual_stats";
const perIndividualDifferences = [
  ...countsThatDiffer(
    `${ofTheIndividuals}: the variants of the pass`,
    perIndividual.num_vars,
    THE_WORKED_EXAMPLE_NUM_VARS,
  ),
  ...countsThatDiffer(
    `${ofTheIndividuals}: the individuals`,
    perIndividual.individuals,
    PER_INDIVIDUAL.individuals,
  ),
  ...numbersThatDiffer(
    `${ofTheIndividuals}: the missing rate`,
    perIndividual.missing_gt_rate,
    PER_INDIVIDUAL.missing_gt_rate,
  ),
  ...numbersThatDiffer(
    `${ofTheIndividuals}: the heterozygosity rate`,
    perIndividual.obs_het_rate,
    PER_INDIVIDUAL.obs_het_rate,
  ),
];
if (perIndividualDifferences.length > 0) {
  failures.push(...perIndividualDifferences);
} else {
  console.log(
    `${ofTheIndividuals}: the two rates of the` +
      ` ${PER_INDIVIDUAL.individuals.length} individuals, as the spec says`,
  );
}

for (const failure of failures) {
  console.error(failure);
}
if (failures.length > 0) {
  process.exit(1);
}
