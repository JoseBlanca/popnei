# popnei under pyodide

pyodide is CPython built for WebAssembly, which runs in a browser tab and
under node. popnei is installed in it as a wheel whose platform tag is
`pyemscripten`, built from the same core crate and the same binding crate
as the native wheel, and what this directory checks is that the wheel
builds, that it installs, and that inside pyodide popnei answers with the
version of the core, reads a VCF and calculates over its variants as the
specs say it does.

From the root of the repository, from a clean checkout:

    scripts/build_pyodide_wheel.sh
    (cd tests/pyodide && npm install)
    node tests/pyodide/smoke.mjs

The first command prints the path of the wheel it left in `dist/`. The
third takes that wheel, installs it in pyodide and checks eight things, and
exits with an error naming each one that differs:

- `popnei.__version__` is the version in `[workspace.package]` of the
  `Cargo.toml` of the repository.
- `popnei.open_vcf` reads `tests/reference/vcf/cases.vcf`, and
  `cases.vcf.gz`, into the four variants that the table of "How it is
  verified" of `docs/specs/io_vcf.md` gives for that file. The position and
  the genotypes of each of the four are in the test as literals, taken from
  that table. The second of the four failed a filter, `q10` in its FILTER
  column, so `open_vcf` with its default leaves it out and the test looks
  for the first, the third and the fourth; with `only_passed=False` it
  looks for the four.
- `popnei.write_vars` writes those four variants into a vars file, the
  arrow file of `docs/specs/io_vars.md`, at a path of the file system of
  emscripten, and the bytes of that file begin and end with `ARROW1`, which
  an arrow IPC file does. It is the crates of arrow that write those bytes,
  so the check says that they linked into the wheel of emscripten and ran
  there; what the file holds is checked natively, by the pytest tests of
  `tests/test_io_vars.py`.
- A VCF of 170000 individuals and no variant, written inside pyodide and
  opened with the ploidy 255, gives its individuals and its ploidy. This is
  the check that belongs here and nowhere else: a count of things in wasm
  is 32 bits and holds 4295 million, and the blocks of the size popnei
  chooses for that file, 100 variants, are 4335 million genotypes. Opening
  a file reads its header and asks for no block, so it answers; blocks of
  10 variants are asked for and the file has no variant to put in one; and
  blocks of 100, and the size popnei chooses, are refused with the
  `ValueError` of a block the machine has not the memory for. The same case
  is in `js/popnei/test/open.test.ts`, under node, where a count of things
  is 64 bits and nothing is refused for its size.
- `popnei.calc_per_var_distribs` gives the means, the histogram counts and
  the counts of the polymorphism ratio of the worked example of "How it is
  verified" of the per variant distributions of `docs/specs/stats.md`: six
  variants of five diploid individuals, written inside pyodide as a VCF,
  with `min_num_individuals` 1 and four bins from 0 to 1, once over the two
  populations of that example and once with no `pops`, which is one
  population of the five individuals. Every number is in the test as a
  literal from that spec, and the means, which the spec prints to six
  digits after the point, are compared within 1e-6.
- `popnei.calc_per_individual_stats` gives the missing rate and the
  heterozygosity rate of the five individuals of the same six variants,
  which the worked example of "How it is verified" of the per individual
  statistics gives. These two checks are what says that the calculations
  reach the same values in a build that has one thread and is not the
  native one; the same two examples run natively in
  `crates/popnei/src/stats.rs`.
- `popnei.calc_pairwise_kosman_dists` gives the three distances of the
  diploid worked example of "How it is verified" of `docs/specs/dists.md`,
  4 variants of 3 individuals, written as a VCF into the file system of
  emscripten and opened with `open_vcf`: 1 over 4, 5 over 6 and 2 over 6,
  compared exactly, since the calculation divides those two whole numbers
  once and a division is rounded the same in JavaScript as in Rust. It
  checks with them the names of the three individuals and the 4 variants
  the pass took, and that with `min_num_snps=3` the pair that was called
  at 2 variants has no distance and the other two keep theirs. The same
  example, with the same VCF, is a pytest test of `tests/test_dists.py`.
- `popnei.calc_pop_diversity` gives, over those same six variants and
  those same two populations with `min_num_individuals` 1, the numbers of
  the worked example of "How it is verified" of each item of
  `docs/specs/diversity.md`, which is the worked example of
  `docs/specs/stats.md`. It is called twice, because the two calls
  exercise different things. The plain call, with no statistic named and
  no `num_called_alleles`, asks for the four statistics that need no draw
  and gives the alleles `pop1` and `pop2` called, 9 and 8, the private
  ones among them, 2 and 1, the variants that vary in each, 3 and 2 of the
  4 that counted, and F_IS, 0 and 0.3478260870, with no spectrum and with
  the three standardized values missing. The second names the five
  statistics and a draw of 4 called alleles, which `pop1` reaches at its
  four variants and `pop2` at three of them, and gives the three
  standardized values and the folded spectrum, 1, 3 and 0 for `pop1` and
  1.0666666667, 1.5333333333 and 0.4 for `pop2`. Every number is in the
  test as a literal from that spec, and the floats are compared within
  1e-6, as the means of the worked example above are. The same example is
  a cargo test of `crates/popnei/src/diversity.rs`.

Neither `dist/` nor `node_modules/` is in git.

`time_pca.mjs`, beside the smoke test, is not a test: it times
`do_pca_from_variants` inside pyodide over a vars file, with the wheel that
`dist/` holds, and it is what task 4.2 of `docs/plans/pca.md` measured the
wheel with. `docs/reports/pca-measurement.md` has its numbers and the file
it read them on.

    node tests/pyodide/time_pca.mjs <path to a vars file> [runs]

## What has to be installed

Outside the repository, on the machine that builds. "What has to be in
place" of `docs/plans/vcf-to-blocks.md` has the commands that installed
them here.

- The Rust target `wasm32-unknown-emscripten`, `rustup target add
  wasm32-unknown-emscripten`. cargo names it when it is not installed.
- emsdk, the installer of emscripten, in `~/devel/emsdk`, with 5.0.3
  activated in it. emscripten is the C compiler and linker that compiles to
  WebAssembly, and the CPython of pyodide 314.0.7 was compiled with that
  version of it. The build script sources `emsdk_env.sh`, which puts `emcc`
  on the PATH together with the node and the Python that emsdk downloaded
  for it: `emcc` refuses to run under a Python older than 3.10, and the
  only `python3` of the PATH on this Mac is Apple's 3.9. Set `EMSDK` to the
  directory of another emsdk to build with that one.
- pyodide-build in the virtual environment `~/devel/pyodide-venv`, with the
  cross build environment of pyodide 314.0.7 installed in it, which is
  where the headers and the libraries of the CPython of pyodide are. Set
  `PYODIDE_VENV` to the directory of another virtual environment.
- node, for the smoke test, which takes pyodide itself from npm.

The build script stops with a message when emsdk or pyodide-build is not
where it looks, when the host Python of pyodide-build is the free threaded
build of CPython, and when the version of `emcc` is not the one pyodide
asks for.

## Moving to another version of pyodide

Three things name the version, and the three have to agree, because a wheel
built against one build of the CPython of pyodide does not load in another:
the cross build environment, `pyodide xbuildenv install <version>`; the
emscripten that was activated in `~/devel/emsdk`; and the `pyodide`
dependency of `package.json` beside this file, which is the runtime the
smoke test loads.

The build script reads the version of emscripten from `pyodide config get`
and checks `emcc` against it, so changing the cross build environment and
the emsdk is enough on that side; the npm dependency is the one line to
change by hand.

A newer pyodide brings a newer numpy and a newer pandas with it, and the
smoke test prints which ones answered. Each has to be one that
`dependencies` of `pyproject.toml` accepts, and the numpy one that the
`numpy` crate of the binding crate supports; when the numpy is not, the
wheel still installs and the failure comes when the genotypes are asked
for, and when the pandas is not, micropip refuses to install the wheel,
since it finds no other pandas to fetch for a wasm build.

## Two traps

Both were found when the trial crate in Rust that `docs/rust_core.md`
reports was built for pyodide.

- The host Python that pyodide-build runs under must not be the free
  threaded build of CPython, the one without the global interpreter lock,
  which `uv venv --python 3.14` picks on this Mac. With it pyodide-build
  writes the sysconfig data of emscripten under `python3.14t`, the `t` of
  which marks that build, where nothing looks for it. The build script
  stops when the host Python of pyodide-build is that one, and
  `.python-version`, the file uv reads to choose the interpreter of the
  `.venv` of popnei, names 3.14.5 and not 3.14 for the same reason.
- `pyodide config get rust_toolchain` answers 1.93.0, older than the
  `rust-version = "1.98"` of the workspace, and that is not a problem:
  pyodide-build installs and switches to that toolchain in
  `pyodide build-recipes`, which builds the packages of the pyodide
  distribution from its own recipes, and not in `pyodide build`, which is
  what this script runs. The wheel here is compiled by whichever `cargo`
  is on the PATH, 1.98.0 on this machine.

## What runs where

The code running inside pyodide sees the file system that emscripten gives
it and not the one of node, and nothing of the repository is in it. So the
smoke test writes the wheel there and installs it from there with micropip,
the package installer of pyodide; the `emfs:` prefix of the path is what
asks micropip for a file of that file system instead of a package of an
index. It writes the two VCFs there as well, under `/vcf`, because
`open_vcf` takes a path and opens it, as it does natively.

Neither micropip, nor numpy, nor pandas is in the npm package of pyodide:
the first run downloads the three of them, and the python-dateutil, the pytz
and the six that pandas brings with it, from the CDN of pyodide and caches
them under `node_modules/`, so that first run needs the network and the
later ones do not. The whole test takes 2.3 s once they are cached,
measured on 22 September 2026 on the owner's Apple M5 Pro, 18 cores, 64 GB,
macOS 27.0.

A block is a run of consecutive variants of a source, held as arrays, and
it is how the genotypes leave popnei: they arrive in a numpy array of
variants x individuals x ploidy. The mean of a statistic and the counts of
its histogram leave it as pandas series and frames, one value per
population, and the square matrix of a `Distances` is a pandas frame, which
`popnei/dists.py` imports when popnei is imported. So popnei cannot be
imported before numpy and pandas are loaded. pyodide 314.0.7 brings numpy
2.4.6 and pandas 3.0.2, which are what `dependencies` of `pyproject.toml`
asks for, the numpy being what the `numpy` crate 0.29 of the binding crate
was compiled against and the pandas being the one pyodide ships, since
micropip refuses a wheel that asks for more than pyodide has; the versions
the test prints are the ones that answered.
