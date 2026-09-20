# Plan: from a VCF to blocks, in Rust, Python and TypeScript

September 2026. Approved by the owner on 20 September 2026, who ordered
it run with the code written by subagents on Opus. Under way since that
day, with its report in `docs/reports/vcf-to-blocks.md`. This plan builds
the first code of popnei: the workspace with its four builds, the variant
record, the VCF reader, and the blocks through which genotypes reach
Python and TypeScript. It is built from three specs, all reviewed and
with no open point:

- `docs/specs/variant.md`, the record and the reader trait;
- `docs/specs/io_vcf.md`, the VCF reader;
- `docs/specs/block.md`, the block and its collector.

When it is done a user can call `open_vcf` on a VCF, in Python natively
and under pyodide, and `openVcf` in TypeScript, and get its genotypes
block by block, the same ones pyNei reads, at the speed the spec of the
reader asks for. The owner answered its breakdown on 20 September 2026:
five work packages, with the second and the third ending at cargo tests
and the comparison with pyNei first made in the fourth.

## In and out

In: what the three specs describe, and the layout of section 8 of
`docs/architecture.md`, with its two binding crates and its two packages.

Out, with where it goes:

- A `File` that the user picked in a page, read inside a web worker.
  TypeScript reads a `Uint8Array` in this plan. Section 11 of the
  architecture has the design, and it goes with the first web
  application.
- The vars file, the filters and the VCF writer: their specs are not
  written, and with them the rest of the walking skeleton of section 10
  of the architecture.
- The row helpers of the `variant` module, and the items of the `block`
  spec that are not written: the view of one variant, `reblock`, the
  blocks of the vars file.
- Publishing to PyPI or to npm, continuous integration and the
  development container.
- Work on the speed of the reader beyond measuring it. If the target of
  the spec is missed, the `performance-review` skill takes it from the
  measurement of task 5.2.

## What has to be in place

Checked by running it, on the owner's Mac, where the whole plan runs.
When one of these is missing the plan does not start and the owner is
told, as the `following-plans` skill says; what installs the wasm target,
wasm-bindgen and the toolchain of pyodide is at the end of this part.

- `rustc --version` gives 1.98 or later, and `rustup target list
  --installed` has `wasm32-unknown-unknown` and
  `wasm32-unknown-emscripten`.
- `wasm-bindgen --version` answers. Task 1.3 pins the `wasm-bindgen`
  crate to that version, because the two have to be the same.
- `uv --version`, `node --version`, 26 or later, and `npm --version`
  answer.
- `bcftools --version` gives 1.24 and `which bgzip plink2` finds both.
- For the wheel of pyodide: `~/devel/emsdk/emsdk_env.sh` exists and,
  after it is sourced, `emcc --version` gives 5.0.3;
  `~/devel/pyodide-venv/bin/pyodide --version` answers, and
  `~/devel/pyodide-venv/bin/python -c "import sys;
  print(sys._is_gil_enabled())"` prints `True`, because the host Python of
  pyodide-build must not be the free threaded one, the trap of section
  3.2 of `docs/rust_core.md`; and `pyodide xbuildenv versions` shows
  314.0.7 installed.
- `/Users/jose/devel/pynei` is a checkout of pyNei at ef0ca6e or later,
  in which `uv run pytest test/test_vcf.py` passes.
- `tests/reference/vcf/` has the three VCFs of `docs/specs/io_vcf.md`,
  the two gzipped ones and the three outputs of bcftools. `python3
  tests/reference/vcf/make_reference.py` writes them again and `git
  status` shows no change.

What installed the three that were missing when the plan was written.
The first `./emsdk install` needs a Python of 3.10 or later, and on this
Mac the only `python3` of the PATH is Apple's 3.9, so it was run with
`EMSDK_PYTHON` set to the 3.14.5 of uv; after it, `emsdk_env.sh` points
at a Python that emsdk downloaded, and nothing else has to be set:

    rustup target add wasm32-unknown-unknown
    cargo install wasm-bindgen-cli
    cd ~/devel && git clone --depth 1 https://github.com/emscripten-core/emsdk.git
    (cd emsdk && ./emsdk install 5.0.3 && ./emsdk activate 5.0.3)
    uv venv --python cpython-3.14.5-macos-aarch64-none ~/devel/pyodide-venv
    VIRTUAL_ENV=~/devel/pyodide-venv uv pip install pyodide-build
    source ~/devel/emsdk/emsdk_env.sh
    ~/devel/pyodide-venv/bin/pyodide xbuildenv install 314.0.7

There is no code, so none of the five commands of "Before the work is
called done" of the `coding` skill runs today: `cargo test` answers
"could not find `Cargo.toml`". The three cargo commands run from task
1.1 on, the two Python ones from task 1.2 on. The checks of the
TypeScript side and of the wheel of pyodide are not in that skill, and
this plan gives them with their deliverables.

One thing is the owner's to settle before task 4.2. `docs/objectives.md`
has pyNei as a path dependency on the sibling checkout, `../pynei`. The
plan runs in a worktree under `.claude/worktrees/`, where that relative
path points nowhere. The options are an absolute path, which works on
this machine alone, or a dependency on pyNei's git repository at a
commit, which also fixes which pyNei is the oracle. Until the owner
answers, task 4.2 uses the absolute path `/Users/jose/devel/pynei`.

## Work package 1: the workspace and its four builds

**What it gives.** Nothing to a user yet. It gives the next work packages
a workspace in which the core crate is built and tested natively, called
from Python, called from TypeScript, and built as a wheel for pyodide,
each proved with one function that returns the version of the core. It is
first because it holds what nobody has tried: the core compiled for
`wasm32-unknown-unknown` with wasm-bindgen, and the TypeScript package
around it.

**Deliverables.**

1. The workspace of section 8 of the architecture with the core crate
   `crates/popnei`, which has `flate2`, the crate that the VCF reader
   will decompress gzip with, with its default backend among its
   dependencies, so that every build of this work package compiles it.
   Check: `cargo fmt --all --check`, `cargo clippy --workspace
   --all-targets -- -D warnings` and `cargo test --workspace` pass, and
   `cargo test -p popnei -- --list` shows 1 test or more.
2. The Python binding crate and the Python package. Check: `uv run
   maturin develop && uv run pytest` passes with 1 test or more, which
   asserts that `popnei.__version__` is the version of the core crate;
   `uv run ruff format --check && uv run ruff check` is clean.
3. The JavaScript binding crate and the TypeScript package. Check: one
   command, written in `js/popnei/README.md`, builds the wasm package
   from a clean checkout, and `npm test` in `js/popnei` runs 1 test or
   more under node, which asserts the same version, and the TypeScript
   compiler finds no error.
4. The wheel of pyodide. Check: a script, `scripts/build_pyodide_wheel.sh`,
   leaves in `dist/` a wheel whose name has the `pyemscripten` tag, and
   `node tests/pyodide/smoke.mjs` loads pyodide 314.0.7, installs that
   wheel with micropip and prints the same version.

**What it stands on.** What has to be in place, above.

**Tasks.**

- [x] 1.1 The cargo workspace and the core crate with its version
  function, the lints and the settings that the `coding` skill asks of
  the core, a `.gitignore`, and `flate2` as a dependency. From section 8
  of the architecture. Serves deliverable 1.
- [x] 1.2 `crates/popnei-python` with pyo3 and its module
  `popnei._core`, `python/popnei`, `pyproject.toml` with maturin, ruff
  and pytest as development dependencies, and the first pytest test.
  From section 8 of the architecture and `.claude/skills/coding/pyo3.md`.
  Serves deliverable 2. Needs 1.1. Side by side with 1.3.
- [x] 1.3 `crates/popnei-js` with wasm-bindgen, pinned to the version of
  the command line that is installed, `js/popnei` with its
  `package.json`, the TypeScript compiler as a development dependency,
  the build command and the first node test. From sections 8 and 11 of
  the architecture. Serves deliverable 3. Needs 1.1. Side by side with
  1.2.
- [ ] 1.4 The script that builds the wheel of pyodide, from the steps of
  `spike/README.md` of the pyNei repository, the trial crate in Rust that
  `docs/rust_core.md` reports and that was built for pyodide, with the versions read from
  `pyodide config get` as `.claude/skills/coding/pyo3.md` asks, and the
  smoke test under node. Serves deliverable 4. Needs 1.2.

**What could go wrong.** The wasm package has to load under node, for
the tests, and in a bundler or a page, for an application, and
wasm-bindgen builds for one of those at a time; task 1.3 says in the
README which one the tests use and how the other is built. The `coding`
skill has a reference for pyo3 and none for wasm-bindgen, so 1.3 has
less to go on than 1.2, and what its subagent learns goes into the work
report for that reference to be written from. If the core does not build
for `wasm32-unknown-unknown`, the plan stops and the owner is told,
because every later work package builds on it.

## Work package 2: the variant module

**What it gives.** To work package 3, the `Variant` that a reader fills,
`Needs`, the `ChromTable`, the `VariantReader` trait and the error of
the crate. It stops inside the core: a user sees nothing of it until
there is a reader, and the comparison with pyNei is made in work package
4.

**Deliverables.**

1. The module `variant` of the core crate with everything of "The Rust
   interface" of `docs/specs/variant.md`. Check: `cargo test -p popnei
   --lib variant:: -- --list` shows 3 tests or more, one for each of the
   three types that "How it is verified" names, and the three cargo
   commands pass.
2. `VariantReader` works as a boxed trait object. Check: a cargo test
   with a reader written in the test, which gives two variants, read
   through a `Box<dyn VariantReader>`.

**What it stands on.** Work package 1.

**Tasks.**

- [ ] 2.1 The module: the two constants, `Needs`, `ChromTable`,
  `Variant`, the trait with its implementation for a box, the error of
  the crate with its first case, and their tests. From "The Rust
  interface", "Fields that were not asked for, and sources that lack
  one", "How it runs" and "How it is verified" of
  `docs/specs/variant.md`. Serves deliverables 1 and 2.

## Work package 3: the VCF reader, on one thread

**What it gives.** To work package 4, a `VcfReader` over any source of
bytes that reads the three reference VCFs as bcftools does, with every
rule of the spec, parsing one line after another. It stops inside the
core, and the comparison with pyNei is made in work package 4. The
threads come in work package 5, when every test of 3 and 4 is there to
guard them.

**Deliverables.**

1. `VcfReader::new` and `from_path`: the gzip bytes, the header, the
   options and their errors. Check: a cargo test for each of the four
   cases of "How it is verified" of `docs/specs/io_vcf.md` that are
   errors of `new`, and one that reads the individuals of `cases.vcf` and
   of `cases.vcf.gz`.
2. `read_variant` with every rule of "What it gives" and of "The cases a
   reader of the rules would not guess". Check: cargo tests with the
   literals of the tables of `cases.vcf` and `differences.vcf`, with the
   default and with every variant given; one for each error of a data
   line and for each case that is not an error in the list of "How it is
   verified"; and the two of what `Needs` does.
3. `many.vcf` and `many.vcf.gz`. Check: cargo tests that compare every
   genotype, chromosome and position with `many.bcftools.tsv` and the
   sixteen counts of the table of "How it is verified" as literals, for
   the plain and the gzipped file.
4. All of it. Check: `cargo test -p popnei --lib io::vcf:: -- --list`
   shows 25 tests or more, and the three cargo commands pass.

**What it stands on.** Work package 2, whose types every task here
uses, so its subagents read `docs/specs/variant.md` too, and the
reference files.

**Tasks.**

- [ ] 3.1 The source and the header: finding the gzip bytes without
  consuming them, `MultiGzDecoder`, `VcfOptions`, the `#CHROM` line, the
  individuals, the cases of the error that this module adds, and
  `from_path`. From "The Rust interface" and the first two paragraphs of
  "The cases a reader of the rules would not guess" of
  `docs/specs/io_vcf.md`. Serves deliverable 1.
- [ ] 3.2 The data line: the columns, the FILTER, the genotypes with
  their separators, the ploidy, the allele numbers against ALT and
  against 127, `Needs`, the line numbers of the errors, and the order of
  the chromosome numbers. From "What it gives", the rest of "The cases a
  reader of the rules would not guess" and "How it is verified". Serves
  deliverable 2. Needs 3.1.
- [ ] 3.3 The tests on `many.vcf` and `many.vcf.gz` against the stored
  output of bcftools, with a helper of the tests that reads that output.
  Serves deliverables 3 and 4. Needs 3.2. A wrong genotype here is
  silent anywhere else, so it is a commit of its own, and deliverable 3
  is what guards it.

**What could go wrong.** A decoder that stops at the first gzip member
gives the header and no variant, which the spec says is not an error, so
only the counts of deliverable 3 on the gzipped file find it. The tests
read the reference files from `tests/reference/vcf/`, at the root of the
repository and not inside the crate, and the path they use has to work
from `cargo test --workspace` and from `cargo test -p popnei`.

## Work package 4: blocks, and the Python and TypeScript sides

**What it gives.** A user calls `open_vcf(path, ploidy=2,
only_passed=True)` in Python, or `openVcf` in TypeScript on the bytes of
a VCF, gets a `Variants` with its individuals and its ploidy, and
iterates over its blocks with `iter_blocks`. It is the first time popnei
and pyNei read the same files side by side.

**Deliverables.**

1. `Block`, `AllelesColumn`, `BlockCollector` and
   `default_num_vars_per_block` in the module `block` of the core. Check:
   the cargo tests of "How it is verified" of `docs/specs/block.md`;
   `cargo test -p popnei --lib block:: -- --list` shows 6 tests or more,
   and the three cargo commands pass.
2. `open_vcf`, `Variants`, `iter_blocks` and `Block` in the Python
   package. Check: `uv run maturin develop && uv run pytest` passes with
   the test against pyNei of `docs/specs/block.md`, on the four files,
   with `num_vars_per_block` of 7 and the default, and the pytest tests
   that "How it is verified" of `docs/specs/io_vcf.md` lists; ruff is
   clean.
3. `openVcf`, `Variants` and `iterBlocks` in the TypeScript package.
   Check: `npm test` in `js/popnei` passes with the node test of each of
   the two specs, and the TypeScript compiler finds no error.
4. The same under pyodide. Check: the script of work package 1 builds
   the wheel, and `node tests/pyodide/smoke.mjs` reads `cases.vcf` with
   `open_vcf` under pyodide and finds the genotypes of the table of
   `docs/specs/io_vcf.md` in its blocks.

**What it stands on.** Work package 3. The owner's answer about how
pyNei is depended on, or the absolute path meanwhile.

**Tasks.**

- [ ] 4.1 The module `block` and its cargo tests. From
  `docs/specs/block.md`, everything but "In Python and in TypeScript".
  Serves deliverable 1.
- [ ] 4.2 The Python side: in the binding crate, a class that holds the
  source and its options and one that gives blocks, with the genotypes
  handed to numpy without a copy and the errors of the core turned into
  `ValueError` and `OSError`; in the package, `open_vcf`, `Variants` and
  the `Block` dataclass; pyNei as a development dependency; the pytest
  tests. From "In Python and in TypeScript" and "How it is verified" of
  `docs/specs/block.md`, "Its Python and TypeScript functions" and "How
  it is verified" of `docs/specs/io_vcf.md`, and "What a Python and a
  TypeScript user see" of `docs/specs/variant.md`. Serves deliverable 2.
  Needs 4.1. Side by side with 4.3.
- [ ] 4.3 The TypeScript side: in the binding crate, the same two
  classes over a `Uint8Array`, with the arrays copied out of the memory
  of wasm, the errors thrown as `Error`, and `free()`; in the package,
  `openVcf`, `Variants` and `iterBlocks` usable in a `for ... of`; the
  node tests. From the same parts of the three specs. Serves deliverable
  3. Needs 4.1. Side by side with 4.2.
- [ ] 4.4 The smoke test under pyodide, extended to read `cases.vcf`.
  Serves deliverable 4. Needs 4.2.

**What could go wrong.** The name of a chromosome number is looked up
after a block is collected, as `docs/specs/block.md` says, and a binding
that looks it up before gets an empty table. pyNei gives `pandas.NA` for
a missing id, where popnei has `None`, and for a missing quality, where
popnei has NaN, and the comparison has to take each pair as equal. Under pyodide numpy has to be loaded before the wheel
is imported.

## Work package 5: the parallel reader, and its speed

**What it gives.** The same reader, natively parsing its lines in
batches with rayon, and the first measurement of popnei: how long it
takes to read the 400 MB VCF of `docs/rust_core.md`.

**Deliverables.**

1. The reader parses in batches on the threads of rayon natively and one
   line after another in wasm, and gives the same. Check: every test of
   work packages 3 and 4 passes untouched; a new cargo test reads
   `many.vcf` with rayon pools of 1 and of 4 threads and gets the same
   variants and the same chromosome numbers; another has a wrong line
   after good ones and gets the good ones first and then the error, with
   4 threads; the wasm package and the wheel of pyodide still build.
2. Nothing is allocated from one variant to the next. Check: made once
   by hand with a counting allocator over a second pass of `many.vcf`,
   as "How it is verified" of `docs/specs/io_vcf.md` says, and written
   into the work report with its numbers. It is not a test that stays.
3. The measurement. Check: `benches/` or a script that the work report
   names reads the VCF of 100000 variants x 1000 individuals with `GTS`
   asked for, on 1 thread and on the 18 cores of the owner's M5 Pro,
   plain and gzipped, and the report
   has the four times beside the targets of "Speed" of
   `docs/specs/io_vcf.md`, the machine, and how the file was made.

**What it stands on.** Work package 4.

**Tasks.**

- [ ] 5.1 The batches: the lines read into a batch, parsed with rayon
  under `cfg(not(target_family = "wasm"))` with the serial version
  beside it, handed out in order by swapping buffers with the lent
  `Variant`, the chromosome numbers given at that moment, and the two
  new tests. From "How it runs" of `docs/specs/io_vcf.md` and section 3
  of the architecture. Serves deliverables 1 and 2.
- [ ] 5.2 The measurement, taken as "What measurement there is" of the
  `performance-review` skill asks: the file written with `write_vcf` of
  `test/gwas_reference/make_reference.py` of pyNei, which
  `docs/rust_core.md` names as its origin, the four times, and plink2 on
  the same file and machine for the comparison. Serves deliverable 3.
  Needs 5.1.

**What could go wrong.** How many lines a batch holds is a number nobody
has measured, and 5.2 is where it is tried. The spike did not check the
ploidy of each genotype, the allele numbers or the FILTER, so the target
may be missed by those checks alone; the plan is done when the
measurement is reported, and a missed target goes to the owner with the
numbers and not into more tasks of this plan.

## How the whole plan is checked

Each work package is checked by its own deliverables. At the end, from a
clean checkout of the branch: the five commands of the `coding` skill,
`npm test` in `js/popnei`, the build of the wheel of pyodide with its
smoke test, and the work report with the times of deliverable 3 of work
package 5.
