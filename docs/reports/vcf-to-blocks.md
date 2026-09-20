# Work report: from a VCF to blocks

The report of the run of `docs/plans/vcf-to-blocks.md`, on the branch
`plan/vcf-to-blocks`, in the worktree `.claude/worktrees/vcf-to-blocks`.
It is written while the work goes, a section after each work package.
The run started on 20 September 2026, on the owner's Apple M5 Pro. The
subagents that write the code run on Opus, as the owner ordered, and the
orchestrator and the reviewers on the model of the session.

## Before the first task

Everything of "What has to be in place" of the plan was checked by
running it, on main at 134c184: `emcc --version` 5.0.3 after sourcing
`~/devel/emsdk/emsdk_env.sh`; `~/devel/pyodide-venv/bin/pyodide` on a
Python 3.14.5 with the GIL enabled, with the cross build env of 314.0.7
installed; `wasm32-unknown-unknown` and `wasm32-unknown-emscripten`
among the installed targets; `wasm-bindgen 0.2.128`; rustc 1.98.0,
uv 0.12.15, node 26.8.2; bcftools 1.24, bgzip and plink2 in
`/opt/homebrew/bin`; pyNei at ef0ca6e with `uv run pytest
test/test_vcf.py` giving `9 passed`; and
`tests/reference/vcf/make_reference.py`, which wrote its files again
with no change for git.

The owner has not yet said how pyNei is depended on, a path or a git
commit. Task 4.2 follows the plan's meanwhile, the absolute path.

## Work package 1: the workspace and its four builds

The tasks, in the order they ran.

Task 1.1, commit 71279e5, one subagent run of 82 thousand tokens and 3
minutes. The orchestrator ran the checks again: `cargo fmt --all
--check` exit 0, clippy finished with no warning, `cargo test
--workspace` `2 passed`, and `cargo test -p popnei --lib -- --list`
`2 tests`. The subagent added a test that gzips a line and reads it
back, so that flate2 is called and not only compiled, and built the core
for both wasm targets: `wasm32-unknown-unknown` and
`wasm32-unknown-emscripten` both finish, with miniz_oxide and no C
compiler. The plan's check `cargo test -p popnei -- --list` ends with the
line of the doc tests, `0 tests`; the count of the library is read with
`--lib`.

Change to the plan: tasks 1.2 and 1.3 were marked as side by side, and
they run one after the other. Both add a member to the root
`Cargo.toml` and both write `Cargo.lock`, and a `cargo clippy
--workspace` of one would meet the half written crate of the other. One
tree has one writer per file.

Task 1.2, commit 63b6f66, one subagent run of 116 thousand tokens and 9
minutes. The orchestrator ran the five checks again: fmt exit 0, clippy
with no warning, `cargo test --workspace` `2 passed`, ruff `2 files
already formatted` and `All checks passed!`, `maturin develop` `Installed
popnei-0.1.0` and pytest `1 passed`. Four things it decided that the task
did not say, each with its reason in the commit message:

- The binding crate has `test = false` and `doctest = false`. `cargo test
  --workspace` linked a test binary of it against the Python of the PATH,
  Apple's 3.9.6, and failed at the link. The crate holds no calculation
  and its tests are the pytest ones; clippy still checks it.
- `.python-version` says 3.14.5 and not 3.14. On this Mac uv answers
  `3.14` with the free threaded 3.14.7, the newest 3.14 it has, the same
  trap `docs/rust_core.md` records for pyodide-build.
- `[tool.uv] package = false`, so that maturin alone installs popnei. An
  explicit `uv sync` takes the module out, and `uv run maturin develop`
  has to follow it.
- ruff is given `python/`, `tests/*.py` and `pyproject.toml`. Over the
  whole repository ruff 0.16.8 also formats the Python inside Markdown
  fences, 45 files, 32 of them under `.claude/` and 9 under `docs/`.

It also confirmed what `.claude/skills/coding/pyo3.md` left to confirm:
maturin builds the module without the `extension-module` feature of
pyo3 and without a warning.

Task 1.3, commit f03f6b6, one subagent run of 111 thousand tokens and 9
minutes. The orchestrator ran again the three cargo checks, which pass
with the new member, `npm run build` in `js/popnei`, which ends at the
TypeScript compiler with no error, and `npm test`, `pass 2`, `fail 0`.
The core built for `wasm32-unknown-unknown` with wasm-bindgen at the
first try, which was the risk this work package was put first for. What
it decided: one build of wasm-bindgen, `--target web`, and two entry
points that the `node` condition of `exports` in `package.json` chooses
between, node reading the `.wasm` file from disk and a page or a bundler
fetching it; `version()` called before `await init()` throws a named
`Error`. What it learned, for the reference on wasm-bindgen that the
`coding` skill lacks: the generated code passes `forbid(unsafe_code)`;
natively the crate needs only `test = false`; `--target web` under node
fails with `fetch failed` on a `file:` URL, which is why node has its own
entry; `--target nodejs` gives CommonJS, which a package with `"type":
"module"` cannot use; `--target bundler` imports the `.wasm` as a module
and needs no `init`; the `.d.ts` carries the Rust doc comments. Not
tried: `--target bundler` and a real page. `@types/node` is 24, because
npm has no types for node 26.

Task 1.4, commit c70aa0b, one subagent run of 114 thousand tokens and 9
minutes. `pyodide config get rust_toolchain` says 1.93.0 and the
workspace asks for 1.98, which looked like a clash: that toolchain is
used by `pyodide build-recipes` alone, and `pyodide build` compiled the
wheel with the rustc 1.98.0 of the PATH. The smoke test downloads
micropip from the CDN of pyodide the first time and keeps it under
`node_modules/`, so it needs the network once.

### The deliverables, run by the orchestrator at c70aa0b

1. `cargo fmt --all --check` exit 0; `cargo clippy --workspace
   --all-targets -- -D warnings` no warning; `cargo test --workspace`
   `2 passed`; `cargo test -p popnei --lib -- --list` `2 tests`.
2. `uv run maturin develop` `Installed popnei-0.1.0`, `uv run pytest`
   `1 passed`, ruff `2 files already formatted` and `All checks passed!`.
3. `npm run build` in `js/popnei` ends at the TypeScript compiler with no
   error, `npm test` `pass 2`, `fail 0`.
4. `bash scripts/build_pyodide_wheel.sh`, with `dist/` removed first,
   left `popnei-0.1.0-cp314-cp314-pyemscripten_2026_0_wasm32.whl`, and
   `node tests/pyodide/smoke.mjs` printed `popnei.__version__ is 0.1.0`
   and exited with 0.

### The review

Four reviewers over 1ec6bc3..c70aa0b: `spec`, against the work package
and sections 5, 8 and 11 of the architecture, since no spec of a module
is behind this work, `tests`, `architecture` and `binding`; 79, 84, 75
and 96 thousand tokens. The orchestrator checked the first three
findings below by running them. Every finding that held went back to
the subagent that wrote that side, in two commits.

What held and was fixed:

- `js/popnei/package.json` said `"license": "MIT"`. No document of the
  repository names a license and the owner has not chosen one; three
  reviewers reported it. An earlier version of this report told the
  owner that no manifest had a license field, from the word of the
  subagent of task 1.1, which was true of the cargo manifests when it
  was said. The field is removed, and the license is the owner's to
  choose before anything is published.
- `cargo build --workspace` failed, `ld: library 'python3.9' not found`,
  on the Python binding crate. clippy never links and `test = false`
  kept `cargo test` away from it, so no check saw it. Fixed in 69b8b90
  without naming an interpreter, which would have broken cargo on a
  checkout with no `.venv`: `.cargo/config.toml` sets
  `PYO3_BUILD_EXTENSION_MODULE`, the variable pyo3 0.29 reads now that
  its `extension-module` feature is deprecated, and a `build.rs` of the
  crate calls `pyo3_build_config::add_extension_module_link_args()`.
  Both are needed.
- Nothing ran a wasm build between this work package and the fourth.
  `cargo wasm-check` is a new alias that checks the core for the two
  wasm targets, and the plan has it among the checks of work packages 2
  and 3. It does not find a thread, which compiles for
  `wasm32-unknown-unknown` and fails when it runs.
- No test checked the version the Python distribution is published
  under: with 0.9.0 in the binding crate's manifest everything passed
  and the wheel was named 0.9.0. There is a test now, seen to fail.
- The doc comment on the `#[pyfunction]` was the docstring of
  `popnei._core.version`, which `.claude/skills/coding/pyo3.md` asks to
  be empty. It is `None` now.
- `rust-version = "1.98"` has a comment: it follows the toolchain the
  project is built with, and no older one was tried.
- In the TypeScript side, the comment of `crates/popnei-js/src/lib.rs`
  that told a writer to silence a lint that does not fire; the two
  documented behaviours without a test, the `Error` before `init()` and
  `init()` loading the wasm once; `src/web.ts`, the entry of every
  consumer that is not node, which nothing ran; the README's claim that
  a bundler copies the `.wasm`, which a reviewer found true of vite
  8.3.0 and false of esbuild 0.28.2; no `main` and `types` beside
  `exports`, so TypeScript with `moduleResolution: node10` could not
  import the package; a reader of the version in the tests that refused
  `version="0.1.0"` without spaces; and `npm test` not type checking the
  tests.

The fixes are 69b8b90, the Python and cargo side, and e82157c, the
TypeScript side, each by the subagent that wrote that side, resumed,
with 31 and 37 thousand tokens more. After them the orchestrator ran
everything again: fmt exit 0, clippy no warning, `cargo test
--workspace` `2 passed`, `cargo build --workspace` finished, `cargo
wasm-check` finished, ruff clean, pytest `2 passed`, `npm run build`
with no error and `npm test` `tests 6`, `pass 6`, `fail 0`, the wheel
built and the smoke test printed 0.1.0.

What held and was not changed:

- `gzip_is_read_back_as_it_was_written` uses flate2 from test code only,
  so the test stays green with flate2 moved to the development
  dependencies or given its C backend. Task 3.1 uses flate2 in the
  reader, and from then on `cargo wasm-check` compiles it for wasm; a
  function with no caller was not added to stand in for that.
- Two numbers of commit messages. 63b6f66 says ruff took 45 files over
  the whole repository, and the reviewer counted 46. c70aa0b says the
  wheel builds in 3.5 s, and the reviewer measured 1.5 s warm and 4 s
  when pyo3 is recompiled for emscripten. The commits are not rewritten.

Seen outside the scope, and left: `pyodide build` also leaves a wheel
under `target/wheels/`, and `package.json` has no `prepack` step, which
matters when the package is published.

### For the owner

The license of popnei is not chosen. Nothing needs it until something
is published, and no manifest names one now.

The `coding` skill has no reference for wasm-bindgen. What task 1.3 and
its review learned is above and in `js/popnei/README.md`, to write one
from.

## Work package 2: the variant module

Task 2.1, the whole work package, commits 813dacb, the spec, and
175163d, the code; one subagent run of 106 thousand tokens and 6
minutes. It finished as planned.

The deliverables, run by the orchestrator at 175163d: `cargo test -p
popnei --lib variant:: -- --list` `4 tests`, the three types and the two
variants read through a `Box<dyn VariantReader>`; fmt exit 0, clippy no
warning, `cargo test --workspace` `7 passed`, `cargo wasm-check`
finished.

What the subagent added to `docs/specs/variant.md`, in a commit of its
own before the code, for the owner to look at. None changes what a user
of Python or TypeScript sees:

- A `ChromTable` that is full, 4294967295 names, gives `u32::MAX` for a
  new name and does not keep it, and `name` gives `None` for that
  number. `intern` returns a `u32` and has no error to give, and the
  rule is what keeps it from panicking.
- The error of a field that was asked for and not filled carries the
  fields as a `Needs`, so a consumer that lacks two reports both.
- `ChromTable::is_empty`, which clippy asks for beside `len`. The code
  also got `Default` for `ChromTable` and `Variant`, beside their `new`;
  an earlier version of this report said the spec had them too, on the
  word of the subagent, and it had not. The review found it, and the
  spec has them since 34f7d30.
- `Needs` is written by hand, and `thiserror` is a new dependency of the
  core, which the `coding` skill asks errors to be written with. It
  builds for both wasm targets.

For what comes next: `Variant::clear` empties `alleles` by dropping its
strings and keeps only the capacity of the vector, so a reader that
fills the alleles through `clear` allocates a string for each allele of
each variant. Section 1 of the architecture wants those buffers reused.
The spec asks only for the state a read leaves, so the VCF reader can
write over the strings that are there; task 3.2 is told, and the check
by hand of task 5.1 is where it shows.

The review of this work package is made together with that of work
package 3, as the `following-plans` skill allows when the two are one
piece of code: the types have no caller until the VCF reader, and what
a reviewer can say of them alone is little.

## Work package 3: the VCF reader, on one thread

Tasks 3.1 and 3.2 went to one subagent, commits fc9200a and 6b32bdd,
each after a commit of the spec, 93996b4 and cdcac50; one run of 216
thousand tokens and 19 minutes. Task 3.3 went to another subagent, so
that the reader was tested by somebody who had not written it: commit
10b8dd8, 127 thousand tokens and 5 minutes. It finished as planned.

### The deliverables, run by the orchestrator at 10b8dd8

`cargo test -p popnei --lib io::vcf:: -- --list` `36 tests`, where the
plan asks for 25; `cargo test --workspace` `43 passed`; fmt exit 0,
clippy no warning, `cargo wasm-check` finished. The tests read the three
reference VCFs, plain and gzipped, with the literals of the spec, and
every genotype, chromosome and position of `many.vcf` against the stored
output of bcftools. The subagent of 3.3 put `GzDecoder` in the place of
`MultiGzDecoder` and saw `many.vcf.gz` give 0 variants and both tests of
the counts fail, which is what the spec says must happen.

### The review, of work packages 2 and 3 together

Six reviewers over 9cfdb70..10b8dd8: `spec`, `tests`, `numbers`,
`errors`, `api` and `architecture`; 136, 145, 103, 100, 96 and 117
thousand tokens. The `tests` reviewer made 40 mutations of the code and
ran the tests on each. Every finding that held went back to the subagent
that wrote the reader, which fixed them in 13 commits from 34f7d30, the
specs first, with 132 thousand tokens more. After them: `cargo test
--workspace` `65 passed`, `io::vcf::` `57 tests`, `variant::` `5 tests`,
clippy no warning, `cargo wasm-check` and `cargo build --workspace`
finished, `cargo doc -p popnei --no-deps` with no warning, pytest `2
passed`, ruff clean.

What held and was fixed, the ones that matter first:

- The gzip bytes were looked for in one fill of the buffer, which may
  hold one byte: `cases.vcf.gz` through a `BufReader` of capacity 1 was
  refused as "not a VCF". A pipe or a source of bytes in the browser can
  do that. Three reviewers reported it.
- A variant with more alleles than the one before it allocated a
  string. A counting allocator gave 54 allocations in the last 400
  variants of a second pass over `many.vcf`, its 54 variants with two
  alternative alleles, and 2 after the fix, the interning of `chr2`.
  The reader now keeps the strings it does not use.
- The reader repeated the list of fields of `Variant::clear`, so a field
  added later would have kept the value of the variant before. There is
  one list now, `Variant::clear_but_the_alleles`.
- A FORMAT with no `GT` and a line with too few columns were errors only
  when the genotypes were asked for, where the spec lists them with no
  condition. The orchestrator took the spec's list: the nine fixed
  columns, the `GT` key and one column after FORMAT are checked whatever
  is asked for, and only the count and the contents of the columns of
  the individuals wait for the genotypes. The spec says it now.
- Accepted with no error, and refused now: an ALT with an empty piece,
  which was counted as an allele that a genotype could carry; an
  individual with no name, from a header that ends in a tab; `+1/0`.
- A ploidy of `usize::MAX` panicked with "capacity overflow" at a
  missing genotype, and 2^40 allocated until it was killed. A ploidy
  above `MAX_PLOIDY`, 255, is refused. The orchestrator chose 255; no
  organism is near it.
- A line that is not UTF-8 was an error of the input with no line
  number, an `OSError` in Python for a malformed line. It is an error of
  the line with its number.
- A file that is not there gave "No such file or directory" with no
  name. The error has a case with the path, so that the Python binding
  can fill `OSError.filename`.
- The message of a genotype of another ploidy did not say that popnei
  does not read a VCF of mixed ploidies, which the spec asks of it, and
  said "1 alleles".
- Six public items were in the code and not in the specs, and
  `ChromTable`, `Variant` and `VcfReader` had no `Debug`.
  `missing_docs` is denied now, in the workspace and in
  `.claude/skills/coding/lints.toml`.
- Two rules of the spec that no test would have noticed broken: that a
  field no longer asked for is emptied, and that the columns of the
  individuals are not looked at without the genotypes. Both have tests,
  and so have a QUAL that is not a number, a gzip member cut in the
  middle, a read after an error, a last line with no end and an empty
  source.
- One sentence that cdcac50 added to the spec claimed more than the code
  does, and is corrected.

What held and was not changed:

- `VcfReader::from_path` is compiled for wasm too. It is needed there:
  under pyodide, which is wasm, the files of emscripten are real to the
  reader.
- The message of fc9200a says that `cases.vcf.gz` and `many.vcf.gz` are
  four and two gzip members. They are three and four, as the spec says.
  The commit is not rewritten.
- `/.` is read as a missing genotype, where bcftools refuses the line.
  The rule of the leading separator and the rule of the single dot both
  apply, and the spec does not say which wins.

Passed on to work package 5: the chromosome is given its number while
the line is parsed, and task 5.1 has to move that to the moment a
variant is handed out, as the spec says, because the workers of rayon
cannot share the table. And a first measurement, by the `architecture`
reviewer, on a VCF it generated of 10000 variants x 1000 individuals,
one thread, release, this machine, genotypes asked for: 0.134 s, which
is 1.34 s for the 100000 variants of the spec's target of 0.55 s, nearly
all of it in the genotypes. The fix of the allele numbers made the parse
of a genotype one pass where it was two, and nothing was measured after.
Task 5.2 measures on the real file.

### For the owner

Three decisions, asked in chat on 20 September 2026, each with what the
code does meanwhile:

1. A bgzipped file cut at the boundary of a gzip member gives fewer
   variants and no error, where bcftools says "no BGZF EOF marker".
   Recommended: detect it. Meanwhile it is as it is.
2. A QUAL of `nan`, `inf` or `1e400` parses as a float and is taken, as
   bcftools does. Recommended: refuse a quality that is not finite.
   Meanwhile it is taken.
3. The specs have one error enum for the crate, and the `coding` skill
   says "one error type for each operation that fails in its own way,
   not one for the crate". The code follows the specs. Recommended: keep
   the one enum and correct the skill.
