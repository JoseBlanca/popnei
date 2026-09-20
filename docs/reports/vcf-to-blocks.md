# Work report: from a VCF to blocks

The plan `docs/plans/vcf-to-blocks.md` is done, on the branch
`plan/vcf-to-blocks`, 20 September 2026. Its five work packages
finished as planned, each was reviewed and its findings fixed, and the
final check passes from a clean clone. The speed targets of the spec are
missed, which the plan foresaw as something to report and not to work
on. The orchestrator, in this report, is the session of the assistant
that ran the plan: it sent each task to a subagent, checked what came
back and had it reviewed.

While the last work package was being reviewed, another session of the
assistant, the one in which the owner is writing the spec of the vars
file, sent word that the owner decided on 20 September 2026 to drop the
single variant as the interface of the core: readers will give blocks,
and the `Variant` that a reader fills, `read_variant` and the
`BlockCollector` go. This session has not heard it from the owner, and
the plan was already at its end, so nothing was stopped or changed for
it. If it stands, what that session expects to survive of this branch
is the parse of a data line, the header, gzip, the errors, the batches
of lines, both binding crates, both packages and every test made at
`iter_blocks` and `iterBlocks`, since Python and TypeScript never saw a
variant; and what changes is how the VCF reader hands out what it
parsed, a row written into the arrays of a block where there is now a
swap with the consumer's `Variant`, with the cargo tests made at
`read_variant`. The numbers of the bench below are what the parse
straight into a block will be compared with. Whether the branch is
merged as it is, before that change, is part of the first thing asked
of the owner.

What exists now that did not: a cargo workspace with the core crate,
the pyo3 binding and the Python package, the wasm-bindgen binding and
the TypeScript package, and the build of the wheel for pyodide; in the
core, the variant record and the reader trait, the VCF reader, plain and
gzipped, parsed in batches on the threads of rayon natively and on one
thread in wasm, and the collector of blocks; `popnei.open_vcf` and
`openVcf`, which give a `Variants` whose `iter_blocks` gives the
genotypes block by block. The reader reads the reference VCFs as
bcftools 1.24 does, 500 variants of 50 individuals among them, and the
blocks of a file, put one after another, equal the chunks that pyNei
gives for it, put one after another, in every column. 92 cargo
tests, 38 pytest tests, 39 node tests and a smoke test under pyodide.

What is asked of the owner:

1. The merge of `plan/vcf-to-blocks` into `main`, which is theirs to
   order.
2. The speed. On the 400 MB VCF of 100000 variants x 1000 individuals,
   this Mac, release: 1.24 s on one thread against a target of 0.55 s,
   and 0.160 s on 18 cores against 0.11 s. The spike, the trial parser
   in Rust of `docs/rust_core.md` that the targets come from, gives
   0.54 s and 0.098 s on the same file the same day, and pyNei takes
   13.5 s. A profile puts 95 in 100 of the one thread time in the parse
   of the columns of the individuals. If nothing is done the reader
   stays 2.3 times slower than the spike on one thread and 1.6 on 18.
   `.claude/skills/performance-review/` is how it would be taken up,
   with `crates/popnei/benches/read_vcf.rs`, the bench that this plan
   added. The gzipped target of the spec, 0.55 s for a gzipped file of
   53 MB, cannot be compared with anything: this file gzips to 38 MB,
   and the spike itself takes 0.84 s on it. That row of the spec and of
   `docs/rust_core.md` needs another number or a note. The details are
   under work package 5.
3. Four decisions that the work went on without, each with what the
   code does meanwhile:
   - A bgzipped file cut at the boundary of a gzip member gives fewer
     variants and no error, where bcftools says "no BGZF EOF marker".
     Recommended: detect it.
   - A QUAL of `nan`, `inf` or `1e400` is taken as a float, as bcftools
     does, and `nan` is then the same as no quality. Four reviewers
     reported it. Recommended: refuse a quality that is not finite.
   - The specs have one error enum for the crate and the `coding` skill
     asks for one error type for each operation. The code follows the
     specs. Recommended: keep the enum and correct the skill.
   - pyNei is a development dependency by the absolute path
     `/Users/jose/devel/pynei`, because the `../pynei` of the objectives
     points nowhere from a worktree. Recommended: pyNei's repository on
     GitHub at a commit.
4. Three things the orchestrator decided that a user sees, which the
   owner can reverse: a ploidy above 255 is refused; the TypeScript
   block has `numIndividuals` and `ploidy`; `many.vcf` has ids and
   qualities, with its genotypes and the counts of the spec unchanged.
5. Smaller: no license is chosen, and no manifest names one. The
   `coding` skill has no reference for wasm-bindgen, and
   `js/popnei/README.md` has what was learned, to write one from.

What the next plan and the skills should take from this one. Every
review found things that were wrong and that no check had seen, several
of them by running the case: a session killed by an argument, a license
nobody chose, a build that plain cargo could not make, a test that
guarded nothing. The reviews cost more than the writing: 2.3 million
tokens for 20 reviewers against 1.9 million for the subagents that wrote
the 13 tasks and 0.6 million for those that fixed the findings, each
figure added from the usage that every run of a subagent reports. Two
pairs of tasks that the plan marked as side by side were
run one after the other, because both of each pair wrote the workspace
manifest and `Cargo.lock`: the next plan should mark tasks as side by
side only when they share no manifest. The plan's check `cargo test -p
popnei -- --list` ends with the count of the doc tests, 0, and needs
`--lib`. Twice this report repeated what a subagent said without
checking it, a manifest with no license and a spec that had two
`Default`, and a reviewer found each: what a subagent says of a file is
checked before it is written here. And a read only reviewer checked out
a commit in the plan's own worktree and left it on a detached HEAD,
which cost nothing but has to be told to every reviewer that reads
there; the reviewers of the last work package each got a worktree of
their own.

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

## Work package 4: blocks, and the Python and TypeScript sides

Four tasks, four subagents, one after another and not side by side as
the plan allowed for 4.2 and 4.3, for the reason given in work package
1: both write the workspace manifest and `Cargo.lock`. Task 4.1, the
collector in the core, 7b04ef9 for the spec and df0618f, 184 thousand
tokens and 11 minutes. Task 4.2, the Python side, c3409ba, 195 thousand
and 14 minutes. Task 4.3, the TypeScript side, d84ce87, 181 thousand and
12 minutes. Task 4.4, the smoke test under pyodide, 315d406, by the
subagent of task 1.4, resumed, 40 thousand more and 5 minutes. It
finished as planned. pyNei is a development dependency by the absolute
path `/Users/jose/devel/pynei`, the plan's meanwhile.

### The deliverables, run by the orchestrator at 315d406 and again after the fixes

1. `cargo test -p popnei --lib block:: -- --list` `14 tests`, 19 after
   the fixes, where the plan asks for 6; `cargo test --workspace`
   `79 passed`, then 84; fmt, clippy, `cargo wasm-check` clean.
2. `uv run maturin develop && uv run pytest` `30 passed`, then 38, 8 of
   them the comparison of popnei's blocks, joined, with pyNei's chunks,
   joined, on `cases.vcf`, `cases.vcf.gz`, `many.vcf` and `many.vcf.gz`,
   with blocks of 7 variants and of the default size; ruff clean.
3. `npm run build` with no TypeScript error and `npm test` `tests 27`,
   then `tests 39`, `fail 0`.
4. The wheel of pyodide, 147 KB where it was 40 KB before numpy, and
   `node tests/pyodide/smoke.mjs`, which reads `cases.vcf` and
   `cases.vcf.gz` under pyodide 314.0.7 with its numpy 2.4.6 and finds
   the genotypes of the spec's table, 3 variants by default and 4 with
   every variant, exit 0, in 1.2 s once micropip and numpy are cached.

### The review

Seven reviewers over c0a9fa7..315d406: `spec`, `tests`, `numbers`,
`errors`, `api`, `architecture` and `binding`; 147, 159, 134, 112, 99,
145 and 123 thousand tokens. The findings that held went back to the
three subagents that wrote each side, the core first: 6 commits from
2e8e31a, 79 thousand tokens more; the Python side, 13 commits from
04a1618, 129 thousand; the TypeScript side, 4 commits from 989d1d9, 97
thousand.

What held and was fixed, the ones that matter first:

- A `num_vars_per_block` that a user gives killed the session.
  `iter_blocks(num_vars_per_block=10**13)` on `cases.vcf` printed "memory
  allocation of 80000000000000 bytes failed" and the interpreter exited
  with 134; in wasm 4e9 genotypes trapped and left the module unusable.
  The collector checked only that the size fits a `usize`, and a failed
  allocation is an abort that no `except` catches. Four reviewers ran
  it. Every column is now reserved with `try_reserve_exact` before
  anything is read, and the same call is a `ValueError` at the first
  block.
- After the error of a reader the collector went on reading. It held
  only because `VcfReader` stops itself: a reader that failed once and
  then gave more variants got a block after its error. The vars file
  reader could have been that reader.
- TypeScript took numbers that Python refuses, because wasm-bindgen
  turns a JavaScript number into a `usize` by truncation and modulo
  2^32: `{ploidy: 2.5}` read with 2, `numVarsPerBlock: 4294967297` gave
  blocks of 1 variant, and -1 an error about 4294967295. The package
  checks every whole number before the call.
- `Block.gts` was read only through its view alone: `gts.base` was
  writeable over the same buffer, and the test tried only direct
  assignment. The genotypes are built with their three dimensions in the
  binding, still with no copy, and every array of a block is read only.
- In the browser the VCF was held twice in the memory of wasm: 185.0 MB
  after `openVcf` for a VCF of 91.9 MB, 93.1 MB now.
- The Python binding copied the whole table of chromosomes for every
  block, asked for or not. 20000 variants of 20 individuals with 20000
  contigs, blocks of 100, no field asked for: 0.279 s, 0.176 s now,
  against 0.166 s for the same file with one contig.
- Which name of the API is which field of the reader was written twice,
  the same 24 lines in the two binding crates. It is in the core,
  `needs_of_the_fields`, and both call it.
- A variant with another number of alleles than individuals x ploidy is
  an error in the collector. A reader that broke that promise would have
  been stopped by numpy in Python and read wrong with no sign in
  TypeScript, where the genotypes cross flat.
- Smaller: a bare string for `fields` was split into letters in both
  languages, and is a `TypeError` in Python and an `Error` in
  TypeScript; a negative ploidy was an `OverflowError` and is the
  `ValueError` the spec promises; an error of the input lost its errno
  and its path, so that `except IsADirectoryError` missed it;
  `repr(block)` printed 300 KB; the Python classes reported `builtins`
  as their module; a position above 2^53 was rounded on its way to
  JavaScript with no error; the count of open passes that the TypeScript
  tests lean on never came back to 0 for an iterator that was never
  started; a trap inside wasm was hidden by the error of the `free()`
  that followed it; `Variants` had no `[Symbol.dispose]`.
- Tests that were missing: the default fields of a block in TypeScript,
  which could be all five with every test passing, and that the arrays
  of a block are copies that survive the growth of the memory of wasm.

Two things the orchestrator decided, which change what a user sees and
which the owner can reverse:

- The TypeScript block has `numIndividuals` and `ploidy`. Its `gts` is
  flat, and without them a block could not be indexed without the
  `Variants` it came from. Three reviewers and the subagent of 4.3 asked
  for it. The Python block needs neither, the shape of `gts` carries
  them, and `docs/specs/block.md` says so.
- `many.vcf` has ids and qualities, worked out from the index of the
  variant with no random draw, so that no genotype changed and the
  sixteen counts of the spec stand. Before, the comparison with pyNei of
  those two columns rested on the 4 variants of `cases.vcf`: an id upper
  cased in the core failed 4 of the 8 tests, and fails the 8 now. The
  file is 117346 bytes and its gzip members hold 617, 65252, 51477 and 0
  bytes of text, which the spec has.

What held and was not changed: the id column of a block is a
`Vec<String>`, one allocation per variant, which section 2 of the
architecture dictates; a reviewer measured no cost, 0.284 s with every
field and with the genotypes alone for 20000 variants x 1000
individuals in release. `maturin develop` builds the core with the dev
profile, so pytest runs an unoptimised core, 7.4 s through `iter_blocks`
for a file the release core reads in 0.285 s; the measurement of work
package 5 uses a release build.

### How the work went

The read only reviewers were sent to read in the plan's own worktree,
and one of them checked out the commit under review there, which left
the tree on a detached HEAD. The 19 commits of fixes that followed went
onto that detached line, and the orchestrator's tick of task 4.4 onto
the branch. The subagent of the Python fixes saw it and said so. Nothing
was lost: the tick was cherry-picked onto the line and the branch moved
to it, 78 commits, clean. From now on a reviewer that reads in the shared
worktree is told to switch no commit, and a fixer to check that it is on
the branch before it commits. `main` was not touched.

### For the owner

Besides the two decisions above, the three of work package 3 are still
open. Two more reviewers, independently, reported the second one here: a
QUAL of `nan` in the file is stored as the NaN that means that the
variant has no quality, and the two cannot be told apart.

What was learned about wasm-bindgen is in `js/popnei/README.md`, for the
reference that the `coding` skill lacks: what crosses as a copy, that
names are not turned to camelCase and a `pub const` cannot be exported,
that a number is truncated at the boundary, that the `finally` of a
generator does not run when it was never started, and what an argument
of bytes costs.

## Work package 5: the parallel reader, and its speed

Task 5.1, the batches, commits a24a4ee, a refactor that changes no
result, and 41812a8; one subagent run of 227 thousand tokens and 20
minutes. Task 5.2, the measurement, commit 159f758, another subagent,
180 thousand tokens and 14 minutes. The work package finished as
planned, and both targets of the spec that can be compared are missed.

### The deliverables

1. The reader reads 1024 lines, or 8 MiB of text if that comes first,
   into slots that it owns, parses them on the threads of the rayon pool
   it is called in, and hands the variants out in the order of the file
   by swapping buffers with the lent `Variant`. A chromosome gets its
   number when its variant is handed out. In wasm the same function
   parses one line after another, and neither wasm build has rayon.
   Check, run by the orchestrator: every test of work packages 3 and 4
   passes untouched; `cargo test --workspace` `87 passed` at 41812a8 and
   `92 passed` after the fixes, with the new tests of the same variants
   and chromosome numbers in pools of 1 and of 4 threads, compared with
   the stored output of bcftools too, and of a wrong line in the fourth
   batch of 64 lines, which gives the good variants first and then its
   error; pytest `38 passed`, `npm test` `tests 39`, `fail 0`, the wheel
   of pyodide built and its smoke test passed.
2. The check by hand of the allocations, with a counting allocator, a
   second pass over `many.vcf` with every field asked for: 2 allocations,
   8 bytes, in the last 400 variants, both the interning of `chr2` at
   variant 251, and the same with the genotypes alone. A reviewer
   confirmed the numbers and followed the allocations over 30 batches of
   lines of changing lengths: new allocations stop after the third batch
   and reallocations after the twenty third, all of them at the first
   read of a batch, the buffers of the slots still growing, and none
   between two variants of a batch.
3. The measurement, below. The bench is
   `crates/popnei/benches/read_vcf.rs`, run with `cargo bench --bench
   read_vcf -- <path> --threads 18 --runs 5`, and `make_big_vcf.py`
   beside it makes the file in 3 s.

### The measurement

The file: `simulate_genotypes` and `write_vcf` of
`test/gwas_reference/make_reference.py` of pyNei, which
`docs/rust_core.md` names as the origin of the panel, with seed 42, 1000
individuals and 100000 variants, 3 in 100 of the genotypes missing, `.`
in every FILTER. 403 MB, and 38 MB with bgzip. The original file cannot
be made again bit for bit; this one has its shape. Apple M5 Pro, 18
cores, release build, the file in the page cache, the genotypes asked
for, the default options, the median of 5 runs; VS Code was open and the
load average was 1.8. The spike is `spike/pynei_spike` of pyNei, built
from a copy outside that repository, on the same file and the same day.

| | popnei | the spike | the target, and a tenth above it | met |
|---|---|---|---|---|
| plain, 1 thread | 1.24 s | 0.54 s | 0.55 s, 0.605 s | no |
| plain, 18 threads | 0.160 s | 0.098 s | 0.11 s, 0.121 s | no |
| bgzipped, 1 thread | 1.58 s | 0.84 s | 0.55 s, 0.605 s | no |
| bgzipped, 18 threads | 0.50 s | 0.40 s | none | |

plink2 v2.0.0-a.7.7, `--vcf <file> --make-pgen --threads 1`: 0.273 s,
the mean of 5 runs. It is the 0.27 s of `docs/rust_core.md`, and the
spike gives its own 0.55 s and 0.11 s, so this file is as hard as the
one those numbers came from, and the difference is popnei's: 2.3 times
the spike on one thread and 1.6 on 18. Against pyNei's 13.5 s it is 11
times faster on one thread and 84 on 18. The orchestrator ran the bench
again: 1.25 to 1.31 s on 1 thread and 0.159 to 0.161 s on 18, and after
the fixes of the review 1.27 s and 0.158 s.

Where the time goes on one thread, from a sampling profile of 15107
samples: 95.5 in 100 in the parse and 4.5 in reading the bytes and
cutting the lines. Of the self time, the filling of the genotypes has
45.7 in 100, the split of a text at a character 25.2, the search of a
byte 15.6 and the comparison of strings 10.0. The nine first columns,
the FILTER and the count of the alleles are under 1 in 100 together. The
spike makes three checks less than the reader, the ploidy of each
genotype, its allele numbers against ALT, and the FILTER, and they are
not where the 0.7 s went.
The columns of the 1000 individuals are.

On 18 threads the batch is a barrier: no line of the next batch is read
while this one is parsed and handed out. A reviewer measured, on 5000
variants of 1000 individuals, 10.8 ms for the read of the lines alone on
one thread and 4.1 ms on eight, against 92.5 ms and 13.8 ms with the
genotypes, and 4.1 + (92.5 - 10.8) / 8 = 14.3 ms predicts the 13.8 ms: the
serial read is the floor, and the read ahead thread of section 3 of the
architecture is what would remove it.

The size of a batch, which the spec leaves to a measurement: 256, 1024
and 4096 lines give 0.206, 0.160 and 0.147 s on 18 threads and 1.28,
1.24 and 1.28 s on one. 1024 stays: 4096 gains 8 in 100 on 18 threads,
loses 3 on one, and no size reaches the target.

The gzipped target of the spec does not reproduce with anything. The
spec and `docs/rust_core.md` give 0.55 s for a gzipped file of 53 MB.
This file bgzips to 38 MB, 37 MB with plain gzip, and the spike itself
takes 0.84 s on it. That row is of a file that nobody can make again.

### The review

Three reviewers over a24a4ee, 41812a8 and 159f758, each in a worktree of
its own: `spec`, `tests`, and `architecture` with the errors; 155, 124
and 112 thousand tokens. The findings that held went back to the
subagent of 5.1, 11 commits from 9cb1d93, 93 thousand tokens more.

- A panic in a worker of rayon left the reader going on with lines
  dropped and no error: a reviewer injected one at line 1501 of 3000,
  caught it, and read 2996 variants, the positions 1501 to 1504 lost. No
  input can make the parse panic today, but pyo3 turns a panic into an
  exception that an `except Exception` swallows. A reader whose parse did
  not come back now gives an error at every later read.
- The batch was bounded in lines and not in bytes, so the memory of a
  reader grew with the individuals: 13.3 MB of resident memory for 1000
  individuals, measured, where the doc comment of the constant said
  6 MB, 82.7 MB for 10000, and near 0.8 GB for 100000. A batch
  is now 8 MiB of text at most, which leaves the 1024 lines of the
  400 MB file as they were, and the bench says it costs nothing.
- Three uses of rayon, in two tests and in the bench, were not behind
  the `cfg` of the wasm targets: `cargo wasm-check` passed because it
  checked the library alone. The alias checks all targets now, and both
  wasm targets pass.
- The test named for the reuse of the buffers guarded nothing since the
  batches: throwing every buffer away at the swap passed all 87 tests.
  The test is now the reviewer's probe, 3 distinct buffers for 8 lines
  in batches of 2, and 8 with that mutation.
- No test ran the serial parse natively, none noticed a batch of one
  line, the errors that the batching moved were tested with one batch
  alone, and the bench could not vary the size of a batch, which the
  subagent of 5.2 had changed by hand. All have tests or arguments now.
  The doc of the constant said that nobody had measured it, in the
  commit that measured it.

What held and was not changed: the message of 41812a8 says the first
batch of a second pass costs 3396 allocations and 403 KB, and a reviewer
counted 3458 and 376 KB with its own harness; the two numbers that the
spec leans on, above, reproduced exactly.

## The final check

From a clean clone of the branch at 59bf624, outside the repository:
`cargo fmt --all --check` exit 0; `cargo clippy --workspace
--all-targets -- -D warnings` no warning; `cargo test --workspace` `92
passed`; `cargo wasm-check`, both wasm targets, every target of the
crate, finished; `uv run ruff format --check` `8 files already
formatted`, `uv run ruff check` `All checks passed!`; `uv run maturin
develop && uv run pytest` `38 passed`; `npm run build` and `npm test` in
`js/popnei` `tests 39`, `pass 39`, `fail 0`;
`scripts/build_pyodide_wheel.sh` built the wheel and
`node tests/pyodide/smoke.mjs`, after `npm install` beside it, read
`cases.vcf` and `cases.vcf.gz` under pyodide and exited with 0.
