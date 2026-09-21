# Work report: readers that give blocks

The plan `docs/plans/block-readers.md` is done, on the branch
`plan/block-readers`, in the worktree `.claude/worktrees/block-readers`,
21 September 2026. Its three work packages finished as planned, each was
reviewed and its findings fixed. The owner then read this report and
answered what it left him, and a fourth work package carried out his
answers the same day; it is at the end of this report, before the final
check, and the list of what is asked of him below is what is left after
it. The final check passes from a clean clone. The four speed targets of the VCF reader spec are met, one of them
at its edge. The orchestrator, in this report, is the session of the
assistant that ran the plan: it sent each task to a subagent, checked
what came back and had it reviewed.

What exists now that did not. The variants flow in blocks, runs of
consecutive variants held as arrays, from the file to whoever reads them,
and nothing is left of the single variant that a reader filled, of its
trait or of the collector that copied such variants into blocks: a search
of the code for their names finds nothing, where it found 98 lines. In
the core there is one trait for everything that gives blocks, the view of
one variant of a block, a block that compacts itself in place and checks
its own arrays, and `reblock`, the reader over a reader that cuts and
joins blocks to a size. The VCF reader parses its lines as bytes straight
into the rows of a block, on the threads natively and one after another
in wasm, and both bindings hold it as a boxed reader with a `reblock` at
the end of every pass. It reads the 403 MB VCF of 100000 variants x 1000
individuals, on the owner's Apple M5 Pro, release, in 0.56 to 0.62 s on
one thread and 0.093 s on 18 cores, where the reader before this plan
took 1.24 s and 0.160 s and the target was 0.594 s and 0.108 s; the
bgzipped file in 0.89 s and 0.394 s, against 1.58 s and 0.50 s and
targets of 0.924 s and 0.44 s. Two files that were read before are
refused: one with a quality that is not finite, and a bgzipped one that
was cut short, wherever it was cut. The repository is under the MIT
license, takes pyNei from its GitHub repository at ef0ca6e, and its
`coding` skill describes the one error enum that the code has. A
bgzipped file is read by the size that each of its gzip members states,
each member checked, so that a file that is corrupted is an error and not
an empty file; in Python a file that was cut short or is corrupted is an
`OSError` with its name, a defect of popnei a `RuntimeError`, and every
error of a file names the file. 173 cargo tests, 72 pytest tests, 8 of
them against pyNei, 46 node tests and a smoke test under pyodide, where
the plan started with 92, 38 and 39.

For a user nothing else changed: `open_vcf`, `openVcf`, `Variants` and
`iter_blocks` are what they were, and no test or file of the two packages
that was there changed, but the smoke test of pyodide, to which one check
was added.
What did change for them is under "For the owner" of work package 2, with
the six decisions that the work went on without.

What is asked of the owner:

1. Nothing stops the merge, which he ordered on 21 September 2026,
   without the timings of work package 4.
2. The speed of the reader as the last fixes left it, which nobody has
   timed: another job had the machine. The bgzipped file is the one that
   changed, and the one set taken before those fixes says it is no
   slower, 0.378 s on 18 threads against 0.394 s. It is for the
   performance reviews to come, with the two things that would make the
   reader faster: the read ahead thread, and the decompression of the
   members on several threads, which the reader of the members is now
   split for.
3. One option that work package 4 leaves: a subclass of `OSError` for a
   file that was cut or corrupted, so that a pipeline that wants to fetch
   the file again catches it by name and not by `errno` being `None`.
4. The commits that changed the specs, which the owner has not read:
   `docs/specs/block.md` in 11de065 and f9dc27d, `docs/specs/io_vcf.md`
   in 9e20eb7, 93b0710, e297e63, d37f5b1 and b00be26. Most correct a
   sentence that said what the code does not do; the last two are his
   decisions of 21 September.
5. The answers he gave are recorded at the head of work package 4. "For
   the owner" of work package 2 is as it was written before them.

What the next plan and the skills should take from this one is under
"How the work went": the reviews of the first three work packages cost 2.2
million tokens against 1.3 million for the writing, and found a loop with no end, a copy that grew
with the square of the cuts, a panic on Ctrl-C and variants read after an
error, none of which a check had seen; a task's prompt should say "test
first, and see it fail" and ask for the report of how each test failed;
and one subagent can carry three tasks of one file with profit.

## Before the first task

The owner approved the plan in chat on 20 September 2026, and it was
committed on `main` as a11aea9, which the branch starts from. Everything
of "What has to be in place" was checked by running it, in the worktree:

- The three specs, the architecture and the glossary are on `main` since
  8f19650 and cae6a4c, and `git status --short` was empty. None of the
  specs has an open point.
- `cargo fmt --all --check` exit 0; `cargo clippy --workspace
  --all-targets -- -D warnings` no warning; `cargo test --workspace` `92
  passed`; `cargo wasm-check` finished; `uv run ruff format --check` `8
  files already formatted`; `uv run ruff check` `All checks passed!`;
  `uv run maturin develop && uv run pytest` `38 passed`; `npm run build`
  and `npm test` in `js/popnei` `tests 39`, `pass 39`, `fail 0`;
  `bash scripts/build_pyodide_wheel.sh` built the wheel and `node
  tests/pyodide/smoke.mjs` exited with 0.
- `which` finds cargo, rustc, uv, node, npm, wasm-bindgen, bcftools, bgzip
  and plink2; `~/devel/emsdk/emsdk_env.sh` and
  `~/devel/pyodide-venv/bin/pyodide` are there; both wasm targets are
  installed.
- `git ls-remote https://github.com/JoseBlanca/pynei main` gives ef0ca6e.
- The checks that say the work is done fail, as they must: 98 mentions of
  `read_variant`, `VariantReader` or `BlockCollector` in seven files of
  `crates`, none of `BlockReader`, `VariantRef` or `reblock`, no `LICENSE`,
  pyNei at `/Users/jose/devel/pynei` in `pyproject.toml`. The second gzip
  member of `tests/reference/vcf/many.vcf.gz` ends at byte 12336 with 280
  variants before it, recomputed from the file.
- The file of the bench was made again with
  `crates/popnei/benches/make_big_vcf.py`, 403572954 bytes, and its
  `bgzip -k`, 37695707 bytes, in `/Users/jose/devel/popnei-bench/`, a
  directory that the orchestrator made outside the repository and that
  the owner can delete when the plan is merged.

The reader as built, timed again before any change, `cargo bench --bench
read_vcf -- big.vcf --threads N --runs 5`, the owner's Apple M5 Pro,
release, the file in the page cache, the genotypes asked for: a median of
1.307 s on one thread and of 0.161 s on 18, where the report of the plan
before this one has 1.24 s and 0.160 s. The one thread number of today is
the one that the new reader is compared with.

What the orchestrator changed in the plan before it was approved, none of
which changes what is built: the opening names a third change that a user
sees, a column that was not asked for is no longer checked, which "How it
runs" of `docs/specs/io_vcf.md` decides; task 1.1 corrects section 8 of
`docs/architecture.md` too, which has the path dependency on pyNei as
`docs/objectives.md` has; tasks 2.2 and 2.3 run one after the other,
because `memchr` is not a dependency today and 2.3 writes the manifest of
the core crate and `Cargo.lock`, and 2.6 runs after 2.5; tasks 2.4 and
2.5 name "The Rust interface" of the VCF reader spec.

## Work package 1: the five decisions

Tasks 1.1 and 1.2 went to one subagent in one prompt, since both are
small: commits 8ca8935 and a16ea06, 113 thousand tokens and 7 minutes.

### The deliverables, run by the orchestrator at a16ea06

1. The license. `LICENSE` at the root starts with `MIT License` and
   `Copyright (c) 2026 Jose Blanca`; `cargo metadata --no-deps
   --format-version 1` gives `MIT` for `popnei`, `popnei-python` and
   `popnei-js`, which take it from the workspace manifest;
   `pyproject.toml` has `license = "MIT"` and `license-files`, and
   `js/popnei/package.json` `"license": "MIT"`; the wheel of pyodide holds
   `popnei-0.1.0.dist-info/licenses/LICENSE`, and `npm pack --dry-run` in
   `js/popnei` lists `LICENSE` among its 21 files. npm packs the files of
   the directory of the package and leaves a symbolic link out, so a
   `prepack` script copies the `LICENSE` of the root beside
   `package.json` before every pack, and `.gitignore` holds that copy:
   git has one `LICENSE`.
2. pyNei by its repository. `pyproject.toml` has pyNei as a git source of
   uv at ef0ca6e177be5a18c2dba1cc78901940ad847cf9; `grep -n "Users/jose"
   pyproject.toml uv.lock` finds nothing; `uv sync`, `uv run maturin
   develop` and `uv run pytest -k pynei` give `8 passed, 30 deselected`,
   and the whole of pytest `38 passed`; `grep -rn "path dependency"
   docs/objectives.md docs/architecture.md` finds nothing.
3. The `coding` skill. "Errors, and no panics" says that the core crate
   has one error enum, `non_exhaustive`, to which each module adds its
   cases, with the owner's two reasons, and how each binding crate turns
   it into what its language throws, in `errors.rs` of each. The
   orchestrator read the paragraph against the code: the two files, the
   enums `PyPopneiError` and `JsPopneiError`, and the cases it names are
   there.

### What was changed in the plan, and what the tasks did beyond it

`.claude/skills/code-review/categories.md` asked a reviewer of errors for
one error type per operation, so a reviewer would have reported the single
enum as a defect; the orchestrator added its correction to task 1.2.
The subagent found two more places that said what the tasks correct, and
corrected them: "Tests" of the `coding` skill, which had pyNei as a path
dependency, and `pyo3.md` beside it, which described a newtype with a
`From` for each error type of the core where the binding crate has an
enum with one case for the one error of the core.

### The review

One reviewer, of the category `spec`, in a worktree of its own at a16ea06,
63 thousand tokens. A review of popnei sends one reviewer for each
category of the `code-review` skill, each of which starts with nothing but
the code and the spec: `spec` reads the code against what it was written
from and runs the cases, `tests` breaks the code to see whether a test
fails, `numbers` reads the arithmetic, `errors` looks for panics and reads
the messages, `api` the names and the doc comments, `architecture` the
allocations and the threads, `binding` the two binding crates. Here one
was enough: the work package has no code to break, so `tests` had nothing to
run, and the numbers its commits claim are few and were checked by the
orchestrator above. It confirmed by running them: the `LICENSE` is the
MIT text of SPDX, whole; the native wheel and the source distribution
that maturin builds, the wheel of pyodide and the tarball of npm from a
fresh checkout all carry it; `uv sync` in a tree that never saw the
owner's path installs pyNei from GitHub at ef0ca6e, and pytest gives `38
passed` there.

Two findings held, both about the new paragraph of the `coding` skill
against the code, and the subagent that wrote it corrected them in
b63bd9b, 10 thousand tokens more:

- It said that no call site of the Python binding crate maps an error by
  hand. `crates/popnei-python/src/vcf.rs` has four, one of them only
  because the functions that pyo3 exports return pyo3's own result. The
  skill now states the rule, every function of a binding crate returns
  the result of that crate and a call site maps an error only to add the
  path of the file, and task 2.2 makes the code follow it.
- It named two exceptions of Python where the binding crate raises three:
  a `RuntimeError` for a defect of the binding. It names the three now,
  and the four cases of each enum.

The orchestrator had read that paragraph against the code and had not
seen either. `docs/rust_core.md` also had pyNei as a path dependency, in
item 2 of its list of decisions, and has the sentence of the architecture
now: `grep -rn "path dependency\|sibling checkout" docs .claude/skills
.claude/agents`, the plans and the reports aside, finds nothing.

### For the owner

- Three files that tell a reviewer where pyNei is still name the checkout
  on this machine, `/Users/jose/devel/pynei`: the `spec` section of
  `.claude/skills/code-review/categories.md`,
  `.claude/agents/spec-reviewer.md` and
  `.claude/skills/performance-review/profiling_environment.md`. That
  checkout is at ef0ca6e today, the commit that the tests now pin, so no
  number is wrong. When the checkout moves on, a reviewer will compare
  with another pyNei than the tests. A reviewer also reads pyNei's tests
  and its spike there, which the installed package does not have, so the
  fix is not a change of path alone. Not changed: it is outside this
  plan.
- `npm pack` in `js/popnei` from a fresh checkout makes a tarball of
  three files, with no code, unless `npm run build` was run first: the
  package has no script that builds before a pack. Publishing is outside
  this plan.

## Work package 2: the VCF reader gives blocks

Notes taken task by task; the deliverables, the review and what the owner
should know are written when the work package ends.

Task 2.1, the trait of the readers, what a block gains, the view of one
variant and `reblock`, in the core: commit 1c9f2e4, one subagent, 207
thousand tokens and 15 minutes. Run by the orchestrator: `cargo test
--workspace` `107 passed`, 15 of them new, pytest `38 passed`, both wasm
targets checked. Three of the six cases of the error that the spec lists
were there from the collector and fit; three are new. What it did that
the task or the spec did not say, which the reviewers of the work package
were told to look at: `retain_vars` runs `check` before it moves a row,
because a block with a defect would make it index outside an array;
`reblock` skips a block of no variants that its source gives, where the
spec says only that a reader gives none; the collector got a `set_needs`
of four lines so that the reader that stands on it until work package 3
can keep the contract of the trait; and the subagent wrote the code before
the tests, against the `coding` skill, and then broke the code in three
places to see the tests fail, which showed a defect of its first version,
a block of no variants given after a cut.

Task 2.2, both binding crates on a boxed reader of blocks: commit
319d435, one subagent, 185 thousand tokens and 10 minutes. Run by the
orchestrator: the checks of the `coding` skill, `107 passed`, `38
passed`, `tests 39`, `fail 0`; `git diff --stat 1c9f2e4 -- python
js/popnei/src tests js/popnei/test` prints nothing, so no test and no
file of the two packages changed; `grep -rn
"VariantReader\|BlockCollector" crates/popnei-python crates/popnei-js`
finds nothing; no manifest changed. Every function of the Python binding
crate returns the result of the crate, as the corrected skill asks, and
the three `map_err` that are left on an error of the core add the path of
the file. Two cases of the error that say a reader has a defect are a
`RuntimeError` in Python. No test sees that choice: neither binding crate
can build cargo tests, and no reader that Python reaches gives those
errors.

Task 2.3, the parser of one data line, as bytes, into one row of a block:
commit 71b0d9d, one subagent, 202 thousand tokens and 15 minutes. Run by
the orchestrator: the checks of the `coding` skill, `132 passed`, 25 of
them new, `38 passed`, `tests 39`, `fail 0`, both wasm targets checked
with `memchr`, which is the one new dependency, pure Rust. Nothing calls
the parser yet but its tests, one of which parses every data line of the
three reference VCFs with the new parser and with the old one and
compares them; it goes with the old parser in work package 3. The
subagent timed the parse alone, the 100000 lines of the file of the bench
held in memory, the genotypes asked for, one thread, release, the median
of 3 runs: 0.576 s with the new parser and 1.229 s with the old one. The
whole read adds the reading of the lines to that, and its target is
0.55 s, or 0.605 s with the tenth that the spec allows.

### The review of tasks 2.1 and 2.2

Tasks 2.1 and 2.2 are one piece of code, the trait with `reblock` and the
two bindings that hold it, so they were reviewed together while task 2.3
was written: seven reviewers at 319d435, one for each category, each in a
worktree of its own, 1.0 million tokens together (`spec` 154 thousand,
`tests` 151, `numbers` 157, `errors` 140, `api` 133, `architecture` 142,
`binding` 131). No reviewer found a wrong result. What held went back to
the two subagents that wrote the code, test first, each new test seen to
fail: the core in six commits from 11de065 to 9755872, 129 thousand
tokens, and the bindings in five from 4e6b0f2 to f6bed12, 79 thousand.
After them, run by the orchestrator: `cargo test --workspace` `143
passed`, pytest `50 passed`, `npm test` `tests 40`, `fail 0`, both wasm
targets checked, `cargo doc -p popnei --no-deps` with no warning.

What was found that mattered, and is fixed:

- A source that gives a block of no variants made `reblock` ask it again
  for ever: two reviewers ran it, one stopped after 30 million calls. From
  Python that loop would run with the interpreter released and the reader
  locked, where Ctrl-C does not reach. `reblock` now refuses such a block,
  with a new case of the error, the seventh of the block spec.
- A cut of `reblock` copied everything after the cut, so the next cut
  copied it again: one block of the 500 variants of `many.vcf` cut into
  blocks of 1 allocated 14.6 MB where no cut allocates 0.6 MB, and each
  block given kept the capacity of the whole source block, 50000
  genotypes for the 10000 it held. Nothing cuts today; the vars file
  reader will. A cut now copies the rows that leave, once, into a block
  of their size.
- A Ctrl-C while the first block of a process was read gave a
  `PanicException`: the signal stayed pending until numpy's C interface
  was first fetched, which panics on it. It is older than this plan. The
  signals are now checked when the read returns, and
  `tests/test_interrupt.py` interrupts a read in a subprocess; it failed
  with that panic before the fix and passed 20 of 20 runs for the
  subagent and 5 of 5 for the orchestrator.
- A pass of a binding that failed on its own account, after the reader
  gave it a block, read on at the next call: at the class that the wasm
  package exports, a VCF of 6 variants with a position above 2^53 in its
  third, in blocks of 2, gave the first block, the error, and then the
  third block, the variant at position 400 lost without a word. Both
  passes now end at any error, with a node test at that class.
- `num_vars_per_block=-1` was answered with "0 or more, and at most
  18446744073709551615", 0 being refused too, and 2^64 with Python's own
  `OverflowError`, which names no argument. Four reviewers reported it; it
  is older than this plan. Both arguments are now converted by the binding
  crate, and `tests/test_counts.py` has five values for each.
- Three reviewers found that no test cut a block that carried ids,
  qualities or alleles, and that every allele of every fixture is one
  byte, so that an offset counted in alleles where bytes were meant gave
  the same answer. The code was right, which two of them showed with
  alleles of 1 to 8 bytes; the tests now use the rows of
  `differences.vcf`, and the test over `many.vcf` asks for every field and
  compares with the file read one variant at a time, where it compared
  the four sizes with each other.
- Eleven guards and choices that could be taken out with every test
  passing have a test each; the case of the spec in which `reblock` loses
  what it kept, 28 blocks of 7 and then the error, has its test; the doc
  comment of the trait, which the vars file reader and the filters will be
  written from, says what `reblock` refuses at run time; the constructors
  of the alleles column are visible inside the crate, which task 2.4
  needs; the six places that go through the columns of a block one by one
  stop compiling when a column is added and forgotten.

`docs/specs/block.md` changed in 11de065, in commits of its own before
the code: what `reblock` does with a block of no variants; what a cut
copies; that `retain_vars` runs `check` before it moves a row; that the
views of a block that does not pass `check` stop early, so a consumer that
did not get its block from `reblock` or a binding calls `check` first; and
seven cases of the error where there were six.

Findings not taken:

- Taking the `reblock`, or the check of a block, out of either binding
  fails no test, because the source gives blocks of the size the user
  asked for and `reblock` has nothing to cut or join. The fix proposed,
  a source that always gives the default size, would copy every block of
  every pass with a size of its own for the sake of the test. The first
  filter, which leaves blocks of uneven sizes, will exercise it.
- The ids of a block are a `String` each, one allocation per variant when
  they are asked for, 672 more over the 500 variants of `many.vcf`. The
  spec fixes that field; task 2.6 says whether it costs.
- numpy's `into_pyarray` has no form that gives an error, so a failed
  allocation there is a panic; the loops over variants inside the binding
  crates and a Python pass with no `close` are older than the plan.

### Tasks 2.4 to 2.6

Task 2.4, `VcfReader` as a reader of blocks on the parser of task 2.3,
went to the subagent that wrote that parser: commits 9e20eb7, a sentence
of the spec, 83349df, the reader, both bindings and the bench, 2bf494d
and 8a887a5; 241 thousand tokens and 28 minutes. Run by the orchestrator:
the checks of the `coding` skill, `147 passed`, pytest `50 passed`, `npm
test` `tests 42`, `fail 0`, both wasm targets checked; between f6bed12 and
8a887a5 no file of `python/`, `js/popnei/src/`, `tests/` or
`js/popnei/test/` changed but one new node test; `grep -rn
"VariantReader\|BlockCollector\|CollectedBlocks" crates/popnei-python
crates/popnei-js` finds nothing. `read_variant` and the old parser of a
line are gone with this task, since the reader could not be both at
once; 53 mentions of the record level are left in `block.rs` and
`variant.rs` for task 3.1. No pytest or node test leaned on the old rule
that a position is checked when it was not asked for.

The bench on the new reader, run by the orchestrator right after the
task, the same command and file as the baseline above, load average 1.3:
the plain file in a median of 0.650, 0.598 and 0.606 s on one thread in
three runs of 5, and of 0.120 s on 18; the bgzipped one in 0.881 s and
0.431 s. The subagent had 0.586, 0.123, 0.874 and 0.434 s. Task 2.6 is
the measurement that says which targets are met.

`docs/specs/io_vcf.md` changed in 9e20eb7: a size of the blocks that the
caller asked for is refused when the reader is built, and the default
size, which nobody asked for, when the first block is built. The reason
is the finding of the review of the bindings: in wasm a header of 170000
individuals opened with the ploidy 255 was refused at `openVcf`, because
the default block of 100 variants holds more genotypes than 32 bits
count, although that file is read in blocks of 10. It is a corner that a
user sees, decided by the subagent and kept by the orchestrator because
it follows from "nothing is read when it is called but the header" of
section 5 of the architecture; the owner can reverse it.

Task 2.5, the two new errors, went to the same subagent: commit b67f652,
48 thousand tokens and 8 minutes. Run by the orchestrator: `153 passed`,
pytest `59 passed`, `npm test` `tests 44`, `fail 0`, both wasm targets
checked. A quality of `1e39` is refused like `1e400`, since a quality of
a block is a float of 32 bits and reads both as infinite. What a user
sees of a file that was cut: the error comes after the variants, and
through `iter_blocks` the variants that `reblock` was keeping for its
next block are lost with it, as the block spec says of any error, so
`many.vcf.gz` cut after its second gzip member gives a Python user, in
blocks of 100, 200 of its 280 variants and then the error, where the
reader itself gives the 280. The message says that the file is cut short
and has to be fetched again.

The order of the rest was changed, and the plan says so: the review of
tasks 2.3 to 2.5 and task 3.1, the removal, run at the same time, in
other trees and other files; then the fixes; then task 2.6, the
measurement, with no build beside it and on the code that the plan
leaves.

## Work package 3: the record level goes

Task 3.1 ran while the reviewers of tasks 2.3 to 2.5 read the VCF reader,
which is in another file, and went to the subagent that wrote the block
module: commits c5a5041, one test written over the VCF reader before the
test that it replaces went, and 56f5251, the removal; 74 thousand tokens
and 10 minutes.

### The deliverables, run by the orchestrator at 56f5251

1. No record level. `grep -rn
   "read_variant\|VariantReader\|BlockCollector\|CollectedBlocks" crates
   --include='*.rs'` finds nothing, where the plan started with 98 lines,
   and `"BlockReader\|VariantRef\|reblock"` gives 137, where it gave 0.
   `cargo doc -p popnei --no-deps` has no warning, and the pages of the
   `variant` module are `ChromTable`, `Needs`, `VariantRef` and the two
   constants: no `Variant`. One case of the error went,
   `VariantOfAnotherSize`; what it guarded holds of a block now, and
   `Block::check` finds it. Neither binding crate named a case that went.
2. Everything still passes: `cargo test --workspace` `143 passed`, 154
   before the removal and 11 tests gone with what they called; pytest `59
   passed`; `npm test` `tests 44`, `fail 0`; both wasm targets checked;
   the wheel of pyodide built and its smoke test exited with 0, by the
   subagent, and again in the final check.

The eleven tests that went are in the second table at the end of this
report, each with what it checked and the test that checks it now. One
rule would have lost its only test, a block refused for the memory of its
columns and not only of its genotypes, and got a test over the VCF reader
first. Two things have no test any more and need none: the reader that
stood on the collector, and the capacity of the buffers of a variant that
was filled again and again, whose rule the spec of the VCF reader leaves
to a check by hand.

The case of the error for a consumer that did not get a field it depends
on was called `FieldsNotFilled`, after the variant that a reader filled.
It is `FieldsNotInTheBlock` now. Nothing gives it yet, and no binding
crate or test outside `error.rs` names it; it is a public name of the
core, so the owner may want another.

### The review of tasks 2.3 to 2.5

Seven reviewers at b67f652, one for each category, each in a worktree of
its own, 1.1 million tokens together (`spec` 155 thousand, `tests` 181,
`numbers` 123, `errors` 181, `api` 128, `architecture` 155, `binding`
184). No reviewer found a wrong genotype, position, id, allele or
quality. What they ran: the reader against bcftools 1.24 and pyNei on
VCFs of their own, with ploidies of 1, 3 and 4, 127 and 128 alleles,
every form of a genotype, `GT` at any place of the FORMAT, 1500 variants
equal in the three programs; every read with the same bytes for five
sizes of block, six bounds of a batch and 1 or 4 threads; the 32
combinations of the fields with the values of all the fields; 640 passes
from Python against pyNei, 4 files by 32 sets of fields by 5 sizes, with
no mismatch, and the same blocks from node; 210 thousand mutated inputs,
a file cut at every length, every byte changed and removed, random bytes
after a header, on 1 to 8 threads, with no panic and no hang; a Ctrl-C
at three moments of a pass of 2 million variants on a release build,
always a `KeyboardInterrupt`; the memory of wasm the same after 100
passes. The sixteen counts of the spec, recomputed from the stored output
of bcftools, are exact.

What held went back to the subagent that wrote the reader, test first,
six commits from 93b0710 to 2ed860b, 102 thousand tokens, and three small
ones of the bindings to the subagent that wrote those, b385e77 to aaa7464,
37 thousand tokens. After them, run by the orchestrator: `153 passed`,
pytest `61 passed`, `npm test` `tests 45`, `fail 0`, both wasm targets
checked, the wheel of pyodide built and its smoke test exited with 0. What was found
that mattered, and is fixed:

- A bgzipped file that was cut got the new error only when the cut fell
  where a gzip member ends: of the 21904 ways to cut `many.vcf.gz` short,
  3 gave it and 21899 gave "incomplete deflate stream" or "unexpected end
  of file", in Python an `OSError` with no file name in it. A bgzipped
  source that ends early is now that one error wherever it was cut, a
  `ValueError` in Python, with the variants before the cut given first by
  the reader; an error of the file system stays the error of the input,
  and has a test with a source that fails after some bytes, which nothing
  drove before.
- Eight rules of the reader could be broken with its 101 tests passing,
  and each has a test now that was seen to fail: no block after an error,
  which every test had checked with the wrong line at the end of its
  file, where the reader stops for lack of lines; the text of a batch
  cleared for the next one, which is what the bound of 8 MiB is for, and
  without which the reader held every line of the file; both halves of
  what tells a bgzip file from a gzip one; the window of the last 28
  compressed bytes; two messages that were asserted by their place alone;
  and an allele of exactly 127, which only refusals tested.
- A default size of the blocks that does not fit, which only wasm meets,
  told a user who gave no size to ask for fewer variants in a block. It
  has words of its own, and the check of a size is one function that
  `reblock` and the VCF reader both call, where it was written in three
  places that had drifted apart.
- The bench took any argument it did not know as the path of the file, so
  `--lines-per-batch=4096` ran with the default and said nothing. It
  refuses what it does not know, and has a `--help`.
- A `ploidy` of -1 was told that up to 18446744073709551615 is allowed;
  an error of the input with no number from the system lost the name of
  its file on the way to Python; and the header of 170000 individuals that
  has to open where a size is 32 bits was tested under node and not under
  pyodide, which is such a build.

`docs/specs/io_vcf.md` changed in 93b0710 and `docs/specs/block.md` in
f9dc27d, in commits of their own. Six sentences said what the code does
not do, and the code was right: that the reader as built parses every
position, which was true of the old one; that `bcftools view -f .,PASS`
makes the same choice as popnei, which it does not for a FILTER that
names `PASS` beside a filter that failed, `PASS;q10`, which bcftools
keeps and popnei skips; that a line whose bytes are not text is an error
of the line, where the nine first columns are checked and a byte that is
not text in the column of an individual is "not an allele number" with
the name of that individual; that the positions and the qualities are
written from the threads, where only the genotypes are, the rest being
appended after each batch, which adds 1 ms to the 127 ms that the read of
the 403 MB file takes on 18 threads; that two
bytes of the source are looked at, where it is sixteen. The seventh is
the decision of the first item above, which extends the owner's decision
on the mark of the end to a cut inside a member.

Findings not taken, each because its fix is the owner's:

- A bgzip file whose second member has the two bytes of the length of its
  extra field changed to one precise value makes the decoder land on the
  last, empty member: the mark of the end is there, and the file gives no
  variant and no error. One reviewer found it, the only silent case among
  the 101745 changes of one byte of `cases.vcf.gz`; the reader before this
  plan had it too. The fix is to read a bgzip file by the size that each
  of its members states, as bcftools does, where the spec names flate2's
  decoder of several members.
- Through `iter_blocks` a file that was cut gives fewer variants than the
  reader does before the error, because `reblock` loses what it was
  keeping, as the block spec says of any error: `many.vcf.gz` without its
  last 28 bytes gives 500 variants in blocks of 1 and of 100, 497 in
  blocks of 7, and none with the default size, the whole file being one
  short block that was waiting.
- No error of a data line names the file, only the line and the column.
- A parse that did not come back is a `ValueError`, as the VCF reader spec
  says, where the `coding` skill keeps `RuntimeError` for a defect of
  popnei; and neither binding reaches that error: a panic of the parse is
  a `PanicException` in Python.
- A position written `+5` is read as 5, as pyNei reads it; the spec says
  that an allele number is a run of digits and says nothing of a
  position.

### The measurement, task 2.6

A fresh subagent, 222 thousand tokens and 32 minutes, with nothing else
building on the machine; commits 3c4bf5d, 52e506f and 5cbf389, one
constant each, and e297e63, the table of "Speed" of the spec. Its whole
account, with every command as it was run, the sweeps, the memory and the
two profiles, is `docs/reports/block-readers-measurement.md`. What the
owner needs of it is here.

The file is the VCF of 100000 variants x 1000 individuals of the bench,
403 MB plain and 38 MB bgzipped; the owner's Apple M5 Pro, 18 cores,
release, the file in the page cache, the genotypes alone asked for, the
default options, `cargo bench --bench read_vcf -- <file> --threads N
--runs 5`, the median of 5 runs, three sets of each and fourteen of the
bgzipped read on one thread, on 21 September 2026, with load averages of
1.4 to 2.4: `mediaanalysisd`, a background process of macOS, ran
throughout. The reader
before this plan and the spike of 20 September are from
`docs/reports/vcf-to-blocks.md`; the target is that spike's time and a
tenth more, as "Speed" of the spec has it.

| | before this plan | the spike, 20 September | the target | the reader now, and its sets | met |
|---|---|---|---|---|---|
| plain, 1 thread | 1.24 s | 0.54 s | 0.594 s | 0.563 s, 0.560 to 0.563 | yes, see below |
| plain, 18 threads | 0.160 s | 0.098 s | 0.108 s | 0.093 s, 0.093 to 0.094 | yes |
| bgzipped, 1 thread | 1.58 s | 0.84 s | 0.924 s | 0.890 s, 0.859 to 0.940 | yes |
| bgzipped, 18 threads | 0.50 s | 0.40 s | 0.44 s | 0.394 s, 0.393 to 0.396 | yes |

The orchestrator ran the bench again after the subagent, with a load
average of 1.7, right after a run of all the checks: 0.093 s and 0.397 s
on 18 threads and 0.875 s bgzipped on one, which agree, and the plain
file on one thread in four medians of 0.622, 0.569, 0.614 and 0.614 s,
which do not: one of the four is where the subagent's three were, and
three are 3 to 5 in 100 above the target. The subagent saw the same two
clusters in the bgzipped read on one thread, near 0.866 and near 0.913 s,
with no load average to tell them apart. So on one thread the reader is
at the target or a twentieth above it, depending on something of the
machine that nobody has named; which core a single thread lands on is a
guess that was not tested. Against the reader before this plan, which the
orchestrator timed at 1.307 s on the same machine the day before, it is
2.1 to 2.3 times faster on one thread, and 1.7 times on 18.

The spike was built again from a copy outside both repositories and
timed the same day: 0.523, 0.103, 0.806 and 0.392 s. On 18 threads the
reader is now faster than the spike on the plain file, 0.093 s against
0.103, and the same on the bgzipped one; on one thread it is 8 in 100
behind it, and it makes three checks that the spike does not: the ploidy
of each genotype, its allele numbers against ALT, and the FILTER.

The three constants that the spec leaves to a measurement, each changed
in a commit of its own with its numbers, plain file, 18 threads: the
buffer of a file that is opened by its path, 8 KiB to 256 KiB, a new
named constant, 0.118 s to 0.106, flat beyond; the lines of a batch, 1024
to 4096, 0.105 s to 0.098, and 0.139 s with 256; the text of a batch, 8
MiB to 16 MiB, 0.098 s to 0.094, and 0.411 s to 0.392 bgzipped. None
costs anything on one thread. Together they cost memory: the most bytes
alive at once in a reader, by a counting allocator, went from 16.6 MB to
35.7 MB for 1000 individuals and from 22.5 MB to 34.6 MB for 10000, a
file of 3000 variants and 120 MB written for this; the maximum resident
size of the bench went from 20.0 to 33.5 MB for the 1000 individuals and
from 24.4 to 32.8 MB for the 10000. The
memory still does not grow with the individuals, which is what the bound
in bytes is for. The last of the three buys 4 in 100 on 18 threads for
about 12 MB, and is the one the owner may want back at 8 MiB. None of it
reaches wasm, which reads one line at a time from bytes in memory.

Where the time goes, from `/usr/bin/sample` over 30 s: on one thread
91.6 in 100 of the 23038 samples are in the columns of the individuals,
6.1 in the read of the lines, 1.3 in the nine first columns. On 18
threads the read of the lines, which is serial, is 4308 of the 9363
samples of the thread that reads, and the 18 workers wait in 61 in 100 of
theirs: it is still the floor, as the report of the plan before this one
found, and the read ahead thread of section 3 of the architecture is what
removes it. The targets are met without it.

### The deliverables of work package 2, run by the orchestrator at e297e63

1. The trait, what a block gains, the view of one variant and `reblock`:
   `cargo test -p popnei --lib block -- --list` gives `53 tests`, among
   them the ones "How it is verified" of the block spec lists, made at
   `retain_vars`, `variants` and `check`, the two of `reblock` over a
   reader written in the test, and
   `a_boxed_reader_of_blocks_is_read_through_reblock`; all pass.
2. The row parser: `cargo test -p popnei --lib io::vcf -- --list` gives
   `110 tests`, where the file had 65 before the plan. The first table at
   the end of this report has each of the 65 with the test that took its
   place; the orchestrator checked with a script that every old test has
   a row and that every test the table names is in the code, and the
   `tests` reviewer read the rows against the two versions of the file.
3. `VcfReader` is a reader of blocks: the tests on `cases.vcf`,
   `differences.vcf`, `many.vcf` and their gzipped forms against the
   stored output of bcftools and the sixteen counts of the spec, the sizes
   of the blocks, the pools of 1 and of 4 threads, the wrong line, and the
   serial parse against the parallel one are among those 110, and pass.
4. The two new errors: cargo tests, `tests/test_truncated_and_infinite.py`
   and `js/popnei/test/cut_short.test.ts`; both are a `ValueError`.
5. Both bindings hold a boxed reader of blocks: pytest `61 passed`, 38
   before the plan; `npm test` `tests 45`, `fail 0`, 39 before; between
   a11aea9 and e297e63 no file of `python/`, `js/popnei/src/`, `tests/` or
   `js/popnei/test/` that was there changed, but the smoke test of
   pyodide, to which one check was added, and seven files of tests are
   new; the wheel
   of pyodide builds and its smoke test exits with 0; a search of the two
   binding crates for `VariantReader` and `BlockCollector` finds nothing.
6. The measurement: above.

### How the work went

Every task was done at its first sending, and no subagent had to be
replaced. The tokens, each figure from the usage that a run of a subagent
reports: 1.3 million for the eight runs that wrote the tasks (tasks 1.1
and 1.2 together 113 thousand, 2.1 207, 2.2 185, 2.3 202, 2.4 241, 2.5
48, 3.1 74, 2.6 222), 0.36 million for the five that fixed what the
reviews found, and 2.2 million for 15 reviewers. The reviews cost more
than the writing again, and again found what no check had seen: a loop
with no end, a copy that grows with the square of the cuts, a panic on
Ctrl-C, variants read after an error, a truncated file told of a deflate
stream, 19 rules that no test held to.

The subagent of task 2.1 wrote the code before the tests, against the
`coding` skill, and said so. Its review found eleven guards that no test
held to, and the fixes went back with "test first and seen to fail" in
the message, which every later task and fix then did and reported. The
prompt of a task should say it from the start.

Tasks 2.3, 2.4 and 2.5, and the fixes of their review, went to one
subagent, which kept what it had learned of the parser: 593 thousand
tokens at its end, and no sign that it had lost the start.

What the lessons of the plan before this one gave: no two subagents wrote
in the tree at once, every reviewer had a worktree of its own, and the
branch was never left; the orchestrator committed the plan and the report
by the name of their files while a subagent worked, with no collision.
Twice the orchestrator's own check of a subagent's text missed what a
reviewer then found, the paragraph of the `coding` skill and the
allocation counts of task 2.4, which a reviewer could not reproduce, 347
allocations for a pass over `many.vcf` with the genotypes alone against
the 292 that the subagent reported; those counts are not in this report, and the ones that are,
one allocation for each block and none for each variant with the
genotypes alone, are the `architecture` reviewer's own.

### For the owner

The decisions that the work went on without, each with what the code does
meanwhile and what the orchestrator recommends:

1. A bgzip file with the length of the extra field of one member
   corrupted to one precise value reads as a file with no variants and
   no error. The options: leave it, since it takes that precise
   corruption and the reader before this plan had it too; or read a bgzip
   file by the size that each member states, as bcftools does, which
   replaces the decoder that the spec names, about a day of work, and is
   also the way to decompress on several threads, where the bgzipped file
   now takes 0.394 s on 18 threads against 0.093 s plain. Recommended:
   the second, as a plan of its own with the read ahead thread.
2. Through `iter_blocks` a file that was cut gives fewer variants before
   its error than the reader does, and none with the default size of the
   blocks when the file is shorter than a block, because `reblock` loses
   what it was keeping, which the block spec decides for any error. The
   options: leave it, since the error always comes and a file that was
   cut is of no use; or have `reblock` give what it kept as a last short
   block before it passes the error on, for every error, which reverses
   that rule of the block spec. Recommended: leave it, and say in "The
   cases a reader of the rules would not guess" of the VCF reader spec
   what a user of `iter_blocks` gets.
3. The errors of a data line name the line and the column and not the
   file. Recommended: the binding crates put the path in front of the
   message of every error of a file, since they have it.
4. A parse that did not come back is a `ValueError`, by the VCF reader
   spec, and the `coding` skill keeps `RuntimeError` for a defect of
   popnei. Neither binding reaches it today: a panic of the parse is a
   `PanicException` in Python. Recommended: `RuntimeError`, one line of
   the spec and one arm of the Python binding crate.
5. A quality that is missing is NaN inside the core, by the block spec,
   and "Floats" of the `coding` skill asks for an `Option` inside the
   core and NaN at the boundary with Python alone. The doc comment of the
   view says that a caller tests for NaN before any arithmetic.
   Recommended: decide it with the first calculation that reads the
   qualities.
6. The constructors of the alleles column of a block are visible inside
   the crate and not outside it, so a reader of blocks written outside
   popnei cannot fill the alleles. Recommended: leave it until somebody
   writes one.

What the orchestrator and the subagents decided that a user sees, which
the owner can reverse: a bgzipped file cut inside a gzip member is the
error of the mark of the end, a `ValueError`, where it was an `OSError`
about a deflate stream; the default size of the blocks is checked when
the first block is built and not when the file is opened; a reader that
gives a block of no variants is a new case of the error,
`ReaderGaveABlockOfNoVariants`, a `RuntimeError` in Python with the two
other defects of a reader; `FieldsNotFilled` is `FieldsNotInTheBlock`; a
quality of `1e39` is refused like `1e400`; a size or a ploidy that Python
cannot fit in 64 bits is a `ValueError` that names the argument, where it
was Python's `OverflowError`; the three constants of the measurement,
above. And what changed in the specs, always in commits of their own:
`docs/specs/block.md` in 11de065 and f9dc27d, `docs/specs/io_vcf.md` in
9e20eb7, 93b0710 and e297e63.

Smaller, and outside this plan: the three files that send a reviewer to
the checkout of pyNei on this machine, and `npm pack` that packs no code
without a build before it, both under work package 1; the spec and the
message of the cut file call a gzip member of a bgzip file a "block",
which is bgzip's word and also the glossary's word for the variants; a
position written `+5` is read as 5, as pyNei reads it; the ids of a block
are a `String` each, one allocation for each variant when they are asked
for, which the measurement did not reach because the bench asks for the
genotypes alone; a pass of Python has no `close`, so a pass that is
abandoned keeps its file open until the collector runs; and the
directory `/Users/jose/devel/popnei-bench/`, with the files of the bench,
the builds of the sweeps and a clean clone, which can be deleted after
the merge.

## Work package 4: the owner's decisions of 21 September 2026

The owner read this report on 21 September 2026 and answered in chat what
it left him, with a rule, an error never passes silently, and a
convention for the exceptions of Python: a `ValueError` is a wrong input
of a function, a `RuntimeError` a defect of popnei, an `OSError` a file
that cannot be read, that was cut short or that is corrupted. He asked
for the work of those decisions before the merge, and the orchestrator
added this work package to the plan, in fa063d0. His answers, by the
numbers of "For the owner" of work package 2: the corrupted bgzip file is
an error, now; a file that was cut gives its error in the iteration as
soon as the cut is found. The orchestrator read that as the rule of the
block spec kept, that when an error comes `reblock`, the reader that
puts the blocks to a size at the end of every `iter_blocks`, loses the
block it had not yet filled and does not give it first, and told him so,
since it is what makes a short file that was cut give no variant at all
with the default size of the blocks; every error of a file names the file;
a parse that did not come back is a `RuntimeError`; a missing quality is
NaN inside the core too; the constructors of the alleles column stay
inside the crate. Of what was decided without him: a bgzipped file that
was cut is an `OSError` and not a `ValueError`; the rest stays, the
16 MiB of a batch among it.

Task 4.1, the reader of a bgzip file by the sizes of its members, went to
a fresh subagent: commits d37f5b1, the spec, and aed0273, the code; 306
thousand tokens and 31 minutes. A file that bgzip wrote is a row of gzip
members, each with at most 64 KB of text, and an empty one of 28 fixed
bytes that marks the end. The header of each member has an extra field,
and in it a subfield called `BC` that holds the size of the member; at
its end a member has the CRC32 of its text, a checksum, and the length of
that text. The decoder that popnei used goes from one member to the next
on its own and never reads that size. Such a file is now cut into its
members by the size that the header of each states, and each
member is decompressed on its own, with flate2's raw deflate, and checked:
that its data ends where the member does, its CRC32, and the length of
its text. The cut and the decompression are two steps, which is what
decompressing on several threads would stand on; that is not built. A
gzip file that bgzip did not write keeps the decoder it had. What watched
the last 28 compressed bytes go by is gone: a source is whole when its
last member holds no text. No dependency was added. Run by the
orchestrator: `cargo test --workspace` `161 passed`, `1 ignored`, pytest
`61 passed`, `npm test` `tests 45`, `fail 0`, both wasm targets checked.

The file of the review, `many.vcf.gz` with its bytes 320 and 321 changed
from `06 00` to `44 54`, read by the orchestrator from Python: before
this task no variant and no error, and now the error that names the
member, where it starts and what is wrong with it, "the size it states,
12026 bytes, leaves no room for its data after the 12 bytes of its
header, the 21572 of its extra field", and says that the file has to be
fetched again. The report had said that reading by the sizes of the
members is what bcftools does. The orchestrator ran `bcftools view -H` of
bcftools 1.24 on that file: 0 variants, no message, exit 0. popnei is
stricter than bcftools here, and the spec says so. pyNei gives
`ValueError: Empty VCF file, it has no variants`.

The test that the owner's rule rests on changes every byte of
`cases.vcf.gz` to each of its 255 other values, 101745 files, and asks of
each either an error or exactly the variants of the whole file: 96562
errors, 5183 read as the whole file, none silent, where the reader before
had one. It runs in 1.2 s and is a test like any other. The same over the
first two members of `many.vcf.gz` and its end mark, 3.15 million files,
none silent, takes 14 s in release on 18 threads and is a test that is
run on request. Both counts are the subagent's; the orchestrator ran the
first as part of `cargo test`.

Task 4.2, the exceptions of Python and the name of the file, went to the
subagent that wrote the bindings: commits 35527ab, 842ea22, 3464932 and
2bc1ae5, 74 thousand tokens and 11 minutes. Run by the orchestrator:
`161 passed`, pytest `65 passed`, `npm test` `tests 46`, `fail 0`. A
bgzipped file that was cut short or that is corrupted is an `OSError`
with the path in `filename` and no number; a parse that did not come
back is a `RuntimeError`; the message of a `ValueError` and of a
`RuntimeError` of a file starts with its path, "/tmp/x/wrong_line.vcf:
line 3 of the VCF, the column POS: `x` is not a position", and an
argument that is refused names no file. Python prints an `OSError` that
has no number as "[Errno None] message: 'path'", which is its own format.
JavaScript reads bytes in memory and has no path to give.

### The review

Seven reviewers at 2bc1ae5, one for each category, each in a worktree of
its own, 0.96 million tokens together (`spec` 133 thousand, `tests` 208,
`numbers` 114, `errors` 164, `api` 108, `architecture` 116, `binding`
114). What held went back to the two subagents, test first: the core in
eight commits from b00be26 to 3ea216a, 147 thousand tokens, and the
bindings in seven from 884f8a6 to fe79f32, 73 thousand. After them, run
by the orchestrator: `cargo test --workspace` `173 passed`, `1 ignored`,
pytest `72 passed`, `npm test` `tests 46`, `fail 0`, `cargo wasm-check`
with no warning, `cargo doc -p popnei --no-deps` with no warning.

What was found that mattered, and is fixed:

- One file was still read with other variants than it holds and no
  error, and four of the seven reviewers found it, each from its own
  side. What said that bgzip wrote a file was its `BC` subfield at the
  bytes 12 and 13, where bgzip puts it; BGZF lets another subfield come
  first, and such a file was taken for a plain gzip, with no check of its
  members and no check of its end: `many.vcf.gz` so changed gave no
  variant with the length of the review, 220 with another, and the
  variants up to the cut when it was cut where a member ends, always with
  no error. htslib asks for `BC` at byte 12 too, so bgzip never writes
  such a file and bcftools reads it the same silent way. The `BC` of the
  first member is now found wherever it is.
- Eight of the twelve checks of a member could be taken out with all 161
  tests passing, the one that catches the file of the review among them,
  because another check then refused that file; with three of them gone a
  corrupted member was read with no error. The sweep of every byte
  catches the CRC32 check and no other: a corruption whose text still
  comes out right is, by its rule, the variants of the whole file. So it
  proves that no change of one byte gives other variants in silence, and
  not that each check is reached; its doc comment says so now, and each
  check has a test that asserts the message it gives. Two checks were
  missing and were added: a subfield after `BC` that runs past the extra
  field, and a text longer than a member holds.
- After a corrupted member the reader of the members handed out, at its
  next call, the text whose CRC32 had failed. The VCF reader stops at its
  first error, so no user saw it; the module is written to be built on.
  It gives no text after an error.
- The new module of the core that reads the members, `io::bgzf`, keeps
  the cut of a member and its decompression as two steps, so that some
  later plan can decompress several members at once on several threads.
  The two steps could not be run apart: they shared one buffer, and
  a second cut before the first decompression gave the text of member 2
  under the CRC32 of member 1. The decompression is a function that
  borrows what it needs and no reader, and a test runs it on two threads
  at once.
- Both messages of a bgzipped file ended with "The variants before it
  were given", which is not so for a user of `iter_blocks`: with the
  default size of the blocks a file that was cut or corrupted gives a
  Python or node user no variant at all, the whole file being one short
  block that `reblock` was keeping, 0 of the 456 that the reader gives of
  `many.vcf.gz` cut at 21000 bytes. The owner decided that `reblock`
  keeps its rule, so the sentence went out of the messages, the spec says
  "the reader gives", the docstrings of `open_vcf` and `iter_blocks` say
  what a user gets, and a pytest test holds the four counts of the spec,
  500, 497, 500 and 0 variants before the error with blocks of 1, 7, 100
  and the default, where only the 100 had one.
- A bgzipped file cut inside its first member was told that its header
  has no `#CHROM` line, a `ValueError`; it is the error of a file that
  was cut short. Bytes after the mark of the end were told that the file
  is cut short; they are a corrupted file. The messages say "member"
  where they said "block", and the glossary has the word.
- In Python: a ploidy out of range named the file, though it is an
  argument that the user typed; a path whose bytes are not UTF-8 came
  back mangled in `filename`, which a reviewer showed with a file that is
  not there, since this file system refuses to make one; the number of
  the system was said twice, "[Errno 2] ... (os error 2)"; a filter that
  gives a mask of another length than its block, a defect of popnei, was
  a `ValueError`, which the orchestrator had asked for before the owner
  gave his convention; two statements of the skills were false, nine
  cases of the VCF reader where there are ten, and that every message has
  the path in front.
- `cargo wasm-check` exited with 0 and warned, for both wasm targets, of
  the imports of a test that is not compiled there. The orchestrator saw
  it in the output of its own check. The command now denies warnings, as
  `cargo clippy` does.

`docs/specs/io_vcf.md` changed in d37f5b1 and b00be26, in commits of
their own, with `docs/glossary.md` and the table of the modules of
`docs/architecture.md`. Besides how a bgzip file is read and which
exception each case is, it now says where the owner's rule cannot be
kept: every member of a bgzip file is valid on its own and none records
its place, so a member that is removed, repeated or moved cannot be seen;
a reviewer ran it on `many.vcf.gz`, 220 variants, 780, and the same 500
in another order, with no error in popnei nor in bcftools.

Findings not taken, and what is left for the owner:

- A file that was cut or corrupted and a file that is not there are both
  an `OSError`, and a pipeline that wants to fetch a file again can tell
  them apart only by `errno` being `None`. A subclass of `OSError` that a
  user catches by name would settle it; pyNei defines no exception of its
  own. Not built: it is a name of the public API.
- The `RuntimeError` of a parse that did not come back has no test: no
  file gives it, and the binding crate is a module that the interpreter
  loads, which no test binary can link without changing how it is built.
- A plain gzip that was cut gives an error of the input and no variant,
  where a bgzip gives the variants first; nothing in the spec says which
  it should be.

### The deliverables of work package 4, run by the orchestrator at fe79f32

1. A bgzip file is read by the sizes of its members: the file of the
   review gives an error in a cargo test, in `tests/test_corrupted_bgzip.py`
   and in a node test, and from Python in the orchestrator's own hands;
   the sweep of every byte of `cases.vcf.gz` is one of the `173 passed`;
   the tests of the gzipped and bgzipped files that were there pass
   untouched; a gzip that bgzip did not write is read; `cargo wasm-check`
   passes with no warning.
2. The exceptions of Python: pytest `72 passed`, among them an `OSError`
   whose `filename` is the path for `many.vcf.gz` cut at 12336 bytes,
   inside a member and without its last 28 bytes, and for the corrupted
   file; a wrong data line a `ValueError` whose message starts with the
   path.
3. The documents agree with the code: the VCF reader spec, the glossary,
   the table of the modules of the architecture, and "Errors, and no
   panics" and "Floats" of the `coding` skill, which has the owner's
   convention and his rule in his words, and the quality as the one
   missing value that is NaN inside the core.
4. The speed: not measured, by the owner's order. While the fixes of the
   review were made another job of the owner's took ten of the 18 cores,
   with load averages of 8 to 13, and he chose to merge without the
   timings. What there is: the subagent of task 4.1, before the fixes of
   the review, on a machine with a load average of 2.4 to 2.8, the
   bgzipped file of the bench, release, the median of 5 runs, two sets:
   0.861 and 0.910 s on one thread, where task 2.6 had 0.890 s, and
   0.378 s twice on 18 threads, where it had 0.394 s. The new reader
   makes 6251 fewer allocations over that file than the decoder it
   replaced, one for each member, and a reviewer's profile of one thread
   puts 65 in 100 of the read in the parse, 31 in the inflate and 1 in the
   CRC32. Nobody has timed the reader as the fixes left it.

### How the work went, in work package 4

Both tasks were done at their first sending. The tokens: 380 thousand for
the two tasks, 220 thousand for the two rounds of fixes, 957 thousand for
the seven reviewers. The review found what the plan was reopened for, a
file read in silence, still there in another form, and four reviewers
found it without being told where to look: the owner's rule was in the
prompt of each, and each went looking for a file that breaks it.

The orchestrator told the owner that reading a bgzip file by the sizes of
its members is what bcftools does to the same effect, without having run
it. The subagent of task 4.1 ran it: bcftools reads the corrupted file as
no variants with no message. The orchestrator ran it too, corrected the
plan, and told the owner.

The report of one round of a subagent never reached the orchestrator,
only the notice that it had ended. The orchestrator read the branch
instead of the report, which is what it should do in any case.

A reviewer wrote a script into the scratch directory of the orchestrator
and wrote over a file of the same name there, a helper that was no longer
needed. The prompt of a reviewer names its own worktree as the only place
to write; it should name the scratch directories too.

## The final check

From a clean clone of the branch at fe79f32, the last commit that changes
code, outside the repository, on 21 September 2026: `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `173 passed`, `1 ignored`, the sweep
of 3.15 million files that is run on request; `cargo wasm-check`, both
wasm targets, every target of the crate, finished with no warning, which
it now denies; `uv sync`, which fetched pyNei from GitHub at ef0ca6e,
then `uv run ruff format --check` `13 files already formatted`, `uv run
ruff check` `All checks passed!`, `uv run maturin develop && uv run
pytest` `72 passed`; `npm run build` and `npm test` in `js/popnei` `tests
46`, `pass 46`, `fail 0`; `bash scripts/build_pyodide_wheel.sh` built the
wheel and `node tests/pyodide/smoke.mjs` exited with 0, the header of
170000 individuals among the lines it printed; the search of `crates` for
`read_variant`, `VariantReader` and `BlockCollector` gives 0 lines, where
the plan started with 98, and the one for `BlockReader`, `VariantRef` and
`reblock` gives 137, where it gave 0; `LICENSE` is there, nothing of
`/Users/jose` is in `pyproject.toml` or `uv.lock`, and `cargo doc -p
popnei --no-deps` has no warning. The same check had passed at e297e63,
the end of the first three work packages, with `153 passed`, `61 passed`
and `tests 45`. The commits after fe79f32 change this report and the plan
alone.

## The 65 tests that the VCF reader had before the plan, and what took their place

Task 2.4 of `docs/plans/block-readers.md` rewrote the VCF reader: it gives
blocks and no longer fills one `Variant` at a time. `crates/popnei/src/io/vcf.rs`
had 65 tests at commit a11aea9, made at `read_variant` and at
`VcfReader::new`, and had 95 when task 2.4 ended; it has 110 at the end of
the plan. This table has, for each of the 65, the
test that took its place, where it is made, and what it covers that the old
one did not. Two of the 65 went with no test of the same name, and the last
rows say which and why.

Where a test is made: `new` is `VcfReader::new`, `from_path` is
`VcfReader::from_path`, `next_block` is the reader through the trait
`BlockReader`, and `parse_row` is the parser of one data line into one row
of a block, which task 2.3 built and whose tests are made at a line with no
reader around it.

| the test at a11aea9 | the test now | made at |
|---|---|---|
| the_individuals_of_a_plain_and_of_a_gzipped_vcf_are_read | the same name | `from_path` |
| the_gzip_bytes_are_found_in_a_source_that_is_not_a_file | the same name | `new` |
| a_path_that_no_file_is_at_gives_an_error_that_carries_the_path | the same name | `from_path` |
| a_source_that_gives_one_byte_at_a_time_is_read_gzipped_and_plain | the same name | `new` and `next_block` |
| a_source_that_starts_with_neither_a_hash_nor_the_gzip_bytes_is_refused | the same name | `new` |
| a_header_with_no_format_column_is_refused | the same name | `new` |
| a_header_with_a_format_column_and_no_individual_is_refused | the same name | `new` |
| two_individuals_with_the_same_name_are_refused | the same name | `new` |
| a_header_that_ends_before_the_chrom_line_is_refused | the same name | `new` |
| an_empty_source_is_refused | the same name | `new` |
| a_gzipped_source_that_is_not_a_vcf_is_refused | the same name | `new` |
| a_gzipped_source_cut_in_the_middle_of_a_member_is_refused | the same name | `next_block` |
| a_last_line_with_no_end_of_line_is_read | the same name, and a_line_that_ends_in_an_end_of_line_is_read_as_one_that_does_not | `next_block`, `parse_row` |
| a_read_after_an_error_and_after_the_last_variant_gives_no_variant | a_block_after_an_error_and_after_the_last_block_is_no_block | `next_block` |
| the_four_variants_of_cases_vcf_are_read_when_every_variant_is_given | the same name, and the_four_lines_of_cases_vcf_are_parsed_into_their_rows | `next_block`, `parse_row` |
| the_default_leaves_out_the_variant_of_cases_vcf_that_failed_its_filter | the same name | `next_block` |
| the_leading_separators_and_the_dot_of_differences_vcf_are_read | the same name, and the_two_lines_of_differences_vcf_are_parsed_into_their_rows | `next_block`, `parse_row` |
| the_chromosomes_are_numbered_in_the_order_of_the_variants_that_are_given | the same name, now over blocks of 1, 2 and 100 variants, since the numbers are given serially after each batch | `next_block` |
| a_tetraploid_genotype_with_a_ploidy_of_two_is_refused | the same name, and a_genotype_of_another_ploidy_is_refused_with_its_individual | `next_block`, `parse_row` |
| a_haploid_genotype_with_a_ploidy_of_two_is_refused_after_the_variants_before_it | a_haploid_genotype_with_a_ploidy_of_two_is_refused_after_the_blocks_before_it | `next_block` |
| a_tetraploid_vcf_read_with_a_ploidy_of_four_is_read | the same name, and a_tetraploid_line_read_with_the_ploidy_four_gives_four_alleles | `next_block`, `parse_row` |
| an_allele_that_the_variant_does_not_declare_is_refused | the same name, and an_allele_the_variant_does_not_declare_is_refused_with_the_genotypes_alone | `next_block`, `parse_row` |
| an_allele_above_the_largest_one_popnei_holds_is_refused | the same name, and an_allele_above_the_largest_one_is_refused | `next_block`, `parse_row` |
| a_line_with_the_genotypes_of_two_individuals_under_a_header_of_three_is_refused | the same name, and a_line_with_another_number_of_columns_of_individuals_is_refused | `next_block`, `parse_row` |
| a_line_with_the_genotypes_of_four_individuals_under_a_header_of_three_is_refused | the same name, and a_line_with_another_number_of_columns_of_individuals_is_refused | `next_block`, `parse_row` |
| a_format_with_no_gt_is_refused | the same name, and a_format_with_no_gt_is_refused_in_a_row | `next_block`, `parse_row` |
| a_position_that_is_not_a_number_is_refused | the same name, and a_position_that_is_not_a_number_is_read_when_the_position_was_not_asked_for, which is the rule that changed | `next_block` |
| a_quality_that_is_not_a_number_is_refused | the same name, and a_quality_that_is_not_a_number_is_refused_only_when_the_quality_was_asked_for | `next_block`, `parse_row` |
| an_allele_number_is_a_run_of_digits_and_nothing_else | the same name | `next_block` |
| an_allele_number_that_no_i8_holds_is_refused_for_being_above_the_largest | the same name | `next_block` |
| a_gt_that_is_not_the_first_key_of_the_format_is_read | the same name, and a_gt_that_is_not_the_first_key_of_the_format_is_read_into_a_row | `next_block`, `parse_row` |
| a_vcf_whose_lines_end_in_a_carriage_return_is_read | the same name | `next_block` |
| an_empty_line_at_the_end_is_skipped | the same name | `next_block` |
| a_vcf_with_a_header_and_no_variant_gives_no_variant_and_no_error | a_vcf_with_a_header_and_no_variant_gives_no_block_and_no_error | `next_block` |
| an_allele_that_alt_declares_and_no_genotype_carries_is_read | the same name, and an_allele_that_no_genotype_carries_is_read | `next_block`, `parse_row` |
| a_line_that_failed_its_filter_is_skipped_before_its_genotypes_are_read | the same name | `next_block` |
| a_line_that_ends_before_its_filter_is_refused_with_the_default | the same name | `next_block` |
| a_line_that_failed_its_filter_is_not_read_before_its_filter_either | the same name | `next_block` |
| with_the_genotypes_alone_the_id_and_the_alleles_are_not_filled | with_the_genotypes_alone_a_block_has_no_column, and only_the_fields_that_were_asked_for_are_parsed | `next_block`, `parse_row` |
| with_the_id_and_the_alleles_alone_the_genotypes_are_not_filled | with_the_id_and_the_alleles_alone_a_block_has_those_two_columns_and_no_genotype, and only_the_fields_that_were_asked_for_are_parsed | `next_block`, `parse_row` |
| a_data_line_whose_bytes_are_not_text_is_refused_with_its_number | the same name, and a_line_whose_bytes_are_not_valid_utf8_is_refused | `next_block`, `parse_row` |
| a_header_line_whose_bytes_are_not_text_is_refused | the same name | `new` |
| an_alt_that_ends_in_a_comma_is_refused | the same name, and an_allele_with_no_letter_in_it_is_refused_in_ref_and_in_alt | `next_block`, `parse_row` |
| a_ref_with_no_letter_in_it_is_refused | the same name, and an_allele_with_no_letter_in_it_is_refused_in_ref_and_in_alt | `next_block`, `parse_row` |
| an_individual_with_no_name_is_refused | the same name | `new` |
| a_line_of_seven_columns_is_refused_with_the_genotypes_not_asked_for | the same name, and a_line_of_seven_columns_is_refused_when_the_genotypes_are_not_asked_for | `next_block`, `parse_row` |
| a_format_with_no_gt_is_refused_with_the_genotypes_not_asked_for | the same name, and a_format_of_dp_is_refused_when_the_genotypes_are_not_asked_for | `next_block`, `parse_row` |
| what_is_in_the_columns_of_the_individuals_is_not_read_without_the_genotypes | the same name, and the_columns_of_the_individuals_are_not_read_when_the_genotypes_are_not_asked_for | `next_block`, `parse_row` |
| a_reader_asked_for_the_genotypes_alone_empties_the_alleles_it_filled_before | a_reader_asked_for_the_genotypes_alone_gives_its_next_block_without_the_columns | `next_block` |
| the_buffers_of_a_variant_are_written_over_from_one_variant_to_the_next | none of that name: see below | |
| the_genotypes_of_a_variant_come_back_in_the_buffers_of_the_variants_before_it | none of that name: see below | |
| every_variant_of_many_vcf_is_read_as_bcftools_read_it | the same name | `next_block` |
| the_default_gives_the_variants_of_many_vcf_whose_filter_passed | the same name | `next_block` |
| the_first_genotypes_of_many_vcf_are_the_ones_of_the_spec | the same name | `next_block` |
| the_counts_of_many_vcf_with_every_variant_given_are_the_ones_of_the_spec | the same name | `next_block` |
| the_counts_of_many_vcf_with_the_default_are_the_ones_of_the_spec | the same name | `next_block` |
| many_vcf_gives_the_same_variants_in_a_pool_of_one_thread_and_in_one_of_four | the same name | `next_block` |
| a_reader_whose_parse_did_not_come_back_gives_an_error_and_no_variant | a_reader_whose_parse_did_not_come_back_gives_an_error_and_no_block | `next_block` |
| a_wrong_line_of_a_later_batch_comes_after_the_variants_that_were_read_before_it | a_wrong_line_of_a_later_batch_comes_after_the_blocks_that_were_read_before_it | `next_block` |
| the_lines_parsed_one_after_another_give_what_the_threads_give | the same name, over `parse_rows` and `parse_rows_one_by_one` and the 500 lines of many.vcf | the parse of a batch |
| a_batch_holds_more_than_one_line_where_there_are_threads | the same name | the constant `LINES_PER_BATCH` |
| a_bound_of_bytes_smaller_than_the_file_cuts_the_batches_and_changes_no_result | the same name, with the file in one block of 1000 variants | `next_block` |
| a_file_read_one_line_at_a_time_gives_what_it_gives_in_one_batch | the same name | `next_block` |
| a_ploidy_of_zero_and_one_above_the_largest_are_refused | the same name | `new` |
| the_largest_ploidy_is_read | the same name | `next_block` |

### The two that went

`the_buffers_of_a_variant_are_written_over_from_one_variant_to_the_next` and
`the_genotypes_of_a_variant_come_back_in_the_buffers_of_the_variants_before_it`
were about the `Variant` that the consumer owned and lent to the reader: that
its buffers were written over instead of allocated again, and that the
buffers of the variants of a batch came back to the reader through the swap
of every read. There is no such variant now. What they were guarding, that a
reader of a file of any length allocates nothing for each variant, is
guarded by `the_rows_of_a_batch_are_the_ones_the_next_batch_is_parsed_into`,
at `next_block`: `many.vcf` read in batches of 8 lines gives its 475
variants, the reader keeps 8 rows and no more, and it filled more than 60
batches into them.

### What else moved

Two tests of `crates/popnei/src/block.rs` were about the VCF reader through
the collector and are in `io::vcf` now, where the reader is:
`a_reader_with_no_variants_gives_no_block`, which is
`a_vcf_with_a_header_and_no_variant_gives_no_block_and_no_error`, and
`the_error_of_the_third_variant_comes_after_the_block_of_the_two_before_it`,
which keeps its name and now reads two wrong lines, so that it also shows
that the error is the one of the first of them. The other tests of
`block.rs` that read a reference VCF through `BlockCollector` are made at
the VCF reader itself, with the same assertions: the sizes of the blocks of
`many.vcf`, what they hold joined, the columns of `cases.vcf`, the blocks of
a tetraploid VCF and the two tests of `reblock` over a real file.

## The eleven tests that went with the record level, task 3.1

21 September 2026. Task 3.1 of `docs/plans/block-readers.md` took the
single variant that a reader filled, `Variant`, the `VariantReader` trait,
the `BlockCollector` that copied its variants into blocks and the
`CollectedBlocks` that gave those blocks as a `BlockReader`, out of the
core crate. Eleven tests called them and went with them: 154 cargo tests
before the removal, 143 after.

This table is what a reviewer of the tests reads against. For each test
that went: what it checked, and the test that checks it now, by its name
and its file, or why nothing of it is left to check. The commit before the
removal, "the memory of the columns of a block is refused over the reader
that allocates", wrote the one test that had to exist first, because the
rule it checks still holds and no other test had it.

The names of files below are relative to the root of the repository.

| The test that went | What it checked | What checks it now |
|---|---|---|
| `block::tests::a_collector_of_blocks_of_no_variant_is_refused` | A reader given a size of 0 is refused where it is built, with `BlockOfNoVariants`, whose message names the 0. | `block::tests::a_reblock_of_blocks_of_no_variant_is_refused` for `Reblock::new`, and `io::vcf::tests::a_reader_of_blocks_of_no_variant_is_refused` for the VCF reader, both in the crate `popnei`. |
| `block::tests::a_block_of_more_genotypes_than_a_usize_holds_is_refused_when_the_collector_is_built` | The variants of a block times the individuals times the ploidy is a multiplication that a size a caller wrote carries beyond a `usize`, and that is an error where the reader is built and not a panic. | `block::tests::a_reblock_of_more_genotypes_than_a_usize_holds_is_refused_when_it_is_built` and `io::vcf::tests::a_block_of_more_genotypes_than_a_usize_holds_is_refused_when_the_reader_is_built`. |
| `block::tests::a_block_of_more_memory_than_the_machine_gives_is_refused_before_a_variant_is_read` | A block is refused for the memory of its columns and not only for genotypes that no `usize` counts: the positions of a variant are 8 bytes whatever the individuals are, so a source of no individual, whose blocks hold no genotype, is still refused. | `block::tests::a_block_of_more_memory_than_the_machine_gives_is_refused_when_it_is_started`, written in the commit before the removal over the VCF reader, with the chromosomes and the positions asked for and not the genotypes. A source of no individual is gone with the collector: a VCF with no individual is refused by its header. |
| `block::tests::a_collector_asks_its_reader_for_its_columns_and_the_genotypes` | What a reader of blocks asks its source for is what its blocks hold, the genotypes among them. | `block::tests::a_reader_asked_for_the_genotypes_alone_gives_a_block_with_no_other_column`, and `io::vcf::tests::with_the_genotypes_alone_a_block_has_no_column`, `io::vcf::tests::with_the_id_and_the_alleles_alone_a_block_has_those_two_columns_and_no_genotype` and `io::vcf::tests::what_is_in_the_columns_of_the_individuals_is_not_read_without_the_genotypes` for the reader that parses only what it was asked for. |
| `block::tests::a_variant_of_a_number_of_alleles_other_than_the_individuals_is_an_error` | A reader that gives a variant of a number of alleles other than its individuals times its ploidy is an error, `VariantOfAnotherSize`, and not a block whose genotypes are each at the place of another. | The rule holds of a block and not of a variant now, and `Block::check` is what finds it: `block::tests::a_block_whose_genotypes_lost_an_allele_fails_check` and `block::tests::reblock_refuses_a_block_whose_arrays_are_not_of_its_size`. The case `VariantOfAnotherSize` went with the test. |
| `block::tests::the_collector_gives_no_block_after_the_error_of_its_reader` | A reader over another reader gives no block after its source failed, and does not lean on the source stopping itself. | `block::tests::reblock_calls_a_source_that_failed_once_and_no_more`, which also counts the calls to the source. |
| `block::tests::the_blocks_of_a_reader_of_single_variants_are_a_reader_of_blocks` | `CollectedBlocks` was a `BlockReader`: its individuals, its ploidy, its chromosomes, its blocks of the size asked for and the `None` at the end. | Nothing: the thing it tested is gone. The VCF reader is the `BlockReader` that its tests read, `io::vcf::tests::the_blocks_of_many_vcf_have_the_sizes_of_the_spec` and the others made at `next_block`. |
| `block::tests::a_field_that_was_asked_for_and_that_the_reader_does_not_fill_is_an_error` | Two things: a reader held as `Box<dyn VariantReader>` works through the box; and a consumer that did not get a field it depends on gets the error that names it. | The box: `block::tests::what_is_asked_of_a_boxed_reader_reaches_the_reader_inside_it` and `block::tests::a_boxed_reader_of_blocks_is_read_through_reblock`, both over `Box<dyn BlockReader>`. The error: it is now `FieldsNotInTheBlock`, which a consumer builds from `asked_for.difference(block.fields())`, and no consumer is written yet; `error::tests::the_message_of_a_field_the_block_does_not_hold_names_the_field` keeps its message tested, and `block::tests::the_views_of_a_block_give_none_for_a_column_it_does_not_hold` checks what `Block::fields` reports of a block that lacks a column. |
| `variant::tests::a_cleared_variant_is_empty_and_keeps_the_capacity_of_its_buffers` | `Variant::clear` empties every field and keeps the capacity of the buffers, so a reader refilling a million variants allocates nothing after the first few. | Nothing: no variant is refilled any more. The rule it served, that a reader allocates a few times for each block and not for each variant, has no test that stays: "How it is verified" of `docs/specs/io_vcf.md` decides that it is checked by hand with a counting allocator when a reader is written. |
| `variant::tests::a_variant_cleared_but_the_alleles_keeps_them_and_empties_the_rest` | The same, for the reader that wrote over the strings of the alleles instead of dropping them. | Nothing: the alleles of a block are one buffer of text, `AllelesColumn`, whose buffers are checked by `block::tests::the_buffers_of_an_alleles_column_hold_a_block_of_biallelic_variants`. |
| `variant::tests::two_variants_are_read_through_a_boxed_reader` | Three things: a reader through `Box<dyn VariantReader>`; the chromosomes numbered in the order in which their names first appear; and a reader asked for other fields between two reads giving the next one without the columns that were dropped. | The box: `block::tests::what_is_asked_of_a_boxed_reader_reaches_the_reader_inside_it`. The chromosomes: `variant::tests::a_chrom_table_numbers_the_names_in_the_order_they_first_appear` and `io::vcf::tests::the_chromosomes_are_numbered_in_the_order_of_the_variants_that_are_given`. The change of fields: `io::vcf::tests::a_reader_asked_for_the_genotypes_alone_gives_its_next_block_without_the_columns`. |
