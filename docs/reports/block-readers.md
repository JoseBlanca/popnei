# Work report: readers that give blocks

The plan `docs/plans/block-readers.md` is under way, on the branch
`plan/block-readers`, in the worktree `.claude/worktrees/block-readers`,
since 20 September 2026. The orchestrator, in this report, is the session
of the assistant that runs the plan: it sends each task to a subagent,
checks what comes back and has it reviewed. What the owner reads first,
whether the plan is done, what exists now and what is asked of them, is
written here when the plan ends.

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

One reviewer, `spec`, in a worktree of its own at a16ea06, 63 thousand
tokens: the work package has no code to break, so `tests` had nothing to
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
appended after each batch at a cost of 1 ms of 127 on 18 threads; that two
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
