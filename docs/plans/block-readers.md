# Plan: readers that give blocks

September 2026. Approved by the owner on 20 September 2026 and done on 21
September 2026, on the branch `plan/block-readers`, with its work report
in `docs/reports/block-readers.md`. The owner decided on 20 September
2026 that the variants flow through popnei in blocks, from the source to
the calculation, and that the single variant that a reader filled goes;
`docs/architecture.md`, as revised that day, has the decision and its
reasons in section 1. This plan takes the code that
`docs/plans/vcf-to-blocks.md` built, which has that single variant at its
centre, to the new shape: a `BlockReader` trait, which everything that
gives blocks implements, a VCF reader that parses its lines straight into
the rows of a block, and none of the three things of the code as built:
`Variant`, the single variant that the caller owns and a reader fills
through `read_variant`; `VariantReader`, the trait with that method; and
`BlockCollector`, which copies the variants of such a reader into
blocks. It also carries out five smaller
decisions that the owner took the same day. The owner answered its
breakdown in chat on 20 September 2026.

It was written before the specs it builds from were revised, on the
owner's order. Another session of the assistant revised them for blocks
the same day, each through the spec reviewer and the first reader, and
the plan was then checked against them: every part of a spec that a task
names is a heading they have, and `reblock` moved from the third work
package to the second, because the revised `docs/specs/block.md` puts a
`reblock` at the end of every `iter_blocks`, so the bindings need it. The revised
specs are on `main` since 8f19650, and the revised architecture, glossary
and skills since cae6a4c. The session that will run the plan checked it
once more against them as committed, at d211355, and what that changed is
small: task 1.1 also corrects section 8 of the architecture, tasks 2.4
and 2.5 name "The Rust interface" of the VCF reader spec, tasks 2.2 and
2.3 run one after the other, and the paragraph below names a third
change that a user sees. The specs:

- `docs/specs/variant.md`, what is left of the variant module: the set of
  wanted fields, the table of chromosomes, the view of one variant of a
  block;
- `docs/specs/block.md`, the block, the `BlockReader` trait and
  `reblock`;
- `docs/specs/io_vcf.md`, the VCF reader.

When it is done three things change for a user, and the speed. Two errors
that were not there. And a column that was not asked for is no longer
checked, as "How it runs" of `docs/specs/io_vcf.md` decides: the reader
as built parses the position of every line, so a position that is not a
number is an error today whatever `fields` is, and after this plan it is
an error only when `fields` has `"chrom"` or `"pos"`, which the default
has. `open_vcf`, `openVcf`, `Variants` and `iter_blocks` are otherwise
what they were, their tests pass untouched, and the VCF reader is
measured again against the targets that it missed.

## In and out

In:

- The `BlockReader` trait and the view of one variant of a block, of
  sections 1 and 2 of the architecture.
- The VCF reader rewritten to give blocks, with the rows parsed on the
  threads of rayon straight into the arrays of the block, and the columns
  of the individuals parsed as bytes and not as text.
- Two decisions of the owner of 20 September 2026 about what the VCF
  reader refuses: a QUAL that is not finite, and a bgzipped source that
  ends without the empty block that marks the end of a bgzip file.
- `reblock`, a reader over a reader that cuts and joins the blocks of its
  source to a size.
- Three more decisions of that day: the license is MIT; pyNei is a
  development dependency on its GitHub repository at a commit and not a
  path; the crate keeps one error enum, and the line of the `coding` skill
  that asks for one error type for each operation is corrected.

Out, with where it goes:

- The reader and the writer of the vars file: `docs/specs/io_vars.md` and
  a plan of its own.
- The filters, the row helpers of one variant, the dosages, the masks and
  the allele counts, and every calculation: their specs are not written.
- The read ahead thread of section 3 of the architecture. The review of
  the batched reader measured that the serial read of the lines is what
  bounds the 18 threads, and that thread is what removes it; it comes
  with the first consumer that has work to do while the next block is
  read. If the 18 thread target is missed for that reason alone, task 2.4
  says so with the numbers.
- A `File` that the user picked in a page, read inside a web worker,
  which was out of the plan before this one and stays out.
- Publishing, continuous integration and the development container.

## What has to be in place

- The three specs above, revised for blocks and reviewed, and the
  revised `docs/architecture.md` and `docs/glossary.md`, all committed on
  `main`: `git status --short docs .claude` shows none of them as
  changed. Checked by reading: none of the three specs names
  `read_variant`, `VariantReader` or `BlockCollector` as something to
  build, `docs/specs/block.md` has `BlockReader` and `Reblock` in "The
  Rust interface", and each has "Open points: None". The revised
  `docs/specs/io_vcf.md` says, in "How it runs", what the
  tasks of work package 2 build from: how a line that gives no variant,
  one that failed its FILTER or an empty one, gets no row while the rows
  are parsed side by side; how many lines are held as text at a time, so
  that the memory of a reader does not grow without limit, and how a
  block still gets the size that was asked for; how the alleles and the
  ids, which are one buffer each, are filled from rows parsed side by
  side; when the chromosomes get their numbers; that an error loses the
  block it is in and the blocks before it are given first; what a reader
  whose parse panicked does at the next call; the two new errors, and
  with the second what tells a file that bgzip wrote from a plain gzip
  one, which has no such end and is still read; and the numbers of
  "Speed".
- `main` at 13fe55f or later, where the five commands of "Before the
  work is called done" of the `coding` skill pass, `cargo wasm-check`
  passes, `npm test` in `js/popnei` gives `tests 39`, `fail 0`, and
  `bash scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs`
  exits with 0. Run on 20 September 2026 after the merge of the plan
  before this one: `92 passed`, `38 passed`, `tests 39`, exit 0.
- The toolchain that `docs/plans/vcf-to-blocks.md` lists under "What has
  to be in place", checked the same way.
- The file of the bench: `uv run --no-project --with numpy python
  crates/popnei/benches/make_big_vcf.py <dir>/big.vcf`, 3 s, 403 MB, and
  `bgzip -k` of it, outside the repository. The numbers that this plan
  compares with, from `docs/reports/vcf-to-blocks.md`, the owner's Apple
  M5 Pro, release, the file in the page cache, the genotypes asked for,
  the median of 5 runs. The reader as built: the plain file in 1.24 s on
  one thread and 0.160 s on 18, the bgzipped one in 1.58 s on one thread
  and 0.50 s on 18. The spike, the trial parser in Rust of
  `docs/rust_core.md`, on the same files: plain 0.54 s and 0.098 s,
  bgzipped 0.84 s and 0.40 s. "Speed" of `docs/specs/io_vcf.md` has the
  table.
- pyNei's commit ef0ca6e is on `origin/main` of
  `https://github.com/JoseBlanca/pynei`, checked on 20 September 2026.

The checks that say the work is done fail today, as they must: `grep -rn
"read_variant\|VariantReader\|BlockCollector" crates --include='*.rs' |
wc -l` gives 98, in seven files, and the same for
`"BlockReader\|VariantRef\|reblock"`, the trait, the view of one variant
of a block and the reader that resizes blocks, gives 0; there is no `LICENSE` file
and no manifest names a license; `pyproject.toml` has pyNei at
`/Users/jose/devel/pynei`.

The plan runs in a worktree, and what `docs/reports/vcf-to-blocks.md`
learned holds: tasks that write the workspace manifest or `Cargo.lock`
run one after another; every reviewer works in a worktree of its own;
what a subagent says of a file is checked before it goes into the
report; the count of the tests of the library is read with `--lib`.

## Work package 1: the five decisions

**What it gives.** A repository that says under which license it is, that
builds its tests on any machine with the network and not on the owner's
alone, and whose `coding` skill agrees with its code about errors. It
comes first because it is small, stands on nothing, and the next work
packages add errors to the enum that the skill will then describe.

**Deliverables.**

1. The license. Check: a `LICENSE` file at the root with the text of the
   MIT license, the year 2026 and the owner's name as `git log` has it;
   `cargo metadata --no-deps --format-version 1` shows `"license":"MIT"`
   for the three crates; `pyproject.toml` and `js/popnei/package.json`
   name it; the wheel of pyodide and `npm pack --dry-run` in `js/popnei`
   carry it.
2. pyNei by its repository. Check: `grep -n "Users/jose" pyproject.toml
   uv.lock` finds nothing; `uv sync` followed by `uv run maturin develop
   && uv run pytest -k pynei` gives `8 passed`; the two documents that
   say path dependency say what the dependency now is, objective 1 of
   `docs/objectives.md` and section 8 of `docs/architecture.md`, and
   `grep -rn "path dependency" docs/objectives.md docs/architecture.md`
   finds nothing, where today it finds those two.
3. The `coding` skill. Check: "Errors, and no panics" of
   `.claude/skills/coding/SKILL.md` says that the core crate has one
   error enum, `non_exhaustive`, to which each module adds its cases, and
   why, that the two binding crates turn it into the exceptions of their
   language in one place each, and no longer asks for one error type for
   each operation.

**What it stands on.** The owner's decisions of 20 September 2026. The
`coding` skill as committed at cae6a4c still asks, in "Errors, and no
panics", for one error type for each operation.

**Tasks.**

- [x] 1.1 The license and the dependency on pyNei: the `LICENSE` file,
  the field in the workspace manifest, in `pyproject.toml` and in
  `package.json`; pyNei as a git source of uv at ef0ca6e, with the
  comment of `pyproject.toml` that explains the absolute path replaced;
  the sentence of objective 1 of `docs/objectives.md` and the one of
  section 8 of `docs/architecture.md`. Serves deliverables 1 and 2.
- [x] 1.2 The paragraph of the `coding` skill, written as the `writing`
  skill asks, with the reasons the owner took: one type is what the two
  bindings map, and `non_exhaustive` lets a module add a case without
  breaking a caller. The line of "errors" in
  `.claude/skills/code-review/categories.md` that asks a reviewer for one
  error type per operation is corrected with it. Serves deliverable 3.

**What could go wrong.** `uv` fetches pyNei from GitHub, so the first
`uv sync` needs the network, and pyNei's own dependencies, pandas and
pyarrow, are built or fetched for Python 3.14.

## Work package 2: the VCF reader gives blocks

**What it gives.** The same `open_vcf` and `openVcf`, on a reader that
parses its lines into blocks and no longer into single variants, faster,
and that refuses two files it took before. It is the work package that
would change the plan if it failed, and everything after it is removal.

**Deliverables.**

1. In the core, as "The Rust interface" of `docs/specs/block.md` and of
   `docs/specs/variant.md` have them: the `BlockReader` trait with its
   contract; what `Block` gains, `fields`, `variants`, `variant`,
   `retain_vars` and `check`; the view of one variant; and `Reblock`.
   Check: the cargo tests that "How it is verified" of
   `docs/specs/block.md` lists for the block, made at `retain_vars`,
   `variants` and `check`, and those of `reblock` that are made over a
   reader written in the test, the blocks of one variant that go through
   with no copy and the source that is called once and no more after its
   error, exist and pass; a test uses a reader as a boxed
   `dyn BlockReader`.
2. The row parser: one data line, as bytes, written into one row of a
   block, the genotypes into a slice of individuals x ploidy alleles.
   Check: every case of a data line that "How it is verified" of
   `docs/specs/io_vcf.md` lists, an error or not, has a cargo test made
   at this function or at the reader with the spec's literals; the cases
   that the plan before this one tested at `read_variant` are all still
   tested, and a table in the work report says, for each of the 65 tests
   that `crates/popnei/src/io/vcf.rs` has today, which test took its
   place or why it went.
3. `VcfReader` is a `BlockReader`. Check: the cargo tests on
   `cases.vcf`, `differences.vcf`, `many.vcf` and their gzipped forms
   against the stored output of bcftools and the sixteen counts of the
   spec pass, made at `next_block`; the blocks of `many.vcf` have the
   sizes that "How it is verified" of `docs/specs/block.md` gives, for
   the sizes asked for there, and the tests of `reblock` made over a
   `VcfReader`, the 72 blocks of 7 and the tetraploid genotype of line
   251, pass; the same blocks come out in rayon pools of
   1 and of 4 threads, with the same chromosome numbers; a wrong line
   gives the blocks before it and then its error and then nothing; the
   serial parse, which is the one wasm runs, is tested natively against
   the parallel one.
4. The two new errors, each with a cargo test and a pytest test that
   sees a `ValueError`: a QUAL of `nan`, of `inf` and of `1e400`; and
   `many.vcf.gz` cut at the end of its second gzip member, its first
   12336 bytes, which today gives the variants of that member and no
   error, as a review of the plan before this one found on the file of
   that day. A gzipped file that bgzip did not write, plain
   gzip, is still read.
5. Both bindings hold a boxed `BlockReader`, with a `Reblock` at the
   end of every `iter_blocks` and `Block::check` before the genotypes
   cross, as `docs/specs/block.md` asks. Check: `uv run maturin
   develop && uv run pytest` `38 passed` and more, with no test of
   `tests/` changed but for the new ones; `npm test` in `js/popnei`
   `tests 39` and more, `fail 0`, none changed; the wheel of pyodide
   builds and its smoke test exits with 0; `grep -rn
   "VariantReader\|BlockCollector" crates/popnei-python crates/popnei-js`
   finds nothing.
6. The measurement, with `crates/popnei/benches/read_vcf.rs` reading
   through `next_block`: the four timings of "What has to be in place",
   taken the same way, in the work report beside the numbers of today,
   those of the spike and the targets of "Speed" of the spec; a sampling
   profile of the one thread run; and the memory of a reader, measured,
   for 1000 and for 10000 individuals. Check: the report has them, and
   says of each target whether it is met.

**What it stands on.** Work package 1 for the enum that the skill
describes, no more. The three revised specs.

**Tasks.**

- [x] 2.1 In the core: `BlockReader`, the methods that `Block` gains,
  the view of one variant, `Reblock`, the cases of the error that they
  add, and their tests; and a reader, marked as one that goes in work
  package 3, that gives the blocks of the `BlockCollector` that exists,
  so that the next task can move the bindings before the VCF reader
  changes. From "What it gives", "What a reader of the rules would not
  guess", "How it is verified" and "The Rust interface" of
  `docs/specs/block.md`, and "The Rust interface" of
  `docs/specs/variant.md` for the view. Serves deliverable 1.
- [x] 2.2 The two binding crates hold a `Box<dyn BlockReader>`, put a
  `Reblock` over it for the size that `iter_blocks` was asked for, and
  call `Block::check` before the genotypes cross. No file
  of `python/`, `js/popnei/src/`, `tests/` or `js/popnei/test/` changes.
  From "In Python and in TypeScript" of `docs/specs/block.md`. Serves
  deliverable 5. Needs 2.1.
- [x] 2.3 The row parser over bytes and its tests, a function with no
  reader around it. From "What it gives" and "The cases a reader of the
  rules would not guess" of `docs/specs/io_vcf.md`. Serves deliverable 2.
  It needs nothing of 2.2 and runs after it all the same: "Speed" of the
  spec has the implementer search the bytes with `memchr`, which is not a
  dependency today, so the task writes the manifest of the core crate and
  `Cargo.lock`, and a core that is half written breaks the build of the
  binding crates in the same tree. A wrong genotype here is silent anywhere else, so it is a
  commit of its own, and the comparison with bcftools of deliverable 3 is
  what guards it.
- [x] 2.4 `VcfReader` as a `BlockReader` on the row parser: the lines of
  a block, the rows given to the lines that have one, the threads, the
  numbers of the chromosomes, the errors in their order, the reader that
  refuses to go on after a parse that did not come back; the tests of
  the reader made at `next_block`; the bench reading through
  `next_block`; and the table of the 65 tests. From "How it runs",
  "How it is verified" and "The Rust interface" of
  `docs/specs/io_vcf.md`, which has the size of the blocks among the
  options of the reader and the error of a parse that did not come back.
  Serves deliverables 2, 3 and 5. Needs 2.2 and 2.3.
- [x] 2.5 The two new errors, in the core and seen from Python. From "The
  cases a reader of the rules would not guess", the last paragraph of
  the blocks in "How it is verified" and the cases of the error in "The
  Rust interface" of `docs/specs/io_vcf.md`. Serves deliverable 4.
  Needs 2.4.
- [x] 2.6 The measurement, as "What measurement there is" of the
  `performance-review` skill asks. It changes no code but the constants
  that the spec leaves to a measurement, each in a commit of its own with
  its numbers. Serves deliverable 6. Needs 2.4, and runs after 2.5: both
  write `crates/popnei/src/io/vcf.rs`, and a build beside a timing
  changes the timing. It also runs after task 3.1 and after the fixes of
  the review of tasks 2.3 to 2.5, an order that the orchestrator chose
  while the plan ran: the reviewers and task 3.1 work at the same time, in
  other trees and other files, no build runs beside the timing, and what
  is timed is the code that the plan leaves.

**What could go wrong.** The profile of the reader as it is puts 95 in
100 of the one thread time in the columns of the individuals read as
text, the split at a character, the search of a byte, the comparison of
strings, and not in how a variant is handed out; rows written into a
block do not close that alone, which is why 2.3 parses bytes. The spike,
`spike/pynei_spike` of pyNei, is the code whose numbers are the target,
and its row parser is there to read. For the end of a bgzip file, two
facts that task 2.5 meets: that end is a fixed block of 28 bytes, and
flate2 does not tell its caller where a gzip member ends. If a
target of "Speed" is still missed, the plan is done when the measurement
is reported, and the owner gets the numbers.

## Work package 3: the record level goes

**What it gives.** A core in which a block is the one shape of the
variants, as the architecture says, with nothing left of the single
variant that a reader filled.

**Deliverables.**

1. No record level. Check: `grep -rn
   "read_variant\|VariantReader\|BlockCollector" crates --include='*.rs'`
   finds nothing; `Variant` is not a public item of the core, which
   `cargo doc -p popnei --no-deps` shows; the reader of task 2.1 that
   stood on the collector is gone; the cases of the error enum that only
   the record level could give are gone, and the two binding crates still
   map every case that is left.
2. Everything still passes. Check: the final check below.

**What it stands on.** Work package 2.

**Tasks.**

- [x] 3.1 The removal: `Variant` and `VariantReader` out of
  `crates/popnei/src/variant.rs`, `BlockCollector` and the reader of
  task 2.1 out of `crates/popnei/src/block.rs`, the cases of the error
  that go, the tests that went with them, each named in the work report
  with the test of work package 2 that covers what it covered. From
  "Not in this spec" and "The Rust interface" of `docs/specs/variant.md`
  and `docs/specs/block.md`. Serves deliverable 1.

**What could go wrong.** A test that is deleted because the thing it
called is gone can take with it the only check of a rule that still
holds, the chromosome numbers in the order of the variants that are
given, the line number of an error: the table of task 2.4 is what the
reviewer of the tests reads against.

## How the whole plan is checked

From a clean clone of the branch: the five commands of the `coding`
skill; `cargo wasm-check`; `npm run build` and `npm test` in
`js/popnei`; the build of the wheel of pyodide and its smoke test; the
two `grep` of "What has to be in place", which now give 0 and more than
0; and the work report, with the measurement of deliverable 6 of work
package 2 and the table of the tests.
