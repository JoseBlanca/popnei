# Plan: readers that give blocks

September 2026. Draft, not approved. The owner decided on 20 September
2026 that the variants flow through popnei in blocks, from the source to
the calculation, and that the single variant that a reader filled goes;
`docs/architecture.md`, as revised that day, has the decision and its
reasons in section 1. This plan takes the code that
`docs/plans/vcf-to-blocks.md` built, which has that single variant at its
centre, to the new shape: a `BlockReader` trait, a VCF reader that parses
its lines straight into the rows of a block, and no `Variant`,
`VariantReader` or `BlockCollector`. It also carries out five smaller
decisions that the owner took the same day. The owner answered its
breakdown in chat on 20 September 2026.

It was written before the specs it builds from were revised, on the
owner's order, and that is its first condition: the three specs below are
revised for blocks and reviewed before the plan is approved, and the
parts of them that its tasks name are then checked against their
headings. The specs, which another session of the assistant revises:

- `docs/specs/variant.md`, what is left of the variant module: the set of
  wanted fields, the table of chromosomes, the view of one variant of a
  block;
- `docs/specs/block.md`, the block, the `BlockReader` trait and
  `reblock`;
- `docs/specs/io_vcf.md`, the VCF reader.

When it is done nothing changes for a user but two errors that were not
there and the speed: `open_vcf`, `openVcf`, `Variants` and `iter_blocks`
are what they were, their tests pass untouched, and the VCF reader is
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

- The three specs above, revised for blocks, each through the first
  reader and the spec reviewer, and the revised `docs/architecture.md`
  and `docs/glossary.md` committed. Checked by reading: none of the three
  names `read_variant`, `VariantReader` or `BlockCollector` as something
  to build, and `docs/specs/block.md` has `BlockReader` in "The Rust
  interface". The revised `docs/specs/io_vcf.md` has to say, because the
  tasks of work package 2 build from it: how a line that gives no variant,
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
  the median of 5 runs: 1.24 s on one thread and 0.160 s on 18 plain,
  1.58 s and 0.50 s bgzipped; the spike of `docs/rust_core.md` on the same
  file 0.54 s and 0.098 s, 0.84 s and 0.40 s.
- pyNei's commit ef0ca6e is on `origin/main` of
  `https://github.com/JoseBlanca/pynei`, checked on 20 September 2026.

The checks that say the work is done fail today, as they must: `grep -rn
"read_variant\|VariantReader\|BlockCollector" crates --include='*.rs' |
wc -l` gives 98, in seven files, and the same for
`"BlockReader\|VariantRef\|reblock"` gives 0; there is no `LICENSE` file
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
   && uv run pytest -k pynei` gives `8 passed`; `docs/objectives.md` says
   what the dependency now is, in the two places where it says path
   dependency.
3. The `coding` skill. Check: "Errors, and no panics" of
   `.claude/skills/coding/SKILL.md` says that the core crate has one
   error enum, `non_exhaustive`, to which each module adds its cases, and
   why, that the two binding crates turn it into the exceptions of their
   language in one place each, and no longer asks for one error type for
   each operation.

**What it stands on.** The owner's decisions of 20 September 2026. The
`coding` skill has edits of the other session that are not committed:
deliverable 3 waits until they are.

**Tasks.**

- [ ] 1.1 The license and the dependency on pyNei: the `LICENSE` file,
  the field in the workspace manifest, in `pyproject.toml` and in
  `package.json`; pyNei as a git source of uv at ef0ca6e, with the
  comment of `pyproject.toml` that explains the absolute path replaced;
  the two sentences of `docs/objectives.md`. Serves deliverables 1 and 2.
- [ ] 1.2 The paragraph of the `coding` skill, written as the `writing`
  skill asks, with the reasons the owner took: one type is what the two
  bindings map, and `non_exhaustive` lets a module add a case without
  breaking a caller. Serves deliverable 3.

**What could go wrong.** `uv` fetches pyNei from GitHub, so the first
`uv sync` needs the network, and pyNei's own dependencies, pandas and
pyarrow, are built or fetched for Python 3.14.

## Work package 2: the VCF reader gives blocks

**What it gives.** The same `open_vcf` and `openVcf`, on a reader that
parses its lines into blocks and no longer into single variants, faster,
and that refuses two files it took before. It is the work package that
would change the plan if it failed, and everything after it is removal.

**Deliverables.**

1. `BlockReader` in the core, and `Block::variant` with its view, as
   "The Rust interface" of `docs/specs/block.md` has them. Check: `cargo
   test -p popnei --lib block:: -- --list` names tests of the view, one
   variant of `cases.vcf` seen through it with the literals of the table
   of `docs/specs/io_vcf.md`, and of a reader used as a boxed
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
   the sizes asked for there; the same blocks come out in rayon pools of
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
5. Both bindings hold a boxed `BlockReader`. Check: `uv run maturin
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

- [ ] 2.1 `BlockReader`, `Block::variant` and its view, in the core,
  with their tests, and a reader, marked as one that goes in work package
  3, that gives the blocks of the `BlockCollector` that exists, so that
  the next task can move the bindings before the VCF reader changes. From
  "The Rust interface" of `docs/specs/block.md`. Serves deliverable 1.
- [ ] 2.2 The two binding crates hold a `Box<dyn BlockReader>`. No file
  of `python/`, `js/popnei/src/`, `tests/` or `js/popnei/test/` changes.
  From "In Python and in TypeScript" of `docs/specs/block.md`. Serves
  deliverable 5. Needs 2.1.
- [ ] 2.3 The row parser over bytes and its tests, a function with no
  reader around it. From "What it gives" and "The cases a reader of the
  rules would not guess" of `docs/specs/io_vcf.md`. Serves deliverable 2.
  Side by side with 2.2: it touches `crates/popnei/src/io/vcf.rs` and no
  manifest. A wrong genotype here is silent anywhere else, so it is a
  commit of its own, and the comparison with bcftools of deliverable 3 is
  what guards it.
- [ ] 2.4 `VcfReader` as a `BlockReader` on the row parser: the lines of
  a block, the rows given to the lines that have one, the threads, the
  numbers of the chromosomes, the errors in their order, the reader that
  refuses to go on after a parse that did not come back; the tests of
  the reader made at `next_block`; the bench reading through
  `next_block`; and the table of the 65 tests. From "How it runs" and
  "How it is verified" of `docs/specs/io_vcf.md`. Serves deliverables 2,
  3 and 5. Needs 2.2 and 2.3.
- [ ] 2.5 The two new errors, in the core and seen from Python. From the
  parts of `docs/specs/io_vcf.md` that state them. Serves deliverable 4.
  Needs 2.4.
- [ ] 2.6 The measurement, as "What measurement there is" of the
  `performance-review` skill asks. It changes no code but the constants
  that the spec leaves to a measurement, each in a commit of its own with
  its numbers. Serves deliverable 6. Needs 2.4.

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

## Work package 3: the record level goes, and reblock comes

**What it gives.** A core in which a block is the one shape of the
variants, as the architecture says, and the reader over a reader that
the vars file, the filters and `iter_blocks` with a size will need.

**Deliverables.**

1. No record level. Check: `grep -rn
   "read_variant\|VariantReader\|BlockCollector" crates --include='*.rs'`
   finds nothing; `Variant` is not a public item of the core, which
   `cargo doc -p popnei --no-deps` shows; the reader of task 2.1 that
   stood on the collector is gone; the cases of the error enum that only
   the record level could give are gone, and the two binding crates still
   map every case that is left.
2. `reblock`, as `docs/specs/block.md` has it. Check: cargo tests with
   a reader written in the test that gives blocks of uneven sizes, 3, 0
   is never given, 5 and 1 variants: asked for 4 it gives 4, 4 and 1, and
   asked for 100 it gives one of 9; every column that is there is cut and
   joined with the genotypes, the alleles among them; the blocks joined
   are those of the source; the error of the source comes after the
   blocks that were whole before it.
3. Everything still passes. Check: the final check below.

**What it stands on.** Work package 2.

**Tasks.**

- [ ] 3.1 The removal: `Variant` and `VariantReader` out of
  `crates/popnei/src/variant.rs`, `BlockCollector` and the reader of
  task 2.1 out of `crates/popnei/src/block.rs`, the cases of the error
  that go, the tests that went with them, each named in the work report
  with the test of work package 2 that covers what it covered. From
  "Not in this spec" and "The Rust interface" of `docs/specs/variant.md`
  and `docs/specs/block.md`. Serves deliverable 1.
- [ ] 3.2 `reblock` and its tests. From the part of `docs/specs/block.md`
  that has it. Serves deliverable 2. Side by side with 3.1 only if the
  spec puts it in a file of its own; otherwise after it.

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
