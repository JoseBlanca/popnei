# Plan: the vars file, its writer and its reader

September 2026. Under way since 21 September 2026 on the branch
`plan/vars-file`, with its work report in `docs/reports/vars-file.md`.
Written on 20 September 2026; the owner answered its breakdown in chat
that day, and ordered it run as soon as `docs/plans/block-readers.md` was
merged into `main`, which it was on 21 September 2026, at 2d99d64. It
builds the module `io::vars` of the core crate and what a user reaches it
through: the vars file is the arrow file in which popnei keeps variants
once a VCF has been read, and there is no code of it. It is built from
three specs, all reviewed and with no open point:

- `docs/specs/io_vars.md`, the format of the file, the writer and the
  reader, which is the spec every task names unless it says another;
- `docs/specs/block.md`, the block, the run of consecutive variants held
  as arrays that the variants flow in, the `BlockReader` trait of
  everything that gives blocks, and `reblock`, the reader over a reader
  that puts blocks back to one size;
- `docs/specs/variant.md`, the `Needs` that say which fields a consumer
  wants, the `ChromTable` of the chromosome names, and the `Variants`
  handle that a Python or a TypeScript user holds.

When it is done a user can call `write_vars` on the `Variants` of a VCF
and `open_vars` on the file it wrote, in Python natively and under
pyodide, and `writeVars` and `openVars` in TypeScript, and any program
with an arrow library opens that file as a table, with its six columns,
`chrom`, `pos`, `id`, `alleles`, `qual` and `gts`. How long a pass over
the genotypes of such a file takes is measured against the 21 ms of
"Speed" of the spec.

It could not start before another plan, `docs/plans/block-readers.md`,
was done and merged into `main`, and "What has to be in place" says how
that was checked.

## In and out

In: what `docs/specs/io_vars.md` describes, in the core crate, in both
binding crates, in the Python package and in the TypeScript package; the
reference file that popnei cannot write itself, a vars file compressed
with zstd; and the measurement of "Speed".

Out, with where it goes:

- A `File` that the user picked in a page as the `source` of `openVars`,
  which "Its Python and TypeScript functions" of the reader gives beside
  the `Uint8Array`. The owner left it out on 20 September 2026, as
  `docs/plans/block-readers.md` leaves it out for `openVcf`: it is read
  through `FileReaderSync`, which exists only inside a web worker, and
  the tests run under node. It goes to a plan that does it for both
  readers. `VarsReader` is generic over `Read + Seek`, which is what that
  plan will need.
- What "Not in this spec" lists: the function that asks for the variants
  of a region, the read ahead thread, the genotypes packed in 2 bits,
  the vars file of pyNei, the VCF writer.
- Any work on speed. Work package 4 measures, and if the number of
  "Speed" is missed the owner gets the measurement.

The spec has no open point.

One sentence of the spec is corrected by task 1.1. "How it is verified"
of the writer names pyarrow 23.0.0 as the reference outside the project.
That is the pyarrow of pyNei's environment, where the trial of that
paragraph was run; popnei's environment has pyarrow 25.0.1, which comes
with pyNei and which no manifest of popnei names. The owner decided on 20
September 2026 that pyarrow becomes a development dependency of popnei
with 23 as its lowest version, and that the spec says that the tests run
on the pyarrow of `uv.lock`, 25.0.1 that day, and that the trial was run
on 23.0.0. The option not taken was to pin 23.0.0.

Two more sentences of the spec are written by task 1.1, each from a rule
the owner gave, and the work report lists both for him to reverse.

The exceptions of Python. "The Rust interface" says that every case this
module adds is a `ValueError` but the error of the input or the output.
On 21 September 2026 the owner gave the convention that "Errors, and no
panics" of the `coding` skill now has, and asked that the spec of a
module list its cases with the exception each one is. For this module: a
`ValueError` is a source that is not a vars file, whatever it lacks, a
format version that is not 1, a column of another type, a `gts` width
that does not match, a null where there can be none, a footer whose
entries do not match the batches, a batch with another number of
variants than its entry, a file compressed with zstd, two individuals
of one name, a path that exists and a `num_vars_per_block` of 0; an
`OSError`, with the path in `filename`, is an error of the input or the
output, a file that starts as an arrow file and was cut short, and a
batch that arrow-rs cannot decode or decompress; a `RuntimeError`, a
defect of popnei, is a block that does not fit the writer, a block whose
columns differ from those of the first one written, and a chromosome
number with no name. Every error of a file names the file, which the
Python binding crate adds with `PyPopneiError::of_the_file`.

A batch of no variants. `docs/specs/block.md` says that a reader never
gives a block of no variants, and the spec of the vars file does not say
what its reader does with a batch of none, which popnei's writer never
makes and another arrow program can. The reader takes the next batch, as
a filter does with a block it emptied, when the entry of the footer says
0 too. The owner was asked on 21 September 2026 and had not answered
when the work started; the option not taken is an error.

## What has to be in place

The first item is checked on `main` and gates the worktree. The others
are checked in the worktree, once it is made, before the first task.

- `docs/plans/block-readers.md` run to its end and merged into `main` on
  the owner's order. The whole plan and not its task 2.1, where the
  things this plan builds on first appear: its work package 3 removes
  `Variant`, `VariantReader` and `BlockCollector`, its task 2.1 adds a
  reader that work package 3 takes out again, and its task 2.2 changes
  how both binding crates hold a reader. Checked on `main`, by running:
  `git worktree list` and `git branch --merged main` show
  `plan/block-readers` merged; `docs/reports/block-readers.md` is there
  and says the plan is done; `grep -rn
  "read_variant\|VariantReader\|BlockCollector" crates --include='*.rs'`
  finds nothing; `grep -n "pub trait BlockReader\|pub struct Reblock\|pub
  fn check\|pub fn fields\|pub fn retain_vars"
  crates/popnei/src/block.rs` finds the five; `grep -n "Users/jose"
  pyproject.toml uv.lock` finds nothing. On 20 September 2026, at a11aea9,
  the branch had no commit of its own, the first grep found 98 lines and
  the second none. The session of the assistant that runs that plan, which
  `ListAgents` showed that day under the name `popnei-d5`, agreed to send
  a message when the merge is on `main`; with no message, the commands
  above are what says it.
- The worktree and the branch of this plan are made after that merge and
  not before. Both plans write the workspace `Cargo.toml`,
  `crates/popnei/Cargo.toml`, `Cargo.lock`, `crates/popnei/src/block.rs`,
  `error.rs` and `io/mod.rs`, the error mapping of both binding crates,
  `pyproject.toml` and `uv.lock`.
- On `main` after the merge: the five commands of "Before the work is
  called done" of the `coding` skill pass, `cargo wasm-check` passes,
  `npm run build` and `npm test` in `js/popnei` pass, and `bash
  scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` exits
  with 0. The orchestrator writes into the work report the three counts
  it gets there, of `cargo test --workspace`, of pytest and of `npm
  test`, and the size in bytes of `js/popnei/wasm/popnei_bg.wasm`: every
  "and more" of this plan is counted from them. Before that merge they
  were `92 passed`, `38 passed`, `tests 39` and 189872 bytes, and on 21
  September 2026, in the worktree at 2d99d64, `173 passed` and 1 ignored,
  `72 passed`, `tests 46` and 225344 bytes.
- The toolchain that `docs/plans/vcf-to-blocks.md` lists under "What has
  to be in place", each program by running the command that plan gives
  for it.
- The three specs on `main` with "Open points: None". Checked on 20
  September 2026 at d211355.

What was run while this plan was written, on 20 September 2026 on the
owner's Apple M5 Pro:

- A trial crate outside the repository with `arrow-ipc` 60.0.0, its
  default features off and `lz4` on, `arrow-array`, `arrow-schema`,
  `arrow-buffer` and `serde_json`. A test that writes a batch compressed
  with lz4 and a key of the footer into a `Vec<u8>` and reads it back
  passes natively, and `cargo check` finishes for
  `wasm32-unknown-unknown` and for `wasm32-unknown-emscripten`. Not
  tried: the link of the wheel of pyodide, which task 1.3 is the first
  to make.
- pyarrow 25.0.1 wrote a file of four variants with the six columns of
  "What it holds", compressed with zstd, with the `popnei` key in the
  schema and the `popnei_batches` key in the footer, which
  `pyarrow.ipc.new_file` takes as its `metadata` argument, and read both
  keys back. So the script of task 2.1 can be written.
- From `tests/reference/vcf/many.bcftools.tsv`: 500 rows, 167 with a dot
  as the id and 100 with a dot as the quality, and the two tables of
  regions of "How it is verified" of the writer, the one of the 500
  variants and the one of the 475 that passed, row for row.

The checks of the deliverables fail today, as they must: `cargo test -p
popnei --lib io::vars -- --list` prints `0 tests`; `grep -rn
"open_vars\|write_vars\|openVars\|writeVars" crates python js/popnei/src
tests/*.py js/popnei/test` finds nothing; `tests/test_io_vars.py`,
`tests/reference/vars/` and a bench of the vars file are not there;
`Cargo.lock` has no line with `arrow`; `pyproject.toml` does not name
pyarrow. The orchestrator runs them again on the commit the branch is cut
from.

What `docs/reports/vcf-to-blocks.md` learned holds here. Every task
writes a manifest, the core or a binding, so the tasks run one after
another, 2.1 aside: a task that writes a manifest or `Cargo.lock` does
not run beside another, and a core that is half written breaks the build
of the binding crates in the same tree. Every reviewer works in a
worktree of its own. What a subagent says of a file is checked before it
goes into the report. The count of the tests of the library is read with
`--lib`. The subagents that write code run on Opus.

## Work package 1: `write_vars` in Python

**What it gives.** A Python user calls `popnei.write_vars(variants, path)`
on the `Variants` of a VCF and gets a vars file that pandas, R and polars
open as a table. The writer comes before the reader because the tests of
the reader need vars files and popnei's writer is what makes them, and a
reader tested on files that pyarrow has opened is tested on files that a
program outside the project vouched for. It is also the work package that
would change the plan if it failed: its third task is the first to link
arrow-rs into the wheel of pyodide.

**Deliverables.**

1. arrow-rs and a json crate in the core, and what the file says about
   itself. Check: `cargo wasm-check` passes with `arrow-ipc`, its `lz4`
   feature on and no `zstd` crate, among the dependencies of the core,
   and `cargo tree -p popnei | grep -i zstd` finds nothing; `cargo test
   -p popnei --lib io::vars -- --list` names tests, 6 or more, that the
   value of each of the two keys is written and parsed back as "What it
   holds" has it, the example of `popnei_batches` of that part among
   them, a batch with no regions, and a value that is not json or lacks
   one of the four keys of `popnei` as the error that says so;
   `docs/specs/io_vars.md` has the sentence about pyarrow as "In and out"
   gives it.
2. `VarsWriter` and the core `write_vars`, as "The Rust interface" has
   them. Check: cargo tests made at `write_block` and at `write_vars`
   into a `Vec<u8>`, which the test opens with the `FileReader` of
   arrow-rs. From the block built by hand of the four variants of the
   `cases.vcf` table of `docs/specs/io_vcf.md`, with `num_vars_per_block`
   3: two batches, of 3 variants and of 1, the six columns with the types
   and the nulls of the table of "What it holds", null ids in the last
   three variants, both keys with the values the first cargo test of "How
   it is verified" of the writer gives. The second and the fourth cargo
   tests of that part, the variants that are not sorted and the source
   with the genotypes alone, read the same way. A source with no variants
   gives a file with both keys, a `gts` column and no batch. And one test
   for each error of `write_block` that "The Rust interface" lists: a
   block of other individuals, of another ploidy, one that fails
   `Block::check`, a chromosome number with no name, and a later block
   with other columns than the first, which names the field; after each,
   nothing of that block is in the sink. `cargo test -p popnei --lib
   io::vars -- --list` names 16 tests or more, the 6 of deliverable 1
   among them.
3. `popnei.write_vars`. Check: `uv run maturin develop && uv run pytest
   tests/test_io_vars.py -k write_vars` runs 6 tests or more and all
   pass: the four pytest tests of "How it is verified" of the writer,
   the first of them in its three forms, `only_passed=False` with
   batches of 100, the default of `only_passed`, and no
   `num_vars_per_block`, with the literals of the spec, the two tables of
   regions and the 167 and the 100 nulls among them, the file opened
   with pyarrow and compared with `many.bcftools.tsv`. `uv run pytest`
   gives the count of "What has to be in place" and those more, with no
   test that was there changed. pyarrow is in the `dev` group of
   `pyproject.toml`.
4. The wheel of pyodide links arrow-rs. Check: `bash
   scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` exits
   with 0, and the smoke test calls `write_vars` on `cases.vcf` and finds
   a file that starts and ends with the six bytes `ARROW1`, which an
   arrow IPC file does.

**What it stands on.** `BlockReader`, `Reblock`, `Block::check` and
`Block::fields` of the core, and the binding crates, as
`docs/plans/block-readers.md` leaves them.

**Tasks.**

- [ ] 1.1 The dependencies and the two keys. In the workspace manifest
  and in `crates/popnei/Cargo.toml`, each with the comment that says why
  it is there and that it builds for both wasm targets, as the manifest
  does for the others: the crates of arrow-rs that the module needs,
  `arrow-ipc` with its default features off and `lz4` on, and a json
  crate in pure Rust. In `crates/popnei/src/io/vars.rs`, new:
  `VarsMetadata`, `BatchInfo`, `Region`, the json of the two keys,
  written and parsed, and the cases that "The Rust interface" adds to
  the error of the crate, all of them, so that the next tasks add none.
  In the spec, in a commit of its own before the code: the sentence
  about pyarrow, the cases of the error with the exception each one is,
  and the batch of no variants, the three as "In and out" gives them.
  From "What it holds" and "The Rust interface". Serves deliverable 1.
- [ ] 1.2 `VarsWriter` and `write_vars` in the core, in
  `crates/popnei/src/io/vars.rs`, and their tests.
  From "What it gives" and "How it runs" of the writer, "How it is
  verified" of the writer for the cases, and "The Rust interface". The
  block goes in by value and the vector of its genotypes becomes the
  buffer of the `gts` column with no copy: a test sees it by the address
  of that buffer, which is the one the vector had, or the task says in
  the report why a test cannot see it. Serves deliverable 2. Needs
  1.1. A wrong null or a wrong region here is silent, so it is a commit
  of its own, and the comparison with `many.bcftools.tsv` of deliverable
  3 is what guards it.
- [ ] 1.3 The Python side: in `crates/popnei-python`, the function that
  opens the file, refuses a path that exists, calls the core
  `write_vars` with the interpreter released, and removes the file when
  the core gives an error; `write_vars` in the Python package, in a
  module `io_vars.py`, exported by the package; pyarrow in the `dev`
  group and in `uv.lock`; `tests/test_io_vars.py` with the four tests;
  the call in the smoke test of pyodide. From "Its Python and TypeScript
  functions" and "How it is verified" of the writer, and
  `.claude/skills/coding/pyo3.md`. Serves deliverables 3 and 4. Needs
  1.2.

**What could go wrong.** The link of arrow-rs under emscripten has not
been tried; `cargo check` for that target passes, and section 5 of
`docs/rust_core.md` has how the link of a C dependency failed there. If
task 1.3 cannot build the wheel, the orchestrator stops and the owner
gets what the linker said: the plan has no way around it. arrow-rs pulls
in crates that the core did not have, `chrono` and `half` among them:
the task reads what `cargo tree` gives for both wasm targets and turns
off the features it does not use. The key of the footer is written with
`write_metadata` of the `FileWriter` of arrow-rs, which the trial used.

## Work package 2: `open_vars` in Python

**What it gives.** A Python user calls `popnei.open_vars(path)` and gets
the same `Variants` handle as from `open_vcf`, over a file that is read
in blocks with only the columns that were asked for decompressed.

**Deliverables.**

1. The file that popnei cannot write. Check: `uv run python
   tests/reference/vars/make_reference.py` writes
   `tests/reference/vars/zstd.vars` again with no change for git; pyarrow
   opens it and finds the four variants of `cases.vcf`, both keys, one
   batch, and buffers compressed with zstd.
2. A vars file is opened. Check: cargo tests made at `VarsReader::new`
   over a `Cursor<Vec<u8>>`, on files built in the test with arrow-rs:
   `metadata`, `batches` and `num_vars` of a file that the writer made
   give what went in before any batch is read; each error of the file as
   a whole of "What it refuses" and of "How it is verified" of the
   reader, a file with no `popnei` key, a `format_version` of `2.0` with
   `2.0` in the message, a `gts` width of 7 with 3 individuals and a
   ploidy of 2, a `pos` column of `Int32`, two batches and one entry in
   `popnei_batches`, bytes that are not an arrow file, and two
   individuals of one name; a `format_version` of `1.7` and a seventh
   column, `depth`, are read; `from_path` on a directory is an error.
   `cargo test -p popnei --lib io::vars -- --list` names 26 tests or
   more, the 16 of work package 1 among them.
3. A vars file gives its blocks. Check: cargo tests made at `next_block`:
   the first and the third cargo tests of "How it is verified" of the
   writer, the round trip of the four variants of `cases.vcf` with every
   field and the one of `many.vcf` through the VCF reader, and the second
   and the fourth read now through `VarsReader`; a null position, with
   the column and the variant in the error; a batch that holds another
   number of variants than its entry of the footer; `zstd.vars`, which
   opens and gives the error of a file compressed with zstd at its first
   `next_block`; a file with no compression, which is read; a file of
   batches of 100, which gives blocks of 100; with `Needs` of the
   genotypes alone the blocks have no other column, and a change of
   `Needs` between two blocks holds from the next; after an error the
   reader gives `None` at every call; a file built in the test with a
   batch of no variants between two others gives the blocks of the two;
   and a test reads a `VarsReader` as
   a boxed `dyn BlockReader`. The same command names 38 tests or more.
4. `popnei.open_vars`. Check: `uv run pytest tests/test_io_vars.py -k
   open_vars` runs 5 tests or more and all pass: the two pytest tests of
   "How it is verified" of the reader, the first with
   `num_vars_per_block` of 7 and with the default and with the six
   counts of the spec as literals; a file that is not a vars file, a VCF,
   is a `ValueError` at the call of `open_vars`; a path that is not
   there is an `OSError` with the path in `filename`. The smoke test of
   pyodide reads back the file it wrote and finds the four variants of
   its table, and exits with 0.

**What it stands on.** Work package 1: the types of the two keys, the
error cases, and the writer, which makes the files of deliverables 3 and
4.

**Tasks.**

- [ ] 2.1 `tests/reference/vars/make_reference.py` and `zstd.vars`: the
  four variants of the `cases.vcf` table of `docs/specs/io_vcf.md`, three
  individuals, the six columns, both keys, one batch, zstd, written with
  pyarrow, the script with a header that says what it writes, with which
  pyarrow and how to run it, as `tests/reference/vcf/make_reference.py`
  does. From "How it is verified" of the reader and "The compression".
  Serves deliverable 1. It writes no code and no manifest, and runs
  beside 2.2.
- [ ] 2.2 `VarsReader::new`, `from_path`, `metadata`, `batches` and
  `num_vars`, in `crates/popnei/src/io/vars.rs`, the checks of the schema
  and of the two keys, and their tests. From "What it refuses", "How it
  runs" of the reader and "The Rust interface". Serves deliverable 2.
- [ ] 2.3 `VarsReader` as a `BlockReader`, in the same file: the
  projection that a `Needs` becomes, a batch into a block column by
  column, the nulls, the chromosome numbers in the order of first
  appearance, the error of zstd wrapped, the contract of a reader after an
  error of "What it gives" of `docs/specs/block.md`, and the tests. From
  "What it gives", "What it refuses" and "How it runs" of the reader, and
  "How it is verified" of the writer and of the reader. Serves deliverable
  3. Needs 2.1 and 2.2. A wrong genotype or a wrong chromosome number here
  is silent, so it is a commit of its own, and the round trip of
  `many.vcf` and the six counts of deliverable 4 are what guard it.
- [ ] 2.4 The Python side: in `crates/popnei-python`, a source that is a
  vars file beside the one that is a VCF, opened again at every pass,
  held as the binding holds a reader after `docs/plans/block-readers.md`,
  with the `Reblock` and the `Block::check` that every `iter_blocks` has;
  `open_vars` in `io_vars.py`; the `Variants` of the package over either
  source, and `write_vars` taking either; the tests; the read in the
  smoke test of pyodide. From "Its Python and TypeScript functions" of
  the reader, "What a Python and a TypeScript user see" of
  `docs/specs/variant.md` and "In Python and in TypeScript" of
  `docs/specs/block.md`. Serves deliverable 4. Needs 2.3.

**What could go wrong.** arrow-rs decompresses a batch when it reads it
and says `zstd IPC decompression requires the zstd feature`; how that
error is told from another one of arrow-rs, by its kind or by its text,
is for task 2.3 to find, and a test on `zstd.vars` holds it when arrow-rs
changes the text. The projection of arrow-rs takes column indices, and a
file may lack a column or have one that the reader does not know, so the
indices come from the schema of the file and not from the table of the
spec. pyNei is not run by the tests of this work package or of the one
before, because it reads another file: what stands in for it is pyarrow
and `many.bcftools.tsv` for the writer, and for the reader the blocks of
the VCF reader, which `docs/specs/block.md` checks against pyNei.

## Work package 3: `writeVars` and `openVars` in TypeScript

**What it gives.** A web application calls `writeVars(variants,
{numVarsPerBlock})` and gets the bytes of a vars file as a `Uint8Array`,
to offer as a download, and `openVars(bytes)` gives it a `Variants`.

**Deliverables.**

1. The two functions. Check: `npm run build` in `js/popnei` ends with no
   error of TypeScript and `npm test` gives the count of "What has to be
   in place" and 4 more or more, `fail 0`, with no test that was there
   changed: the TypeScript test of "How it is verified" of the writer,
   `cases.vcf` from a `Uint8Array` with `onlyPassed` false, written with
   `numVarsPerBlock` 3, read back with `openVars` and compared with the
   table of `cases.vcf`, two blocks of 3 variants and of 1 with every
   field; bytes that are not a vars file are an `Error` at `openVars`; a
   `source` that is not a `Uint8Array` and a `numVarsPerBlock` that is
   not a whole number of 1 or more are an `Error` that names the value;
   and `numberOfOpenPasses` of `js/popnei/src/variant.ts`, which counts
   the passes that still hold memory of wasm, gives 0 after a
   `writeVars`.
2. The size of `js/popnei/wasm/popnei_bg.wasm`, in the work report
   beside the one of "What has to be in place".

**What it stands on.** Work packages 1 and 2.

**Tasks.**

- [ ] 3.1 In `crates/popnei-js`, the source that is a vars file in
  memory and the function that writes one into a `Vec<u8>`; in
  `js/popnei`, `io_vars.ts` with the two functions, exported by both
  entry points, `Variants` over either source; the tests, in
  `js/popnei/test/vars.test.ts`. From "Its Python and TypeScript
  functions" of the writer and of the reader, section 11 of
  `docs/architecture.md`, and `js/popnei/README.md` for what was learned
  about wasm-bindgen. Serves deliverables 1 and 2.

**What could go wrong.** The review of `openVcf` found the VCF held twice
in the memory of wasm, 185.0 MB for a file of 91.9 MB, because the bytes
were copied for every pass. A vars file is opened again at every pass
too, so the passes share the bytes that `openVars` took and none copies
them; the reviewer of the binding measures it.

## Work package 4: the speed

**What it gives.** The owner knows how long popnei takes to read the
genotypes of a vars file and to write one, on the panel that "The
compression" and "Speed" measured arrow-rs on, before any work on speed.
It stops inside the core: a timing of the reader has no Python side,
and pyNei is not timed, because it reads another file.

**Deliverables.**

1. The bench and its file. Check: `cargo bench --bench vars_file --
   <path> --runs 5` prints one wall time for each run, of a pass over a
   vars file held in memory, through `next_block`, with the genotypes
   alone asked for and the genotypes summed, on one thread, and of the
   writing of its blocks, held in memory, into a `Vec<u8>`; the header of
   the bench says how the file is made: the VCF of 1000 individuals and
   20000 variants from `crates/popnei/benches/make_big_vcf.py`, which
   takes the number of variants as an argument, written as a vars file
   by popnei with the default size of block.
2. The measurement, in the work report. Check: it has the size of the file
   and its number of batches; the best of 5 runs of the read, as "The
   compression" took its numbers, and their median, beside the 21 ms of
   "Speed", which are a tenth above what arrow-rs took there to decompress
   the batches and copy the genotypes of each into a vector, 19.1 to 19.3
   ms, and says whether the number is met; the same two figures for the
   write, which has no number to reach; the machine, the build and the
   load; and, if the 21 ms are missed, a sampling profile of the read.

**What it stands on.** Work packages 1 and 2. The orchestrator runs it
with nothing else building on the machine.

**Tasks.**

- [ ] 4.1 The bench, with no harness, as `read_vcf` has none, the
  argument of `make_big_vcf.py`, and the measurement, as "What
  measurement there is" of the `performance-review` skill asks. It
  changes no code of the library. From "Speed" and "The compression".
  Serves deliverables 1 and 2.

**What could go wrong.** The file of the bench is not the file of the
spec: that one had the `gts` column alone and was written by pyarrow,
and this one has six columns and is written by popnei, from genotypes
that `make_big_vcf.py` draws the missing ones of in another order than
pyNei's script. The genotypes are of the same simulation and the read
asks for them alone, so the 21 ms stand as the number to reach, and the
report says what differs. If it is missed the plan is done when the
measurement is reported.

## How the whole plan is checked

From a clean clone of the branch: the five commands of the `coding`
skill; `cargo wasm-check`; `cargo tree -p popnei | grep -i zstd`, which
finds nothing; `npm run build` and `npm test` in `js/popnei`; the build
of the wheel of pyodide and its smoke test; and the work report, with
the counts of "What has to be in place" beside the ones at the end, the
two sizes of the wasm file and the measurement of work package 4.
