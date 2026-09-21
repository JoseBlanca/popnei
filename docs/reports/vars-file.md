# Work report: the vars file, its writer and its reader

The plan `docs/plans/vars-file.md` is under way, on the branch
`plan/vars-file`, in the worktree `.claude/worktrees/vars-file`, since 21
September 2026. The orchestrator, in this report, is the session of the
assistant that runs the plan: it sends each task to a subagent on Opus,
checks what comes back and has each work package reviewed.

## Before the first task

The owner answered the breakdown of the plan on 20 September 2026 and
ordered the plan run as soon as `docs/plans/block-readers.md` was merged
into `main`. The session that ran that plan sent word of the merge on 21
September 2026, and the orchestrator checked it on `main` at 2d99d64 with
the commands of "What has to be in place": `git branch --merged main`
lists `plan/block-readers`; `docs/reports/block-readers.md` says that
plan is done; the search for `read_variant`, `VariantReader` and
`BlockCollector` in the crates finds nothing; `BlockReader`, `Reblock`,
`check`, `fields` and `retain_vars` are in `crates/popnei/src/block.rs`;
and no manifest holds the owner's absolute path to pyNei. The worktree
was made after that, from 2d99d64.

In the worktree, run by the orchestrator: `cargo fmt --all --check` exit
0; `cargo clippy --workspace --all-targets -- -D warnings` no warning;
`cargo test --workspace` `173 passed`, 1 ignored; `cargo wasm-check`
finished; ruff `13 files already formatted` and `All checks passed!`;
`uv run maturin develop && uv run pytest` `72 passed`; `npm run build`
with no error and `npm test` `tests 46`, `fail 0`;
`js/popnei/wasm/popnei_bg.wasm` is 225344 bytes; `bash
scripts/build_pyodide_wheel.sh` built the wheel and `node
tests/pyodide/smoke.mjs` exited with 0. These are the counts that every
"and more" of the plan is counted from.

Changes to the plan before the first task, from what the other plan left
on `main`:

- The owner gave on 21 September 2026 a convention for the exceptions of
  Python, which "Errors, and no panics" of the `coding` skill has, and it
  asks the spec of a module to list its cases with the exception each
  one is. `docs/specs/io_vars.md` says that every case is a `ValueError`
  but the error of the input or the output. Task 1.1 writes the list
  into the spec, and "In and out" of the plan has how the orchestrator
  sorted the cases: the content of a file that is not what the format
  holds is a `ValueError`, a file cut short or a batch that cannot be
  decoded an `OSError`, and a block that does not fit the writer, a
  defect of a reader, a `RuntimeError`.
- A batch of no variants, which the spec does not mention, is skipped by
  the reader. The owner was asked on 21 September 2026 and had not
  answered when the work started.

### For the owner

- The two changes above are his to reverse.
- His rule of 21 September 2026 that an error never passes silently, and
  that a file that was cut or damaged is an error however improbable the
  damage, asks more of the vars file than its spec gives. An arrow IPC
  file has no checksum, and arrow-rs writes its lz4 frames without the
  optional checksums of lz4, so a byte changed inside the compressed
  genotypes can decompress into other genotypes with no error. What
  would close it is a checksum of each batch in `popnei_batches`, which
  is a change of the format and so of the spec. The plan goes on without
  it: the key is json, and a field added to its entries does not change
  what is built before. The reviews of work package 2 will say what a
  sweep of changed bytes finds.

## Work package 1: `write_vars` in Python

Task 1.1, commits 5905926, the spec, and 43deb76, the code; one subagent
run of 203 thousand tokens and 16 minutes. The orchestrator ran the
checks again: fmt exit 0, clippy no warning, `cargo test --workspace`
`183 passed`, 1 ignored, `cargo test -p popnei --lib io::vars -- --list`
`10 tests`, where the plan asks for 6, `cargo wasm-check` finished, and
`cargo tree -p popnei | grep -i zstd` finds nothing. The core has four
crates of arrow-rs 60, `arrow-array`, `arrow-buffer`, `arrow-schema` and
`arrow-ipc`, all with their default features off and `arrow-ipc` with
`lz4`, and `serde_json` with `std` alone; the subagent reports 59 crates
more in `Cargo.lock`. The spec lists the cases of the error as nine
`ValueError`, three `OSError` and three `RuntimeError`, two of them new,
a file cut short and a batch that arrow-rs cannot decode. All 14 cases
are in the error enum, so the later tasks add none, and the Python
binding crate names the five that are not a `ValueError`. The subagent
also built the wheel of pyodide, which linked; nothing of the binding
calls arrow-rs yet, so task 1.3 is still the first link that counts.

Task 1.2, commits a81a379, the spec, and 1965bb8, the code; one subagent
run of 281 thousand tokens and 18 minutes. The orchestrator ran the
checks again: fmt exit 0, clippy no warning, `cargo test --workspace`
`195 passed`, 1 ignored, `cargo test -p popnei --lib io::vars -- --list`
`22 tests`, where the plan asks for 16, `cargo wasm-check` finished. What
the subagent decided that the task did not say:

- The spec did not say what the values inside the two list columns,
  `alleles` and `gts`, are called, which arrow keeps as a field of its
  own. "What it holds" now says that popnei writes what pyarrow writes
  for any list, a field named `item` that can be null, so that the file
  of task 2.1, which pyarrow writes, has the same columns.
- A `gts` width above 2147483647, the most a fixed size list of arrow
  holds, 2 GB of genotypes for one variant, is the error of a block that
  is too large, which the crate had. No case was added for it.
- The test of the genotypes that are not copied compares addresses at
  the private function that builds the `gts` column, because once the
  batch is written the buffer is inside arrow-rs.

Task 1.3, commits 352f950, the spec, and 5ec27ee, the code; one subagent
run of 214 thousand tokens and 13 minutes. `popnei.write_vars` makes its
file with a call that fails when the path is taken, runs the core with
the interpreter released, and removes the file on an error and on a
Ctrl-C. pyarrow is in the `dev` group with 23 as its lowest version, and
`uv.lock` changed in the two lines that make it a direct dependency.

### The deliverables, run by the orchestrator at 5ec27ee

1. `cargo wasm-check` finished, `cargo tree -p popnei | grep -i zstd`
   finds nothing, and the spec has the sentence about pyarrow.
2. `cargo test -p popnei --lib io::vars -- --list` `22 tests`, where the
   plan asks for 16; `cargo test --workspace` `195 passed`, 1 ignored.
3. `uv run maturin develop && uv run pytest` `80 passed`, 72 before, and
   `uv run pytest tests/test_io_vars.py -k write_vars` `8 passed`, where
   the plan asks for 6; no test that was there changed.
4. `bash scripts/build_pyodide_wheel.sh` built the wheel, and `node
   tests/pyodide/smoke.mjs` printed `cases.vcf written as a vars file of
   3338 bytes, with "ARROW1" at both ends` and exited with 0. arrow-rs
   links under emscripten, which was the one thing of the plan that had
   not been tried.

### The review

Six reviewers over the six commits of the work package, at 5ec27ee:
`spec` and `tests`, each in a worktree of its own, and `errors`, `api`,
`architecture` and `binding`; 149, 130, 115, 95, 149 and 106 thousand
tokens. The `tests` reviewer made 14 mutations of the code and ran the
tests on each. The `spec` reviewer recomputed every literal of the pytest
tests from `many.bcftools.tsv` and opened a file of popnei with pyarrow
25.0.1, pandas and `feather.read_table`: the six columns in the order of
the spec, five batches, both keys, 167 and 100 nulls, lz4.

What held and was fixed, the ones that matter first:

- A write that failed was told as a read of the VCF. With the size of
  the output limited, so that only the write could fail, a user got
  `OSError: [Errno 27] the file could not be read: File too large:
  'many.vcf'`. Four reviewers ran it. The subagent of task 1.3 had
  written into the spec that every error names the VCF, because the core
  did not say which file it was at, and the orchestrator reversed that:
  the owner's rule is that an error of a file names that file. The core
  has a case of its own for a vars file that could not be written, and
  in Python it is an `OSError` with the path of the vars file.
- One batch with more than 2 GiB of text in its alleles, its ids or its
  chromosome names panicked inside arrow-rs, whose text columns count
  their bytes in 32 bits, and the panic left a half written file. A
  sequence resolved VCF of structural variants could reach it, and a
  `num_vars_per_block` that a user gives has no bound. It is an error
  now, before anything of the block is written, and it says to write
  with fewer variants in a block.
- Two errors passed in silence: a partial file that could not be
  removed, and a close of the file that failed, which is where a network
  file system says that it is full.
- `VarsWriter::new` took a ploidy and a size of block of 0 and wrote a
  file that popnei's own parser of the key refuses, and a source with no
  individual gave a file with no `gts` column, which the spec says every
  file has. The three are refused.
- A directory as the path was "a file is already there", and a VCF path
  given where the `Variants` goes was an `AttributeError` about
  `_source`.
- Tests that guarded nothing: sorting the regions by the name of the
  chromosome passed every test, because every fixture had chr1 first;
  the two checks of a Ctrl-C could be deleted; so could the request for
  every field, because the VCF reader asks for every field by default;
  the test of the genotypes that are not copied did not see a copy made
  one function above; no test gave the writer a sink that fails, and
  none the path in a directory that is not there.
- Smaller: a public constant that the spec has not; the comment of the
  manifest said that the tree of wasm is pure Rust, and under emscripten
  chrono brings `iana-time-zone`, which calls the runtime of emscripten,
  though nothing in the path of a vars file asks for a clock; section 6
  of the architecture had `qual` before `alleles`.

What held and was not changed:

- `num_vars_per_block=True` is read as 1, as in `iter_blocks`: a bool of
  Python is an integer.
- The json of the two keys is built by joining texts, with `serde_json`
  for the escapes. Names with quotes and backslashes are written and
  parsed back in the tests.
- The message of 1965bb8 says 6 tests of the two keys and 16 of the
  writer; they are 10 and 12. The commit is not rewritten.

The findings went back to the two subagents that wrote each side, the
core first, since a core that is half written breaks the build of the
binding: 8 commits from 9efa6e0, the spec first, 99 thousand tokens
more, and 8 commits from df741a9, the spec first, 85 thousand more. All
18 findings that were sent held for the subagents too. After them, run
by the orchestrator: fmt exit 0, clippy no warning, `cargo test
--workspace` `200 passed`, 1 ignored, `io::vars` `27 tests`, `cargo
wasm-check` finished, no zstd crate, ruff clean, `uv run pytest` `85
passed` and `-k write_vars` `13 passed`, the wheel of pyodide built and
its smoke test exited with 0.

Three cases were added to the error of the crate, which task 1.1 was
meant to leave complete: a vars file that could not be written, an
`OSError` with the path of that file; a block with more text than a
column of an arrow file takes, and a source with no individual or a
ploidy of 0, both a `ValueError`. The spec lists them.

### For the owner

Decisions of this work package that a user sees, each written in the
spec with the option not taken, for the owner to reverse:

- A Ctrl-C during `write_vars` is raised when the pass over the source
  is over, since the whole file is written inside one call of the core,
  and it removes the file, also a file that the call had finished. The
  option not taken keeps a finished file. One that is pending when the
  call starts is raised before any file is made.
- `write_vars` returns when the bytes are on the disc: it calls
  `sync_all`, so that a file system that reports its error at the close
  does not give a file that looks whole. What that costs in time is
  measured in work package 4.
- When the file of a call that failed cannot be removed, the exception
  keeps the first error and gets a note, in `__notes__`, that says a
  file was left at the path.
- A directory at the path is an `IsADirectoryError`, and a first
  argument that is not a `Variants` a `TypeError`.

### How the work went

Each task went to one subagent and came back right the first time. The
reviews cost 744 thousand tokens, the three tasks 698 thousand and the
fixes 184 thousand. The subagent of task 1.3 wrote a decision into the
spec, that an error of the write names the VCF, which the review
reversed: a subagent that finds the core lacking something a rule of the
owner asks for should say so and not write the way around it into the
spec. The prompts of the next tasks say it.

## Work package 2: `open_vars` in Python

Tasks 2.1 and 2.2 ran side by side, as the plan allows, with two
subagents that shared no file.

Task 2.1, commit a0a37b7, 92 thousand tokens and 3 minutes:
`tests/reference/vars/make_reference.py` and `zstd.vars`, 2418 bytes. The
orchestrator ran the script again, which changed no file, and opened
the file with pyarrow: one batch, four variants, both keys, and the mark
of a zstd frame 12 times, once for each buffer. The subagent compared it
with the file popnei writes from `cases.vcf`: the same columns, the same
values and the same `popnei_batches`, and `num_vars_per_block` 4 where
popnei's default says 10000.

Task 2.2, commit fcc585a, 254 thousand tokens and 25 minutes:
`VarsReader::new`, `from_path`, `metadata`, `batches` and `num_vars`. The
orchestrator ran the checks again: fmt exit 0, clippy no warning, `cargo
test --workspace` `217 passed`, 1 ignored, `io::vars` `44 tests`, where
the plan asks for 26, `cargo wasm-check` finished. The reader parses the
footer itself and will decode each batch with the `FileDecoder` of
arrow-rs, because its `FileReader` fixes the columns it reads when it is
built and a `Needs` can change between two blocks. A file whose footer
is not all there is a file cut short, an `OSError` in Python; bytes that
were never an arrow file are not a vars file, a `ValueError`. It left a
file with no `gts` column as one that is read; the orchestrator decided
that such a file is refused as not a vars file, since the spec says that
`gts` is always there, and task 2.3 does it.

Task 2.3, commits 1763dad, the spec, and 6fb7007, the code; one subagent
run of 297 thousand tokens and 18 minutes. `VarsReader` gives each batch
as a block. The orchestrator ran the checks again: fmt exit 0, clippy no
warning, `cargo test --workspace` `238 passed`, 1 ignored, `io::vars`
`65 tests`, where the plan asks for 38, `cargo wasm-check` finished. What
the subagent found and decided:

- arrow-rs gives the error of a zstd file as an invalid argument whose
  text holds `zstd`, and the reader takes that kind with that word for
  the case of a file compressed with zstd; every other error of arrow-rs
  while a batch is read is a batch that could not be read. The test on
  `zstd.vars` fails if arrow-rs changes either.
- That only the columns asked for reach arrow-rs is seen on `zstd.vars`:
  asked for no field it gives its 4 variants and no error, and asked for
  the genotypes it gives the error of zstd.
- arrow-rs panics on a batch whose message is shorter than 8 bytes,
  which it indexes, so the reader refuses such a batch before arrow-rs
  sees it. The subagent said no test could reach it; the review showed
  one that does, by changing the length of a batch in the footer of a
  good file, and it is there now.
- The spec got three sentences: a file with no `gts` column is not a
  vars file, a null inside a list of alleles or of genotypes is a null,
  and the block that the machine has no memory for, the case of
  `docs/specs/block.md`.

Task 2.4, commit 699a6a7; one subagent run of 214 thousand tokens and 12
minutes. `popnei.open_vars` gives the same `Variants` as `open_vcf`. The
two sources share one trait in the binding crate, the file a source
reads and the reader of one pass over it, and everything after that is
one copy: the fields, the size, the `Reblock`, `Block::check` and the
columns that go to numpy. `write_vars` takes a `Variants` of either
source.

### The deliverables, run by the orchestrator at 699a6a7

1. `uv run python tests/reference/vars/make_reference.py` changed no
   file; pyarrow opens `zstd.vars` with one batch, four variants, both
   keys and 12 zstd frames.
2. and 3. `cargo test -p popnei --lib io::vars -- --list` `65 tests`,
   where the plan asks for 38; `cargo test --workspace` `238 passed`, 1
   ignored; `cargo wasm-check` finished.
4. `uv run maturin develop && uv run pytest` `94 passed`, and `-k
   open_vars` `9 passed`, where the plan asks for 5; the wheel of pyodide
   built, and its smoke test printed `cases.vars read back: 4 variants,
   as the spec says` and exited with 0.

### The review

Six reviewers at 699a6a7: `spec`, `tests` and `errors`, each in a
worktree of its own, and `numbers`, `architecture` and `binding`; 196,
199, 183, 151, 219 and 183 thousand tokens. `numbers` also read the
writer of work package 1, which had not had that category. The `tests`
reviewer made 40 mutations of the reader and the binding, of which 30
were caught. The `spec` reviewer recomputed the six counts of the pytest
test and the fifteen regions of the spec from `many.bcftools.tsv`, and
wrote with pyarrow some thirty files that popnei did not write, among
them columns in another order, a column more, a column less, another
version, an empty batch and nulls in every place, and read each.

What a damaged file does, which the owner's rule of 21 September 2026
asks about. The `errors` reviewer changed single bytes of two vars files
that popnei wrote and read each with every field, 208 s and 54 s of
work, and both counts came out the same on a second run:

| | the 4 variants of `cases.vcf`, 3338 bytes, each byte to all 255 other values | the 500 of `many.vcf`, 49426 bytes, each byte to 4 other values |
|---|---|---|
| files | 851190 | 176030 |
| an error | 295748 | 104420 |
| the same blocks | 458465 | 18250 |
| other blocks and no error | 23880 | 48849 |
| a panic inside arrow-rs | 70243 | 3916 |
| the process killed | 2854 | 595 |

Cut at every length, both files gave an error every time. The files
that gave other blocks with no error are 28 in 100 of the changes of
the larger file, whose buffers are compressed: lz4 as arrow-rs writes it
has no checksum, and one byte changed in the genotypes of the first
batch turned 25 genotypes from 1 into -127. One byte of the footer's
copy of the schema, the `p` of `pos`, made a file give no chromosomes
and no positions, since a column whose name the reader does not know is
ignored. The `binding` reviewer found the same from Python, with 510
files. The panics and the killed processes are arrow-rs allocating and
slicing by numbers it read from the file.

What held and was fixed, the ones that matter first:

- A damaged file reached panics and allocations inside arrow-rs, which
  slices and allocates by numbers it reads from the file. The reader now
  parses the message of a batch before arrow-rs sees it and checks the
  rows against the entry of the footer, every buffer against the body,
  and the size that every lz4 buffer says it decompresses to against
  255 bytes for each of its bytes; it asks for the memory of a batch
  with `try_reserve`; and natively it catches a panic of arrow-rs and
  gives the error of a batch that could not be read. After it, the sweep
  of every byte of the file of `cases.vcf` to every other value, which
  is a test that is run on request, 1299990 files in 8.2 s: 550055 an
  error, 726033 the same blocks, 23902 other blocks and no error, 2783
  panics caught, and no process killed. The 2783 are two asserts of
  arrow-rs that only its own walk of the schema would see, and in wasm
  nothing catches them. A sweep of four values of each byte, 16296
  files in 0.13 s, runs with every `cargo test`.
- A `popnei` key that names no individual was read, and gave blocks of
  variants with no genotypes that passed `Block::check`. Refused.
- An infinite quality, and a NaN that is a value, got into a block from
  a file of another program. The VCF reader refuses both, and the rule
  that NaN is a variant with no quality rests on that. Refused, by the
  orchestrator's decision, written in the spec.
- The writer marked the values inside the `gts` and `alleles` lists as
  ones that can be null, and arrow-rs then writes a bitmap of all ones,
  a bit for each allele, that the reader decompresses and throws away.
  On the panel of the spec, the genotypes alone, the best of 10 passes
  on one thread: 20.52 ms and 15691298 bytes with it, 19.49 ms and
  15680146 bytes without. The writer no longer marks them so; the
  reader takes either, so the file of pyarrow still reads.
- A null in a column that the file says has none is refused by arrow-rs
  before popnei looks, as a batch that could not be read, an `OSError`,
  and not the `ValueError` with the variant that the spec promised. The
  spec says it now.
- The message of a batch too large for the machine said to ask for
  fewer variants in a block, which does nothing for a file: it says to
  write the file again with fewer.
- Tests that were missing, each shown by a mutation that left every
  test passing: the batch shorter than 8 bytes; a null list of alleles
  and of genotypes, without whose check the reader gave genotypes read
  out of padding; the number of the variant and of the batch of an
  error beyond the first batch; a second pass over one `open_vars`
  handle; a damaged batch as an `OSError`; columns in another order.
- In the binding: the errors about the file being written that named
  the file that was read; the docstring that said a `Variants` cannot be
  copied, which `copy.copy` can.

What held and was not changed:

- A `num_vars_per_block` of 2^40 is refused for a VCF, whose reader
  would reserve that memory, and taken for a vars file, which never
  needs it.
- A pass for the genotypes alone reads every byte of the file, 16272664
  of 16280410, and decompresses only the genotypes, which is what the
  spec promises. In a tab, with a `File` as the source, the whole file
  would cross for every pass; that is for the plan that does the `File`.
- Nothing checks that an allele of `gts` is -1 or an allele of its
  variant, so a damaged file can give an allele 7 in a variant of two.
  It costs a pass over the genotypes, and the checksum comes first.

After the fixes, 18 commits by the subagent of the reader, the spec
before the code, 155 thousand tokens more, and 7 by the one of the
Python side, 77 thousand: `cargo test --workspace` `249 passed`, 2
ignored, `io::vars` `77 tests`, `uv run pytest` `99 passed`, fmt, clippy,
ruff and `cargo wasm-check` clean, the wheel of pyodide built and its
smoke test exited with 0.

A first timing, by the `architecture` reviewer, before the fix of the
bitmap: a pass with the genotypes alone over the panel of the spec, four
batches of 5000, release, one thread, the best of 10 runs: 21.0 to 21.2
ms, against the 21 ms of the spec. Work package 4 measures it.

### For the owner

- The checksum. 23902 of 1299990 files with one byte changed, and 28 in
  100 of those of a file whose buffers are compressed, give other
  variants with no error. The owner was asked on 21 September 2026 with
  three options: a CRC32 of each batch in `popnei_batches` and one over
  the two keys and the columns, checked before arrow-rs sees the bytes,
  which the orchestrator recommends; the same but optional, so that a
  file of another program is read without it; or none, with the spec
  saying that a vars file is not protected. It changes the format, so
  nothing was built.
- Two decisions of this work package written in the spec for him to
  reverse: a quality that is not finite is refused, and the values
  inside the two lists are written as ones that cannot be null.

### How the work went

The reviews cost 1131 thousand tokens, the four tasks 857 thousand and
the fixes 232 thousand. No task was sent twice.

## Work package 3: `writeVars` and `openVars` in TypeScript

Task 3.1, commit 1ed9a5a; one subagent run of 231 thousand tokens and 16
minutes. `openVars(source)` and `writeVars(variants, {numVarsPerBlock})`
are exported by both entry points of the package. In the JavaScript
binding crate the two sources share a trait, the reader of one pass, and
the blocks code that was in the VCF file is one copy for both. No
manifest changed.

### The deliverables, run by the orchestrator at 1ed9a5a

1. `npm run build` in `js/popnei` ended with no error of TypeScript, and
   `npm test` gave `tests 62`, `pass 62`, `fail 0`, 46 before, none
   changed. fmt, clippy, `cargo test --workspace` `249 passed`, `cargo
   wasm-check` and `uv run pytest` `99 passed` as before.
2. `js/popnei/wasm/popnei_bg.wasm` is 1679488 bytes, and was 225344
   before arrow-rs: 1.45 MB more, release, with no work on its size. The
   trial crate of the spec gave 2.64 MB for arrow-rs with lz4 alone.

The memory of wasm, measured by the subagent on a vars file of 26.9 MB,
20000 variants of 1000 individuals in batches of 1000: `openVars` takes
27.1 MB, one copy of the bytes; a pass over the genotypes 12.2 MB more,
the batch being decoded; a second pass nothing. A test opens twelve
passes at once and fails when a pass copies the bytes. `writeVars` of
that file from its VCF of 76.9 MB took 106.8 MB beyond the source with
batches of 1000 and 158.5 MB with batches of 10000.

### The review

Four reviewers at 1ed9a5a: `spec` and `tests`, each in a worktree of its
own, `errors`, in one too, and `binding` with the memory; 174, 116, 162
and 151 thousand tokens. The owner stopped the last two while they ran,
on 21 September 2026, and ordered them sent again. The `spec` reviewer
wrote `many.vcf` with `writeVars` under node and opened the file with
pyarrow 25: the six columns, five batches of 100, 167 and 100 nulls and
the five entries of the spec's table of regions; the file of `cases.vcf`
is byte for byte the one that Python's `write_vars` makes; and a file
that Python wrote, read with `openVars`, gives the 500 variants of
`openVcf` in every field.

What a damaged file does in wasm, where a panic is a trap, a
`RuntimeError: unreachable` with no message of popnei, and nothing
catches it. The `errors` reviewer wrote `cases.vcf` in batches of 3,
5098 bytes, set every byte in turn to four values and read each of the
16296 damaged files with every field, 2 s under node, starting the
process again after each trap:

| what it gave | files |
|---|---|
| an `Error` at `openVars` | 2775 |
| an `Error` at a block | 4136 |
| the same four variants | 8988 |
| other variants, or a column gone, and no error | 366 |
| a trap | 31 |

The 31 are the files that natively reach the two asserts of arrow-rs,
where the core catches them and gives the error of a batch that could
not be read. After a trap the module still works, and the memory the
pass held never comes back: about two copies of the file for each trap.
The binding crate cannot catch a trap, so only keeping damaged bytes
from arrow-rs helps, which is what the checksum would do.

What held and was fixed, the ones that matter first:

- A size that a damaged batch declared could take the memory of a tab.
  The bound the core had put on the decompressed size of an lz4 buffer,
  255 bytes for each compressed byte, is above what wasm32 addresses
  for any batch over 16 MB: on a file of one batch of 25 MB, the eight
  bytes that say 40000000 for the genotypes, changed to 2000000000,
  left the memory of wasm at 2.07 GB for the life of the tab, and
  changed to 3000000000 trapped. The reader now walks the buffers of a
  batch in the order of the schema and holds each to what its column
  can hold, the genotypes to rows x individuals x ploidy bytes exactly.
- A `Uint8Array` that did not fit in the memory left trapped, in
  `openVars` and in `openVcf`, because the glue of wasm-bindgen asks for
  the memory before any code of popnei runs. The package now asks the
  binding whether there is room, and gives an `Error`.
- `writeVars` took 3.4 times the file it wrote in the memory of wasm,
  which never shrinks, because its buffer doubled: 62.4 MB for a file of
  18.3 MB. It collects the file in pieces of 1 MiB that cross one at a
  time and are put together in JavaScript: 34.2 MB. And
  `numVarsPerBlock` did not reach the VCF reader, which built blocks of
  the default size whatever was asked, so the advice of the README to
  lower it did nothing: written from the VCF of 77 MB in batches of 100,
  51.3 MB before and 20.8 MB after, 1.12 times the file.
- Four tests that could not fail: the one that guards against a pass
  copying the bytes of the file, which measured 0 bytes of growth with
  and without the copy because the memory of wasm had room left by the
  tests before it, and now runs in a process of its own, 0 bytes against
  36765696; that `writeVars` gives a copy and not a view into the memory
  of wasm; the ploidy, only ever read from a diploid file; the size of
  the batches written. `npm test` ran a `dist/` that could be older than
  `src/`, and builds the TypeScript first now.
- A detached `Uint8Array`, whose buffer a page had given to a worker,
  was a `TypeError` of the generated code.
- The wasm file carried 443 KB of names of functions for a debugger. It
  is built without them, 1242562 bytes where it was 1687938; a trap then
  shows no names in its stack.
- The numbers of memory that the subagent had put in the README did not
  reproduce for two reviewers, and are measured again, each with its
  size of block and how it was measured: on a VCF of 80692954 bytes from
  `make_big_vcf.py` with 20000 variants, whose vars file in batches of
  1000 is 19185674 bytes, `openVars` takes 18.4 MB, the first pass 11.7
  MB with blocks of 1000, 39.2 MB with no size given, which is 5000, and
  62.6 MB with 10000, and a second pass nothing.

What held and was not changed:

- The 31 damaged files that trap. The binding cannot catch a trap.
- A file written in wasm is larger than the same one written natively,
  53650 bytes against 49426 for `many.vcf` in batches of 100, because
  `lz4_flex` hashes four bytes on a 32 bit target and five on a 64 bit
  one. Both are valid and hold the same table. The README says it.
- The message of 1ed9a5a has the numbers that did not reproduce. The
  commit is not rewritten.

After the fixes, 2 commits of the core, 40 thousand tokens more, and 11
of the TypeScript side, 128 thousand: `cargo test --workspace` `250
passed`, 2 ignored, `io::vars` `78 tests`, `npm test` `tests 70`, `fail
0`, `uv run pytest` `99 passed`, fmt, clippy and `cargo wasm-check`
clean.

### For the owner

- The 31 traps and the 366 files read wrong with no error, of 16296, are
  one more reason for the checksum: in wasm nothing else keeps a damaged
  batch from arrow-rs.
- The wasm file is built without the names of its functions, 445 KB
  less. With them the stack of a trap names the function it was in. It
  is one flag of `build:wasm` in `js/popnei/package.json` to reverse.

### How the work went

The reviews cost 603 thousand tokens, the task 231 thousand and the
fixes 168 thousand. The owner stopped two reviewers while they ran, and
they were sent again when he said to go on. The subagent of the task
wrote its tests after the code, and four of them could not fail; its
numbers of memory were measured in a way it did not write down and did
not reproduce. The prompt of a task that measures should ask for the
procedure beside each number.
