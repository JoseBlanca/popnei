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
