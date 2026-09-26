# Report: the VCF writer, the filter by regions, the missing rate and the density of the variants

26 September 2026. It records how `docs/plans/writer-regions-density.md`
was carried out, on the branch `plan/writer-regions-density` in the
worktree `.claude/worktrees/plan-writer-regions-density`, which branches
from `spec/writer-regions-density` at 838459f and not from `main`, because
neither the specs nor the plan are on `main` yet.

The plan is under way.

## Before the first task

The checks of the `coding` skill, run at 838459f on 26 September 2026, gave
what the plan says they gave on 7331638. `cargo fmt`, `cargo clippy`,
`cargo wasm-check`, `cargo wasm-check-js` and `ruff` passed.
`cargo test --workspace` gave 1022 passed with 2 ignored and 150 passed;
`cargo test -p popnei --no-default-features` gave 1022 passed.
`uv run maturin develop --release && uv run pytest` gave 550 passed and 6
skipped at the first run and 551 passed and 5 skipped at the second: the
tests of `tests/test_interrupt.py` skip themselves when the interrupt
arrives after the work it was meant to stop, which depends on timing, so
the count of skipped tests can move by one or more between runs.
`npm run build && npm test` in `js/popnei` gave 444 tests, 443 passed and
the one failure the plan names, in `test/gwas.test.ts`.

bcftools 1.24, tabix 1.24, bgzip 1.24 and plink2 v2.0.0-a.7.7 are in the
PATH, and `/Users/jose/devel/popnei-bench/big.vcf`, 403 MB, and `big.vars`,
81 MB, are there. The machine is an Apple M5 Pro with 18 cores.

The board held no message of another branch that is still at work: the
three of `main` announce merges, and the two of
`spec/writer-regions-density` are this plan's own spec and the failure of
the node test above. No file of this plan is being changed on another
branch.

## 1. The reader trait and the header of a source

Task 1.1, commit edcac3d. `SourceHeader` and the three methods of
`BlockReader`, `header`, `skip_outside` and `num_skipped`, are in every one
of the 24 implementations of the trait, where the plan counted 21: the
crate has 19, among them the ones for `Box` and `&mut`, and the benches 5.
The threshold filters and the filter by linkage disequilibrium refuse the
offer of regions; the filter of individuals, `reblock`, the reader one
block ahead, `Box` and `&mut` hand it to their source. The reader one
block ahead sends the offer to its reading thread and waits for the
answer, and gives first a block the thread built before it saw the offer.
`skip_outside` takes a `RegionSelection`, which the plan builds in task
4.1, so that type and an empty `Regions` came forward into
`crates/popnei/src/filters.rs`, with no constructor outside the tests;
task 4.1 fills them in. Four tests under `filters::tests::source_header`;
the subagent broke the code three ways and each broke them.

Tasks 1.2 and 1.3, commits f6393b4 and 3ba0021, 2c46528 and 343bc4f, and
f60c4e3. The VCF reader keeps the lines of its header before `#CHROM` and
the length of each `##contig` line, in the order of the file.
`tests/reference/vcf/write.vcf` is the file of "How it is verified" of the
writer, and gives chr1 2000 and chr2 1500 and its nine meta lines. A vars
file keeps the lengths in `chrom_lengths` of its `popnei` key and says
`format_version` 1.1; a file of 1.0 reads with no lengths.

Six cases the spec did not settle were settled with the code and written
into `docs/specs/io_vcf.md` and `docs/specs/io_vars.md`. A length is
written in digits alone, so `+5` is the wrong header, as `0` and `abc`
are. Two lines of one ID with the same length give it once. A comma
inside quotes, in a `Description="a, b"`, does not end a field. A
`##contig` line with a length and no ID is the wrong header. The
subagent had first made this last one give no length and no error; it was
sent back, because the VCF format requires the ID and an error in popnei
never passes silently. In a vars file, a `chrom_lengths` that is not a
list of names with whole numbers above 0, or that gives one chromosome
twice, makes the file not a vars file.

Deliverable 1 said that no test of before would change but for the
version. The longer key made a few more change, none of them an
assertion made weaker: the literals that hold the key whole gained
`"chrom_lengths": []`, and the sizes of the vars files in
`js/popnei/test/progress.test.ts` grew by 16 and 24 bytes. The plan says
so now.

Deliverables, run at f60c4e3:

- 1: `cargo test --workspace` gave 1035 passed, 2 ignored, and 150
  passed; `--no-default-features` 1035 passed; pytest 551 passed and 6
  skipped, which is the 551 of before, one test added, and one more skip
  of the interrupt tests by timing; `npm test` 443 of 444, the old
  failure. fmt, clippy, both wasm checks and ruff clean.
- 2 and 3: `cargo test -p popnei --lib source_header -- --list` counts 13
  tests, where it counted 0 before, and all 13 pass.

The review ran in six categories, spec, tests, errors, api, architecture
and numbers, over 838459f..f60c4e3; the binding category was left out
because neither binding crate changed. What it found that mattered, all of
it fixed in bba6e43 to 1d3a003 and one commit after them:

- The `##contig` lines were read by a parser stricter than htslib in some
  places and looser in others, and four reviewers found it from four
  sides. A space after a comma, a blank after the `>`, or an escaped quote
  `\"` inside a Description made a length go missing with no error, or
  gave a wrong one. When one line gave `length=` twice the last one won,
  where bcftools 1.24 keeps the first. The reader now trims the blanks as
  htslib does and takes `\"` as text. It refuses, as the wrong header with
  the line, a line with two `length=` or two `ID=` fields, an empty ID, no
  closing `>`, an unclosed quote, or no `<`. `docs/specs/io_vcf.md` records
  these as decisions of 26 September 2026.
- The test that a comma inside quotes does not end a field could not
  fail: its fixture had a space before the second `length=` and the real
  length last, and a reviewer who broke the quote handling saw all 1035
  tests pass. With the fixture changed, 4 tests fail on that breakage.
- Each new chromosome name was looked up in a list, so the time to open a
  file grew with the square of its `##contig` lines. A VCF with 200000 of
  them opened in 13.37 s; it opens in 0.030 s now, with the names in a
  hash, where `bcftools view -H` takes 0.086 s (release build, Apple M5
  Pro, the fastest of 3 runs). Draft plant assemblies have 1e5 to 1e6
  scaffolds.
- The reader one block ahead stopped serving the questions of its caller
  once its source had given its last block, so an offer of regions made
  at that moment was answered true or false by the timing of the threads:
  a reviewer counted 161 of 600 rounds one way and 439 the other. The
  thread now answers until its caller lets it go, and the test that its
  held block is not lost no longer depends on the scheduler.
- An offer of regions a source accepts holds for the life of that source.
  The trait says so now, and `docs/specs/block.md` says that the variants
  a source may skip are those the selection does not keep, which with
  `exclude` are the ones inside the regions. Every pass of popnei builds
  its chain with `chain_of`, which owns its source, so no pass can meet a
  source that skips after its filter is gone; a Rust caller that builds
  the filter by regions over a borrowed source and reads the source again
  can, and the doc is what warns them.
- Smaller: an answer arriving where the reader one block ahead waits for
  a block is now an error of popnei's own, a `RuntimeError` in Python; a
  wrong `chrom_lengths` names its pair and not the whole list; the error
  for a length past 18446744073709551615 says so; and three doc comments
  were put right.

Not taken: the test sources that give a header with no individuals, which
no test reads.

## 2. The missing rate

Task 2.1, commit ecfc6cf. The missing rate is the sixth statistic of
`calc_per_var_distribs`: the missing genotypes of a variant over its
called and missing ones, a half called genotype counted as missing, from
the counts the observed heterozygosity already makes. Six tests have
`per_var_missing_rate` in their names: the worked example of the spec with
and without populations and at `min_num_individuals` 20, the rates 0, 0.25,
0.5, 0.75 and 1 on the edges of 4 bins, 3 of 20 in bin 5 of 40, and the
table of plink2 for `many.vcf` over all and in popA and popB, whose
histograms match exactly and whose means match within 1e-12 relative.
Both binding crates only gained the field in their destructuring, so that
they compile; until task 2.2 `stats: ["missing_rate"]` in TypeScript is
taken and gives nothing. The spec's paragraph on `MAX_NUM_BINS` counted
four statistics with a histogram and 3.2 MB for each chunk of rows; with
the missing rate they are five and 4 MB, and the spec says so now.

Task 2.2, commit 833da72. `PerVarStat.MISSING_RATE` in Python and
`"missing_rate"` in TypeScript ask for the rate, and the results give it as
`missing_rate` and `missingRate`; it is among the six statistics computed
when none are named. `tests/reference/stats/make_reference.py` keeps
`many.vmiss`, `many.popA.vmiss` and `many.popB.vmiss` of plink2
v2.0.0-a.7.7 and checks the table of the spec against them.

Deliverables, run at 78a214f:

- 1: `cargo test -p popnei --lib per_var_missing_rate` gave 6 passed,
  where there were none.
- 2: the reference script ran and its check passed; the pytest test
  against the three files, over all and over popA and popB, passed.
- 3: the pytest comparison with pyNei leaves the missing rate out and
  passes; the test of one pass against one pass each covers six.
- 4: the TypeScript test of the means over all and in popA passed.
- All the checks: `cargo test --workspace` 1050 passed, 2 ignored, and
  150; pytest 555 passed, 5 skipped; `npm test` 444 of 445, the old
  failure; `npm run test:browser` 8 passed; fmt, clippy, both wasm checks
  and ruff clean.

The review ran in five categories, spec, tests, numbers, errors with api,
and binding, over ecfc6cf and 833da72; architecture was left out, since
nothing there touches the readers or the threads. No reviewer found a
wrong number: a sweep against `numpy.histogram` over populations of 1 to
61 individuals, 1 to 1000 bins and three ranges gave no difference, and
the mutants each reviewer tried in the core and in Python were caught.
What they found, all fixed in 752d12b to 9b0038e:

- The TypeScript tests could not see two mistakes a reviewer made on
  purpose: the missing rate dropped from the statistics computed when
  none are named, and the means of two populations swapped in the binding
  crate. Each now makes a test fail; the test of plink2's table asks for
  popA and popB in one call and asserts the whole histogram of each.
- A population with no individual left its variants out of the mean of
  the missing rate with no error. No reader gives such a population, but
  the Rust `Pops::all(0)` builds one; it is now an error of popnei's own,
  a `RuntimeError` in Python.
- A missing rate asked for with the string `"missing_rate"` got a hint to
  write `PerVarStat.MAF`; it names `PerVarStat.MISSING_RATE` now.
- Two Python tests, of the type of the counts and of the order of the
  populations, left the missing rate out; the README of the TypeScript
  package, the help of the benchmark of the pass and a docstring still
  spoke of five statistics, or of a mean that can be NaN; two lines of
  `docs/specs/stats.md` asked for the rate with a string the code refuses
  and said it had no code.
- `tests/reference/stats/make_reference.py` rewrote
  `tests/reference/stats/panel.vcf.gz` at every run with other bytes and
  the same content, which left the checkout changed; it leaves the file
  alone now. This was so before this plan.

Seen outside this plan and not fixed: `stats: [5]` in TypeScript is
refused with a message that says an `Array` was given, where a number was.
It gives no wrong result, and it can become a GitHub issue if the owner
wants one.

After the fixes, `cargo test --workspace` gave 1067 passed, 2 ignored, and
150; pytest 555 passed, 6 skipped; `npm test` 444 of 445, the old failure;
`npm run test:browser` 8 passed.

What the owner should know: the pass that `calc_per_var_distribs` makes
when no statistic is named now computes six, so a time taken of it from
now on is of one more statistic than the times in "Speed" of
`docs/specs/stats.md`. It was not measured again.

## 3. The VCF writer

Task 3.1, commits 9af8d7a and ab9666f. The VCF reader keeps the text of
each line when it is asked for `Needs::VCF_TEXT`, and `retain_vars`,
`retain_individuals` and `reblock` keep each line's text with its
variant. Sixteen tests have `vcf_text` in their names, on `write.vcf` and
on `many.vcf` in blocks of 7, and each checks every line's text against
the file and against the position and genotypes at its place in the
block. With the text, the reader also counts the columns of the
individuals and checks that they are UTF-8, and refuses a line longer
than 4294967295 bytes; the spec says so now. Reading all of `big.vcf`
takes 0.575 to 0.604 s on one thread and 0.083 to 0.085 s on 18, with and
without the change, at a load average of 5.

Task 3.2, commits 13bf74b and 19a7f80. `write_vcf` writes each line as the
source VCF had it, or builds it from the columns and the header when the
source is a vars file, and takes AC and AN out of INFO when the pass has
fewer individuals than its source; it finds them by the INFO key and by
the `ID` of the `##INFO` lines, as bcftools does, which the spec says now.
`tests/reference/vcf/make_reference.py` stores the outputs of bcftools
1.24 for the writer's cases: the filter of individuals, the rows of the
file written from the vars file, and the 215 lines of the missing data
filter at 0.04. Twelve tests compare whole files; a change of one byte in
four places of the writer made 7, 2, 2 and 1 of them fail.

Task 3.3, commits 7e31380 and b9c3a2b. `write_vcf` bgzips in members of
65280 bytes of text compressed on the threads of rayon, one after another
in wasm, with the empty member of 28 bytes at the end; a member whose
deflate does not fit in 65536 bytes is stored. Every case of the writer is
also written bgzipped and decompresses to the bytes of the plain case;
`bgzip -t`, `tabix -p vcf` and the query at chr1:900-1100 pass on each
file, and a test prints that it was skipped when a program is not in the
PATH. `vcf_text_num_vars_per_block` is a fifth of the genotypes divided
by the individuals, from 100 to 10000 variants.

### Waiting for the owner: the deflate of the bgzip

The plan asks the owner to choose when the bgzipped write of `big.vcf` on
one thread is more than a tenth above bcftools' 8.9 s, that is above
9.79 s. With flate2 over miniz_oxide at level 6, the deflate the core
had, it takes 10.6 s. The times below are of `write_vcf` on one thread,
from the VCF with its text, to a file, release build, Apple M5 Pro, three
runs in one process, on 26 September 2026, at a load average that fell
from 4.9 to 2.1 over the session; they are on the branch
`exp/writer-deflate`, commit 9dbd912, in its own worktree, with the
program in the session's scratch directory.

| deflate | level | seconds, three runs | MB | builds for both wasm targets |
|---|---|---|---|---|
| miniz_oxide, now | 5 | 5.02, 5.01, 5.00 | 41.6 | yes |
| miniz_oxide, now | 6 | 10.68, 10.63, 10.60 | 38.4 | yes |
| zlib-rs | 5 | 3.12, 3.11, 3.11 | 39.8 | yes |
| zlib-rs | 6 | 5.56, 5.55, 5.54 | 37.4 | yes |
| libdeflate | 6 | 4.49, 4.50, 4.48 | 39.8 | no |

In the same session `bcftools view -Oz` took 8.87 to 8.93 s and `bgzip
-@1` 7.41 to 7.44 s, both for 37.7 MB, and the plain write 0.98 to 1.00 s.
zlib-rs is pure Rust, `zlib-rs` 0.6.8 behind the `zlib-rs` feature of
flate2 1.1.10, and with it on flate2 decompresses with it too, so the
readers of bgzipped files change their inflate as well, which was not
timed. libdeflate, `libdeflater` 1.26.1, is C: it failed to build for
`wasm32-unknown-emscripten` without `emcc` and for
`wasm32-unknown-unknown` without C headers.

Task 3.4, commits fe5c9cb and b72ce55. `write_vcf(variants, path)` in
Python writes the file bgzipped when the path ends in `.gz`, and
`writeVcf(variants, {bgzip})` in TypeScript gives the bytes, bgzipped
unless `bgzip` is false. Python shares with `write_vars` one way of
handling the path: a path where a file already is is refused, a failed
write is an `OSError` naming the file, and the file is taken away on an
error or a Ctrl-C. `writeVcf` is the thirteenth consumer of the TypeScript
package, which `docs/specs/js_sources.md` counts now. At 59c7716 the checks
gave 1086 passed, 2 ignored, and 150 in cargo, 563 passed and 6 skipped in
pytest, 449 of 450 in node with the old failure, and 8 in the browser.

The review ran in all seven categories over the eight commits of tasks 3.1
to 3.4. No reviewer found a wrong line in a file the tests write: one
compared a filter of individuals followed by the missing data filter on a
file of 3500 variants and 1200 individuals with the same steps in bcftools
and got the same 857 lines byte for byte, in Python and in TypeScript.
What they found, and what is being fixed:

- Two ways to write a file bcftools would not. When AC and AN come out,
  empty INFO values were kept, so `AC=1;AN=4;` became an empty INFO
  column where bcftools writes `.`. And a vars file without the `id`
  column, which `docs/specs/io_vars.md` calls valid and `write_vars`
  writes, was refused by the writer as a defect of popnei; the specs
  contradicted each other there, and the owner's convention for errors
  settles it: a missing ID or QUAL is written `.`, and a missing
  position, chromosome, alleles or genotypes is a `ValueError` naming
  the file and the column.
- A panic during a write left a partial file at the path, for `write_vars`
  as well, which was so before this plan.
- Tests that could not fail, each shown by a mutation that passed the
  whole suite: the order of the bgzip members, since no test wrote more
  than one member at once; ploidies other than 2 written from a vars
  file; the empty member at the end, compared with the code's own
  constant; and the member stored when deflate cannot shrink its text,
  a branch no test reached because miniz_oxide stores such text itself.
- Speed and memory. The writer asked a VCF source for every column and
  then read only the text: the plain write of `big.vcf` took 0.95 s on
  one thread where the text alone takes 0.52 s. The filter of
  individuals rewrote the text on one thread while the genotypes beside
  it were gathered on all, which doubled a write with that filter on 18
  threads, 0.17 s to 0.33 s. A bgzipped write held each block's text
  three times, 147 MB at its peak on a file whose reading with its text
  took 56 MB.
- Smaller: an error named `VarsFileNotWritten` that the VCF writer gives
  too; `.gz` matched only in lower case, so `a.VCF.GZ` was written plain;
  the choice of the writer's block size made in both binding crates
  rather than in the core; paths that would give an empty text or a
  default value instead of an error if a check before them moved.

The fixes are the 33 commits from 5cd7d40 to 586828e. Each finding that
changed a behaviour has a test that failed before its fix. Not fixed as
asked:

- The removal of a partial file after a panic, and the mapping of the
  writer's two defects to `RuntimeError`, have no automated test: no
  input makes the core panic or reach those defects, and a test would
  need an entry point of `_core` for tests only. The first was checked by
  hand with a panic put into the code: two calls in a row each raised a
  `PanicException` and left no file.
- The block size of the writer is now chosen by one function of the core,
  tested there; no test of either binding counts the blocks, because the
  progress of a pass counts bytes.
- The comment that said one deflate per thread says one per job of rayon
  now; nothing else changed there.

Peak memory of a bgzipped write fell from 366 MB to 341 MB on 18 threads,
on a VCF of 245 MB, 2000 lines of 5000 individuals with GT, AD, DP, GQ and
PL. The plain write of that file peaks at 314 MB, which nobody has looked
into yet. The filter of individuals now writes its text on the threads:
keeping 500 of the 1000 individuals of `big.vcf` on 18 threads, the
write went from 0.336 to 0.379 s to 0.268 to 0.286 s.

The subagent also made the writer ask a VCF source for the text of its
lines alone. That took the plain write of `big.vcf` from 0.96 to 0.98 s
down to 0.52 to 0.54 s on one thread, but a line whose POS is not a
number, or whose genotype does not parse, was then copied into the output
without an error. The owner's rule is that a corrupt input is reported
even when the check costs time, so that change was undone in db0d83b and
d03d89a. Three tests now hold it: `0/x`, three alleles in a diploid VCF
and a POS of `nine` are each refused. After the undo, on one thread, the
plain write takes 0.98 to 1.00 s and the bgzipped one 10.8 to 11.0 s. The
owner may still prefer the faster write; it is among the questions at the
end.
