# Report: the VCF writer, the filter by regions, the missing rate and the density of the variants

26 September 2026, finished on 27 September 2026. It records how
`docs/plans/writer-regions-density.md` was carried out, on the branch
`plan/writer-regions-density` in the worktree
`.claude/worktrees/plan-writer-regions-density`, which branches from
`spec/writer-regions-density` at 838459f and not from `main`, because
neither the specs nor the plan are on `main` yet.

**The plan is done but for one choice that is the owner's: the deflate
of the bgzipped VCF.** Every task is ticked and every deliverable was
checked. One target is missed: the bgzipped write of `big.vcf` on one
thread takes 10.65 to 10.96 s, where the plan allows 9.79 s, bcftools'
8.9 s and a tenth. The options and their times are under "Waiting for
the owner" in work package 3, and the first of the questions at the end.
Nothing was merged into `main` and nothing pushed; the last commit of the
work is 116671d.

**What exists now that did not.** Through the core, both binding crates,
the Python package and the TypeScript package:

- `write_vcf(variants, "kept.vcf.gz")` and `writeVcf(variants)` write the
  variants of a pass as a VCF, bgzipped or plain. Each line is as the
  source VCF had it, with AC and AN taken out when the pass has fewer
  individuals than its source; from a vars file, each line is what the
  file holds. On the cases of the spec the output is byte for byte that
  of bcftools 1.24, and `bgzip -t` and `tabix -p vcf` accept it.
- `variants.filter_by_regions("genes.bed")`, and `exclude=True` for the
  other side, keep the variants inside or outside the regions of a BED
  file, as bcftools `view -T` and plink2 `--extract bed0` keep them. When
  it is the first filter of variants, the VCF reader does not parse the
  lines outside the regions and the vars file reader does not read the
  batches outside them: 1000 of the 100000 variants of `big.vcf` come in
  0.039 s from the plain file on one thread, where reading all of it
  takes 0.54 s.
- The missing rate of each variant, a sixth statistic of
  `calc_per_var_distribs`, checked against plink2.
- `calc_var_density(variants, 100000)`, the number of variants in each
  window along each chromosome, checked against tabix.
- What they rest on: every reader gives the header of its source (its
  individuals, the length of each chromosome and, for a VCF, the lines of
  its header), and a vars file keeps the lengths, at version 1.1 of its
  format.

The suites went from 1022 cargo tests, 551 pytest and 443 of 444 node to
1185, 612 and 473 of 474; the one node test that fails is the same one,
in `test/gwas.test.ts`, which fails on `main` and is not of this plan.

**What the reviews found.** Each of the six work packages was reviewed in
four to seven categories. In five of the six a review found at least one
defect that gave a wrong result with no error, ten in all, each fixed
with a test that fails without the fix; the missing rate had none: `##contig` lengths read wrong or dropped, an empty INFO
column written, the regions of a chromosome named `tracks1` dropped as a
header, a misspelt option of TypeScript that flipped the side of the
regions kept, a cut VCF read as whole when the cut fell outside the
regions. The owner should know of three changes that reach beyond this
plan: every options object of the TypeScript package now refuses a key it
does not know; `write_vars` in Python, like `write_vcf`, now takes away
its file after a panic; and a VCF with 200000 `##contig` lines opens in
0.030 s where it took 13.37 s.

**What is asked of the owner.** The choice of the deflate, five smaller
choices, and the merge. They are under "Questions for the owner" at the
end, each with the options and a recommendation.

This page was written while the work went, one work package at a time.
Its last section, "How the work went", is for whoever next revises a
skill or writes an implementation plan, and not for the owner; it says so
in its first line.

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

### Task 3.4, the review and the fixes

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

The measurement of task 3.5, commit f4e0635, is in "### The writer" of
"Speed" of `docs/specs/io_vcf.md`, taken on 27 September 2026 at load
averages of 3.5 to 4.0, three runs each, in seconds:

| write of big.vcf | bcftools, same session | target | from big.vcf | from big.vars |
|---|---|---|---|---|
| plain, 1 thread | 1.78, 1.69, 1.70 | 1.65 | 1.001, 0.954, 0.954 | 1.499, 1.485, 1.496 |
| bgzipped, 1 thread | 8.84, 8.90, 9.22 | 8.9 | 10.649, 10.664, 10.961 | 10.859, 10.800, 10.842 |
| bgzipped, 18 threads | 1.59, 1.66, 1.62 | 1.63 | 0.897, 0.890, 0.894 | 0.859, 0.854, 0.863 |

Two of the three targets are met; the bgzipped write on one thread is
not, and waits for the owner's choice of the deflate. The bgzipped file
is 38437792 bytes, 2 in 100 more than bcftools' 37695742.

## 4. The filter by regions

Task 4.1, commits 556e47a and c19cb56. `Regions` reads a BED, plain or
gzipped, and `RegionsReader` keeps the variants inside its regions, or
outside them with `exclude`, looking each position up with a binary
search over the joined regions of its chromosome. It offers the regions
to its source, which refuses them until work package 5, and adds to its
counts what the source skipped when the source took them. The code is in
`crates/popnei/src/filters/regions.rs`, a module of its own beside the
5600 lines of `filters.rs`. `tests/reference/filters/make_reference.py`
runs bcftools 1.24 `view -T` and `-T ^` and plink2 v2.0.0-a.7.7
`--extract bed0` and `--exclude bed0` with the BEDs of the spec; the two
programs agree, and give the 45 and the 455 of the spec on `many.vcf`.
The spec gained how a BED line is read (a carriage return dropped, `+99`
refused, names compared byte for byte) and two error cases the spec had
not named, a second filter by regions and a block whose chromosome number
has no name, a defect.

Task 4.2, commit b5b001c. `variants.filter_by_regions(bed_path,
exclude=False)` in Python and `variants.filterByRegions(bed, {exclude})`
in TypeScript, with the BED as bytes; the steps are `"regions"` and
`"excluded_regions"`. Nineteen pytest tests and five TypeScript ones.

Deliverables, run at 497d712:

- 1: the reference script ran, its check of the 45 and the 455 passed,
  and it left the tree clean.
- 2: `cargo test -p popnei --lib by_regions` gave 21 passed, where there
  were none.
- 3: pytest 583 passed, 6 skipped; `npm test` 455 of 456, the old
  failure; `npm run test:browser` 8 passed. `cargo test --workspace`
  gave 1126 passed, 2 ignored, and 150.

The review ran in five reports: spec, tests, numbers, errors with api,
and binding with architecture. A reviewer checked the filter against
bcftools and plink2 on a BED with unsorted lines, regions that overlap,
touch and nest, a region at a chromosome's last position, and a
chromosome of the BED the VCF lacks and the reverse: the three kept the
same 14 variants and excluded to the same 31. What they found:

- Three ways to get the wrong variants with no error. In TypeScript an
  option with a misspelt name, `{excluded: true}`, was dropped without a
  word, so the filter kept the regions it was asked to take out; every
  options object of the package did so, before this plan too. A BED line
  whose chromosome name begins with `track` or `browser`, `tracks1`, was
  skipped as a header line, where bcftools and plink2 keep it. And a BED
  starting with a byte order mark gave its first chromosome a name no VCF
  has.
- The test that the counts are the same whether the source skips or not
  could not fail: its source skipped exactly what the filter would drop,
  so a filter that trusted the source and dropped nothing passed. Work
  package 5 rests on it.
- A BED of 89.8 MB grew the memory of wasm by 314 MB, several copies of
  it, which a tab never gets back.
- The docstrings said the counts are the same wherever the filter stands
  among the steps, which is not so after a threshold filter.

The fixes are the 16 commits from fc05524 to f879804, each finding that
changed a behaviour with a test that failed first. Each mutation a
reviewer had shown to pass now fails a test. What the owner should know:

- Every function and method of the TypeScript package that takes options,
  fifteen of them, now refuses a key it does not know and names it. User
  code that passed an extra key and relied on its being ignored gets an
  `Error` now.
- A BED line is a header only when it is `track` or `browser` followed by
  a blank or the end of the line, or starts with `#`; an empty chromosome
  name is refused; a leading byte order mark is skipped. The spec says so.
- The memory of wasm that a BED of 89.8 MB takes went from 315 MB to
  225 MB, of which about 90 MB is the copy wasm-bindgen makes of the
  bytes; a test holds it below three times the size of the BED. Reading
  the BED twice, once to count the regions of each chromosome and once to
  fill them, would bring it to about 160 MB, and was not done.
- The regions of one chromosome are now reached by the bytes of its name
  through one handle, which the VCF reader of work package 5 keeps from
  line to line.

Not taken: the error of a second filter by regions names neither BED
file. That is what the spec chose; naming the one already set would help
the user find it, and it is among the questions at the end.

After the fixes: `cargo test --workspace` 1165 passed, 2 ignored, and 150;
pytest 586 passed, 6 skipped; `npm test` 460 of 461, the old failure;
`npm run test:browser` 8 passed.

## 5. The skip of the sources

Task 5.1, commits e543c6b and 7b50d43. When the filter by regions is the
first filter of variants, the VCF reader reads CHROM and POS of each line
in its serial pass and gives no row to a line the regions keep none of,
so its individuals' columns are not parsed; the vars file reader does not
seek to nor decompress a batch whose footer says each of its chromosomes
lies outside what the regions keep. Both count what they pass over, and
the filter adds it to its counts. The vars file of `many.vcf` in batches
of 100 reads its batches 0, 1, 2 and 4 with the BED of the spec, the
fourth batch skipped, and all five with `exclude`. The tests of work
package 4 run with the skip and with the offer refused and give the same
positions, compared by chromosome name, and the same counts.

The review ran in four reports, spec, tests, numbers with errors, and
architecture. It found the skip less careful than the spec in two ways,
both of them a wrong result with no error, and both fixed in the ten
commits from 4d6925a to 8063d1d:

- A line outside the regions was passed over before any check of its
  shape, while `docs/specs/io_vcf.md` checks the nine first columns of
  every line and the skip gives up only the checks of the genotypes. A
  plain VCF cut inside such a line, whose cut last line is the only sign
  it was cut, read as whole. A skipped line now keeps those checks: the
  nine columns, UTF-8 in them, a `GT` in FORMAT and, when the genotypes
  or the text are asked for, one column per individual of the header; a
  line that fails gets a row and its error. Counting the columns costs
  time: the read of the plain `big.vcf` with a BED that keeps 1000 of
  its variants went from 0.036 to 0.037 s to 0.044 s on one thread, at
  load averages of 6 to 9, under the 0.15 s of "Speed" of
  `docs/specs/filters.md`.
- The vars reader trusted its footer. Editing the smallest position of
  one batch in the footer from 1000 to 2100 gave 17 variants in place of
  45. The footer is now checked against every batch that is read, a
  batch that is skipped has the count of its rows checked against the
  footer without being decompressed, and a range whose smallest position
  is above its largest is refused at open. What cannot be caught without
  reading the batch, a footer that lies about the positions of a batch
  the regions skip, is written into `docs/specs/io_vars.md` as the one
  place the footer is trusted, and is among the questions at the end.
- Five tests could not fail, each shown by a change of the code that
  passed all 1173 tests: a vars file without positions under the skip,
  which gave no variant and no error; a batch whose largest position
  touches a region at one position; a batch counted twice after the
  columns asked for change in the middle of a pass; an empty POS; and the
  skip with `only_passed`, the default in Python, which gives 42 and 433
  kept of 475, as `bcftools view -f PASS,. -T` does.
- In wasm the reader reads one line at a time, so the skip allocated and
  hashed the chromosome name for every line; the reader now keeps both
  from line to line.

The read of the whole plain `big.vcf` with no regions is unchanged by the
skip: 0.572 to 0.600 s on one thread against the target of 0.594 s, the
same as before the change under the same load of 3 to 6, and 0.088 to
0.089 s on 18 threads against 0.108 s.

The measurement of task 5.2, commit 3f9e185, is in "### The filter by
regions" of "Speed" of `docs/specs/filters.md`, taken on 27 September 2026
at load averages of 3.2 to 3.5, with a BED of `chr1 0 1000000` that keeps
1000 of the 100000 variants, five runs each; every target is met:

| case | median of five runs | target |
|---|---|---|
| plain `big.vcf` read with the regions, 1 thread | 0.039 s | 0.15 s |
| bgzipped `big.vcf` read with the regions, 1 thread | 0.312 s | 0.45 s |
| `big.vars` through the filter, 1 thread | 0.005 s | 0.02 s |
| all of plain `big.vcf`, no regions, 1 thread | 0.542 s | 0.594 s |
| the same, 18 threads | 0.082 s | 0.108 s |

## 6. The density of the variants

Task 6.1, commits 317ec23 and f7f7a44. `calc_var_density` counts the
variants in windows of a size along each chromosome, up to its length
when the source gives one and to the last variant otherwise, in
`crates/popnei/src/stats/density.rs`. Twenty-six tests have `var_density`
in their names; with `pos` for `pos - 1` twelve of them fail. The counts
stay exact for lengths up to 18446744073709551615. The spec said the
error for a variant past a length of 10000 comes at chr1 10213; it comes
at 10028, the first of the six past it, and the spec says so now. Four
cases the spec left open were settled with the code and written into it:
a variant at position 0 is a `ValueError`, where tabix 1.24 counts it in
the first window with a warning; a window of more than 4294967295
variants is refused; a chromosome named twice in the lengths is refused;
and the last window of a chromosome with no length ends at
18446744073709551615 when its full width would pass it.

Task 6.2, commits 4713c52, eedb0f3 and 3316fc0. `calc_var_density(variants,
window_size, chrom_lengths=None)` in Python gives a frame of the
chromosome, the start and the end of each window and its count, and
`calcVarDensity(variants, windowSize, {chromLengths})` in TypeScript the
same; `calcVarDensity` is the fourteenth consumer of the TypeScript
package. `tests/reference/stats/make_reference.py` keeps the counts of
tabix 1.24 for each window of 1000 of `many.vcf.gz`. Twenty pytest tests
and five TypeScript ones. In TypeScript a window that ends past 2^53 is an
error, since the ends are numbers of JavaScript, which are exact only up
to 2^53; the core is exact to 18446744073709551615.

The review ran in four reports, spec, tests, numbers, and errors with api
and binding. No reviewer found a wrong count: 20000 random cases agreed
with a count made window by window, and the 31 windows of `many.vcf`
agree with tabix. What they found:

- In TypeScript a `Map` given as `chromLengths` was read as no lengths,
  with no error; and names that are whole numbers, "1" to "22", come out
  first and in ascending order, where Python keeps the order they were
  given, because that is the order in which JavaScript lists the keys of
  an object.
- No test at any layer held the order of `chrom_lengths`: a sort put into
  the core, into Python or into TypeScript passed every test.
- Two errors that only a reader with a defect can give reached Python as
  a `ValueError`, where the owner's convention makes a defect of popnei a
  `RuntimeError`.
- The spec did not say that a variant is counted in the window of its
  POS whatever the length of its REF, where tabix counts a deletion in
  every window it overlaps; the two agree on `many.vcf` only because all
  its variants are one base long.

The fixes are the commits from 66f2637 to 22803af. A `Map`, or any object
that is not a plain one, given as `chromLengths` is refused, naming it;
the TSDoc and the spec say that JavaScript lists the names that are whole
numbers first, and a test holds it. A bad length in Python names its
chromosome. The two errors only a reader with a defect can give are a
`RuntimeError`. The spec says a variant is counted in the window of its
POS, and that tabix counts a record in every window its REF overlaps.
Each layer has a test with the lengths out of the order of their names.

The first measurement of the density, commit d9a058c, is in "Speed" of
`docs/specs/stats.md`, taken on 27 September 2026 at a load average of
3.5: windows of 100000 over `big.vars` take 0.005 s on one thread and on
18, and over `big.vcf` 0.040 to 0.042 s on one thread and 0.039 s on 18,
since the density asks the reader for the positions alone. Both give 1000
windows that add up to 100000 variants. The spec sets no target yet.
The measurer saw the read of `big.vars` on 18 threads take 0.026 s, where
`docs/specs/stats.md` says 0.102 s on both thread counts; nobody has
looked into why.

## Questions for the owner

1. **The deflate of the bgzipped VCF.** With flate2 over miniz_oxide at
   level 6, the deflate the core has, the bgzipped write of `big.vcf` on
   one thread takes 10.65 to 10.96 s, above the 9.79 s the plan allows.
   The options, measured in one session on 26 September 2026 (the table
   under "Waiting for the owner" in work package 3):
   - flate2 over zlib-rs at level 6: 5.55 s and 37.4 MB, smaller than
     bcftools' 37.7 MB. It builds for both wasm targets, is pure Rust,
     and adds `zlib-rs` 0.6.8 behind a feature of flate2; flate2 then
     inflates with it too, so the readers of bgzipped files change their
     inflate, which was not timed.
   - miniz_oxide at level 5: 5.0 s and 41.6 MB, 10 in 100 larger than
     today's file. No dependency changes.
   - Keep level 6 of miniz_oxide and change the bound.

   Recommended: zlib-rs at level 6. It is the only option that is both
   faster than bcftools and as small as bgzip. If you agree, a task
   switches the feature, times the readers of bgzipped files before and
   after, and measures the writer again.
2. **The VCF writer parses every line.** It asks the reader for every
   column of a VCF, so a line with a genotype or a POS that does not parse
   is refused, as your rule that a corrupt input is reported asks. Reading
   only the text of the lines would take the plain write from 0.98 s to
   0.52 s on one thread, and would copy such a line into the output
   without an error. Recommended: keep the full parse.
3. **Position 0 in the density.** `calc_var_density` refuses a variant at
   position 0; tabix 1.24 counts it in the first window, with a warning.
   The VCF format allows position 0 for a telomere. Recommended: count it
   in the first window, as tabix does, and say so in the spec.
4. **The order of `chromLengths` in TypeScript.** JavaScript lists the
   keys of an object that are whole numbers first, in ascending order, so
   chromosomes named "1" to "22" come out in that order whatever order the
   user wrote, where Python keeps it. A `Map` or a list of pairs would keep
   it; today a `Map` is refused. Recommended: take a `Map` as well.
5. **The footer of a vars file is trusted for the batches the regions
   skip.** Every batch that is read is checked against its footer, and a
   skipped batch has its count of rows checked. A footer that lies about
   the positions of a batch the regions skip cannot be caught without
   reading that batch, and would change the variants kept. Recommended:
   accept it, since popnei writes the footer from the rows it writes; the
   spec says so now.
6. **Smaller.** The error of a second filter by regions names neither
   BED file, as the spec chose; naming the one already set would help the
   user find it. Recommended: name it.

Three things were seen outside this plan and not fixed; each can become a
GitHub issue if you want one. The Python binding turns any error of the
core it does not name into a `ValueError`, through a last arm that
catches the rest, so a new defect of popnei would read as a wrong input.
In TypeScript, `stats: [5]` is refused with a message that says an
`Array` was given. And the plain write of a VCF of 5000 individuals with
GT, AD, DP, GQ and PL peaks at 314 MB, which nobody has looked into.

And the merge of `plan/writer-regions-density` into `main`, when you
decide it. The branch also carries the specs of
`spec/writer-regions-density`, which are not on `main` either. The
branch `exp/writer-deflate`, in `.claude/worktrees/exp-writer-deflate`,
holds the timing of the three deflates; it is for the first question and
is not to be merged.

## How the work went

This section is for whoever next revises a skill or writes an
implementation plan, not for the owner, who can stop here.

- **The reviews found what the tasks missed, in every work package.**
  Five of the six had at least one finding that gave a wrong result with
  no error, ten in all, and in every one a reviewer's mutation showed a
  test the writer had meant to guard a case that could not fail. The review cost as much as the
  writing: the subagents that wrote a work package ended with 350000 to
  630000 tokens of context each, across their tasks and their fixes, and
  its reviewers used 60000 to 130000 tokens each, four to seven of them,
  that is 330000 to 750000 a work package.
- **Where the spec is silent, a subagent decides, and twice it decided
  against the owner's rule on errors.** One made a `##contig` line
  without an ID give no length and no error; another made the writer read
  only the text of a VCF for speed, which let a corrupt line through. Both
  were caught by the orchestrator and not by the subagent. The prompt of a
  task should carry the rule in a sentence: a case the spec leaves open is
  decided toward an error, and a check of the input is not traded for
  speed without the owner.
- **The plan's counts were wrong twice in the small.** It counted 21
  implementations of the trait where there are 24, and it said no test of
  before would change but for the version of the vars file, where a
  longer key changes every test that holds the key or a file's size.
  Counting with `grep` when the plan is written would have caught both.
- **One tree, one writer at a time, held.** Tasks that the plan let run
  side by side were run one after another in the tree, or beside reviewers
  that only read; a writer never compiled against another's half-made
  edit. The reviewers that read the shared tree saw the uncommitted files
  of the task running beside them, and said so, which cost a line in each
  prompt and nothing else.
- **Two hand-backs of reviewers had to be asked for again**: one handed
  back "No findings yet." before it had checked anything, and the report
  of the same reviewer did not arrive after its second run. A usage limit
  stopped the measuring subagent for about four hours; it went on from
  its context when resumed.
