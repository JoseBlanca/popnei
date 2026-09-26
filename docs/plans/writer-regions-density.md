# Plan: the VCF writer, the filter by regions, the missing rate and the density of the variants

26 September 2026. State: **under way**, approved by the owner on 26 September 2026. It builds
four items that the specs gained on 26 September 2026, on the branch
`spec/writer-regions-density`, which is not on `main`:

- the VCF writer and what the reader keeps for it, "The VCF writer" of
  `docs/specs/io_vcf.md`;
- the filter by regions and the skip of its sources, "The filter by
  regions" of `docs/specs/filters.md`;
- the missing rate, a sixth statistic of the per variant distributions,
  and the density of the variants, "The missing rate" and "The density of
  the variants along the chromosomes" of `docs/specs/stats.md`;
- what those four add to `docs/specs/block.md` (the three new methods of
  `BlockReader`, `SourceHeader`, `VcfText`), to `docs/specs/variant.md`
  (`Needs::VCF_TEXT`) and to `docs/specs/io_vars.md` (`chrom_lengths` and
  the version 1.1, the batches the reader skips).

Each runs through the core crate, both binding crates, the Python package
and the TypeScript package. The specs have no open point: the owner
decided the two they had on 26 September 2026, and the breakdown below on
the same day.

## In and out

Built: the six work packages below, the reference data they add under
`tests/reference/` with the commands that make it, and the three
measurements the specs ask for, the writer against bcftools, the filter
by regions against its targets, and the density's first number.

Not built, with where each goes: INFO, FILTER or the phase in the vars
file, which the owner will weigh later against the size of the file
("Not in this spec" of `docs/specs/io_vcf.md`); windows that slide by a
step, and matching `chr1` with `1` ("Not in this spec" of
`docs/specs/stats.md` and of `docs/specs/filters.md`); the regions of a
tabix index.

## What has to be in place

- The branch `plan/writer-regions-density`, from `spec/writer-regions-density`
  at the commit that holds this plan, in a worktree of its own under
  `.claude/worktrees/`, with a message on the board, `.claude/board/`,
  that names the shared files of the next paragraph.
- The files other branches may also change: `crates/popnei/src/block.rs`
  (the trait, which every reader of the workspace implements, benches
  included), `crates/popnei/src/variant.rs` (`Needs`), the `popnei` key
  of `crates/popnei/src/io/vars.rs`, `crates/popnei/src/stats.rs`
  (`PerVarStat`), and the entry points of both binding crates and both
  packages.
- The checks of the `coding` skill, run on `spec/writer-regions-density`
  at 7331638 on 26 September 2026, which has the code of `main` at
  cc11fc1: `cargo test --workspace` gave `1022 passed`, 2 ignored, and
  `150 passed`; `uv run maturin develop --release && uv run pytest` gave
  `551 passed, 5 skipped`; `npm test` in `js/popnei`, after `npm run
  build`, gave `tests 444`, `pass 443`, `fail 1`. The one that fails is
  `a kinship that does not tell the two variances apart gives none of
  them` of `test/gwas.test.ts`, an effect of 0.31250000000000355 against
  the literal 0.3125, 1.1e-14 of it against the 7e-15 its tolerance
  allows. It fails on `main`, is not of this plan, and is reported to the
  owner and not fixed here; "no fewer" below counts from 443 passed. `npm run build` in
  `js/popnei` comes before any `npm test` there, since a fresh worktree has
  no `js/popnei/wasm/`, as `docs/plans/diversity.md` found.
- bcftools 1.24, tabix 1.24, bgzip 1.24 and plink2 v2.0.0-a.7.7 in the
  PATH, checked with `which` and `--version` on 26 September 2026; each
  command of the specs' "How it is verified" was run with them that day.
  bedtools is not there and no check needs it.
- For the measurements, `/Users/jose/devel/popnei-bench/big.vcf`, 403 MB,
  and `big.vars`, which the specs' "Speed" sections measure on.

## 1. The reader trait and the header of a source

What it gives: every reader of popnei says what its source said of itself,
its individuals, the length of each chromosome, and for a VCF the lines of
its header, and can be offered regions to skip, which none takes yet. A
vars file keeps the lengths. No user sees it yet; work packages 3, 5 and 6
do.

Deliverables:

1. Every existing test passes untouched: the counts above, no fewer, and no test file of before changed but for the
   `format_version` `"1.1"` that `docs/specs/io_vars.md` now asserts, and
   for what the longer `popnei` key forces: `"chrom_lengths": []` in the
   literals that hold the key whole, and the byte sizes of vars files in
   `js/popnei/test/progress.test.ts`, 16 and 24 bytes more.
2. Cargo tests whose names contain `source_header`, in `io::vcf` and `io::vars`:
   `write.vcf` of "How it is verified" of the writer gives the lengths
   chr1 2000 and chr2 1500 and its nine meta lines; a `##contig` length of
   `0`, of `abc` and two lengths of one ID are the wrong header; the vars
   file written from it has `chrom_lengths` `[["chr1", 2000], ["chr2",
   1500]]` and `format_version` `"1.1"`; a file of 1.0 reads with no
   lengths. `cargo test -p popnei --lib source_header -- --list` counts
   them, and counted 0 on 7331638.
3. A cargo test whose name contains `source_header`, counted with those of
   deliverable 2, that each reader over a reader, the three threshold filters,
   the filter by linkage disequilibrium, the filter of individuals,
   `reblock` and the reader one block ahead, gives its source's header,
   and that the filters answer `skip_outside` with false and `num_skipped`
   with 0.

Stands on: nothing of this plan.

Tasks:

- [x] 1.1 `SourceHeader` and the three methods of `BlockReader` in
  `crates/popnei/src/block.rs`, from "The Rust interface" of
  `docs/specs/block.md`, and their implementation in every reader of the
  crate and of `crates/popnei/benches/`: the readers over a reader pass the
  header on, and the filters of variants refuse the offer. The sources
  answer false for now. Deliverables 1 and 3.
- [x] 1.2 The header of the VCF reader, from "What the reader keeps for
  the writer" of `docs/specs/io_vcf.md`: the meta lines, the lengths of
  the `##contig` lines and their refusals. Deliverables 1 and 2.
- [x] 1.3 `chrom_lengths` in the `popnei` key and the version 1.1, from
  the table of that key in `docs/specs/io_vars.md`: the writer takes the
  lengths from the header of its source, the reader gives them in
  `header()` and in `VarsMetadata`. Deliverables 1 and 2. Needs 1.2 for the
  test that writes `write.vcf` to a vars file.

What could go wrong: the trait has 24 implementations, 5 of them in the
benches (21 were counted when the plan was written), and a reader over a reader that answers with its own empty header
instead of its source's compiles; deliverable 3 is what finds it.

## 2. The missing rate

What it gives: `calc_per_var_distribs(variants, stats=[PerVarStat.MISSING_RATE])`
and `calcPerVarDistribs(variants, {stats: ["missing_rate"]})` give the mean
and the histogram of the missing rate of each population.

Deliverables:

1. Cargo tests whose names contain `per_var_missing_rate`, at
   `calc_per_var_distribs`, with the three means and the three histograms of
   the worked example of "How it is verified" of "The missing rate" as
   literals. 0 on 7331638.
2. The script `tests/reference/stats/make_reference.py` keeps `many.vmiss`,
   `many.popA.vmiss` and `many.popB.vmiss` of plink2, whose commands are in
   that part of the spec, and checks the table of that part against them. A
   pytest test at `calc_per_var_distribs` on `many.vcf` with and without the
   two populations compares with the three files: histograms exactly, means
   within 1e-12 relative, the 51 variants of popA in bin 5.
3. The pytest test of the pass against pyNei still leaves the missing rate
   out and passes; the test that all the statistics of one pass are those of
   one pass each covers six.
4. The TypeScript test of `calcPerVarDistribs` asserts the means of the
   table over all and in popA.

Stands on: nothing of this plan. It can run side by side with work
package 1, since neither changes the files of the other but for
`crates/popnei/src/stats.rs`, whose `BlockReader` of the tests 1.1 touches
and whose statistics 2.1 does; the orchestrator runs 2.1 after 1.1 is
committed.

Tasks:

- [x] 2.1 The statistic in the core, `crates/popnei/src/stats.rs`, from
  "The missing rate" and the six names of "The Rust interface" of
  `docs/specs/stats.md`. Deliverable 1. A variant put in the wrong bin
  gives a wrong count and no error, so this task is a commit of its own,
  and deliverable 1 is what guards it.
- [x] 2.2 Both binding crates, the Python and the TypeScript result, the
  reference files and the tests. Deliverables 2, 3 and 4.

## 3. The VCF writer

What it gives: `write_vcf(variants, "kept.vcf.gz")` and
`writeVcf(variants)` write the variants of a `Variants` after its steps as
a VCF, each line as the source VCF had it, or from a vars file what the
file holds, bgzipped or plain.

Deliverables:

1. `tests/reference/vcf/write.vcf`, the six lines of "How it is verified"
   of the writer, and the commands of that part, run by
   `tests/reference/vcf/make_reference.py`, whose outputs are stored beside
   it: the file of the filter of individuals from bcftools, and the rows of
   `bcftools query` of the file written from the vars file.
2. Cargo tests whose names contain `vcf_text`: the text of the lines that
   the VCF reader gives with `VCF_TEXT`, and that `retain_vars`,
   `retain_individuals` and `reblock` compact, cut and join with the rest
   of the block, on `write.vcf` and on `many.vcf` in blocks of 7 variants.
   0 on 7331638.
3. Cargo tests whose names contain `write_vcf`, at `write_vcf`, for each
   case of "How it is verified" of the writer: the bytes of `write.vcf`
   and of `many.vcf` written back, plain and bgzipped; the default
   `only_passed`; the filter of individuals against the file of bcftools
   of deliverable 1, and the order `c`, `b`, `a` that keeps AC and AN; the
   five lines from the vars file; the 215 lines of the missing data filter
   at 0.04. `bgzip -t`, `tabix -p vcf` and the query at `chr1:900-1100` on
   each bgzipped file, run by the tests when the programs are there and
   reported as skipped when not. 0 on 7331638.
4. The pytest tests of that part at `write_vcf`, and the TypeScript test of
   `writeVcf`.
5. The measurement of "### The writer" of "Speed" of `docs/specs/io_vcf.md`
   on `big.vcf` and `big.vars`, written into that section with the machine,
   the load average and every run, against the three numbers of bcftools there.

Stands on: work package 1, for the header.

Tasks:

- [x] 3.1 `Needs::VCF_TEXT` and `VcfText`, the text the VCF reader keeps
  and the compaction of it by `retain_vars`, `retain_individuals` and
  `reblock`, from "What the reader keeps for the writer" of
  `docs/specs/io_vcf.md` and the paragraph of `docs/specs/block.md` that
  begins "A block can hold one more column". Deliverable 2.
- [ ] 3.2 `write_vcf` in `crates/popnei/src/io/vcf.rs`, the lines from the
  text and from the columns and the header, with AC and AN taken out when
  the pass has fewer individuals than its source, from "What it gives" and
  "How it runs" of the writer. `write.vcf`, the script and its outputs.
  Deliverables 1 and 3. A line that differs from bcftools by one value is
  silent: guarded by deliverable 3.
- [ ] 3.3 The bgzip of the writer, members of 65280 bytes compressed on the
  threads of rayon and the empty member at the end, and
  `vcf_text_num_vars_per_block`. Deliverable 3, bgzipped.
- [ ] 3.4 `write_vcf` and `writeVcf` in both binding crates and both
  packages, with the handling of the path of `write_vars`. Deliverable 4.
- [ ] 3.5 The measurement. Deliverable 5.

What could go wrong: the compression. bcftools bgzips `big.vcf` in 8.9 s
on one thread, and flate2 over miniz_oxide, the deflate the core has, has
not been timed against it. When 3.5 misses that number by more than a
tenth, the orchestrator stops and brings the owner the times of the
compression levels of miniz_oxide and of another deflate crate that builds
for both wasm targets, since the choice changes a dependency.

## 4. The filter by regions

What it gives: `variants.filter_by_regions("genes.bed")`, with
`exclude=True` for the other side, and `filterByRegions(bed, {exclude})`,
which keep the variants inside or outside the regions of a BED file.

Deliverables:

1. `tests/reference/filters/make_reference.py` runs `bcftools view -T` and
   `-T ^` and `plink2 --extract bed0` and `--exclude bed0` with the BED of
   "How it is verified" of the filter on `many.vcf`, and on `write.vcf`
   with the BED of the worked example, and stores the positions each keeps
   beside it, checking the 45 and the 455 of the spec. The BED files are
   stored there too.
2. Cargo tests whose names contain `by_regions`: the worked example at
   `RegionFilter`; the 45 and the 455 at `next_block` of a `RegionsReader`
   over `many.vcf`, in blocks of 7 and of the default size, with the counts
   500 and 45; the four refusals of a BED, each with its line; the lines
   the BED reader skips; a gzipped BED; the error of a source with no
   positions. 0 on 7331638.
3. The pytest tests at `filter_by_regions` of that part, with the `steps`
   of both kinds, and the TypeScript test of the 45 and the 455.

Stands on: work package 1, for `skip_outside` and `num_skipped`, which the
filter calls and whose answer is false from every source until work
package 5.

Tasks:

- [ ] 4.1 `Regions`, `RegionSelection`, `RegionFilter` and `RegionsReader`
  in `crates/popnei/src/filters.rs`, and `PassStep::Regions` in
  `chain_of` and `refuse_a_second_filter_of_a_kind`, from "The filter by
  regions" and its part of "The Rust interface" of `docs/specs/filters.md`.
  The BEDs and the script. Deliverables 1 and 2.
- [ ] 4.2 The step in both binding crates and both packages. Deliverable 3.

## 5. The skip of the sources

What it gives: the same variants and the same counts as work package 4,
faster, when the filter by regions is the first filter of variants: the
VCF reader does not parse the lines outside the regions and the vars file
reader does not read the batches outside them.

Deliverables:

1. The tests of deliverable 2 of work package 4 run with the skip and
   without it and give the same positions and the same counts, compared by
   the names of the chromosomes. A cargo test whose name contains
   `skip_outside` asserts that the vars file of `many.vcf` in batches of 100 has its
   fourth batch skipped with the BED of the spec and none with `exclude`,
   by the count of batches the reader decompressed.
2. Cargo tests whose names contain `skip_outside` for the counts of "How it runs"
   of the filter: a threshold filter before the filter by regions, which
   leaves the source without the regions; a `"regions"` and an
   `"excluded_regions"` step together; the filter of individuals between
   the filter and the source, which hands the regions on; and a POS that
   does not parse, which is the error of that column with the skip.
3. The measurement of "### The filter by regions" of "Speed" of
   `docs/specs/filters.md`, written into that section, against its three
   numbers.

Stands on: work package 4.

Tasks:

- [ ] 5.1 The skip of the VCF reader in its serial pass, and of the vars
  file reader by the footer, with `num_skipped`, from "How it runs" of the
  filter by regions and the paragraph of the reader of
  `docs/specs/io_vars.md` that begins "The reader answers". Deliverables 1
  and 2. A count off by the skipped variants is silent: guarded by
  deliverables 1 and 2.
- [ ] 5.2 The measurement. Deliverable 3.

What could go wrong: the serial pass of the VCF reader is what bounds it
on 18 threads, by "Speed" of `docs/specs/io_vcf.md`, and reading POS there
adds to it. The read of the whole plain file without regions has to stay
within the targets of that section, 0.594 s on one thread and 0.108 s on
18, which 5.2 measures too.

## 6. The density of the variants

What it gives: `calc_var_density(variants, 100000)` and
`calcVarDensity(variants, 100000)` give the number of variants in each
window of 100000 base pairs along each chromosome.

Deliverables:

1. Cargo tests whose names contain `var_density`, at `calc_var_density`: the
   worked example of "How it is verified" of the density, in windows of 500
   and of 600, with the lengths of the header and without them; the
   `ValueError` of a variant past a length; a `window_size` of 0; more
   windows than `MAX_NUM_WINDOWS`. 0 on 7331638.
2. The script `tests/reference/stats/make_reference.py` keeps the counts of
   `tabix` for each window of 1000 of `many.vcf.gz`, with the command of
   that part, and the pytest tests at `calc_var_density` assert the two
   tables of that part, the lengths of `chrom_lengths` and its error.
3. The TypeScript test asserts the first table.
4. A first measurement of the density on `big.vcf` and `big.vars`,
   written into "Speed" of `docs/specs/stats.md`, which gives no number to
   reach until then.

Stands on: work package 1, for the lengths. Task 6.1 edits
`crates/popnei/src/stats.rs`, as 2.1 does, and runs after 2.1 is
committed; otherwise work package 6 can run side by side with 3, 4 and 5.

Tasks:

- [ ] 6.1 `calc_var_density` in `crates/popnei/src/stats.rs`, from its
  item and its part of "The Rust interface". Deliverable 1.
- [ ] 6.2 Both binding crates, both packages, the reference counts, the
  tests and the measurement. Deliverables 2, 3 and 4.

## The whole plan

The sum of its work packages, and the checks of the `coding` skill on the
last commit, the wasm build of `crates/popnei-js` among them.
