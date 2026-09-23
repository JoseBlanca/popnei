# The work on the kinship

23 September 2026. What was done while `docs/plans/kinship.md` was carried
out, written as the work went. The plan builds the genomic relationship
matrix of a set of variants, which says for every pair of individuals how
much of their genome they share beyond what two individuals drawn at random
from the panel share, from `docs/specs/kinship.md`. The work is on the
branch `plan/kinship` in the worktree `.claude/worktrees/kinship`, and
nothing of it is on `main`.

State: under way, work package 1 of 3.

## The starting commit

`85d27a2`, the local `main` at `75c5fd3` with `spec/gwas` merged into it.
The merge had no conflict. What the plan's "What has to be in place" asks
for was there: `wasm-bindgen` at `/Users/jose/.cargo/bin/wasm-bindgen`; the
reference data under `tests/reference/kinship/`, both panels and both
plink2 matrices with their `.id` files; `tests/reference/dists/panel.vcf.gz`
and `tests/reference/gwas/phenotypes.csv`; and the `wasm-check` alias in
`.cargo/config.toml`.

The eight checks the plan lists were run on that commit and each gave the
number the plan gives, so the work started from the state the plan was
written against.

| command | result |
| --- | --- |
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo test --workspace` | 604 passed, 2 ignored in the core crate; 149 in the linear algebra crate |
| `cargo test -p popnei-linalg --no-default-features` | 136 passed |
| `cargo test -p popnei --lib kinship -- --list` | `0 tests, 0 benchmarks` |
| `cargo test -p popnei --lib pca -- --list` | `47 tests, 0 benchmarks` |
| `uv run maturin develop && uv run pytest` | 347 passed |
| `npm run build && npm test` in `js/popnei` | 242 pass, 0 fail |

## Work package 1: one row pass for the PCA and the kinship

Under way. One change to the plan was made before the task started.

### What was changed in the plan, and why

The check of deliverable 1 was `grep -c "fn the_standardized_row"
crates/popnei/src/*.rs` giving 1 in `variant.rs` and 0 in `pca.rs`, "where
today it is 0 and 1". Run on the starting commit it gives 5 in `pca.rs`,
because without the open parenthesis the pattern also counts
`the_standardized_rows`, which drives a whole block and has one definition
for the threads and one for WebAssembly, and
`the_standardized_rows_one_by_one`. It counts the wrapper the
`bench-internals` feature exposes as well. A check whose starting number is
wrong cannot say whether the move happened, so it was replaced by `grep -rn
"fn the_standardized_row(" crates/popnei/src/*.rs`, which names `pca.rs`
twice today, at the function and at that wrapper, and must name only
`variant.rs` when the task is done. The new check is the stronger of the
two: it distinguishes the row from the block driver, which the old one
could not.

### What was built

Task 1.1, in commit `23e82e9`: six files, 768 lines added and 472 removed.
The pass that turns one variant into its standardized dosages,
`the_standardized_row` with its buffers `RowScratch` and the helpers it
calls, left `crates/popnei/src/pca.rs` for `crates/popnei/src/variant.rs`,
and now takes a `DosageOptions` holding `transform_to_biallelic` and a
`DosageScale`, which is either the standard deviation of the dosages, what
the principal components of the variants divide by, or
`sqrt(ploidy * p * (1 - p))`, what the kinship will divide by in work
package 2. The three functions that drive a whole block of variants stayed
in `pca.rs`, one for the threads, one for WebAssembly, which has none, and
one that reads the rows one after another.

`Error::PcaVariantWithMoreThanTwoAlleles` became
`Error::VariantWithMoreThanTwoAlleles`, since the kinship raises it too.
Its message is unchanged character for character.

### The deliverables

Each check was run by the orchestrator on `eae2e59`, the commit of the code
with the correction to the spec on top of it.

| deliverable | command | result |
| --- | --- | --- |
| 1, the pass lives in `variant.rs` | `grep -rn "fn the_standardized_row(" crates/popnei/src/*.rs` | two hits, both `variant.rs`; `pca.rs` named in neither |
| 2, the PCA computes what it computed | `cargo test -p popnei --lib pca -- --list` | `47 tests, 0 benchmarks` |
| 2 | `cargo test --workspace` | 607 passed, 0 failed, 2 ignored in the core crate; 149 in the linear algebra crate |
| 2 | `uv run pytest tests/test_pca.py` | 33 passed |
| 3, one error for both callers | `cargo test -p popnei --lib more_than_two_alleles -- --list` | `2 tests`, one named for the PCA and one for the row pass called with either divisor |

The other checks of the `coding` skill on the same commit: `cargo fmt --all
--check` and `cargo clippy --workspace --all-targets -- -D warnings` clean,
`cargo test -p popnei-linalg --no-default-features` 136 passed, `cargo
wasm-check` clean, `uv run pytest` 347 passed, and `npm run build && npm
test` in `js/popnei` 242 pass with 0 fail.

Deliverable 2 asks for 604 tests in the core crate and there are 607. The
three are new in `variant.rs`, which deliverable 3 asks for; none was
removed. Every line of `pca.rs` holding an assertion is the same before and
after the move, checked by sorting them and comparing the two lists, except
one `const _: () = assert!(GENOTYPES_PER_RUN <= 255)`, which is read when
the code is compiled and not when it is run, and which moved to `variant.rs`
with the constant it is about. One test was renamed, from
`a_variant_of_three_alleles_is_the_error_that_names_its_position` to
`a_variant_of_more_than_two_alleles_is_the_error_that_names_its_position`,
so that the check of deliverable 3 finds a test for each of the two callers.
The 45 test names the two versions of the file hold differ in that one name
and nothing else.

### What the review found

Seven reviewers read the commit, one per category, each with a fresh
context. Sixteen findings held and were fixed in seven commits, `d6b00af`
to `e445cc5`. The checks were run again after them: the core crate passes
610 tests with 2 ignored, the linear algebra crate 149 and 136 without its
default features, Python 347, JavaScript 242, and `cargo fmt`, `cargo
clippy --workspace --all-targets`, `cargo clippy -p popnei --all-targets
--features bench-internals` and `cargo wasm-check` are clean. The three
deliverables were checked again and still hold, with the 47 PCA tests
carrying the same names they carried before task 1.1.

The one that mattered most is the divisor the kinship will use. It was
tested only on diploid variants whose mean dosage is 1, where `p` is 0.5
and `p * (1 - p)` and `p * p` are the same number, so two wrong formulas
passed all 607 tests: `p * p` in place of `p * (1 - p)`, and the ploidy
read as 2 instead of taken from the variant. Both were run and both left
`607 passed; 0 failed`. Work package 2 would have built the kinship on a
divisor that nothing checked. The tests added use a diploid variant of
`p = 0.125`, a tetraploid and a haploid one, with the values computed by
hand and checked against pyNei's `_calc_dosages` with
`sqrt(ploidy * freqs * (1 - freqs))`. Both mutations were applied again
after the fix and each now fails
`variant::tests::the_divisor_under_hardy_weinberg_reads_the_ploidy_and_the_frequency`
and nothing else.

The rest, in the order of what they would have cost:

- The pass refused a ploidy above 254 with an error whose message opens
  "the principal components of the variants cannot be taken on this
  dataset", from the pass the kinship is to share. No user could reach it,
  because `pca_of_variants` refuses the ploidy at its own entry before a
  block is read, but the kinship will. It now has a case of its own,
  `Error::VariantPloidyTooLarge`, in the group of the errors of a variant,
  a `ValueError` in Python like its neighbour. `MAX_PLOIDY_OF_THE_VARIANTS`
  moved to `variant.rs` with it, so `variant.rs` no longer imports from
  `pca.rs`: it had been importing from the module above it, which would
  have pulled `pca` in behind every caller of the row.
- The two buffers of a row, 255 dosage counters and 256 values, are sized
  for a ploidy of 254 and nothing said so. Raising the limit to 255 gave a
  genotype with an allele missing a value of 16.03 where 0.0 is right, with
  no error and no panic. Three relations are now asserted when the crate
  compiles; raising the limit stops the build with "the dosages of the
  largest ploidy, which are the ploidy and one more, each need a counter".
- `codes.resize` could be deleted with all 607 tests passing, because every
  fixture built its buffers at the size of its row. A row longer than the
  buffers would have been left half written.
- The test of a variant with more than two alleles used a variant of four,
  so tightening the refusal from more than two to more than three left it
  green.
- `DosageOptions` and `DosageScale` were public, so a crate outside the
  workspace could build one and pass it nowhere; they are the crate's own
  now. Two helpers had been widened to the crate for no caller.
- The divisor read the ploidy out of a loop bound, so a caller that passed
  anything else got a wrong number and no error; it takes the ploidy.
- Nothing asserted the middle clause of the message deliverable 3 asks to
  keep unchanged, so it could have been rewritten with nothing going red.
- The tail of the row pass is copied into the benchmark and nothing
  compared the two; a test now asserts they give the same row bit for bit.
- Seven doc comments still said PCA or private, one in the Python binding
  claimed the row pass is walked by every calculation that turns a variant
  into dosages, which the r² and the filters disprove by having their own.

### What was not taken

- **The six tests of the moved code stay in `pca.rs`.** Two reviewers asked
  for them to move to `variant.rs`, where the code they test now lives.
  Moving them would drop `cargo test -p popnei --lib pca -- --list` below
  47 and so destroy the check deliverable 2 rests on, which is the only
  evidence that the move changed no number. It is worth doing once that
  check has served, which is when work package 2 is reviewed.
- **The JavaScript message names `transform_to_biallelic` where the
  TypeScript option is `transformToBiallelic`**, so a user greps for a name
  that is not in their code. It is older than this plan and it is what a
  user of the principal components sees today, so it is the owner's and it
  is asked of them below. The crate already rewrites `maxNumVars` and
  `maxAllowedMaf` this way, so there is a pattern to follow.
- **The Python binding maps 29 of the 99 cases of `Error` through a
  wildcard**, so a case added later becomes a `ValueError` without anyone
  choosing that. `Error` is `#[non_exhaustive]`, which makes the wildcard
  compulsory, and every case that falls through it is right today. It is
  older than this plan and nothing of this plan rests on it.

### What was changed in the plan

- Deliverable 2 asked for 604 tests in the core crate. That was the count
  before the work package, and deliverable 3 asks for tests that did not
  exist, so the two could not both be met. The check now reads the list of
  test names of `pca.rs`, which is what says that no test was dropped, and
  asks for no fewer than 604 in the crate.
- Task 2.0 is new. Since task 1.1 the three functions that drive a whole
  block of variants hold nothing of the principal components and take a
  `DosageOptions`, and the kinship needs the same drive over a block, so
  task 2.1 would have copied about 150 lines including the arm for the
  threads and the arm for WebAssembly, which is a second place where the
  two can fall out of step.
- The check of deliverable 1 was replaced before the work started, for the
  reason above.

### What the owner should know

- The spec said the kinship refuses two datasets "both for reasons the
  PCA's row pass already refuses them for". Only the ploidy is refused
  there; the limit of 46340 individuals is checked in `pca_of_variants`,
  which the kinship will never call. Whoever wrote work package 2 from that
  sentence would have got the ploidy for free and lost the other in
  silence. The spec now says which check is where.
- `DosageScale::OfHardyWeinberg` has no caller in the library until the
  kinship exists, so it carries an expectation that it is dead code. The
  commit that builds the kinship has to delete that expectation or the
  build fails, and the text of the expectation says so.
- The worked example of "How it is verified" in the spec, and both of its
  kept variants, have `p = 0.5`. The cargo test work package 2 builds from
  it inherits the blind spot this review found, so the eleven plink2
  literals of that section are what will close it there.

### How the work went

Task 1.1 took one subagent 197000 tokens over 112 tool calls and 15
minutes. The seven reviewers took 756000 tokens between them, from 79000
for the binding to 122000 for the numbers, and ran in parallel in about 13
minutes. The fixes took the subagent that wrote the code another 119000
tokens over 20 minutes; sending them back to it rather than to a fresh one
cost nothing in re-reading.

Sending all seven categories rather than choosing among them was worth it
here: the finding that mattered most came from `tests`, which was the
slowest of the seven and the only one that mutates the code, and three of
the others found the PCA-named ploidy error from three different sides.

## Work package 2: the matrix

Done. `calc_kinship` gives the genomic relationship matrix of a set of
variants in the core crate, in Python and in TypeScript, with its per pair
denominators, `num_vars`, the pass stats and `filter_individuals`.

### The deliverables

Each check was run by the orchestrator after the review's fixes.

| deliverable | command | result |
| --- | --- | --- |
| 1, both panels are plink2's | `uv run pytest tests/test_kinship.py` | 37 passed; every entry of both 200 x 200 matrices within 1e-12 relative of plink2, `num_vars` 1200 for both |
| 2, the worked example and the plink2 literals | `cargo test -p popnei --lib kinship -- --list` | `35 tests`, where the plan asks for at least 12 and the work began at 0 |
| 3, one test per case | `cargo test -p popnei --lib kinship::cases -- --list` | `6 tests`, one for each of the four cases the spec names and two the review added |
| 4, the individuals argument | `uv run pytest tests/test_kinship.py` | 1195 variants and a largest difference of 0.12917 from the same 40 rows and columns of the whole panel |
| 5, popnei and pyNei agree | `uv run pytest tests/test_kinship.py` | every entry of both panels within 1e-12 relative, `num_vars` equal |
| 6, `calcKinship` under node | `npm run build && npm test` in `js/popnei` | 260 pass, 0 fail |

The whole suite on the same commit: `cargo fmt --all --check`, `cargo
clippy --workspace --all-targets -- -D warnings`, `cargo clippy -p popnei
--all-targets --features bench-internals -- -D warnings` and `cargo
wasm-check` clean; `cargo test --workspace` 645 passed with 2 ignored in the
core crate and 149 in the linear algebra crate; `cargo test -p popnei-linalg
--no-default-features` 136; `uv run pytest` 384 passed; `npm run build &&
npm test` 260 pass.

Deliverable 4 settled a reading of the spec: the 0.129 it quotes is the
fully called panel, which gives 0.12917. The panel with genotypes missing
gives 0.1387.

### What the tolerance against plink2 really was

The plan and the spec asked for every entry to be within 1e-5 absolute of
the stored `tests/reference/kinship/*.plink2.rel.gz`, and the first run gave
a largest difference of 4.95e-06 on `panel_called` and 4.93e-06 on `panel`.
That looked like agreement with a factor of two to spare. It was not
agreement at all: those files hold six significant digits of text, and the
whole difference is plink2's rounding. Every one of the 40000 differences is
below the half-ulp of six digits.

Re-run against plink2's binary output, `--make-rel square bin`, popnei is
4.44e-16 from plink2 on `panel_called` and 5.55e-16 on `panel`, absolute.

So the check had no headroom: an arithmetic error of up to 5e-6 would have
passed it. And plink2's rounding is relative while the tolerance was
absolute, so the rule held only while the entries stayed below about 2. On a
panel of 60 individuals and 600 variants where most alleles are private,
which is ordinary in sequence data of a small panel, the diagonal reaches
19.66: there a correct popnei is 1.14e-13 from plink2's binary matrix and
1.83e-05 from its printed text, and would have failed the test.

`make_reference.py` now writes `<name>.plink2.rel.bin.gz` as well, and both
the cargo and the pytest tests compare against those within 1e-12 relative.
The eleven text literals stay for the spec's table, which is what a reader
checks by eye.

The first bound written was 1e-12 of each entry, and it was the wrong shape,
not merely the wrong number. An entry of the matrix is a sum of products
that cancel, so it can come out as near 0 as the data makes it while the
rounding of its sum stays where it was. Asking each entry to be within a
share of *itself* therefore asks the smallest entries for an accuracy no
arithmetic has. It broke on the backend nobody was testing: the faer run
failed at an entry of 1.29e-05 whose difference from plink2 was 1.9e-17,
smaller in absolute terms than the differences at entries a hundred times
larger, which passed. Each entry is now asked to be within 1e-13 of the
**largest** entry of the matrix. Measured over all 40000 entries of each
panel, the worst difference as a share of that largest entry, 1.23 on both
panels, is 3.6e-16 and 4.5e-16 with Accelerate and 3.3e-15 and 2.3e-15 with
faer, so the bound is thirty times the worst of the four and ten times
tighter at the largest entry than the rule it replaced.

faer is seven times further from plink2 than Accelerate on the same data.
A sum of 1200 products can carry 2.7e-13 of the scale, so both are well
inside what the order and the blocking of the sums allow, and the spec
records it as a range rather than a defect.

The agreement with pyNei is bit for bit, a largest difference of 0.0 on both
panels, but that is an artifact of both libraries calling the same `dsyrk`
on a single chunk of 1200 variants. At 30000 variants they differ by
1.6e-14, and on faer they would differ further. The bit equality is not a
property to lean on.

### What the review found

Seven reviewers, one per category. Twenty-two findings held and were fixed
in fourteen commits, `ba94a0b` to `7e58df6`.

**Four wrong numbers that every test passed.** The pass over the blocks
carries, into the denominator of every pair, how many variants the blocks
before it used. Passing how many they *gave* instead, one word apart, left
all 637 tests green; on a dataset of 10002 variants with a variant dropped
before the first missing genotype, every entry of the matrix moves. The same
blindness hid three more: the position handed to the row pass, which names
the variant in the error of more than two alleles; the two variant counts in
the error of a pair never called together, which could be swapped; and the
ploidy, which no kinship test used at anything but 2.

All four are invisible for one reason. `calc_kinship` calls
`Reblock::new(reader, None)`, and `Reblock` joins any dataset under 10000
variants into a single block, so no fixture in the repository has a second
block, and no caller can ask for a smaller one. The plan's own warning said
a failure would show as "a failure on `panel` with a pass on
`panel_called`". It cannot: both panels are one block. The tests are now on
hand-built multi-block readers, and each of the four mutations fails its own
test and no other.

**Two silent wrong results a user could reach.** A `Kinship` built by hand
in Python took a matrix that named one individual twice, and
`filter_individuals` then returned two rows for the one name asked, which
would travel on to `calc_gwas`. And `filter_individuals("a")`, one name
where a sequence is meant, returned a kinship of the individual `a` rather
than refusing: `filter_individuals("ab")` gave two individuals. Both are
refused now, and both are in the spec.

**The two packages refused different matrices, in both directions.** Python
took an individual named twice, which TypeScript refused; TypeScript took a
`NaN` on the diagonal, which its symmetry check never looked at, while
Python refused a `NaN` by calling the matrix asymmetric, which is the wrong
reason, since a matrix can hold a `NaN` and be symmetric and every
comparison with a `NaN` is false. Both now refuse a value that is not a
number first, naming its cell, and then test symmetry.

**Messages that said what was not true.** The error of a pair with no
variant called in both gave positions where the spec twice promises names,
and with `individuals` given those positions were not even the file's; three
reviewers found it. An individual whose sequencing failed, with no called
genotype anywhere, read "the individuals at the positions 0 and 0 ... leave
one of the two out", and is now named on its own. And "every variant has the
same genotype in every individual" is false: the rule is one dosage, not one
genotype, and both an all-missing dataset and a collapsed multiallelic one
contradict it.

**The count the core threw away.** `Kinship` carried only the variants that
were used, while `pass_stats` means the variants the pass gave, and the core
computed the second and dropped it. Every other calculation hands it out:
`Pca::num_cols`, `Stats.num_vars`, `KosmanSums::num_vars`. So each binding
wrapped the reader chain in its own `BlockReader` to count it again, about
60 lines each, written independently and without sight of each other; they
were the only two `impl BlockReader` outside the core. `Kinship` now carries
`num_vars_given` and both wrappers are gone.

The part of this worth keeping is not the duplication, which both subagents
reported themselves. It is that the core counts with `checked_add` and
raises, and both layers above independently chose `saturating_add` and
wrapped in silence. Two layers recomputing a number the layer below had
already computed, and both landing on a weaker overflow rule than it, is a
shape that will recur wherever a binding wants a number the core discards.
The evidence is `crates/popnei-python/src/kinship.rs` and
`crates/popnei-js/src/kinship.rs` as they stood at `8392778`.

**Memory.** `Kinship.__post_init__` checked symmetry with `values -
values.T` and `numpy.abs`, two more matrices of individuals by individuals.
Measured with `tracemalloc` on a frame of 3000 individuals, 72 MB: a peak of
144 MB, which at the 10000 individuals the docstring names gives back about
1.6 GB of the 800 MB that handing the array to numpy without copying it had
saved. It reads one row at a time now: 0.4 MB and 0.043 s on the same frame.

**Smaller.** `num_vars` was a `usize`, so the same count was 32 bits in a
browser and 64 natively; it is a `u64`. `kinship.rs` took the limit on the
individuals from `pca.rs`, the wrong way through the layers, and it sits in
`variant.rs` now with `pca` re-exporting it. Four spec statements had code
and no test. Two test comments explained guarantees the tests did not give,
and one of them, the 37-variant test, now asserts bit equality against the
default reading and so becomes the first real guard on re-blocking.

### What was not taken

- **The six tests of the moved row pass stay in `pca.rs`**, as in work
  package 1: moving them drops the PCA count below 47 and removes the
  evidence that no number of the principal components changed. Worth doing
  once the plan is merged.
- **The Python error mapping sends 29 of the 99 cases through a wildcard**,
  so a case added later becomes a `ValueError` with nobody choosing it.
  `Error` is `#[non_exhaustive]`, which makes the wildcard compulsory, and
  every case falling through it is right today.
- **The JavaScript message names `transform_to_biallelic` where the option
  is `transformToBiallelic`.** Older than this plan, and it is what a user
  of the principal components reads today, so it is the owner's and is asked
  of them. A later fix touches `crates/popnei-js/src/kinship.rs:144`, the
  PCA's own line, and a case in `errors.rs`.
- **`PcaNoVariantWithVariance` carries the same false wording** about every
  variant having the same genotype. Fixing it changes what a user of the
  principal components reads, so it is the owner's.

### What the owner should know

- The kinship walks each block's genotypes twice more, serially, to build
  the mask of called genotypes that the parallel row pass had already
  decided per genotype. No measurement was made; it is for the performance
  session.
- `cargo wasm-check` covers `popnei` and `popnei-linalg` only, not
  `popnei-js`, and `npm test` does not rebuild the wasm. The plan's final
  check says `npm run build && npm test`, which does, so the plan is right;
  anyone running `npm test` alone is testing the previous wasm.
- `cargo doc -p popnei` fails on a link in `ld.rs` that predates this plan,
  so rustdoc is not among the checks that run.
- pyNei's VCF reader raises `IndexError` above a ploidy of 2, so the
  tetraploid test takes its numbers from `Variants.from_gt_array` instead.

### How the work went

Building the work package took four subagents about 940000 tokens: 144000
for the move of the block pass, 427000 for the core and its two fix rounds,
315000 for Python, 251000 for TypeScript. The seven reviewers took 1035000
between them, from 127000 for the architecture to 166000 for the spec.

The review cost more than the building and found four wrong numbers that
every one of the 637 tests passed. The single most useful thing it did was
mutate the code: the `tests` reviewer, the slowest of the seven at thirteen
minutes, found three of the four, and the orchestrator confirmed each by
rerunning the mutation before and after the fix.

Two things about the shape of the work are worth carrying to the next plan.
Splitting the fixes into rounds by file ownership kept three subagents from
editing one file, and left the branch tip uncompilable between the rounds,
which cost another session its baseline; a round that another branch may
take should end green. And running the two binding tasks in parallel, with
neither able to see the other's files, is what produced two copies of the
same counting reader and two copies of the same overflow mistake.

## The five the owner decided on 23 September 2026

The review of work package 2 left five questions that were the owner's,
because each changes something a user reads or a document that governs how
popnei is written. All five were answered the same day and all five are
done.

**A message in TypeScript names its options in TypeScript.** The error of a
variant with more than two alleles told a JavaScript user to pass
`transform_to_biallelic`, while the option is `transformToBiallelic`, so
they would grep their code for a name that is not in it. The owner made it
a rule rather than a fix, so the whole surface was swept. Eight names of
the form `a_b_c` appear in the core's messages; seven have a TypeScript
option of the same meaning. Three needed rewriting, `transformToBiallelic`,
`numBins` and `numVarsPerBlock`; two, `polyThreshold` and `binType`, were
already rewritten before the error left `stats.rs`; and three cannot be
reached from JavaScript at all, `numPrinComps` because the second reader is
opened exactly when it is above 0, and the two of `max_num_vars` because
one is already rewritten and the other is refused by `arguments.ts` against
the same bound the core checks. The eighth, `popnei_batches`, is not an
option: it is the key popnei writes into the footer of a vars file, spelled
that way in the file itself, so both languages keep it.

The rewrite is one helper matching on the error, not a case for each. The
three cases that existed carry something the binding knows and the core
does not, the name of an argument; a spelling does not. `BlockTooLarge`
alone can be raised from seven calls of the wasm crate, and a `map_err` at
each is one that a later call forgets.

**An error said what was not true.** `PcaNoVariantWithVariance` read "every
variant has the same genotype in every individual". The rule is
`dosages_seen < 2`, one dosage and not one genotype, and two datasets reach
that message and contradict it: every genotype uncalled, and `0/1 0/2 0/1`
read with `transform_to_biallelic`. It now reads "no variant has more than
one dosage among its called genotypes, so none of them varies", which is
word for word what the kinship's own case already said. `docs/specs/pca.md`
quoted the old sentence and now states the rule, names pyNei's
`std(mat012) > 0` and says the old wording was a mistake of the message and
not a difference of behaviour. Nothing else in popnei had the same mistake.

**A kinship matrix is labelled with names.** Python took a frame labelled
with integers, which `pandas.DataFrame(matrix)` gives by default, while
TypeScript required strings; it was the last of the five places where the
two packages disagreed about what to accept. Python refuses it now, with a
message that names the mistake a user will actually make and the line to
write instead. The guards of a hand built `Kinship` are now sorted by
kind: a `TypeError` for a wrong type, a matrix that is no frame, a label
that is no name, a bare string where the names are meant, an `individuals`
that cannot be walked; a `ValueError` for a wrong value, a matrix that is
not square, two sides naming different individuals, an individual named
twice, a cell that is no number, a cell that is not finite, a matrix
outside the symmetry tolerance, `individuals` naming nobody, and a name
that is not in the matrix.

**The `coding` skill gained a paragraph.** A binding that works out a
number for itself is a sign the core threw one away. The evidence is in
this plan: both bindings wrapped the reader chain in a `BlockReader` of
their own to count the variants a pass gave, about 60 lines each, written
at the same time without sight of each other, because the core computed
that count, read it once and left it out of its result. The duplication was
the half that showed. The other half is that the core counted with
`checked_add` and raised, and both wrappers used `saturating_add` and
stopped counting in silence, so the three layers refused different datasets
and nobody had decided that.

**Ctrl-C is to interrupt a pass, and does not.** The owner decided it
should be fixed even if the core needs a way to send a message back. It is
not fixed here, and the scope is why. `crates/popnei-python/src/source.rs`
checks for a signal between blocks, because there Python drives the loop.
The five functions that run over a whole dataset, the kinship, the
principal components, the distances, the linkage disequilibrium and writing
a vars file, check only before and after, because the core owns the loop
and returns when the file is done. Writing a vars file looks like a
precedent and is not: it checks after the write and deletes the half
written file.

The design that works needs no change to the core's signatures. The binding
wraps the reader chain in a `BlockReader` whose `next_block` re-attaches to
the interpreter, calls `check_signals`, and turns a `KeyboardInterrupt`
into an error that travels out through the core's `Result`. Two things make
it more than a review fix. The interpreter is released for the whole pass
because rayon's threads deadlock on a caller that holds it, and the claim
that re-attaching between blocks is safe, since no rayon work is in flight
at that moment, has to be checked and not assumed. And it touches five
calculations in two binding crates, one new error case in the core, and
tests that send a real signal in the middle of a pass. It is recommended as
a plan of its own after this one merges, and it wants a sentence in a spec
about what a user sees when they stop a pass.
