# Report: how much variety each population holds

24 September 2026. It records how `docs/plans/diversity.md` was carried
out, on the branch `plan/diversity` in the worktree
`.claude/worktrees/diversity`, which branches from `spec/diversity` and
not from `main`, because neither the spec nor the plan is on `main` yet.

**The plan is under way.** This page is written while the work goes, so
what follows grows one work package at a time.

## Before the first task

The plan's "What has to be in place" was checked by running it, on the
owner's Apple M5 Pro under macOS 27.0, on 24 September 2026, at commit
664485e.

| what | command | what it gave |
|---|---|---|
| the format | `cargo fmt --all --check` | exit 0 |
| the core and the linalg crate | `cargo test --workspace` | `787 passed`, 2 ignored, and `149 passed` |
| the lints | `cargo clippy --workspace --all-targets -- -D warnings` | no warning |
| the two wasm targets | `cargo wasm-check` | exit 0 |
| the faer backend | `cargo test -p popnei --no-default-features` | `787 passed`, 2 ignored |
| the Python lint | `uv run ruff format --check` and `uv run ruff check` | `33 files already formatted`, `All checks passed!` |
| the Python suite | `uv run maturin develop && uv run pytest` | `499 passed` |
| R and its three packages | `Rscript -e 'packageVersion(...)'` | R 4.6.1, `vegan` 2.7.6, `adegenet` 2.1.11, `poppr` 2.9.8 |
| the panel | `ls tests/reference/stats/panel.vcf.gz` | there, with `panel_pops_bcftools.txt` |
| the file of the measurement | `ls /Users/jose/devel/popnei-bench/big.vars` | there, so task 4.1 does not write it |
| the TypeScript package | `npm run build` then `npm test` in `js/popnei` | exit 0, then `tests 325`, `pass 325`, `fail 0` |

Every count matches the one the plan wrote down, so the "more than" of
the deliverables is counted from these.

`npm run build` was run once at the start, as the plan asks, before any
task. Without it `npm test` fails in a fresh worktree on a file nobody
touched.

**`dadi` installs reproducibly**, which is the one thing that could have
sent the spec back to the owner. `uv venv --python 3.12` and `uv pip
install 'dadi==2.4.4'` succeeded, and the `dadi` they installed folded a
projected spectrum of a made up variant, 0.1333333333, 0.5333333333,
0.3333333333 and 0 in a draw of 6 called alleles. So the spec's check of
the folded spectrum stands as written and task 1.2 can be built on it.

The checks the plan says fail today do fail: `cargo test -p popnei --lib
-- diversity:: --list` prints `0 tests`, and `tests/test_diversity.py`,
`js/popnei/test/diversity.test.ts` and `tests/reference/diversity/` do
not exist.

## Work package 1: the reference numbers and the programs that give them

It finished as planned. `tests/reference/diversity/` now holds three
scripts and the six files they write, and every number of the spec's "How
it is verified" that a program outside popnei can give is in the
repository. The four tasks are five commits: `a527da8`, `48b5f5d`,
`42430ee`, `6eae424` and `bc8b7e2`.

### The deliverables, each with the command the orchestrator ran

| deliverable | command | what it gave |
|---|---|---|
| 1, the R script | `Rscript tests/reference/diversity/make_reference.R` | exit 0; `git status --short tests/reference/diversity` empty after a second run |
| 2, the Python script | `uv run python tests/reference/diversity/make_reference.py` | `done`, exit 0; the same emptiness after a second run |
| 3, the enumeration | `uv run python tests/reference/diversity/enumerate_private.py` | exit 0; of its 23 lines, 22 have a difference of `0` and the shared individual has `1/2`, which is the pair the closed form gets wrong |
| 4, the page | `grep` for each program and version in `README.md` | all five named, `vegan` 2.7.6, `adegenet` 2.1.11, `poppr` 2.9.8, `scikit-allel` 1.3.13 and `dadi` 2.4.4, with both Python 3.12 environments |

The stored numbers are the spec's. `panel_num_alleles.tsv` holds 2373,
2377 and 2384 alleles called with the standardized 1.9283948650,
1.9219209943 and 1.9197370844; `panel_private_alleles.tsv` holds 0, 0 and
1; `panel_variable_vars.tsv` holds 1173, 1177 and 1184 with the three
standardized ratios one less than the standardized allele counts, which is
the identity the spec's "How it is verified" of that item gives;
`panel_folded_sfs_dadi.tsv` holds 11 bins by 3 populations, the three
columns summing to 1200 short by 1.3e-11, 7.0e-11 and 4.0e-11, which is
eleven roundings of `dadi`'s own values and 5.8e-14 of 1200 at worst, so
the test of deliverable 4 of work package 3 compares that sum within a
tolerance and not exactly; `panel_fis_plain_allel.tsv` holds
-0.0237536998, -0.0258924700 and -0.0247472838.

The plan's own final check was run for this work package's data: both
`make_reference` scripts again from a clean tree, and `git status --short
tests/reference/diversity` showed nothing. Moving the two Python 3.12
environments aside and re-running rebuilt them in 5.5 s and reproduced
both stored files byte for byte; with both present the script takes 0.7 s.

### What was changed in the plan, and why

**The standardized private alleles of the panel get no file.** The work
package said that every number of the spec's "How it is verified" becomes
a file. Three cannot: 0.0112196177, 0.0099715392 and 0.0089014974, the
standardized private alleles of the three populations at a draw of 20. No
program outside popnei computes any standardized value, which is what the
plan's own final check already said in another way. They stay literals of
a pytest test, which is what deliverable 3 of work package 3 asks for
anyway, so no check got weaker. Task 1.1 found this.

**Two sentences about the project's Python were false**, one in the plan
and one in the spec, and task 1.4 found them. They are in "What the owner
should know" below, because the choice they leave open is the owner's.

### What the review found

Four reviewers were sent over `bc8b7e2`, the categories `spec`, `tests`,
`numbers` and `errors`. `api`, `architecture` and `binding` were not sent:
this work package adds no type, no signature, no reader, no thread, no
dependency of a crate and nothing in either binding crate.

**Every stored number survived.** Two reviewers recomputed all of them
from `panel.vcf.gz` in exact rational arithmetic, independently of the
scripts: `vegan`'s standardized allele counts agree within 1.7e-16
relative, the standardized ratios within 3.2e-16, `scikit-allel`'s F_IS
within 7.3e-15 and `dadi`'s 33 spectrum entries within 6.7e-14, and the
counts 2373/2377/2384, 1173/1177/1184 and 0/0/1 were reproduced exactly.
All 18 enumerated rationals were reproduced twice more, once by a brute
force over labelled gene copies. One reviewer ran 25 mutations across the
three scripts and every one stopped the script it belonged to.

**The enumeration of the private alleles assumes the independence it said
it did not.** Three of the four reviewers found this separately, with the
same evidence. `by_enumeration` lists each population's draws from that
population's own allele counts and weights a combination by the product
across populations, which is the independence the closed form assumes. On
the spec's own overlap case, two populations that are the one diploid
individual `0/1` at a draw of one allele, it gives 1/2, the same as the
closed form, where "The cases" says the truth is 0. So the spec's "It
assumes nothing" was wrong, and so was the sentence crediting the
enumeration with having found the overlap, which was worked out from the
formula.

The check keeps its value, and the spec now says what that value is: the
enumeration never writes the closed form down, so it catches any error in
the algebra, and changing one factor of the closed form made 16 of the 18
pairs differ then and makes 19 of the 22 differ now that four pairs have
been added. For those pairs the product measure is sound rather than
circular, because their populations share no individual and draws from
disjoint sets of gene copies are independent as a fact of the sampling.
The spec was corrected in `e1ff8e9` and a pair enumerated over labelled
gene copies was added, which gives 0 against the closed form's 0.5 and is
the one pair of the file whose difference is not 0.

**Four guards were weaker than they read**, none of them able to store a
wrong number today, all of them able to let one through on the next
dataset or to report the wrong cause:

- Reusing a reference environment is decided on the installed package
  version alone. A reviewer asked for Python 3.13 with the 3.12
  environments present: the script ran on 3.12, printed nothing and
  exited 0, so the interpreter that these numbers' provenance rests on can
  go stale with no message.
- The check that each spectrum column sums to 1200 cannot fire, because
  the eleven per bin comparisons above it already pin the sum to 6e-10 of
  1200 against its threshold of 1e-6. Dropping bin 0 was caught only once
  the bin check was disabled.
- The Python panel reader accepts three inputs the R one refuses: no
  `#CHROM` line, a `FORMAT` column that is not `GT`, and a populations
  file that does not name the VCF's individuals. Dropping one line of the
  populations file silently drops that individual from its population, and
  only a spec literal catches it, naming a spectrum bin and not the cause.
- The 18 pairs are a thin fixture for what they become, the literals of
  cargo tests in work package 3: seven of the eighteen values are exactly
  0 and two exactly 1, and fourteen of the eighteen population slots have
  exactly 4 called alleles, so only two slots could catch an
  implementation that confused one population's called alleles for
  another's. Two cases were added, one of a single population and one of
  three populations with three different counts.

**Everything above was fixed**, in five commits: `0b5c614` for the R
script, `dfc8872` for the Python one, `9babbd0` for the enumeration, and
`87a9b7b` and `5bd51ee` for the page, which took two because the second
waited for the enumeration's new shape to be committed before describing
it. The refusals were provoked to see
them fire, and the orchestrator provoked one of them itself: a populations
file one line short now exits 1 naming the file, the count and the missing
individual, where before it produced a wrong spectrum bin. After the fixes
every deliverable's check was run again and all four pass, the plan's final
check leaves `git status --short tests/reference/diversity` empty, and the
seven commands of the `coding` skill give what they gave before the work
package, no Rust having been written: 787 and 149 cargo tests, 787 on
faer, 499 pytest, clippy and both wasm targets clean, ruff clean.

**Not taken, one finding.** That a failure between two `write` calls could
leave one new file and one old. It exits non zero and `git status` shows
it. The Python script took the fix, a temporary name and a rename, in four
lines of a shared writer. The R script did not: making it safe there needs
the temporary files cleaned up on every exit path, which its writer judged
more than a few lines, and the orchestrator agreed rather than spend a
second round on a case that cannot arise without a full disk or a
permission change.

**One thing the fixes turned up that the review had not.** Closing the
enumeration finding meant writing a guard, and testing that guard found
that it could pass for the wrong reason: a mutation that drops one gene
copy of a diploid individual leaves the shared individual's pair reading 0
against 0.5, which is the right answer arrived at wrongly, because both
populations then hold a single copy of the same allele. The guard now also
checks that the labelled copies of each population hold the called alleles
the closed form was given.

**One number was wrong in the work's own account of itself**, and it is
recorded because it was repeated to the owner before it was checked.
Writing 1/12 to 16 decimal places was said to make a cargo test that could
only fail. It is two units in the last place from the nearest float64, not
one, and that is 3.3e-16 relative, far inside the 1e-12 the spec compares
these values within, so such a test would pass. The rule of 17 places
stays, for the reason the page beside the scripts already gave: the file
then holds each value at the float64 nearest its exact rational.

### What the owner should know

**`dadi` builds on the project's Python, and both the spec and the plan
said it does not.** The spec's "How it is verified" of the folded spectrum
said `dadi` does not build on the project's Python, "3.14 with the free
threading build". Measured in this worktree on 24 September 2026:
`.python-version` holds 3.14.5, that interpreter reports
`Py_GIL_DISABLED` 0 and `sys._is_gil_enabled()` `True`, and
`pyproject.toml` says the patch version is written out precisely because
`uv` reads a bare `3.14` as the newest interpreter of that minor version,
which on this machine is the free threaded 3.14.7. `dadi` 2.4.4 installed
into a 3.14.5 environment folds a projected spectrum, giving the same four
values as on 3.12. `scikit-allel` 1.3.13 installs there too, and the lock
it re-enables was never off on that interpreter. What does fail is 3.14.7,
where `nlopt` compiles from source and stops for want of `cmake`, which is
not installed on this machine; the same missing `cmake` is why `hierfstat`
does not install, through `RcppParallel` and then `gaston`. The spec was
corrected in `f6fe8f9` and the plan in `766c9ee`.

**So one choice is open, and it is the owner's.** The spec's opening
records their decision of 24 September 2026 that the five reference
programs become development dependencies of popnei. The plan departed from
that for the two Python ones, on the two premises above, both false. The
environments the reference script makes were kept, on the one reason that
survives: the stored numbers then come from `scikit-allel` 1.3.13 and
`dadi` 2.4.4 whatever popnei's development dependencies later hold, and no
test of popnei imports either package. The alternative is to put both in
`pyproject.toml` and take about 60 lines of environment building out of
`make_reference.py`, which costs a `scikit-allel` import in the project's
environment that no test needs. Keeping the environments is the
recommendation, and nothing later in the plan turns on it, so the plan went
on to work package 2 without waiting.

### How the work went

Not for the owner, who can stop here. It is for whoever next revises a
skill or writes a plan.

**What the review cost against what the work cost.** The four tasks took
494000 tokens of subagent context to write. The four reviewers took 413000
and the four rounds of fixes 197000, so checking and repairing the work
cost 610000, 1.24 times what writing it cost. That ratio is the number to
size the tasks of the next plan with, and it is high because this work
package is reference data: almost every finding was about a guard that
could not fire or a claim that was not true, which is what a reviewer is
good at and what a writer cannot see in its own work.

**Four reviewers were the right number and three of them found the same
thing.** The independence of the enumeration was found by `spec`, `tests`
and `numbers` separately, each with its own evidence, and the `code-review`
skill is right that this is evidence rather than repetition. But it also
means that a work package of three scripts could have been reviewed with
two categories and lost little: `errors` was the only one whose findings no
other reviewer reached, and `spec` and `tests` between them found
everything else.

**The plan sent 1.1, 1.2 and 1.3 side by side and that worked, at one
cost.** The three scripts had to be made self contained for it, since the
throwaway scripts they grew from read each other's output. That was worth
saying in the prompts and it was: each script derives what it needs from
the panel. The cost was that `git status --short tests/reference/diversity`
is not a usable check while three agents write in one directory, and each
writer had to be told to look only at its own files. A plan that wants
three tasks in one directory should give the check per file.

**A premise of the plan was false and the work found it, not the review.**
The plan and the spec both said the project's Python is a free threading
build. Task 1.4 checked it while writing a page about it, which is what
checking every claim of a document by running it is for. The `writing`
skill's rule that a number goes on the page only if it was seen printed is
what caught a wrong sentence in two other documents.

**Three sentences of the spec described checks that were stronger than the
checks were.** The enumeration's independence, the 14 pair identity that no
script ran, and the panel's standardized private alleles that no file held.
All three were written while the spec was being written, from work done in
a session, and none survived somebody trying to run it. A spec that says "X
is checked by Y" is worth a task of its own that runs Y.

## Work package 2: the counts with no draw, through the four layers

It finished as planned. All six tasks are done, all eight deliverables
pass, and the seven reviewers found twenty findings that held, every one of
them fixed and every deliverable run again afterwards. Five of the seven
had to be run twice, having been cut off by a rate limit the first time,
which cost the work nothing but the wait.

A user can now call `calc_pop_diversity` in Python and `calcPopDiversity`
in TypeScript and get, per population, the alleles it called, the private
ones among them, the variants that vary in it and F_IS, each as a total
with its mean or its ratio, with the counts of the pass beside them. The
draw and the folded spectrum are work package 3.

### The deliverables, each with the command the orchestrator ran

| deliverable | command | what it gave |
|---|---|---|
| 1, the module | `cargo test -p popnei --lib -- diversity:: --list` | `35 tests`, where the plan asks for 14 or more and where it printed `0 tests` |
| 2, the worked example | the five named cargo tests | the alleles 9 and 8 with means 2.25 and 2, the private 2 and 1 with means 0.5 and 0.25, the variable 3 and 2 with ratios 0.75 and 0.5, F_IS 0 and 0.3478260870 |
| 3, which variants count | `four_of_the_six_variants_...` and `a_min_num_individuals_of_zero_...` | both pass |
| 4, the four errors | the four named cargo tests | all four pass |
| 5, the Python function | `uv run pytest tests/test_diversity.py` | `18 passed`, where the plan asks for 8 or more |
| 6, one population | `-k one_population` | `1 passed` |
| 7, the TypeScript function | `npm run build && npm test` in `js/popnei` | `tests 336`, `pass 336`, `fail 0`, against 325 before |
| 8, the agreement with `stats` | `-k stats_module_counts` | `4 passed`, one per `poly_threshold` |

The counts across the plan so far: 787 cargo tests, 499 pytest and 325
node when it started, and 822, 517 and 336 now, with `cargo test -p popnei
--no-default-features` giving the same 822 on the faer backend.

Deliverable 5 reads the panel's counts and means from the four stored
files of work package 1 rather than from literals typed into the test, and
asserts the three unbiased F_IS of the spec within 1e-12 of the value: the
differences are 5.2e-14, 1.2e-14 and 0.

### What was changed in the plan, and why

**The plan's final check gained the public interface.** Task 2.1 gave
`PopDiversity` only the three accessors whose numbers the pass computed at
that point, rather than writing the whole impl block of "The Rust
interface" with ten of them returning 0 or `None`. That is the safer order,
since no caller can read a number that is not yet what the spec promises,
and each later task adds its own accessor with its computation. But nothing
fails if one is forgotten, so the plan's final check now names the ten and
requires them with the signatures the spec gives.


### What the review found

Seven reviewers over `053562f`: all seven categories apply, because this
work package builds a module through the core crate, both binding crates,
the Python package and the TypeScript package. They cost 1.15 million
tokens of subagent context, against 1.24 million for the six tasks that
wrote the work, and they found two ways of breaking the code that pass all
1176 tests of the three suites.

**Everything the four statistics compute was confirmed against
independent arithmetic.** One reviewer recomputed the panel from
`panel.vcf.gz` in exact rational arithmetic: the four counts of all three
populations exactly, and F_IS within 3.9e-14 relative of exact. Two went
outside the panel's regime, which has no half called genotype and no
variant of three alleles: a random VCF of 300 variants with four alleles
and a quarter of the genotypes missing, `many.vcf` with its 54 triallelic
variants under three population layouts including overlapping populations
and a population of one, and a hand-built tetraploid file. popnei matched
an exact oracle on every count, every mean to the bit, and F_IS to
3.4e-14. A third measured the memory of the pass and found the trap the
work package warned of avoided: 4.56 MB of allocation at a block of 10000
rows with 50 populations is 157 chunks of 25.6 KB and not 10000 of them,
which would be 256 MB, with one live at a time.

**Two mutations passed every test.** They are the findings that justify
the review.

- **A variant with only one of the two heterozygosities.** The rule that
  such a variant enters neither of F_IS's means can be replaced by
  counting it as zero in both, and all 822 cargo, 18 pytest and 336 node
  tests still pass. It needs half called genotypes of different alleles,
  which no fixture had. On two diploid individuals with `0/1 0/1` then
  `0/. 1/.` the pass gives -0.4999999999999998, which the orchestrator
  confirmed, and the broken version 0.4.
- **The allele count on the per-population path.** Bounded so that it does
  not trip the whole-row path's 128, exactly one test fails, and nothing
  at all guarded the variable-variants side: a population whose only
  genotype is `3/3` would have been called variable.

**A user who asked for a draw got a wrong number and nothing said so.**
With a valid `num_called_alleles` the three standardized values came back
NaN, the two counts of the draw 0 and the spectrum absent, which are the
values the documented contract reserves for a draw above every
population's called alleles and for a statistic nobody asked for. So the
answer was wrong and read as an answer, against the rule the owner gave on
21 September 2026. It came from the orchestrator's brief to tasks 2.5 and
2.6, which said to pass the argument through and leave those fields as the
spec describes them. Both packages now refuse the draw and the spectrum
outright, with a `NotImplementedError` in Python and a thrown `Error` in
TypeScript, naming work package 3. Three strict `xfail` tests in pytest
and one guard test in node hold the finished assertions and fail the day
the pass fills the values, so the refusal cannot outlive its reason.

**Nothing exercised the parallel reduction or the wasm path.** A chunk is
64 rows and every cargo fixture had at most 6, so `par_chunks` yielded one
chunk in every cargo test and no test varied the thread count, where eight
other modules each have one that does. The panel tests of pytest and node
do span 19 chunks, but reversing the order of the partials leaves the
panel F_IS bit-identical, so they would pass over an unordered reduction.
And `add_the_chunks_one_by_one`, which is the whole wasm path, was run by
no cargo test at all: deleting a line from it left all 822 cargo tests
passing while four node tests failed.

**A table and the thing it described had drifted apart, twice.**
`make_reference.R` divided the mean private alleles by the variants each
population counted while writing that column under a header naming the
variants every population counted, two different numbers that happen to
coincide at 1200 on this panel, so the reference file encoded a rule that
is not the spec's. And `DiversityStats::of_name` swallowed any name not in
its five match arms, which a reviewer demonstrated by adding a sixth name
to the table: the refusal then listed the very name it had just refused.

**Eleven smaller findings held**, each fixed: the list of reasons F_IS is
NaN missed the commonest one in all three layers; `of_name` returned an
`Option` so both bindings wrote the same sentence out; the core took a
`stats` naming no statistic and made a whole pass for nothing, with the
refusal living in both packages instead; the message for a spectrum
without a draw told a user to change an argument they had not written;
TypeScript gave "a whole number of 0 or more" for a `num_called_alleles`
where 0 and 1 are refused; the mean private alleles' divisor was guarded
by no test in either package; two conversions carried a dead
`unwrap_or(MAX)` that would have inflated a count rather than raised;
nothing tested F_IS at a ploidy other than 1 or 2, which is the `coding`
skill's own example of a suite that cannot tell which of two things the
code reads; two cases of the spec had code and no test at the package
layer; the module doc said one statistic was left where four things were;
and a comment in both bindings said every statistic comes out the same
whatever the blocks are, which is false for F_IS, whose panel value moves
6.7e-16, 5.2e-14 relative, between block sizes.

**Two findings were not taken.** The `saturating_add` calls on the per
variant counters: the bounds are true, the fastest growing needing 1.44e17
variants where the objectives' largest dataset reaches 1.28e8, and `stats`
uses the same shape for the same counters, so changing it here alone would
make the two inconsistent; if the convention changes it is crate-wide.
And the `Option<f64>` the `coding` skill asks for inside the core where
this module returns `Some(NaN)`, which is under "What the owner should
know" below.

**Three of the fixers refused the shape the orchestrator asked for, and
were right each time.** Comparing F_IS bit for bit across thread counts
cannot fail, because the division absorbs a one-ulp difference in the
sums, so the test compares the two sums through a test-only accessor. A
fixture of repeating rows makes every 64-row chunk identical, so no order
of them differs. And 200 variants in four chunks cannot tell one thread
from eight, because rayon splits them the same way; the test uses 2000
variants and 32 chunks, where a broken reduction does show.

### After the fixes

The fixes are five commits: `66f2fff` the reference script, `abde1f0` the
core and both binding crates' call sites, `3992bb5` and one follow-up for
the Python layer, and `41fcb68` the TypeScript layer. Every deliverable
and every check was run again.

| what | command | what it gave |
|---|---|---|
| deliverable 1 | `cargo test -p popnei --lib -- diversity:: --list` | `44 tests`, from 35 |
| deliverables 2, 3, 4 | the named cargo tests | all pass |
| deliverable 5 | `uv run pytest tests/test_diversity.py` | `20 passed, 3 xfailed` |
| deliverable 6 | `-k one_population` | `1 passed` |
| deliverable 7 | `npm run build && npm test` in `js/popnei` | `tests 337`, `pass 337`, `fail 0` |
| deliverable 8 | `-k stats_module_counts` | `4 passed` |
| the seven checks | as the `coding` skill gives them | `831 passed` with 2 ignored and `149 passed`, `831` again on faer, `519 passed, 3 xfailed`, the other four clean |
| the plan's final check | both reference scripts again | `git status --short tests/reference/diversity` empty |

Across the plan: 787 cargo tests, 499 pytest and 325 node when it started,
and 831, 519 and 337 now.

### What the owner should know


The owner answered all three on 24 September 2026, after work package 2
was done and before work package 3 began. Two of the answers are built and
the third is written into the `coding` skill.

**The plain call works, and gives the four statistics that need no draw.**
The default of `stats` was every statistic, the folded spectrum needs a
`num_called_alleles`, and asking for the spectrum without one is a
`ValueError`, so `calc_pop_diversity(variants)` refused itself and its
message named an argument the user had not written. The default is now
`PopDiversityStat.WITHOUT_A_DRAW`, the four, with
`DiversityStats::WITHOUT_A_DRAW` beside `ALL` in the core and both binding
crates reading it from there so the four are named once. Asking for the
spectrum by name without a draw is still refused, which is a user asking
for something popnei cannot give. Spec in `71c63cb`, code in `d24083f`.

**The choice of what to compute is given by name.** Three calculations of
popnei take both a set of populations and a choice of what to compute, and
they disagreed about which came second: `calc_per_var_distribs` takes
`stats` there and `calc_pop_dists` takes `pops`, as this one did, so the
same position meant two things. `calc_pop_diversity` now takes `stats`,
`num_called_alleles` and `min_num_individuals` by name. TypeScript needed
nothing, taking everything but the variants in an options object already.

**`calc_per_var_distribs` and `calc_pop_dists` are to follow, and that is
asked of the owner rather than done.** The owner approved the convention
for all three; only this module's function was changed, because it has no
users yet and the change is free, where theirs is a change to two settled
specs and two built modules and so is outside what this plan said it was
building. Until it is done the inconsistency is smaller and not gone: a
caller of those two can still pass the choice positionally and a caller of
this one cannot. `docs/specs/diversity.md` says so under "Its Python
function", so nobody reads the half-done state as finished.

**The rule against NaN in the core now says which values it is about.**
`fis` returns `Some(NaN)` for a population with no F_IS, where the `coding`
skill asked for an `Option<f64>` inside the core. The skill's own reason is
that NaN travelling through Rust arithmetic hides where it was born, which
is about values that go on to be computed with, and it was written as
though it covered every value. The rule now says so, and a value the core
has finished with may be NaN when its doc comment gives every reason it can
be one. That was the cheaper side by a wide margin: the alternative was
`Option<Option<f64>>` or a new enum here, about thirty lines, and the same
again in `gwas` and `kinship` for the rule to mean anything, and those two
have a second reason, an `Option<f64>` being sixteen bytes against eight in
columns as long as the dataset, which the orchestrator had missed when it
first put the question. No code changed; `bd1514c` changed the skill.

**One question of work package 1 is still open**: whether `scikit-allel`
and `dadi` become development dependencies of popnei, as the spec's opening
records the owner deciding, or stay in the environments
`tests/reference/diversity/make_reference.py` builds, which is what the
plan does and what the work kept.

## Work package 3: the draw, through the same four layers

It finished as planned, with one task added. The module is now what the
spec describes: a user gives `num_called_alleles` and gets the alleles
called, the private ones and the variable variants as a draw of that many
called alleles would show them, and the folded site frequency spectrum
projected to it. Six tasks became seven, and the review in all seven
categories found twenty-one findings that held, every one fixed.

Across the plan: 787 cargo tests, 499 pytest and 325 node when it started,
and **868, 530 and 345** now, the 868 the same on the faer backend.

### The deliverables, each with the command the orchestrator ran

| deliverable | command | what it gave |
|---|---|---|
| 1, the draw arithmetic | `cargo test -p popnei --lib -- chance_a_draw` | `8 passed`, `vegan`'s eight per variant values and their two means |
| 2, 3, 4, 5 | `uv run pytest tests/test_diversity.py` | `31 passed` |
| 6, the worked example at a draw of 4 | the five named cargo tests | all pass |
| the seven checks | as the `coding` skill gives them | `868 passed` with 2 ignored and `149 passed`, `868` again on faer, `530 passed`, `35 files already formatted`, the rest clean |
| the TypeScript package | `npm run build && npm test` in `js/popnei` | `tests 345`, `pass 345`, `fail 0` |
| the plan's final check | both reference scripts and the enumeration again | `git status --short tests/reference/diversity` empty |

The numbers the module now gives, against what is stored: the panel's
standardized allele counts agree with `vegan` within 2.3e-16 of the value
and the standardized ratios of variable variants within 1.2e-16; the
standardized private alleles give the spec's ten decimals; all 33 spectrum
values agree with `dadi` within 6.674682e-14 and with exact rational
arithmetic within 4.0e-16.

### What was changed in the plan, and why

**A task was added, 3.4b.** Task 3.4 found that nothing bounded
`num_called_alleles` above 2. The spectrum has one bin per count of the
rarer allele up to half the draw, so a draw at the top of what a `u32`
holds asks for 2147483648 bins for each population, 51 GB of result over
three of them, and popnei would have died allocating it rather than say
what was wrong. The owner decided the refusal on 24 September 2026: a draw
above the individuals of the dataset times the ploidy, which is every gene
copy it holds, is refused naming the largest it allows. The case below that
bound which a dataset's missing data leaves unfillable is untouched and
still gives NaN and zeros. On the panel's 200 diploid individuals, 400 is
taken and 401 refused, and a test holds both sides, because a refusal
creeping one step further would eat the case the spec protects.

### What the review found

Six reviewers over `df6c872`, in all seven categories, paired where the
work was small. **Four mutations passed every one of the 1729 tests of the
three suites, and one of them was a wrong number a user could reach.**

**The folded spectrum lost whole variants at a large draw.** This is the
worst defect the plan found. The recurrence for the bins started at the
smallest term of its range, which underflows to a subnormal and then to 0
long before any bin's own value would, and once it is 0 every later step
stays 0, so the variant contributed nothing to any bin. The failure is not
monotone, which is why nothing caught it. Measured on 600 diploid
individuals, every genotype `0/1`, three variants, where the column must
sum to 3:

| draw | the column summed to |
|---|---|
| 400 | 3.000000000000002 |
| 560 | 3.000000000000001 |
| 580 | 3.043348989886391 |
| 584 to 616 | 0.0 |
| 620 | 3.043348989886391 |
| 1000 | 2.9999999999999973 |

The reviewer measured the threshold with an `f64` replica against exact
rationals: about 525 diploid individuals in a population is where a draw
first loses more than 1e-12 of the mass, and about 560 where the variant
vanishes. `docs/objectives.md` puts popnei's largest dataset at ten
thousand individuals, so this was well inside the range popnei claims, in
the statistic the spec calls the input of the programs that fit a
demographic model. And the result contradicted itself in silence:
`num_vars.in_draw` said three beside a column of zeros, with no error,
against the rule the owner gave on 21 September 2026. All 858 cargo tests
passed, because the panel's largest case is 168 called alleles at a draw of
20.

The chance is now carried as a mantissa with the power of two it is worth,
kept in range by exact multiplications by 2^512, which reuses the one
product the module already has rather than needing a second. It is
bit-identical wherever the old arithmetic stayed normal, and it was checked
against exact rationals over 30633 combinations of the called alleles, the
rarer allele's count and the draw, worst bin 1.2e-14. Every draw of the
fixture above now sums to 3 within 5.3e-15, which the orchestrator
confirmed.

**Three more mutations passed every test.** The spectrum could read the
rarest allele instead of the major one, which changes `[0.1, 0.9]` into
`[0.5, 0.5]` on a variant of counts 3, 2 and 1 at a draw of 3 and was
invisible because every asserted fixture is biallelic or has equal counts,
and for two alleles the fold makes the two choices the same. The
standardized private alleles could be summed over rows not in the draw for
every population, 50 per cent wrong on a probe and invisible because the
panel has every variant in every population's draw. And the largest
allowed draw could hard-code a ploidy of 2, which refuses a legitimate
draw on tetraploid data and accepts an impossible one on haploid, the line
task 3.4b was built to hold, undone by a fixture set that is entirely
diploid. Each now fails exactly one test.

**The memory of the spectrum was measured, not read, and the earlier fix
had solved half of it.** The review of work package 2 predicted about 7 MB
of partials for a particular block; task 3.4 replied with one flat vector
per chunk instead of one per population and reported the allocation count
falling from 7850 to 157. The `architecture` reviewer counted bytes and
found the prediction had come true anyway, because the reduction collects
one set of partials per chunk and holds them all alive:

| individuals | block | populations | draw | before | after |
|---|---|---|---|---|---|
| 200 | default | 3 | 20 | 0.293 MB | 0.045 MB |
| 500 | 10000 rows | 50 | 180 | 7.041 MB | 2.342 MB |
| 1000 | 5000 rows | 50 | 2000 | 33.011 MB | 15.871 MB |
| 10000 | 500 rows | 50 | 20000 | 36.481 MB | 36.481 MB |

The chunks are now read in groups of a few per thread and the partials
added in index order, which keeps the order the float sums depend on. The
last row did not move and the fixer said so rather than claiming it: that
block holds eight chunks, fewer than one group, so the grouping has nothing
to bite on. Whether to read fewer chunks than the pool has threads is a
trade nobody has measured.

**Seventeen smaller findings held**, each fixed. A ploidy the reader states
was refused with the `stats` module's message, naming the `ploidy` and
`exponent` arguments that `calc_pop_diversity` does not have, where
`docs/specs/dists.md` had met the same thing and given `calc_pop_dists` a
case of its own. `Totals::of` bounded the bins per population but nothing
bounded the populations, so what kept its multiplication from saturating
was the machine's memory and the failure was the standard library's panic
rather than an error of popnei. The TypeScript package refused a draw below
2 itself, where the spec puts that refusal in the core and Python leaves it
there, so the one rule was written twice and the two languages gave
different sentences for the same mistake. The TypeScript suite sorted its
populations, so it could not catch a reordering that would hand every
population another's counts and spectrum; a sort added to the package left
all 17 of its diversity tests passing. The four statistics that need no
draw had no TypeScript name, so an application would copy the list the core
owns. The node suite had no fixture where the four counts of variants
differ. The spectrum crossed to Python through three copies. A `///`
comment in the private `_core` module became a `__doc__` that the same file
says belongs to the package. The `# Errors` of `calc_pop_diversity` and the
spec's list disagreed in both directions. Two public items both binding
crates call were in no spec item. `chance_a_draw_misses_an_allele`
contradicted its own doc comment in a corner no caller reaches. And four
test names and three comments claimed more than the thing they sat on
proved, which the fixer found by auditing its own work after one was
pointed out.

**Nothing was set aside.** Every finding of this work package held and
every one is fixed.

### What the owner should know

**The spec was corrected nine times during this work package**, every time
because a task or a reviewer tried to assert what it said and found it
false, incomplete or ambiguous. The three worth naming: a population
compared against a copy of itself was said to have no standardized private
alleles at any draw size, which is true of the totals and false of the
standardized value, because the closed form treats a population's draw and
its copy's draw as independent; nothing bounded the draw; and one
paragraph gave popnei's agreement with `vegan` measured against exact
arithmetic as though it were measured against the stored values, a number
of one quantity given as a number of the other, which the orchestrator had
put into a task's brief and which reached a doc comment.

**Three things outside this plan are waiting**, none urgent and none this
plan's to do.

`crates/popnei/src/io/vars.rs` gives the ploidy its metadata states with no
ceiling, which is what lets a ploidy popnei cannot read reach any pass.
That is `docs/specs/io_vars.md`'s question.

`docs/reports/diversity-method/panel.py` hardcodes paths into a worktree
that no longer exists, so the script this spec names as the source of its
unbiased F_IS and its standardized private alleles cannot be rerun as it
stands. Two reviewers and the orchestrator each had to copy it to reproduce
those numbers. It is the provenance of values that tests assert, so it is
worth an hour.

`calc_per_var_distribs` and `calc_pop_dists` still take their choice of
what to compute positionally, where `calc_pop_diversity` now takes it by
name. The owner approved the convention for all three on 24 September 2026
and only this module's function was changed, its having no users yet;
theirs is a change to two settled specs and two built modules.
