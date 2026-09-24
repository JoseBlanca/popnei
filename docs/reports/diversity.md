# Report: how much variety each population holds

24 September 2026. It records how `docs/plans/diversity.md` was carried
out, on the branch `plan/diversity` in the worktree
`.claude/worktrees/diversity`, which branches from `spec/diversity` and
not from `main`, because neither the spec nor the plan is on `main` yet.

**The plan is done, and what is asked of the owner is the merge.** Nothing
was merged into `main` and nothing pushed; the branch is `plan/diversity`,
whose last commit of the work is `4b21868`, with this page and the plan's
ticks after it.

**What exists now that did not.** The `diversity` module of
`docs/specs/diversity.md`. A user calls `calc_pop_diversity(variants,
pops=...)` in Python or `calcPopDiversity` in TypeScript and gets, for each
population, the alleles it called and the private ones among them, the
variants that vary in it and F_IS, each as a total with its mean or its
ratio. Giving `num_called_alleles` adds the first three as a draw of that
many called alleles would show them, which is what lets populations of
different sizes be compared, and the folded site frequency spectrum
projected to that number, which is the input of the programs that fit a
demographic model. It runs natively, in a browser tab through the wasm
package, and under pyodide, the Python of a tab. The numbers are checked
against `vegan`, `adegenet`, `poppr`, `scikit-allel` and `dadi`, whose
outputs are stored under `tests/reference/diversity/` with the scripts that
made them, and against exact rational arithmetic. The suites went from 787
cargo tests, 499 pytest and 325 node to 868, 530 and 345.

**What the work found that changes what the spec said.** The folded
spectrum, which "Speed" called the one part that could dominate the pass,
is 20 to 23 in 100 of it at a draw of 200 called alleles. "Speed" now holds
a measurement and a speed target where it held a sentence saying no
measurement had been made.

**What is open, none of it blocking.** Five things, each with its options
and a recommendation under "What the owner should know" of the work package
it belongs to: whether `scikit-allel` and `dadi` become development
dependencies of popnei or stay in the environments the reference script
builds, in work package 1; `crates/popnei/src/io/vars.rs` giving the ploidy
its metadata states with no ceiling, `docs/reports/diversity-method/panel.py`
hardcoding paths into a worktree that no longer exists so that the script
the spec names for the provenance of two of its numbers cannot be rerun,
and `calc_per_var_distribs` and `calc_pop_dists` still taking their choice
of what to compute positionally where `calc_pop_diversity` now takes it by
name, all three in work package 3; and whether the thirteen Python
harnesses under `crates/popnei/benches/` come under `ruff`, in work package
4. One thing is left for a performance review of its own and is worth more
than anything the plan optimized: three of the statistics compute the same
products over and over on a variant of two alleles, 0.208 s of a 0.505 s
pass.

This page was written while the work went, so what follows is one work
package at a time, and the last section of each is for whoever next
revises a skill or writes a plan and not for the owner.

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

## Work package 4: the measurement and the browser

It finished as planned, with one change to how the two tasks were run.
`docs/specs/diversity.md` now has a speed target where it had a sentence
saying no measurement had been made, and `calc_pop_diversity` runs under
pyodide, the Python of a browser tab, with the numbers of the spec's worked
example. The two tasks are ten commits, four of the work and six of the
fixes the review asked for.

**The measurement refutes what the spec feared.** "Speed" said the folded
site frequency spectrum, the one statistic whose work grows with the draw
rather than with the alleles a population called, was "the one part that
can dominate". It is 20 to 23 in 100 of the pass at a draw of 200 called
alleles, and the three standardized values together are 1.9 times it. The
reason is that the code already does the thing the spec left open: the
weights of a variant are computed once and shared by the bins of its
spectrum.

### The deliverables, each with the command the orchestrator ran

| deliverable | command | what it gave |
|---|---|---|
| 1, the measurement | the two harnesses over `big.vars` and the panel, and `docs/reports/perf-diversity-2026-09-24.md` | the report, with the target written into "Speed" in `86c1775`, before the code of task 4.2 |
| 2, the shared weights | `cargo bench --bench diversity_pass -- --draw 200 --pops 3`, one thread, shared and then unshared | 0.243 s against 11.042 s for the spectrum alone and 0.500 s against 11.447 s for all five; the sharing is kept |
| 3, the wasm build and node | `cargo wasm-check`, then `npm run build && npm test` in `js/popnei` | exit 0; `tests 345`, `pass 345`, `fail 0`, of which `js/popnei/test/diversity.test.ts` holds 20 |
| 4, the wheel under pyodide | `bash scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` | exit 0, printing the worked example's 9 and 8 alleles, 2 and 1 private, 3 and 2 variable, F_IS 0 and 0.3478260870, and in a draw of 4 the spectra 1, 3, 0 and 1.0666666667, 1.5333333333, 0.4 |

The plan's own final check, all three parts. `Rscript
tests/reference/diversity/make_reference.R`, `uv run python
tests/reference/diversity/make_reference.py` and `uv run python
tests/reference/diversity/enumerate_private.py` each exit 0 and leave `git
status --short tests/reference/diversity` empty. All fourteen methods of
the impl block of "The Rust interface" are in
`crates/popnei/src/diversity.rs` with the signature the spec gives, checked
by a script that reads the block out of the spec and looks for each one.

The counts across the plan: 787 cargo tests, 499 pytest and 325 node when
it started, and **868, 530 and 345** now, the 868 the same on the faer
backend, with 2 ignored.

### The target that is now in the spec

`calc_pop_diversity` is measured against `calc_per_var_distribs` of
`docs/specs/stats.md` over the same file with the same populations at the
same thread count, which the plan asked for: that pass makes the same read
and the same counts of how often each population called each allele, and
then does far less arithmetic on them. With the four statistics that need
no draw, and with all five at a draw of 20, at most 1.1 times it; with all
five at a draw of 200, at most 2 times it. Today it is 0.78 to 1.02 times
and 1.16 to 1.61 times, so the target has headroom at every thread count.

It is a ratio and not a number of seconds because both sides are timed on
one machine on one day, and a ratio does not go stale when the machine
does.

### What was changed in the plan, and why

**The two tasks ran one after the other where the plan said they could run
side by side.** Task 4.1 takes wall times and task 4.2 runs the release
build of `crates/popnei-js` and the build of the pyodide wheel. The
`performance-review` skill takes a wall time with nothing else building on
the machine, because two builds at once measure each other. 4.1 went first,
which deliverable 1 asks for anyway: the commit that writes the target into
the spec comes before the code of any later task.

**Task 4.1 was carried out in two parts**, the harnesses and every number
first, then the target in the spec followed by the experiment of
deliverable 2, so that the target was set against the code as it was
measured.

**No performance reviewers were sent.** The `performance-review` skill
sends one reviewer per category to look for candidates. Work package 4
asked for one measurement and one experiment and named them both, and its
"What could go wrong" says that anything beyond deliverable 2 is a
performance review of its own. The report says what that leaves out.

### What the review found

Six reviewers over `f9064bc`, in the categories `spec`, `tests`, `numbers`,
`errors`, `architecture` and `binding`. `api` was not sent: this work
package adds no type, no signature and no public name, and changes no line
of the library but one doc comment. They found twenty-six findings that
held, every one fixed. Not one was in code that runs: the only line of the
library this work package changed is a doc comment, and what the review
found is that the measurement was right and the documents reporting it were
not.

**Fourteen numbers and sentences of "Speed" and of the performance report
were not what was measured.** The three that mattered:

- **A row of the differencing table is labelled with the wrong four
  statistics.** `the four that need a draw | 200 | 0.371, 0.393, 0.396 s`
  is the four that need *no* draw, run at a draw of 200. The orchestrator
  settled it by running both on one thread: `without_a_draw` gives 0.386 s
  and the four that genuinely need a draw give 0.456 s, against all five at
  0.497 s. The conclusion built on the row survives, because all five less
  the four that need no draw is still the spectrum.
- **The share of the spectrum was wrong in three ways at once.** Its low
  end is 0.097 s and not 0.098, by the pairwise subtraction the table
  supports. Its denominator was the pass over blocks already in memory,
  whose floor carries one copy of the 10 MB block for every block and reads
  no file; over the pass a user runs it is 16 to 18 in 100 and not 20 to
  23. And "8 to 10 in 100 at a draw of 20" rested on nothing: the
  orchestrator measured it in three alternating rounds, `without_a_draw`
  0.205, 0.206 and 0.207 s against all five 0.218, 0.218 and 0.218 s, so
  the spectrum adds 0.011 to 0.013 s to a pass of 0.218, which is 5 to 6 in
  100. The documents now separate the two quantities everywhere: the
  spectrum asked for alone costs more than the spectrum added to the other
  four, because a pass with the spectrum alone still counts the alleles of
  every population.
- **The description of the spectrum's arithmetic had its two branches
  swapped.** The product is of `num_called_alleles` factors where the
  copies of the major allele *can* fill the draw on their own, not where
  they cannot. Two reviewers found it separately; the function's own doc
  comment in the code was right all along.

**A self-check of the new benchmark could be switched off by a legal
refactor, and a reviewer switched it off.** `a_statistic_with_no_value`
zipped the public table `DiversityStats::NAMES_AND_STATS` against a
hand-written array in the same order and compared the two entry by entry.
Moving `("fis", FIS)` to the front of that table, which `NAMES` and
`of_name` both follow so nothing else breaks, made the guard never match,
so the check returned nothing for every statistic and stopped existing: a
pass that gave no value for F_IS went from exit 1 to exit 0. It now looks
each statistic up by its own entry, and the fixer broke it the same way
again to see it fail.

**Three of the four memory shapes are passes in which no variant reaches
the draw**, and both documents presented them as the memory of the pass.
The second has 10 individuals a population against the benchmark's
threshold of 20, so no variant counts at all; the third and the fourth ask
for draws of 2000 and 20000 that populations of 20 and 200 individuals
cannot fill. The figures stand, because the bins are allocated whatever the
data, and a reviewer produced the control: the second shape with 2000
individuals instead of 500, where the variants do count, gives 2.019 MB
against 2.007. Both documents now say so.

**The memory figure is not the same on every run above one thread.** The
counting allocator counts the whole process, and at 18 threads on a busy
machine a reviewer got 0.029 to 0.048 MB in 27 runs of the shape whose
figure is 0.053 MB, reaching it in none of them; with the machine quiet the
same command gave 0.053 MB twelve times out of twelve. At one thread it is
exact every time. That the extra is rayon's own splitting and stealing is a
suspicion and is written down as one.

**Three memory figures in a doc comment of the library contradicted the
spec this work package wrote.** `add_the_block` of
`crates/popnei/src/diversity.rs` carried 2.3, 15.9 and 40 MB where the
measurement gives 2.007, 15.479 and 36.343, and its 40 contradicted its own
next clause, "4 MB of bins for each chunk and 4 MB more for the pass
itself", which over 8 chunks is 36.

**Two gaps in the pyodide check**, which is the only check popnei has of
the wheel that runs in a browser tab. The plain call exists to prove that a
user who names no statistic gets the four that need no draw, and it never
asserted that: a wheel whose default differed died inside Python with
`'NoneType' object is not subscriptable`, preceded by node printing 1.2 MB
of minified pyodide, where the file promises to name each value that
differs. And the variants of the pass were never read, though the three
other calculation checks of that file all read them.

**Eleven smaller findings held**, each fixed: a populations file of three
columns was taken and the population named with a tab inside it, so a run
timed populations the file did not mean and exited 0; a `runs` of 0 died on
`min() iterable argument is empty` after a whole untimed pass; four
messages named neither the argument nor the value; the benchmark's
`default_vars` was a second copy of the public
`block::default_num_vars_per_block`; the spec pointed at load averages the
report does not hold; "the clock resolves 0.001 s" was the harness printing
three decimals and not the clock, whose resolution is 4.17e-08 s; the floor
of 0.126 s had no command anywhere, though every "over the floor" figure
rests on it; the three standardized values were given as 0.198 s in one
place and 0.199 in another where the pass with the three together is 0.208;
the growth of the memory with the threads was stated without its ceiling,
which the fourth shape reaches at 4 threads; the report's account of that
shape left 0.339 MB unexplained; and a doc comment said a pass whose blocks
carry no genotype where the run shows the genotypes are there.

**Nothing was set aside.** Every finding of this work package held and
every one is fixed. One was found by a fixer and not by a reviewer: the
0.339 MB is not a remainder but the counts each population keeps for the
chunk in hand, its private-allele pair, and the per-chunk state of the
pass, which together close all eight figures of the table to within 0.004
MB.

### What the owner should know

**The experiment of deliverable 2 was run and kept nothing, which is the
outcome it was for.** The spectrum's weights are already shared by the
bins. To measure what that sharing is worth, an unshared version was built
as a probe, one that computes each bin's weight from its own product of up
to `num_called_alleles` factors. It passed every correctness check before
any timing was read, so it was a slower version of the same calculation and
not a wrong one, and it took 11.042 s against 0.243 s on the pass with the
spectrum alone and 11.447 s against 0.500 s on the pass with all five, 45.4
and 22.9 times. The probe was reverted and is not committed.

The recurrence would have stayed whatever the timing said, for a reason
that is not about speed: the review of work package 3 found that computing
the bins without it underflowed to 0 and lost whole variants at a draw
above about 560 called alleles in a population. The measurement decides
what the sharing is worth; rightness decides that it stays.

**One thing is left for a performance review of its own, and it is worth
more than the spectrum.** Three of the statistics compute the same products
over and over on a variant of two alleles. Every standardized value rests
on one chance, that a draw of `g` of the `c` copies a population called
holds no copy of an allele it called `n` times; on a biallelic variant
`alleles_a_draw_shows`, `chance_a_draw_varies` and
`the_chance_each_draw_misses_the_allele` all ask for the same two products,
with the same arguments in the same order, so six products of `g` factors
are taken where two would do. The three are 0.208 s of a 0.505 s pass. An
upper bound on what sharing them could save is 0.14 s of that pass, 28 in
100, and the real figure is smaller because the sums and the loops are not
free. It was not implemented: work package 4's "What could go wrong" says
that anything beyond the shared weights of deliverable 2 is a review of its
own. `docs/reports/perf-diversity-2026-09-24.md` holds it under "Seen
outside the scope" with the measurement that would settle it, and with the
one thing to settle first, that the identity fails on a variant of three
alleles.

**No memory target was set**, and the report says why. There is nothing
measured to set one against: no other library's figure for this pass, no
run in a browser, and no dataset whose memory anybody has found too large.
A bound invented there would be a standard stricter than the real one,
which a reader would then hold the code to. The measurement is written down
instead: at 18 threads the pass holds 0.053 MB beside a block at the
panel's shape and 36.343 MB at 10000 individuals with 50 populations and a
draw of 20000, against a block of about 10 MB, and every figure of this
module is stated with the threads it was taken at, because the reduction
holds two chunks of partial sums for each thread of the pool.

**Thirteen Python harnesses under `crates/popnei/benches/` are covered by
no check of the project, and they stay that way for now.** `[tool.ruff]` of
`pyproject.toml` includes `python/**/*.py` and `tests/*.py` only, which a
reviewer proved by appending a duplicate import and a bare assignment to a
harness and watching both `uv run ruff check` and `uv run ruff format
--check` pass. Adding `crates/popnei/benches/*.py` to `include` was tried
and reverted: it wants reformatting in `make_big_vcf.py` in four places and
in `make_pops.py` in one, and `ruff check` raises two real complaints, both
in `make_big_vcf.py`, a `write` in a loop that should be `writelines` and
an `encode()` of a literal that should be a bytes literal. Fixing eleven
harnesses this plan did not write is not this plan's work, so the reason is
now in the comment beside `include` and the decision is the owner's. The
two harnesses this work package wrote pass both commands on their own.

**`crates/popnei/benches/stats_pass.rs` has two of the gaps that were fixed
here**, and was left alone: a numeric argument that is not a number is
refused without echoing the value, and a run in which no variant counts for
any population prints a timing and exits 0. Both are one line each. They
belong to `docs/specs/stats.md`'s module and not to this one.

**The four things that were already waiting are still waiting**, and this
work package added nothing to them. They are under "What the owner should
know" of work packages 1, 2 and 3 above: whether `scikit-allel` and `dadi`
become development dependencies of popnei or stay in the environments
`tests/reference/diversity/make_reference.py` builds;
`crates/popnei/src/io/vars.rs` giving the ploidy its metadata states with
no ceiling, which is what lets a ploidy popnei cannot read reach any pass;
`docs/reports/diversity-method/panel.py` hardcoding paths into a worktree
that no longer exists, so the script this spec names as the source of its
unbiased F_IS and its standardized private alleles cannot be rerun as it
stands; and `calc_per_var_distribs` and `calc_pop_dists` still taking their
choice of what to compute positionally where `calc_pop_diversity` now takes
it by name, a convention the owner approved for all three.

### How the work went

Not for the owner, who can stop here. It is for whoever next revises a
skill or writes a plan.

**A measurement is a claim like any other, and this work package is the
case for reviewing one.** Its three tasks took 631000 tokens of subagent
context to write, its six reviewers 716000 and its two fixers 426000, so
checking and repairing cost 1.14 million, 1.8 times what the work cost.
Work package 1, the only other one whose fixes were counted, cost 1.24
times; work package 2's fixes were not recorded, so its reviewers can be
compared and its total cannot, and they cost 0.93 times its tasks against
1.13 here. The 1.8 bought twenty-six findings, of which fourteen were
numbers or sentences that were not what had been measured, in documents
whose whole purpose is to carry numbers. A work package that produces prose
about measurements needs the same review as one that produces code, and the
categories that found most were `spec` and `numbers`, which are the two
that recompute.

**Three reviewers found the same mislabelled row from different sides**,
and the `code-review` skill is right that this is evidence and not
repetition. What settled it was neither of them: it was the orchestrator
running the two passes the label disagreed about and reading which range
the row fell in. A finding about a number is settled by taking the number.

**One reviewer overstated a finding and the orchestrator caught it by
running the case.** The `tests` reviewer reported that all four shapes of
the memory table count no variant; the first shape counts every one of
them, which one run prints. The reviewer's underlying finding held for the
other three, and the fix would have been wrong had the overstatement been
taken as given. Every finding whose fix is more than a word gets the
command run.

**The plan said tasks 4.1 and 4.2 could run side by side, and they could
not.** Nothing in the plan is wrong about the files they touch; what the
plan did not hold is that one of them takes wall times and the other runs
two long builds. A plan that marks two tasks as parallel should say
whether either of them measures time, because that is a resource of the
machine that no list of files shows.

**A fix that runs in the same tree as another fix needs its files named in
the brief, and that worked.** Two fixers ran at once over disjoint sets of
files, each told to report rather than edit anything belonging to the
other, and one of them handed a single line back through the orchestrator.
Nothing was lost and neither commit swallowed the other's paths.

**The benchmark's self-check was the finding worth the whole review.** It
was written to fail when a statistic has no value, it did fail, and a
reordering of a public table that breaks nothing else turned it off in
silence. The rule that a test must be broken to be believed is not enough
on its own: what found this was asking what else could change and leave the
test passing. A check that pairs two lists by position, where one of them
is public and may be reordered, is the shape to look for.
