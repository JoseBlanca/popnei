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

All six tasks are done and all eight deliverables pass. **The review is not
finished**: five of its seven reviewers were cut off by a rate limit on 24
September 2026 and have to be run again, so this work package is not yet
reported as done. What the two that finished found is below, and one of
their findings is a wrong number a user could see.

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

### What the review found so far

Seven reviewers were sent over `053562f`: all seven categories apply,
because this work package builds a module through the core crate, both
binding crates, the Python package and the TypeScript package. **Two
finished and five were cut off by a rate limit**, so this section will grow
when they are run again.

**A user who asks for a draw gets a wrong number and nothing says so.**
The `api` reviewer found it and the orchestrator reproduced it: with
`num_called_alleles` of 20 on the panel, `num_alleles.in_draw` is NaN for
all three populations, where "How it is verified" of "The number of
alleles" gives 1.9283948650, 1.9219209943 and 1.9197370844, and
`folded_sfs` is `None`. Both values are what the documented contract
reserves for something else: NaN for a draw larger than every population's
called alleles, and `None` for a statistic nobody asked for. So the answer
is wrong and reads as an answer.

It is the orchestrator's doing. Tasks 2.5 and 2.6 were told to pass
`num_called_alleles` through and let the core refuse what it refuses,
leaving the standardized fields as the spec describes them until work
package 3 fills them, and that instruction produced a silent wrong answer,
against the rule the owner gave on 21 September 2026. Both packages will
refuse `num_called_alleles` and the folded spectrum outright, saying that
the draw is not built yet, until task 3.5 fills them. Deliverable 5 of work
package 3 requires tests that a draw above every population's called
alleles gives NaN and that the spectrum works with a draw, so that task
cannot pass its own checks while the refusal is still there: the temporary
guard cannot be forgotten.

**The parallel reduction is never exercised by a cargo test.** The
`architecture` reviewer found that a chunk is 64 rows and every fixture in
the module has at most 6, so `par_chunks` yields one chunk in every cargo
test, the sum of the partials in index order is never run over more than
one, and no test varies the thread count, where eight other modules each
have a `the_number_of_threads_does_not_change_the_measures`. It checked the
property itself outside the tree, 5000 variants of 40 individuals in pools
of 1, 2, 4 and 8 threads, and got bit-identical F_IS, so the code is right
today; what is missing is the test that would catch the next change to it.

**What the pass keeps was measured and is right.** The trap the work
package warned of, an array of the populations by the alleles kept per
block instead of per row, is avoided: 4.56 MB of total allocation at a
block of 10000 rows with 50 populations is 157 chunks of 25.6 KB and not
10000 of them, which would be 256 MB, and one is live at a time. Nothing
accumulates from block to block.

**Two smaller findings hold.** The list of reasons F_IS is NaN, in all
three layers, misses the commonest one, a population whose counted
variants hold no whole called genotype. And `DiversityStats::of_name`
returns an `Option`, so the sentence a user reads for a statistic of no
such name is written out byte for byte in both binding crates, where
`stats::PerVarStat::of_name`, which it says it mirrors, returns a `Result`
and lets each binding write `?`.

**One warning for work package 3**, which the orchestrator had asked for.
When the spectrum arrives, `Totals` gains a per population vector of bins,
and the shape as built would then hold one per chunk at once: about 7 MB
and 7850 allocations for a block of 10000 rows with 50 populations at a
draw of 180. The reviewer suggests one flat buffer of populations by bins,
reused, and that goes into the brief of the task that builds the spectrum.
