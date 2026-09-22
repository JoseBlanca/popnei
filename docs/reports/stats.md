# Work report: the stats module and the filter of individuals

The plan `docs/plans/stats.md` is under way, on the branch `plan/stats`,
in the worktree `.claude/worktrees/stats`, since 22 September 2026. The
orchestrator, in this report, is the session of the assistant that runs
the plan: it sends each task to a subagent on Opus, checks what comes
back and has each work package reviewed.

## Before the first task

The owner approved the plan in chat on 22 September 2026. The branch
stands on `spec/stats` at 7bea42b, which holds the two specs, the plan
and `tests/reference/stats/`, and which is not in `main`; `spec/stats`
is from `main` at 7d8366f.

Run by the orchestrator in the worktree at 7bea42b: `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `306 passed`, 2 ignored; `cargo
wasm-check` finished; ruff `18 files already formatted` and `All checks
passed!`; `uv run maturin develop && uv run pytest` `174 passed`; `npm
install`, `npm run build` and `npm test` in `js/popnei` `tests 126`,
`fail 0`; after `npm install` in `tests/pyodide`, `bash
scripts/build_pyodide_wheel.sh` built the wheel and `node
tests/pyodide/smoke.mjs` exited with 0. `which plink2` gives
`/opt/homebrew/bin/plink2`, `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, and
`which bcftools` `/opt/homebrew/bin/bcftools`, version 1.24. pyNei's
`calc_per_var_distribs`, `calc_per_sample_stats`, `filter_samples` and
`load_vars` import under `uv run`, and
`/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars` is there.
`uv run python tests/reference/stats/make_reference.py` printed `done`
and exited with 0; after it `git status --short tests/reference/stats`
showed `panel.vcf.gz` modified, and only in the two bytes of the
modification time of its gzip header, bytes 5 and 6, the text of the
file being the same: the committed file was put back with `git
checkout`. The bcftools commands of "How it is verified" of the filter
of individuals on `many.vcf` give `ind05 ind00 ind49`, 423 variants with
the positions 1000, 1037, 1074, 1111 and 1148 first, and 26 with the
filter over the 50 individuals. `/Users/jose/devel/popnei-bench/big.vars`
is there, 81356714 bytes, and `big.vcf` beside it, 403572954 bytes.

The start message of the branch is on the board,
`.claude/board/2026-09-22T1020-plan-stats.md`.

## Work package 1: the steps of a pass as one enum

Task 1.1, commit c902c09: `PassStep` in `crates/popnei/src/filters.rs`,
`#[non_exhaustive]` as the spec has it, with `VarFilter` and
`KeepIndividuals`, and `chain_of` and `refuse_a_second_filter_of_a_kind`
over a list of it; the `Step` of each binding crate is a struct that
holds one `PassStep` and the names and values of its arguments as the
user sees them, and the two crates' own loops over criteria are gone.
The error enum of the core gets one case, `PassStepNotBuilt`, a
`RuntimeError` in Python, which both functions give for a
`KeepIndividuals` step until task 2.2 builds it; the subagent chose it
over `Ok(())` in `refuse_a_second_filter_of_a_kind` so that a second
filter of individuals could not pass in silence between the two tasks.
Task 2.2 removes the case, its arm in the Python crate and its test.

Both deliverables are met. `cargo test -p popnei --lib -- PassStep
--list` prints `3 tests`: the kind of each variant, the chain of the
three threshold filters over the worked example of the spec at 0.4, 0.88
and 0.25 keeping variant 5 with the counts 6 and 4, 4 and 3, 3 and 1,
and the error of a step not built. `cargo test --workspace` `309
passed`, 2 ignored, the 306 that were there and the three. `grep -rn
"VarFilteringCriterion" crates/popnei-python/src crates/popnei-js/src`
gives five lines per crate: the import, the three criteria built from
the user's argument, and the signature of the function that adds a
threshold filter, which takes the criterion just built. `uv run pytest`
`174 passed` and `npm test` `tests 126`, `fail 0`, both untouched.
`cargo fmt`, `cargo clippy`, `cargo wasm-check` and ruff clean.

Nothing was changed in the plan. The review of this work package is
made together with work package 2, as the `following-plans` skill allows
when two are one piece of code: it changes nothing a user sees, and task
2.2 rewrites the two functions it changed.

Two things for the next tasks. The filter of `cargo test` is a case
sensitive substring of the test path, so the tests the plan's checks
list by a type name, `PassStep`, `IndividualsReader`, sit in a module of
that name inside `mod tests`, with an `#[expect(non_snake_case)]`. And
since `PassStep` is `#[non_exhaustive]`, a `match` on it outside the
core needs a wildcard arm, so a step added to the core no longer stops a
binding crate from building: neither crate matches on it today, and the
tasks of the Python and TypeScript sides of a new step have to add it by
hand.

The subagent of task 1.1 used 179631 tokens in 55 tool calls.

## Work package 2: the filter of individuals

Tasks 2.1 and 2.2 went to one subagent, in that order, with a commit
each, since the reader of 2.2 is written over the method of 2.1.

Task 2.1, commit fe237ec: `Block::retain_individuals` in
`crates/popnei/src/block.rs`, the gather of each row on rayon off wasm
through a buffer per thread and the pack on one thread, as "How it
runs" of the filter has it. Before it, commit e8efa81 adds a paragraph
to `docs/specs/filters.md`, a choice the spec left to the code: which
error each of its three refusals is. They are three new cases of the
error of the crate, defects of popnei and a `RuntimeError` in Python,
as the `keep` of `retain_vars` with the wrong number of values is,
because `resolve_individuals` refuses the name behind each of them
before a call of a user reaches them.

Task 2.2, commit ef17cc5: `resolve_individuals`, `IndividualsReader`,
the `KeepIndividuals` arm of `chain_of` and of
`refuse_a_second_filter_of_a_kind`, and the four cases of the spec in
the error, each a `ValueError` in Python; the temporary case of work
package 1 is gone with its arm and its test. The JavaScript crate
needed no arm, since it turns every error of the core into one `Error`.

Task 2.3, commit 91b74fc: `filter_individuals` in
`crates/popnei-python/src/steps.rs`, which resolves the names against
the individuals of the source at the call and refuses a second filter
of the kind there; `Variants.filter_individuals` with its docstring,
`individuals` and `num_individuals` read through the steps, and the
`repr`; and `tests/test_filter_individuals.py`, 13 tests. Three things
the spec does not say, decided by the subagent: the `Steps` of the
binding crate takes the individuals of the source when it is built, so
the refusal cannot be given another list; `filter_individuals("ind05")`,
a bare string, is a `TypeError` that says to write `("ind05",)`, as
`iter_blocks(fields="pos")` already is, where otherwise the user would
read that `i` is not an individual; and the `repr` prints every kept
name, `individuals(individuals=('ind05', 'ind00', 'ind49'))`, so a
filter of 500 individuals prints 500 names. The last two are the
owner's to reverse.

Task 2.4, commit 4c5e5ae: the same in `crates/popnei-js/src/steps.rs`
and `js/popnei/src/variant.ts`, `filterIndividuals`, with
`js/popnei/test/filter_individuals.test.ts`, 13 tests. The arguments of
a step used to cross to JavaScript as one number each; now a threshold
crosses in one flat array and the names of the kept individuals in
another, with a count of names per argument that says which is which,
0 names meaning a threshold, which is unambiguous because a filter of
no individual is refused. The `Default` of the `Steps` of both crates
is gone, since `new` takes the individuals.

The five deliverables are met, run by the orchestrator at 4c5e5ae.
`cargo test -p popnei --lib -- retain_individuals --list` prints `8
tests`, where the plan asks 4 or more; `cargo test -p popnei --lib --
resolve_individuals IndividualsReader --list` `15 tests`, where it asks
8 or more; `uv run pytest tests/test_filter_individuals.py` `13 passed`
and `uv run pytest` `187 passed`, from 174; `npm test` `tests 139`, `fail
0`, from 126; `bash scripts/build_pyodide_wheel.sh` built the wheel and
`node tests/pyodide/smoke.mjs` exited with 0. `cargo test --workspace`
`331 passed`, 2 ignored, from 309. `cargo fmt`, `cargo clippy`, `cargo
wasm-check` and ruff clean.

Nothing was changed in the plan.

The subagents used: tasks 2.1 and 2.2 together 203867 tokens in 84 tool
calls; task 2.3 147130 in 63; task 2.4 182877 in 60. None had to be
sent back.

### The review of work packages 1 and 2

Seven reviewers, one per category, at 4c5e5ae. None found a wrong
number: the spec reviewer reproduced every value of the item against
bcftools 1.24 and pyNei, and the gather of 18 random subsets of
`many.vcf`, of haploid and triploid files and of the 50 individuals
reversed against direct indexing; the numbers reviewer compared
`retain_individuals` with a naive gather over 625 cases of sizes,
ploidies and kept sets, all equal. What they found, fixed in nine
commits from 15281d2 to 2d360b6 by one subagent, was in what happens
when a bound is broken, in the tests and in the texts:

- The compaction of a block skipped, in silence, a genotype it could not
  reach and a row it could not pack, when the arrays of the block were
  shorter than it states: unreachable today, since `check` runs first,
  but a wrong genotype with no error the day a caller reaches it. Three
  reviewers found it from three sides. Now the two helpers give the
  error of `check`, and a test pops one allele and asserts the error
  and the block as it was. The reader of the filter dropped in silence
  an index beyond the individuals of its source when it built its
  names, unreachable for the same reason; it gives the error now.
- The buffer of the gather was allocated once per rayon job, about one
  per three rows, and not once per thread as its comment said: 1861 to
  2245 jobs for a block of 5000 variants of 1000 individuals, 500 kept,
  on 18 threads, measured by the fixer on the M5 Pro; with a floor of 64
  rows per job, 63. What the gather costs with either was not measured.
- Two tests could not fail: the one thread against several test had 300
  identical rows, since its 20 alleles per row cycled 5 values, and a
  mutation that copies row 0 into every row passed it; and no test
  covered the reader's own error, so `finished` could be dropped from
  that branch with every test green. Both have their test now.
- The kind of an argument of a step crossed to TypeScript as a count of
  names, 0 meaning a threshold, so a later argument of neither kind
  would have taken the next threshold in silence; two reviewers found
  it. A kind per argument crosses now, and an unknown one is an
  `Error`.
- The rule that gives the individuals of the next pass was written twice,
  by two mechanisms, in Python and in TypeScript. It is one function of
  the core now, `individuals_of`, added to "The Rust interface" of the
  spec in 2d360b6, and each binding crate answers with it.
- In Python, a non-string name and a non-iterable argument gave the
  `TypeError` of Python or of pyo3, naming no argument; both name it now,
  as the TypeScript does. The message of an unknown name points at where
  the names are. A docstring that counted three filter methods, and a
  cut sentence in the comment of the exception table, are mended.

Not taken, with the reason: that the `reblock` at the end of
`iter_blocks` sizes the blocks after the filter for the kept individuals
and not the source's, above 500 individuals, because `docs/specs/block.md`
says its default is for the individuals of the reader it is given, and
the sentence of the filter spec is about the filter's own reader; the
spec of the filter says so now, in 15281d2, which also says that the 500,
423 and 26 of `many.vcf` come out with `only_passed` false. That a chain
built over a reader that already holds a filter of individuals is not
refused, because the spec limits that refusal to the threshold filters,
whose `FilteredReader::new` asks the reader. That `filter_individuals`
accepts a generator under a `Sequence[str]` hint, which harms nobody.

Seen outside the scope, for the owner: `docs/glossary.md` line 154 ends
in the middle of a sentence, from before this branch and in a shared
file, so it is not touched here; the Python package ships no `py.typed`,
so its type hints never reach a user's type checker; and the `Steps` of
the JavaScript crate clones the names of the source at every pass.

After the fixes, at 2d360b6: `cargo test --workspace` `334 passed`, 2
ignored; `retain_individuals --list` `9 tests`; `resolve_individuals
IndividualsReader --list` `16 tests`; `uv run pytest` `189 passed`, of
`tests/test_filter_individuals.py` 15; `npm test` `tests 140`, `fail 0`;
the wheel of pyodide built and its smoke test exited with 0; `cargo
fmt`, `cargo clippy`, `cargo wasm-check` and ruff clean.

The reviewers used, in tokens: spec 148820, tests 151940, numbers
134789, errors 127376, api 113582, architecture 115202, binding
126337. The fixer used 231504 in 129 tool calls.

## Work package 3: the populations and the counts over one

Task 3.1, commit 51e98e9: `Pops` in the new module
`crates/popnei/src/stats.rs`, with `all`, `from_names`, `len`,
`is_empty`, `name`, `individuals`, `is_all` and the named constant of
the population of every individual, `pop`; four cases of the error for
its refusals, each a `ValueError` in Python, and one for an index
beyond the variant, a `RuntimeError`; and `count_gts_of` and
`count_alleles_of` in `crates/popnei/src/variant.rs`, beside the two
counts of the whole row, with what one genotype and one allele add to
the counts in one private function each that both pairs share. Four
things the subagent decided: the map from a name to its index is built
once per population and not once per pass, since "The Rust interface"
has `from_names` call `resolve_individuals` for each, which builds its
own, so 3 populations of 200 individuals build 3 maps; `name`,
`individuals` and `is_all` keep the signatures of the spec, a bare
`usize` and no `Option`, and a number at or beyond `len()` gives an
empty name, no individual and false, written in each doc comment; and
both counts refuse a population whose individuals times the ploidy is
above what a `u32` holds, with the case `count_gts` already has, which
only a caller that names an individual twice can reach.

Both deliverables are met, run by the orchestrator at 51e98e9. `cargo
test -p popnei --lib -- stats::pops --list` prints `9 tests`, where the
plan asks 5 or more; `cargo test -p popnei --lib -- count_gts_of
count_alleles_of --list` `12 tests`, 9 new and 3 that were there whose
names begin with `count_gts_of_the_`, where it asks 4 or more. `cargo
test --workspace` `352 passed`, 2 ignored, from 334; `uv run pytest`
`189 passed`, untouched; `cargo fmt`, `cargo clippy`, `cargo
wasm-check` and ruff clean.

Nothing was changed in the plan. The review of this work package is
made together with work package 4, which calls its two counts from the
statistics and compares them with pyNei through them.

The subagent of task 3.1 used 156366 tokens in 63 tool calls.

## Work package 4: the per variant distributions, through the three layers

Task 4.1, commit 6a0dbdc: `HistBins` with the bins of equal width and of
equal ratio, and `ObsHet`, `Maf` and `ExpHet` with the value of one
variant in one population, in `crates/popnei/src/stats.rs`; four cases
of the error, each a `ValueError` in Python. Two choices of arithmetic
the spec leaves open, both made so that the value is the one numpy
computes or does not become NaN: a frequency is raised to the exponent
by that many multiplications and not with `powf`, so at ploidy 2 it is
numpy's `p * p` on every platform; and the unbiased value multiplies
the k factors of each term one over another, because the product
`c (c - 1) ... (c - k + 1)` for a million called alleles and an
exponent of 255 is above the largest float64 and the value would be NaN
with nothing to say so. The writer checked two readings against numpy
before choosing: `linspace` writes the end of the range into the last
edge instead of `start + num_bins * step`, and `logspace` raises 10 to
the logarithm of the end; popnei does the same.

Task 4.2, commits 1232779 and 1fbbdd1: `calc_per_var_distribs` with
`PerVarStat`, `PerVarDistribsConfig`, `StatsDistrib`, `PolyVarsStats`
and `PerVarDistribs`, rayon over chunks of 64 rows merged in the order
of the block, with the serial version beside it for wasm, and two cases
of the error, a pass that gave no variant and a `poly_threshold` out of
range. Before the code, 1232779 put into the spec the numbers the
worked example lacked: the plain expected heterozygosity and the
polymorphism counts of the case with no populations, and the counts
under a `poly_threshold` of 0.5, all from pyNei at ef0ca6e, so that a
test can tell whether the threshold is read at all. Every number the
spec already had came out of the same run unchanged.

Task 4.3, commits 923b36b, 0fa1c74 and a108621: the Python side. The
function of the binding crate builds the chain with `chain_of`, the
populations with `Pops::from_names` against the individuals of that
chain, runs the pass with the interpreter released and reads the counts
of the filters from the chain afterwards; `python/popnei/stats.py` holds
`calc_per_var_distribs` and the four result types with their pandas
series and frames; `tests/test_stats.py` has 34 tests. Two spec
commits came first. 923b36b makes a key of `hist_kwargs` that popnei
does not know a `ValueError` that names the three it takes, where pyNei
ignores it: a user who writes `nbins` gets the 40 bins of the default
with nothing said, which is not the result they asked for. This
follows the owner's rule that an error never passes silently, and it is
his to reverse. 0fa1c74 writes down what the comparison with pyNei can
and cannot assert, which is the finding of this work package worth the
owner's attention.

Task 4.4, commit 631e281: the TypeScript side, `calcPerVarDistribs` as
a method of each source class, as `write_vars` and `blocks` are, with
the populations crossing flat and the result as a `Float64Array` of
means with NaN and a `Uint32Array` of histogram counts; 14 tests in
`js/popnei/test/stats.test.ts`.

The four deliverables are met, run by the orchestrator at 631e281.
`cargo test -p popnei --lib -- stats::hist stats::obs_het stats::maf
stats::exp_het --list` prints `31 tests`, where the plan asks 16 or
more; `cargo test -p popnei --lib -- stats::distribs --list` `13
tests`, where it asks 8 or more; `uv run pytest tests/test_stats.py -k
per_var` `34 passed`; `npm test` `tests 154`, `fail 0`, from 140.
`cargo test --workspace` `396 passed`, 2 ignored, from 352; `uv run
pytest` `223 passed`, from 189. `cargo fmt`, `cargo clippy`, `cargo
wasm-check` and ruff clean.

What the owner should know from this work package.

The unbiased expected heterozygosity does not agree with pyNei to the
last bit, and cannot: popnei multiplies the factors of each term one
over another and pyNei multiplies the plain value by `c / (c - 1)`. The
two agree to the last bit or the one before it, so the means agree
within 1e-12, but a variant whose value falls on an edge of the
histogram is counted on either side of it. Two of them do, both
measured on 22 September 2026 with the default histogram of 40 bins
from 0 to 1: the variant at position 9695 of `many.vcf` has the allele
counts 55 and 45 of 100, whose value is exactly 0.5, and popnei gives
0.5 where pyNei gives 0.4999999999999999; the variant `var0978` of the
panel in the population p2 has 106 and 54 of 160, value exactly 0.45,
and popnei gives 0.45000000000000007 where pyNei gives
0.44999999999999996. So the test compares the counts of that one
statistic allowing a variant within 1e-9 of an edge to fall on either
side, and the counts of the other four exactly. It is in the spec, in
commit 0fa1c74.

The Python package now depends on pandas 3.0.2, pyNei's version,
because the spec says the results are its series and frames; `uv.lock`
moved with it. In TypeScript the populations of a result come in the
iteration order of the keys of the object the user gave, and JavaScript
puts keys that are whole numbers first, in numeric order, so
populations named "1" and "2" would come before one named "north"
whatever the user wrote. It is documented on the argument. A `Map`
would keep the insertion order of every name, and the spec does not
say which of the two it wants.

The subagents used: task 4.1 213303 tokens in 68 tool calls; task 4.2
221213 in 98; task 4.3 282844 in 135; task 4.4 302305 in 81. None had
to be sent back.

### The review of work packages 3 and 4

Seven reviewers, one per category, at 631e281. This is the review that
found the most, and two of its findings were wrong results a user would
have got with nothing said. Twenty-four findings were fixed in
twenty-five commits, from f23aeef to 27ba61a, by two subagents, the
first on the core and the spec and the second on the binding crates and
the packages.

The two wrong results, both in the histogram, both found by three
reviewers from three sides and confirmed by the orchestrator:

- A range whose two ends are so far apart that their distance overflows,
  `(-1e308, 1e308)`, gave the edges NaN, inf, inf, inf, 1e308. The bin of
  a value is found by a binary search, which needs the edges sorted, so
  every variant fell in the first bin: all 500 of `many.vcf`, with a mean
  computed over all of them and a NaN printed among the edges. numpy
  refuses the same input, and so does pyNei through it. popnei now
  refuses it, and the histogram has a largest number of bins as well,
  100000.
- That largest number is the second one. Nothing bounded the number of
  bins above, and the allocation of the edges is what the user's number
  reaches: `num_bins` of 2^60 ended a Python session with an exception
  that derives from `BaseException`, which `except Exception` does not
  catch, so a notebook dies; 10^12 aborted the interpreter with
  `memory allocation of 8000000000000008 bytes failed` and no exception
  at all; in the browser it is the trap that ends the module. A user who
  writes a number that large has made a mistake, and popnei now says so.

The third was in Python alone: every count a user got was an unsigned
64-bit integer, where pyNei's are signed, so subtracting two of them
wrapped. On `many.vcf`, `num_poly - num_variable` gave
18446744073709551600 where pyNei gives -16. The counts are signed now,
and a test subtracts two of them.

What the reviewers found in the tests is the other half of the review.
Six behaviours the code has, and the spec states, were guarded by no
test at all: a value outside the range of the histogram is in the mean
and in no bin; the last edge is the end of the range, which is what
numpy writes there; the pass refuses a `Pops` built against another
reader; it asks its reader for the genotypes alone; the answers of
`Pops` for a population that is not there; and the serial path that wasm
uses, which no test compared with the rayon one. In each case a reviewer
broke the code and all 396 tests stayed green. They are tests now, each
confirmed by making the same change and seeing the new test fail.

Six sentences of the spec said something the code does not do, and each
was measured before it was corrected. The comparison with pyNei of the
unbiased expected heterozygosity crosses a bin edge at seven variant and
population pairs and not at the two the spec named, and the two
arithmetics differ by at most 2.84e-16, not by a bit or two: two
reviewers recomputed that independently, over every biallelic split up
to 500 called alleles. The sums of a pass are bit-identical across
thread counts, where the spec said they agree to about 1e-15 and not to
the bit, so the test now compares the bits. The lookup of the
individuals of a population is one hash map per population, where one
sentence said once per pass and the interface prescribed per population.
The test of the block sizes compares them against each other now, and
not each against a printed literal. At ploidy 1 the plain expected
heterozygosity is not 0 but -2.2e-16 at a variant with a called allele,
so it falls below the histogram and into the mean; numpy gives the same
value, so popnei was left as it is and the spec now says what floating
point does there. That one is the owner's to reverse: rounding it to 0
would count such a variant in the first bin where pyNei drops it.

The rest were messages a user reads and texts a maintainer reads. A
`TypeError` said that what it refused "is one of them", the opposite of
what it meant, in three places. Four refusals in Python named neither
the argument nor the value, where TypeScript named both. A user who
wrote `ploidy=0` was told about "the exponent", a word they never
wrote. A truth value written for `min_num_individuals` was taken as 1,
where TypeScript refuses it. The four distributions of one result shared
one array of bin edges that a user could write into although the result
is frozen; it is read-only now. The five names of the statistics were
written out in four places and the bin types in two, with the same
sentence copied into both binding crates; the core turns a name into a
statistic and into a histogram now, and both crates call it.

Not taken, with the reason: that the per-chunk allocation of the pass
should be one flat array instead of one per population per statistic,
which is true of what it allocates, 136 allocations per population, but
the change touches five types, the loop over the rows and the result,
which is more than a fix; the comment that said the pass allocates
nothing per variant now says what it does allocate. That the pass
should honour Ctrl-C while it runs, which is the choice the writer of
the vars file made before it and which the report of the filters already
put before the owner. That `IndividualBeyondTheVariant` reaches Python
as a `ValueError`: it is in the arm that gives a `RuntimeError`, as it
should be, and the reviewer had misread the arm.

One finding is for work package 6 and not for a fix. The architecture
reviewer measured the pass on `big.vars`, 100000 variants of 1000
individuals, five statistics and no populations: 0.511 s on one thread
and 0.149 s on 18 cores, where "Speed" of the spec asks 0.25 s and
0.15 s, with the read alone at 0.109 s against the spec's 0.102 s. Four
populations of 250 add 0.13 s on one thread. It was measured on a
machine with a load average of 5.45, which is not the quiet machine work
package 6 measures on, so the number there is the one that counts.

After the fixes, at 27ba61a: `cargo test --workspace` `412 passed`, 2
ignored, from 396; `stats::pops --list` `10 tests`, `count_gts_of
count_alleles_of --list` `12 tests`, `stats::hist stats::obs_het
stats::maf stats::exp_het --list` `38 tests`, `stats::distribs --list`
`19 tests`, where the plan asks 5, 4, 16 and 8; `uv run pytest` `233
passed`, from 223, of `tests/test_stats.py -k per_var` 44; `npm test`
`tests 156`, `fail 0`; the wheel of pyodide built and its smoke test
exited with 0; `cargo fmt`, `cargo clippy`, `cargo wasm-check` and ruff
clean.

The reviewers used, in tokens: spec 208043, tests 242450, numbers
200040, errors 158749, api 180412, architecture 154966, binding 147875.
The two fixers used 291812 in 214 tool calls and 282016 in 186.

## Work package 5: the per individual statistics, through the three layers

Task 5.1, commit 37f5be2: `calc_per_individual_stats` and
`PerIndividualStats` in `crates/popnei/src/stats.rs`. The pass keeps two
counts per individual, the missing genotypes and the heterozygous ones,
reads the rows in chunks of 64 on rayon and adds the chunks in the order
of the block, with a serial version beside it for wasm, and makes the
two divisions once at the end. A chunk allocates one array of two
counts per individual and nothing per variant. Two things the subagent
did beyond the task, both inside the tests: the test reader and the
three helpers that open a reference VCF moved into a module of fixtures
that both test modules read, and what a missing and a heterozygous
genotype are is now one function of the `variant` module that both
passes call, rather than written twice.

Task 5.2, commit 8f40a83: the Python side, which mirrors the per variant
pass: it builds the chain, takes the names of the individuals from it,
runs the pass with the interpreter released, turns the individual with
no called genotype into NaN and reads the counts of the filters from the
chain afterwards. `python/popnei/stats.py` has the frozen dataclass and
the function, exported from `popnei`, and `tests/test_stats.py` eight
more tests. The tests of the five individuals of pyNei's own test write
their own VCF, because the fixture of `tests/conftest.py` has a fixed
header of three.

Task 5.3, commit c93a6e2: the TypeScript side, `calcPerIndividualStats`
as a method of each source class, with six tests. Its tests read
`many.vcf` with `onlyPassed` false, because the literals of the spec are
over its 500 variants and the default of the reader keeps 475, which
the reference script confirms, since its plink2 command has no filter of
the variants that passed.

The three deliverables are met, run by the orchestrator at c93a6e2.
`cargo test -p popnei --lib -- stats::per_individual --list` prints `8
tests`, where the plan asks 5 or more; `uv run pytest tests/test_stats.py
-k per_individual` `8 passed`, with the 44 of the per variant pass
untouched; `npm test` `tests 162`, `fail 0`, from 156. `cargo test
--workspace` `420 passed`, 2 ignored, from 412; `uv run pytest` `241
passed`, from 233. `cargo fmt`, `cargo clippy`, `cargo wasm-check` and
ruff clean.

Nothing was changed in the plan and nothing in either spec: the
subagent of task 5.1 recomputed the literals of the panel and of
`many.vcf` from the files and they are the spec's.

The subagents used: task 5.1 187351 tokens in 97 tool calls; task 5.2
158161 in 76; task 5.3 207184 in 82. None had to be sent back.
