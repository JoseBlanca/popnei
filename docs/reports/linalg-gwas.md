# Work report: the linear algebra the association study needs

23 September 2026. This is the work report of the plan
`docs/plans/linalg-gwas.md`, carried out on the branch
`plan/linalg-gwas` in the worktree `.claude/worktrees/linalg-gwas`,
which built the linear algebra that a genome wide association study
needs into `crates/popnei-linalg` from the spec `docs/specs/linalg.md`.
It was written while the work went, so the sections below are in the
order the work happened, each work package with the commands that
checked it, what its review found and what was done about that. What the
owner reads first is the section below; the decisions that are theirs
are under "What is waiting on the owner"; what the work taught about
running the next plan is the last section.

## What the owner reads first

**The plan is done.** All four work packages are built, reviewed and
fixed, every deliverable checks out, and the plan's final check was run
again on the last commit. The branch is `plan/linalg-gwas`, in the
worktree `.claude/worktrees/linalg-gwas`. Nothing is merged into `main`
and nothing is pushed.

**What exists that did not this morning.** `crates/popnei-linalg` gives
the seven operations a genome wide association study needs, each on both
of its backends: the Cholesky factorization of a symmetric matrix that no
vector makes negative, and the solve, the log of the determinant and the
inverse that come off it; the thin QR of a design, which is the
factorization that fits a linear model to more individuals than
coefficients; the solve against a triangular matrix, either half; and the
rank, how many of a design's columns are independent. Its product of two
matrices now computes all four ways the two can be laid out, and its error
type has one case more, for a matrix that could not be factored. The crate
went from 42 tests to 149, and from 37 to 136 in the build that uses faer,
the linear algebra library that runs in the browser.

**How you know it works.** Every number a test asserts is a literal taken
from numpy 2.5.3, the same library pyNei computes with, and every test
runs twice, once against each backend. Where a reviewer doubted that a
test could fail, it broke the code on purpose and showed which test
caught it. Nothing of the crate is reached from Python or TypeScript yet,
so the check that popnei and pyNei agree end to end is not this plan's to
make; the association study, whose spec another session is writing now,
will be the caller.

**What is asked of you.** The merge, which is yours to order. The branch
merges into `main` at `bed9031`, and the merge was made in a throwaway
worktree, built, tested and thrown away, so it is known to build and pass
and not only to be free of conflicts: "The merge into `main`, tried and
not made" has the numbers. The session writing `docs/specs/gwas.md` and
`docs/specs/kinship.md` asked to be told when this branch is ready, so
that the three can be merged together if you want that.

And four decisions, none of which stops anything and each of which
"What is waiting on the owner" gives with its options and its cost.

## Before the first task

The branch starts from `spec/linalg-gwas` at `7ccc709`, which holds the
spec and merges into `main` with no conflict.

Everything the plan asks to be in place is there, checked by running it
in the worktree on 23 September 2026:

| What the plan asks | Command | What it gave |
| --- | --- | --- |
| the workspace is clean | `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings` | both clean |
| the workspace passes | `cargo test --workspace` | `472 passed` in the core crate, `42 passed` in the linalg crate, 2 ignored |
| both backends pass | `cargo test -p popnei-linalg` and the same `--no-default-features` | `42 passed` and `37 passed` |
| both wasm targets build | `cargo wasm-check`, an alias of `.cargo/config.toml` that compiles the core crate and the linear algebra crate for both of popnei's WebAssembly targets with the warnings denied | clean |
| the Python layer passes | `uv run ruff format --check && uv run ruff check`, `uv run maturin develop && uv run pytest` | clean, `257 passed` |
| the TypeScript layer passes | `npm run build && npm test` in `js/popnei` | `pass 180` |
| none of the seven exists yet | `cargo test -p popnei-linalg --lib <filter> -- --list` for `cholesky`, `qr`, `rank`, `triangular`, `determinant`, `singular`, `four_ways`, `first_operand`, `solve_with`, `invert` | `0 tests` for every one of the ten |
| `product` has no typed first operand | `grep -rl TheFirstOperand crates/` | no file |
| the Cargo.toml comment still names an open point | `grep -c "Open 1" crates/popnei-linalg/Cargo.toml` | `1` |
| the reference the tasks build from works | `cargo run --release` and the same `--no-default-features` in `tmp/linalg_gwas_trial/` | every check of both backends printed `ok` |

One trap for anybody who runs the checks here. `uv run maturin develop
--release` makes the pytest test
`test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file` fail:
it sends the process a ctrl-C 0.1 s into a write that a release build has
already finished. The check of the `coding` skill is the debug build, and
with it all 257 pass.

## Work package 1: the product with its first operand turned

### Task 1.1, the typed first operand and its four callers

`8d7fb27`. `product` takes a `TheFirstOperand` beside its
`TheSecondOperand`, the four calls of it in the core crate name
`ByTheRowsOfTheResult`, and the comment of
`crates/popnei-linalg/Cargo.toml` says the owner's decision about the
rustc flag instead of naming an open point.

Checked by the orchestrator, each command run in the worktree: `cargo fmt
--all --check` clean, `cargo clippy --workspace --all-targets -- -D
warnings` clean, `cargo test --workspace` `472 passed` and `42 passed`
with 2 ignored, `cargo test -p popnei-linalg --no-default-features` `37
passed`, `cargo wasm-check` clean, ruff clean, `uv run maturin develop &&
uv run pytest` `257 passed`. Every one is what it was before the task, so
nothing the principal component analysis or the r² computes moved.
`grep -c ByTheRowsOfTheResult` gives 3 in `pca.rs` and 1 in `ld.rs`, and
`grep -c "Open 1" crates/popnei-linalg/Cargo.toml` gives 0, which are
deliverables 1 and 4. The diff of `pca.rs` and of `ld.rs` is the four
calls and one `use` line and nothing else; no test of either was touched.

One thing the task had to settle. The two combinations whose first
operand is `ByTheValuesSummedOver` have no backend behind them until task
1.2, and the enum lets a caller ask for them, so the task made those two
arms give `Error::Dimension` with a message that says the product is not
built yet, rather than an `unreachable` that would be a wrong matrix or a
panic waiting to happen. Task 1.2 replaces the two arm bodies with one
backend call each. Nothing outside the crate could reach it in the
meantime: no caller in popnei names that case.

### Task 1.2, the two combinations that turn the first operand

`d8655ba`. `c = a' b` and `c = a' b'` on both backends: one `dgemm` each
in `crates/popnei-linalg/src/blas.rs`, with a `trans` flag set and no
copy, and one `matmul` over a transposed matrix reference each in
`crates/popnei-linalg/src/faer.rs`, which is another reference over the
same values. Three tests,
`the_same_matrix_comes_out_of_the_product_four_ways` and two named for
the first operand, and the crate's doc comment now says four products and
not two.

Checked by the orchestrator: fmt, clippy with the warnings denied, `cargo
wasm-check` and ruff clean; `cargo test --workspace` `472 passed` with 2
ignored in the core crate and `45 passed` in the linear algebra crate;
`cargo test -p popnei-linalg --no-default-features` `40 passed`; `uv run
maturin develop && uv run pytest` `257 passed`. The linear algebra crate
went from 42 to 45 tests and from 37 to 40 on faer, which is the three
the task added, and every other number is what it was.

The four pairs of operands of the test are the spec's, read from "How it
is verified" of "The product with its first operand turned" and not
worked out again: A of 2 x 3 with rows (1, 2, 0) and (0, 1, 3), the same
matrix the other way round, the B with rows (1, 1), (2, 0) and (0, 3),
and that one the other way round. All four combinations write the 2 x 2
with rows (5, 1) and (2, 9), asserted exactly.

### The deliverables of work package 1

| Deliverable | Command | What it gave |
| --- | --- | --- |
| 1, the typed first operand and its four callers | `grep -c ByTheRowsOfTheResult crates/popnei/src/pca.rs crates/popnei/src/ld.rs`; `cargo test --workspace` | `3` and `1`; `472 passed` in the core crate, the same tests as before and none of them changed |
| 2, the four combinations give one matrix | `cargo test -p popnei-linalg --lib four_ways -- --list`; `cargo test -p popnei-linalg` and the same `--no-default-features` | `the_same_matrix_comes_out_of_the_product_four_ways`, `1 test`, where the filter gave `0 tests` before; `45 passed` and `40 passed` |
| 3, what the new combinations refuse | `cargo test -p popnei-linalg --lib first_operand -- --list` | the two tests of the dimensions and of the value that is not finite, `2 tests`, where the filter gave `0 tests` before |
| 4, the comment about the vector instructions | `grep -c "Open 1" crates/popnei-linalg/Cargo.toml` | `0`, where it was `1` before |

### The review of work package 1

Six reviewers, one per category, over `8d7fb27` and `d8655ba`: spec,
tests, numbers, errors, api and architecture. `binding` was not sent,
because no function of this crate is reached from Python or from
TypeScript and neither binding crate names it. They used 746000 tokens
between them.

**What they found that mattered, and what was done.**

Five of the six found the same gap, and three showed it by mutation: every
case of the test of the four combinations has a 2 x 2 result, so the rows
and the columns of the result are the same number, and a backend that
exchanged the two passed all 45 tests of the crate and all 40 of its faer
build. The spec's older section asks for "A test of each product on
matrices that are not square" and its new section had dropped it, so the
spec was corrected first, in `90e9f08`, with the two matrices that section
already gives for the same A: `a' b` is the 2 x 1 with rows (1) and (6)
and `a' b'` the 2 x 1 with rows (2) and (3), checked against numpy 2.5.3
before they were written. Then the test, in `839e11d`. The orchestrator
made the mutation itself on both backends afterwards:
`the_product_that_turns_its_first_operand_writes_a_result_that_is_not_square`
is the only test of the crate that fails, and it fails on both.

The `tests` reviewer found that
`the_product_that_turns_its_first_operand_refuses_a_value_that_is_not_finite`
passes word for word against the code of `8d7fb27`, where neither new
combination existed, because `lib.rs` refuses a value that is not finite
before the dispatch. It now also asserts the matrix of a call with nothing
wrong, which is what its sister test already did.

The `spec` and `architecture` reviewers found that the spec, both
`Cargo.toml` files and `docs/reports/pca-measurement.md` said the rustc
flag `-C target-feature=+simd128` is set nowhere and changes no byte, and
that neither is true. That is `68b29fd`, and what is left of it for the
owner is below.

The `api` reviewer found four backend functions of identical signature
told apart by a preposition, `product_by_transpose` being `a b'` and
`product_of_the_transpose` being `a' b`. They are now `product`,
`product_with_the_second_turned`, `product_with_the_first_turned` and
`product_with_both_turned`. The `errors` reviewer found that nothing tells
a caller that `rows`, `inner` and `cols` are the least a buffer may hold,
so that a wrong one gives another matrix and no error, which is what
"Errors" of the spec decides; the `# Errors` of `product` now says it.
Three smaller gaps were closed with them: the module comment of `blas.rs`
called `a b'` the call whose `transa` is `T`, which `a' b'` is too; the
doc of `Error::Dimension` did not name an `inner` of 0; and the test of
the cap of 2147483647 and the test of an `a` of no rows reached only the
unturned layout.

**What was not taken, and why.** An `inner` above 2147483647 with `rows`
at 0 gives an error that names `b`: the sentence is true of `b`, whose
dimensions really would be out of range, and it prints the number at
fault. The timings of "What it costs" were measured at the routine and not
through the crate's checks, which "Speed" of the spec already says of
every number it gives. `eigh_lower` in the faer backend indexes faer's
matrices, which is older than this work package and whose indices are
below `n`.

**The checks after the fixes**, each run by the orchestrator: fmt, clippy
with the warnings denied, `cargo wasm-check` and ruff clean; `cargo test
--workspace` `472 passed` with 2 ignored in the core crate and `47
passed` in the linear algebra crate; `cargo test -p popnei-linalg
--no-default-features` `42 passed`; `uv run maturin develop && uv run
pytest` `257 passed`. The crate went from 42 tests to 47 and from 37 to
42 over the whole work package, and nothing else moved.

## Work package 2: the Cholesky factorization and the three things off it

### Task 2.1, the factorization and the `Singular` case

`eb9c162`. `cholesky_lower` above both backends, `dpotrf` with `uplo` `U`
in the BLAS one, since popnei's lower half is the routine's upper half,
and `cholesky_in_place` of faer's `llt::factor` with faer's default
regularization, which is the one that refuses a pivot that is not
positive; and the `Singular` case of the error enum with the fields, the
doc comment and the message of the spec. Eight tests, all on both
backends.

Checked by the orchestrator: fmt, clippy with the warnings denied, `cargo
wasm-check` and ruff clean; `cargo test --workspace` `472 passed` with 2
ignored in the core crate and `55 passed` in the linear algebra crate;
`cargo test -p popnei-linalg --no-default-features` `50 passed`; `uv run
maturin develop && uv run pytest` `257 passed`. The crate went from 47
tests to 55 and from 42 to 50 on faer, which is the eight the task added.
`cargo test -p popnei-linalg --lib cholesky -- --list` names all eight,
where the filter gave `0 tests` before, and `singular` names the one for
the error.

The row of a matrix that cannot be factored comes out the same on both
backends, counting from 0, although `dpotrf` gives the leading corner
counting from 1 and faer an index from 0. The test's matrix stops at the
middle row of three, so a backend that counted from 1 or from the other
end gives another number.

The task did two things differently from the plan. It checked the buffer
with the helper that names the argument `a` and not with the one the plan
named, which writes `g` into the message and is for the two operations
whose argument is called `g`. And it built the 1000 x 1000 of the
generator with a fixture of its own rather than share the
eigendecomposition's, which would have meant editing a test of
`eigh_lower` that this work package must leave alone.

Worth knowing: faer leaves the upper half of the buffer untouched, which
the spec asserts and nothing had checked until now, and a value that is
not finite there is let through by both backends and stays where it was,
since only the lower half is read.

### Task 2.2, the solve and the log of the determinant

`a5d431c`. `solve_with_cholesky`, `dpotrs` with `uplo` `U` in the BLAS
backend and faer's `solve_in_place_with_conj` with the matrix reference of
`b` turned the other way round, `b` holding one row for each right hand
side; and `log_determinant_with_cholesky`, twice the sum of the logs of
the diagonal, which calls neither library and which reads that diagonal
for a value that is not finite and for the `Singular` a `cholesky_lower`
would have given first. 19 tests, every number a literal of the spec.

Checked by the orchestrator: fmt, clippy with the warnings denied, `cargo
wasm-check` and ruff clean; `cargo test --workspace` `472 passed` with 2
ignored in the core crate and `74 passed` in the linear algebra crate;
`cargo test -p popnei-linalg --no-default-features` `70 passed`; `uv run
maturin develop && uv run pytest` `257 passed`. `--lib solve_with --
--list` names 12 tests and `determinant` 7, where both filters gave `0
tests` before.

**The measurement the plan asked for, and it closes the risk.** "What
could go wrong" of this work package said that nobody had measured how
large a scratch faer wants for the solve at the size the association study
works at, and that if it grew with the right hand sides then whether the
solve should also give `Memory` would be a point for the owner. It does
not grow: `solve_in_place_scratch` of faer 0.24.4 asks for **0 bytes**, at
`n` of 5 with 10000 right hand sides and at `n` of 1000 alike, because the
solve walks the two triangles in place. So the spec is right to give
`Memory` to the inverse alone, nothing is asked of the owner, and the task
added a test in the faer backend that asserts the scratch is 0, so that a
later faer which wanted memory fails a test instead of ending the process.
That test is why the faer build has 70 tests and not 69.

### Task 2.3, the inverse

`8af207e`. `invert_with_cholesky`, `dpotri` on a copy of the lower half of
the factorization in the buffer the caller gave in the BLAS backend, and
faer's `llt::inverse::inverse` in the other, whose scratch of n x n is
asked for with `try_new` of `dyn_stack` so that a machine without the
memory gets `Memory` and not the end of the process. 11 tests. faer asks
8000000 bytes at n = 1000, which is the number the spec records, and a
test of the faer backend now asserts it.

**What the task found, which is a change the spec did not have.** An `l`
whose diagonal holds an entry that is not above 0 is not a factorization
`cholesky_lower` ever gives, since that is where it stops, but a caller
holds `l` and `inverse` as two buffers and can pass one that never was a
factorization. Measured on 23 September 2026 on an `l` with a 0 at its row
1: `dpotri` gives an `info` of 2, and faer's inverse gives **no error at
all** and writes infinities and NaN into the caller's buffer. The
orchestrator checked it by taking the guard out and running the test
again: faer returns `Ok`. That is a silent wrong result, which the owner
ruled out on 21 September 2026, so the diagonal is read in the crate above
both backends and the caller gets the `Singular` that
`log_determinant_with_cholesky` already gives for the same diagonal and
that the spec already has the triangular solve give for the same faer
behaviour. No case was added to the error enum. `4d5cab5` writes it into
"The errors the seven add"; that commit comes after the code and not
before, which is not the order the `coding` skill asks for, because the
case was found while the task was written.

### The deliverables of work package 2

Every check is `cargo test -p popnei-linalg --lib <filter> -- --list`, and
every filter printed `0 tests` before the work package.

A filter matches a test by its name, and the name of a test of the solve,
the log of the determinant or the inverse says which factorization it
takes, so `cholesky` matches the tests of all four operations and not the
eight of the factorization alone.

| Deliverable | Filter | What it gave |
| --- | --- | --- |
| 1, the factorization and the row it stops at | `cholesky` | 34 tests, of which 8 are the factorization's own |
| 2, the solve with its right hand sides as rows | `solve_with` | 12 tests |
| 3, the log of the determinant | `determinant` | 7 tests |
| 4, the inverse | `invert` | 11 tests |
| 5, the `Singular` case and its message | `singular` | 3 tests |

Both backends pass: `cargo test -p popnei-linalg` `85 passed` and
`--no-default-features` `82 passed`, against 47 and 42 before the work
package. The rest, run by the orchestrator: fmt, clippy with the warnings
denied, `cargo wasm-check` and ruff clean; `cargo test --workspace` `472
passed` with 2 ignored in the core crate; `uv run maturin develop && uv
run pytest` `257 passed`.

### The review of work package 2

Six reviewers over `eb9c162`, `a5d431c`, `8af207e` and the spec commit
`4d5cab5`: spec, tests, numbers, errors, and api and architecture
together. `binding` was not sent, for the reason it was not sent for work
package 1. They used 682000 tokens between them.

**The one that would have reached a user.** Three reviewers found, each on
its own, that `solve_with_cholesky` gave `Ok` and a buffer of NaN for an
`l` whose diagonal holds an entry that is not above 0. The orchestrator
ran it on both backends: `dpotrs` wrote NaN, an infinity and an infinity
with its sign turned round, and faer wrote three NaN, both saying nothing
went wrong. numpy refuses the same system with `LinAlgError`. Task 2.3 had
given that check to the inverse because there the two backends disagree,
which is the test the spec used to decide what to check; the solve is the
case where they agree and both are wrong, which is what makes the rule
general. All three operations that take an `l` now read its diagonal
first, and `9a658ab` says so in the spec with the reason for each.

**The spec paragraph that would have misled the GWAS.** "The solve of a
factorized matrix" justified the layout of `b` with its two dimensions the
wrong way round, saying the caller passes "c right hand sides of n numbers
each" where lines 483 and 693 of pyNei's `gwas.py` have one right hand
side for each individual of one number for each coefficient. That is what
the spec's own table of those lines says and what the code does, so only
the explanation was wrong; checked against numpy 2.5.3 on 3 individuals
and 2 coefficients, where `solve` of a 2 x 2 against a 2 x 3 gives a 2 x
3. It matters because `docs/specs/gwas.md` is written from that paragraph.

**What else `9a658ab` corrected in the spec.** The log determinant of the
3 x 3 was to be asserted to the bit, which the `coding` skill does not
allow for a value that went through `ln`: five of the nine ways of moving
`ln 2` and `ln 3` by one unit in the last place give another `f64`, and
the libm Rust uses for `wasm32-unknown-unknown` already differs from this
machine on `ln 3`. It is now 1e-15 relative. The 1e-13 of the large log
determinant had to be said to be relative, the two backends being 2.4e-12
absolute away from it, so an absolute tolerance would have failed both.
And the spec now says that none of the seven checks its own result: a
matrix that is positive definite and nearly not factors, and its inverse
holds an infinity with no error, as numpy's does, which is the fit running
away that the caller detects.

**The row a `Singular` names was pinned by nothing.** Every fixture in the
crate stopped at row 1, so replacing the computed row with the constant 1
left all 85 tests passing on both backends. There are now fixtures that
stop at the first row and at the last. The orchestrator made the mutation
itself in the two places that report a row, the shared check and the BLAS
factorization, and each is caught by one test and by no other.

**Seven smaller findings, all fixed.** The diagonal walk was written twice
and work package 3 needs a third with another test on the entry, so it is
one helper now. `invert_with_cholesky` read the data before it checked the
length of `inverse`, so a caller with a short buffer was told about its
matrix. Three doc comments claimed something was the only one of its kind
when the eigendecomposition already did the same. The doc comments of
`Memory` and of `Singular` did not name their new producers. Nothing told
a caller what happens when the slice given for `l` is not a factorization,
which `product` carries a paragraph for: passing the matrix `a` there
gives `Ok` and a plausible wrong answer, (0.3848, 0.2304, 0.216) where
(1, 2, 3) is right. The workspace manifest named three routines of BLAS
and LAPACK where six are now called. And the log determinant had no test
for its dimension cap, none for a buffer longer than its dimensions, and
none of the four was tested at `n` of 1.

**What was not taken, and why.** The overflow sentinel of `dyn_stack`
would make `Error::Memory` say 0 values, but it needs a 32 bit `usize` and
an `n` of 23171, where `l` alone is 4.3 GB and cannot exist in a 4 GB
address space. The `Singular` arm of the BLAS inverse is unreachable now
that the crate checks the diagonal first; it stays, because a backend that
reports what its routine said is right whether or not anything reaches it.
And `NoConvergence` covering a negative `info` from routines that cannot
fail to converge is what the spec defines that case to be.

**The checks after the fixes**, each run by the orchestrator: fmt, clippy
with the warnings denied, `cargo wasm-check` and ruff clean; `cargo test
--workspace` `472 passed` with 2 ignored in the core crate and `93
passed` in the linear algebra crate; `cargo test -p popnei-linalg
--no-default-features` `90 passed`; `uv run maturin develop && uv run
pytest` `257 passed`. Over the whole work package the crate went from 47
tests to 93 and from 42 to 90.

## Work package 3: the thin QR, the triangular solve and the rank

### Task 3.1, the thin QR

`75855dc`. `thin_qr` above both backends: `dgeqrf` and then `dorgqr` on a
column major copy of the design in the BLAS one, which is the route the
spec chose, and faer's `qr` with `compute_thin_Q` and `thin_R` on the
buffer as it lies in the other. Nine tests, which fix the sign of each
column of `q` by making the diagonal of `r` positive before comparing,
since the sign is the backend's. The crate went from 93 tests to 102 and
from 90 to 99 on faer, and every other check is what it was.

**Deliverable 4 wanted a measurement, and it says the fast route is the
one that runs.** On a design of 10000 x 5, the best of 20 runs, `thin_qr`
takes 0.207 ms on Accelerate and 0.288 ms on faer, against the 0.165 ms
and 0.256 ms that "What the seven of the GWAS cost" measured of the
routines alone. The route this task did not take, `dgelqf` on the buffer
as it lies, is 2.98 ms, so there is no mistaking one for the other. What
the crate adds to the spec's numbers is the walk over the 50000 values of
the design for one that is not finite and the two buffers it allocates for
the caller, which is the same kind of difference "Speed" of the spec
already records for `add_self_product_lower`. `grep -c dgelqf` of the BLAS
backend is 0, which is the other half of that deliverable.

The task made the mutation the orchestrator asked for: with `rows` and
`cols` exchanged in the call to a backend, four of the nine tests fail on
each, LAPACK refusing an argument of `dorgqr` and faer giving other
numbers; the other five are the ones the crate refuses before a backend
runs.

### Tasks 3.2 and 3.3, the triangular solve and the rank

`a2dafbb` and `29ba1c7`, each after a commit of the spec, `7437913` and
`0d1669b`, as the `writing-specs` skill asks when a task needs a number
the spec has not.

`solve_upper_triangular`, `dtrtrs` with the halves and the transpose
turned in the BLAS backend and faer's
`solve_upper_triangular_in_place` in the other, with the diagonal read for
a 0 in `lib.rs` above both, since faer divides by it and gives an
infinity while `dtrtrs` gives an `info`. That check is
`refuse_a_diagonal_entry`, which the review of work package 2 had already
made general for exactly this. 12 tests. `rank`, `dgesdd` with `jobz` `N`
on a column major copy and faer's `singular_values`, counting the values
strictly above numpy's tolerance. 11 tests.

**What the two tasks added to the spec, each before its code.** The
fixture of the triangular solve had one right hand side against an `n` of
2, which cannot tell `sides` from `n`, so the task added two more traits
with their coefficients, from numpy 2.5.3, and the fixture is now three
right hand sides against an `n` of 2. It chose three and not two for the
same reason: two against an `n` of 2 would be square. And the five
matrices of the rank all have a singular value either well above the
tolerance or at it, so the task added a sixth, the 3 x 2 whose values are
all 0, of rank 0, which is the only one that tells a count above the
tolerance from a count at it.

**Deliverable 4's second measurement.** On a design of 10000 x 5, the best
of 20 runs, `rank` takes 0.1881 ms on Accelerate and 0.2124 ms on faer,
against the 0.145 ms and 0.159 ms the spec measured of the routines alone
and the 1.54 ms of the route this task did not take. `thin_qr` measured
again beside it gave 0.1855 ms, where task 3.1 had measured 0.207 ms and
a reviewer later 0.1573 ms, all on this machine and the same design:
three runs of one thing spanning 0.157 to 0.207 ms, which is the spread
between runs on a machine doing other work. Nothing is drawn from that
spread, and nothing needs to be: what the deliverable asks is whether
these are the fast route, and 1.54 ms and 2.98 ms are what the routes the
spec did not take cost.

The mutations the orchestrator asked for: `n` and `sides` exchanged in the
triangular solve fails 4 of its 12 tests on each backend, and the two
dimensions exchanged in the rank makes the 4 x 3 whose third column is the
sum of the first two give 3 instead of 2, on each backend.

### The deliverables of work package 3

Every check is `cargo test -p popnei-linalg --lib <filter> -- --list`, and
every filter printed `0 tests` before the work package.

| Deliverable | Check | What it gave |
| --- | --- | --- |
| 1, the thin QR | filter `qr` | 9 tests |
| 2, the triangular solve and its `Singular` | filter `triangular` | 12 tests |
| 3, the rank and the tolerance | filter `rank` | 11 tests |
| 4, the BLAS backend on a column major copy | `grep -c dgeqrf` 17, `grep -c dgelqf` 0, and the timings above | the fast route on both operations |

Both backends pass: `cargo test -p popnei-linalg` `125 passed` and
`--no-default-features` `122 passed`, against 93 and 90 before the work
package. The rest, run by the orchestrator: fmt, clippy with the warnings
denied, `cargo wasm-check` and ruff clean; `cargo test --workspace` `472
passed` with 2 ignored in the core crate; `uv run maturin develop && uv
run pytest` `257 passed`.

### The review of work package 3

Four reviewers over `75855dc`, `7437913`, `a2dafbb`, `0d1669b` and
`29ba1c7`: spec, tests, numbers, and errors with api and architecture
together. They used 453000 tokens between them, and seventeen findings
held.

**The one that was a wrong number.** `rank` formed numpy's tolerance as
the largest singular value times the larger dimension, and that product
times the distance from 1 to the next `f64`; numpy multiplies the last two
first. The two are the same `f64` for every matrix a study holds, that
distance being a power of two, and they part when the first product
overflows: the 2 x 2 with 1e308 and 1 on its diagonal came back rank 0 on
both backends where numpy gives 2, and numpy warns when asked to compute
it in popnei's order. "The rank" takes numpy's tolerance so that a design
popnei refuses is a design pyNei refuses, which is what this broke. Two
reviewers found it by different routes and the orchestrator ran it against
numpy before the fix and after, and put the order back afterwards to see
the new test fail.

**Two fixtures that did not pin what they were written to pin**, both
shown by mutation. Every matrix pinning the rank's tolerance is 2 x 2,
where the larger dimension and the smaller are one number, so reading the
tolerance off the smaller left all 11 rank tests passing while giving a
wrong rank for any matrix that is not square. And every `thin_qr` fixture
has two columns, so the half of `r` below its diagonal is a single entry,
and zeroing too few of them left all 9 tests passing while leaving a
reflector of LAPACK where the doc comment of `ThinQr` promises 0, which
gives a caller a matrix that is not upper triangular with no error. This
is the fourth time a review of this plan has found a fixture that could
not tell two dimensions apart. `dd11536` put the numbers for both new
fixtures into the spec, from numpy 2.5.3, before the code that uses them.

**An allocation that would end the process.** The BLAS backend writes a
column major copy of the matrix for the thin QR and for the rank, and took
it with `vec!`. Two hundred lines above, `eigh_lower` refuses to do that
for its workspace and says in a comment why. Both operations are public
and bounded only by the 2147483647 values of the crate, which is 17 GB.
`c0d8418` gives them `Memory` in the spec and the code asks for the copy
with `try_reserve_exact`.

**Thirteen smaller findings, all fixed**: `rank`'s `# Errors` gave
`NoConvergence` one of its two meanings; three `// SAFETY:` comments of
the workspace queries said the size goes into a slice it does not go into
and that `info` is not written when it is; `write_the_rows_of` would
truncate in silence for a shape the crate refuses before it, and did not
say so; four messages named an argument that was not the one at fault or
explained a limit by an output the caller never passed; the crate's own
doc comment counted eleven operations and left out two of the seven; the
module comment of `blas.rs` called `dtrtrs` one of the routines numpy
calls, which `nm` on numpy's shared library shows it is not; a doc comment
gave a singular value one digit from the spec's; the two new workspace
queries had no test where the one they follow has five; a doc comment
carried a false bound; a doc comment quoted a timing of the routines as if
it were the crate's; and four cases were untested — a wide matrix for the
rank, a square design for the QR, and the two ends of the `rows < cols`
boundary.

**What was not taken, and why.** faer and LAPACK give different ranks for
values near the ends of what an `f64` holds: faer does not converge for a
4 x 2 of 1e154 and loses subnormal values where LAPACK and numpy do not.
That is faer's own arithmetic, so `dd11536` records that the two backends
agree between about 1e-300 and 1e154, which every design of dosages and
covariates is inside, and the question of scaling the matrix to widen it
is for the owner below. And the fourth reviewer reported the rank's
tolerance as already right: every matrix it tried has a largest singular
value of 1, where both orders give the same `f64`, so it never reached the
case the other two found.

**One thing the fix could not pin, reported by the subagent that made
it.** The zeroing of the half of `r` below its diagonal is guarded on the
BLAS backend alone. faer's `thin_R` already gives 0 there, so the crate's
write changes nothing on that backend and no test can see it; the line is
right and is dead, and nothing will catch it if a later faer stops
zeroing.

**The checks after the fixes**, each run by the orchestrator: fmt, clippy
with the warnings denied, `cargo wasm-check` and ruff clean; `cargo test
--workspace` `472 passed` with 2 ignored in the core crate and `141
passed` in the linear algebra crate; `cargo test -p popnei-linalg
--no-default-features` `128 passed`; `uv run maturin develop && uv run
pytest` `257 passed`. Over the whole work package the crate went from 93
tests to 141 and from 90 to 128; the gap between the two backends is the
13 tests of `blas.rs`, which test helpers only that backend has.

## Work package 4: the solve against a lower triangular matrix

Added by the owner on 23 September 2026, after work packages 1 to 3 were
done, at the request of the session writing `docs/specs/gwas.md`. Its fit
of the null model of the logistic mixed model needs the Cholesky factor of
a covariance solved against with one right hand side for each individual,
and a Cholesky factor fills the lower half of its buffer, so the solve
that existed, which reads the upper half, could not serve it.

### Task 4.1, the spec item and the shape of the interface

`852cb7b`. The orchestrator wrote this one rather than send it out,
because what it decides is an interface: whether the lower half is a
second function beside the upper one, or an argument of the one that
exists. It is an argument. Two functions would take the same four
arguments, and either would accept the other's call and answer with a
different matrix and no error, no length telling them apart, which is the
reason `docs/specs/linalg.md` already gives for the product of two
matrices being one function with typed operands. It is not a hypothetical:
the symmetric buffer with rows (2, 1, 0), (1, 3, 2) and (0, 2, 1), against
the right hand side (8, 40, 27), gives (4, 12, 3) read as the lower half
and (6.333333333333333, -4.666666666666666, 27) read as the upper, and
both are answers a caller could believe. `cholesky_lower` and `eigh_lower`
keep their names, neither having another half in the crate for a call of
them to mean.

The timings that justify the operation are the other session's, measured
with numpy and not with this crate, and the spec quotes them as theirs and
points at their report. This plan owns the operation and not the fit.

### Task 4.2, the code of both halves, and the review of work package 4

`28925db`, and `5eb4e7b` after the review. `solve_upper_triangular` is
renamed `solve_triangular(a, n, half, b, sides)` and gains the public
`TheHalfThatHoldsTheMatrix`, whose two cases are the upper half and the
lower; there is no `solve_upper_triangular` any more: the half is the
`uplo` of `dtrtrs`, turned by the layout, and one of faer's two entry
points. Twenty tests name `triangular` where twelve did, and the twelve
assert what they asserted.

Two reviewers, spec with tests and errors with api and architecture. They
used 243000 tokens, and fourteen findings held, of which **five were in
the spec item the orchestrator had written**: the buffer the whole
one-function decision leans on was named wrongly, and the `l` the section
gives, read as its upper half, is the diagonal 2, 3, 1 and solves to (4,
13.333333333333332, 27) and not to the number quoted; "`dtrtrs`, whose
`uplo` is that half" was the opposite of what the code must pass, since
popnei's buffer read column after column is the transpose; the explanation
of faer's last place was not the mechanism, multiplying by the reciprocal
giving exactly 12 and the bits coming from distributing it; "fourteen
operations" was left beside "the last seven"; and the `Singular` case
still said "upper triangular" in both code blocks. `5385750` corrects all
five, each read back against numpy or the code first.

**The finding that mattered most for the code was a trap and not a live
defect.** The diagonal is checked with `entry == 0.0`, where the three
Cholesky operations two screens above check `entry <= 0.0`. A reviewer
changed the triangular one to match its neighbours and all 148 tests
passed — yet `thin_qr` really gives a negative diagonal, so that change
refuses every least squares fit popnei makes. Every fixture used the `r`
with its sign fixed. `7b419fa` puts the case into "How the seven are
verified" and `5eb4e7b` has the test, which is the only one that fails
under that mutation, on both backends.

Eight smaller findings were fixed with it: the doc comments said faer
answers a 0 diagonal with an infinity where it answers with a NaN as
often, by the half and the row; the enum derived no `PartialEq`, which
every other public fieldless enum of popnei derives; twelve call sites
still called the argument `r`; eleven `// SAFETY:` comments of `blas.rs`
credited `the_i32_of` with refusing a dimension of 0, which `lib.rs` does;
a fixture's doc comment carried three wrong claims, one of them a relative
error that was an absolute one; the lower half's `Singular` test is
guarded by the faer run alone, since `dtrtrs` gives that error itself, and
now says so; the lower half's `NotFinite` test reached three of six
places; and one doc comment named a function that no longer exists.

**What was not taken, and is for the owner.** A diagonal entry that is
subnormal passes the check, and then `dtrtrs` gives an infinity where faer
gives a NaN, each with no error, on an input whose exact answer an `f64`
holds. It is the one place in the crate where the two backends answer
differently and neither says so. Both answers are not finite, so the test
a caller must make on its own result catches either, and numpy gives the
infinity, so refusing such an entry would diverge from the oracle on an
input it answers. `e076edc` records it, and the three Cholesky operations
have the same hole. And that the `Singular` message says a factorization
stopped is true of one of its five producers; it is the spec's message and
the owner's to reword.

## How the whole plan was checked

Run by the orchestrator on the last commit of the branch:

| What | Command | What it gave |
| --- | --- | --- |
| the two backends | `cargo test -p popnei-linalg` and the same `--no-default-features` | `149 passed` and `136 passed`, against 42 and 37 when the plan started |
| the browser target | `cargo check -p popnei-linalg --target wasm32-unknown-unknown --no-default-features` | clean |
| the pyodide target | the same for `wasm32-unknown-emscripten` | clean |
| the whole workspace | `cargo test --workspace` | `472 passed`, 2 ignored, in the core crate. That is the core crate on this branch; the merged tree below gives 604, because `main` gained 132 tests of other sessions' work while this plan ran |
| the Python layer | `uv run maturin develop && uv run pytest` | `257 passed` |
| the TypeScript layer | `npm run build && npm test` in `js/popnei` | `pass 180` |
| the rest of the coding skill | fmt, clippy with the warnings denied, `cargo wasm-check`, ruff | all clean |

The principal component analysis and the r² reach Python and TypeScript
through the `product` that work package 1 changed, and the last three rows
are where a change in what they compute would have been seen by a user.
None of them moved.

### The machine every timing here was measured on

The owner's Apple M5 Pro, rustc 1.98, release builds, unless a number is
said to be someone else's. The BLAS and LAPACK of this machine are
Accelerate, the framework numpy also computes through, and it takes the
threads it finds; a number of faer is of faer built natively, which runs
on the same pool of threads popnei's own loops use.

Two timings quoted here were not measured by this plan. The bit count of
the `dists` module takes 1.8 ms with the rustc flag below and 5.2 ms
without it, which is why that flag is set at all; both are from
`docs/specs/pca.md`, where the input and the machine are. The fit of the
logistic mixed model at 4000 individuals, 5.75 s against pyNei's 9.959
s, was measured with numpy by the session writing `docs/specs/gwas.md`
and is in its own report.

### The merge into `main`, tried and not made

`main` moved while this plan ran: it was `bce303c` when the plan started
and `bed9031` on 23 September 2026, other sessions having merged the
performance review of the linkage disequilibrium and work of the vars
file. So "merges into `main` with no conflict" was checked again at the
end, and then checked further, because a merge that has no conflict can
still not build: both sides changed `Cargo.toml` of the workspace and
`crates/popnei/src/ld.rs`, and that file is where this plan rewrote a call
of `product`.

The merge was made in a throwaway worktree at a detached `main`, built and
thrown away; `main` itself was not touched and is still `bed9031`. It has
no conflict, the one call site both sides changed keeps this branch's form,
and the merged tree passes. Run again on the last commit of the branch,
so that the linear algebra counts are the branch's own: `cargo test
--workspace` `604 passed` with 2 ignored in the core crate, and `149
passed` in the linear algebra crate. The 604 is the merged tree and the
472 of the table above is this branch alone: `main` gained 132 tests of
other sessions' work while this plan ran, and the merge brings them in.
The linear algebra crate is at 149 either way, being this plan's alone;
`cargo test -p popnei-linalg
--no-default-features` `136 passed`; `cargo fmt --all --check`, `cargo
clippy --workspace --all-targets -- -D warnings` and `cargo wasm-check`
all clean.

The merge keeps this branch's documents and not `main`'s, which the tests
do not check, since none of them reads a document. `main` had changed
`docs/specs/linalg.md` after this branch's base, in `06a4ca2`, but that
commit is an ancestor of this branch, so there is nothing to resolve: the
spec in the merged tree is this branch's file, compared by its hash.

And the order of the merges, from the session writing
`docs/specs/gwas.md`. Its branch `spec/gwas` started from `main` with this
spec merged in, so it carries `174d795` and `7ccc709`, and it never
changed `docs/specs/linalg.md`, its copy being the older one: merging it
leaves whatever `main` holds, so the two branches merge in either order
with nothing to resolve. It should still go second, since its text names
`solve_triangular` and `TheHalfThatHoldsTheMatrix`, which `main` does not
have until this branch is in.

The trial crate is left in `tmp/`, not committed, as the plan said. What
replaces it is the cargo tests of the four work packages, which assert the
same literals through the crate's own checks and error type.

## What is waiting on the owner

Four decisions. None of them stops the merge or anything this plan built,
and each is written here with what it would take.

### 1. How your decision about the vector instructions of WebAssembly is worded

WebAssembly has instructions that work on sixteen bytes at a time, and two
different switches ask for them: a cargo feature of the library that does
faer's products, and a flag to the Rust compiler. You decided on 23
September 2026 that the feature stays on and the flag stays off, because
on this compiler the flag changed no byte of what popnei builds. The
sentence that records it is in the "Open points" section of
`docs/specs/linalg.md`.

The flag is on. `.cargo/config.toml` sets it for both WebAssembly targets,
because the `dists` module counts the bits of a pair of individuals
sixteen bytes at a time behind a compile time test for those instructions,
and without the flag the scalar loop beside it is what compiles. That code
went in on 22 September 2026, after the byte comparison your decision
rests on was measured and before the spec quoted it. Built again on 23
September 2026, the WebAssembly of the JavaScript crate in release is
2249734 bytes with the flag and 2246645 without, with different md5 sums.

Your choice is untouched: for the linear algebra the flag still changes no
file. What is wrong is only the sentence that records it. `68b29fd`
corrected every statement of fact around it and left the record itself to
you, since it is your decision.

- **Reword it**, to say that nothing is added for the sake of the linear
  algebra and that the flag is on for the `dists` module. Costs one
  sentence; the spec already says the facts in the section the bullet
  points at.
- **Leave it**, and a reader of "Open points" learns something untrue
  about the build. Costs nothing now.

Recommended: reword it. I have not, because what you decided is yours to
put in your own words.

### 2. Whether a new product replaces a buffer and a loop in the PCA

A reviewer found that `the_components_of_the_product_of_the_rows` of
`crates/popnei/src/pca.rs` writes a matrix of the traits by the
components and then copies it entry by entry into its transpose, and
that one of the two products this plan added, the one that reads both of
its operands the other way round, writes that transpose directly. It ran
both on both backends and got the same numbers.

- **Do it**, as a task on a branch of its own. Costs an afternoon with its
  review, and saves an allocation and a copy of the traits times the
  components once per call of the analysis, not once per variant, so no
  time of a whole analysis is known to change.
- **File it**, so that it is not lost. This repository has no issue open
  and none has ever been filed, so I did not open the first one without
  you.
- **Leave it.** Costs nothing; the code is correct as it stands.

Recommended: file it. It is a real simplification and nothing waits on it.

### 3. Whether a diagonal entry that is subnormal should be refused

The three operations that read the diagonal of a Cholesky factorization
refuse an entry at most 0, and the triangular solve refuses one that is 0.
An entry that is neither, but so small that one divided by it overflows,
passes both. Measured on 23 September 2026 on the 2 x 2 with 4e-309 and 1
on its diagonal against the right hand side (0.5, 1), whose exact answer
is (1.2500000000000008e308, 1) and both of whose entries an `f64` holds:
LAPACK gave an infinity for the first entry and faer a NaN, each reporting
success, and numpy 2.5.3 gave the infinity.

This is the one place in the crate where the two backends answer
differently and neither says so, which the crate's own doc comment
promises does not happen. What makes it bearable is that both answers are
not finite, so the one test a caller has to make on what came back catches
either, and the spec already requires that test of the association study.

- **Refuse it**, by reading the diagonal for an entry whose reciprocal is
  not finite. Costs one line in four places and makes the two backends
  agree; popnei then refuses an input numpy answers, which is a departure
  from the oracle the spec elsewhere works hard to match.
- **Leave it and keep the record.** `e076edc` writes the case, the
  measurement and both answers into the spec. Costs nothing, and the
  divergence stays.

Recommended: leave it. A caller that does not test its own result for
being finite is broken whichever answer it gets, and the spec makes that
test the association study's.

### 4. Whether the message of `Singular` should stop naming a factorization

A matrix that cannot be factored gives an error whose message is "the
matrix a is singular: the factorization stopped at its row 3, counting
from 0". Five operations now raise it and only one factors anything: for
the other four the matrix was handed in as a factorization, or, for the
triangular solve of a QR's `r`, was never a factorization at all.

The condition my earlier recommendation put on this has been met.
`docs/specs/gwas.md` has decided what a user sees for each case, and the
session writing it asked on 23 September 2026 for the wording to change
now, with a reason from the caller's side that is worth more than mine.
It catches this error in two places that mean different things. In the fit
of a mixed model null, where missing genotypes have made the kinship
indefinite, it raises a `ValueError` naming the kinship, and it does not
want the word "singular" reaching the user, because the matrix that failed
is not the one the user gave: it is one built from theirs. In the per
variant fit of the logistic Wald test, the message is never shown at all,
the variant getting NaN for its effect, its standard error and its
p-value, so there the error only has to be cheap to make and to match on.
Both places want the argument's name and the row, which the two fields
already carry, and neither wants a conclusion about the caller's data.

- **Reword it to drop the conclusion.** That session proposes "the
  factorization of a stopped at its row 3, counting from 0". It is true of
  the four operations that take a factorization and not of the fifth,
  where `a` is a triangular matrix and nothing factored it. A wording true
  of all five says only where it failed: "the matrix a failed at its row
  3, counting from 0". Costs a sentence of the spec, one line of the code
  and one assertion of a test, and the two fields do not change, so
  nothing that matches on the error moves.
- **Leave it.** Costs nothing. The other session says it will wrap the
  error in both of its places anyway, so no user sees the wording either
  way; what it costs is a message that is wrong about four of its five
  producers for whoever reads the crate next.

Recommended: reword it, to "the matrix a failed at its row 3, counting
from 0". I have not, because a message is a value a user sees and this is
the decision I put to the owner; a session writing another spec cannot
make it, however good its reason.

## How the work went

Every task and every fixing round went to a subagent on the model the
owner uses for coding, and every review to one reviewer per category. What
they used, which is what tells the right size of a task for the next plan:

| The work | Tokens |
| --- | --- |
| task 1.1, the typed first operand and its callers | 123k |
| task 1.2, the two new products | 144k, and 181k for the fixes of its review |
| task 2.1, the Cholesky and `Singular` | 177k |
| task 2.2, the solve and the log determinant | 156k |
| task 2.3, the inverse | 175k, and 230k for the fixes of the work package's review |
| task 3.1, the thin QR | 195k |
| tasks 3.2 and 3.3, the triangular solve and the rank | 231k, and 291k for the fixes |
| task 4.2, both halves of the triangular solve | 154k, and 219k for the fixes |
| the reviewers, 17 over the four work packages, covering 21 categories, some of them carrying two or three | 92k to 156k each |

Three things worth carrying into the next plan.

**One operation is the size of a task.** The tasks of a single operation
on two backends with its tests came to between 123k and 195k tokens. The
one that ran to 231k was two operations in one prompt, and its round of
fixes to 291k, the largest of the plan.

**One reviewer per category is what made the serious findings certain.**
Of everything the four reviews found, the defects that would have given a
user a wrong number — the Cholesky solve answering with NaN, the rank's
tolerance overflowing, and the two fixtures that could not tell two
dimensions apart — were each found by two or more reviewers independently,
by different routes. What a single reviewer saw alone was, without
exception, documentation, a message or a gap in coverage. The overlap is
what let those be fixed without asking the owner.

**The same defect kept coming back in new clothes.** Four times over, a
test fixture could not tell two things apart, because every case it used
made them equal: a result whose rows and columns were both 2, a `Singular`
that always stopped at the same row, a tolerance pinned only on square
matrices, and a triangular `r` always given with its sign fixed. Three of
the four were found by mutation and by nothing else. A plan that asks for
a fixture at each size is not enough; what caught these was reviewers
breaking the code to see which test noticed.

**And one thing the orchestrator did worst.** Of the fourteen findings of
work package 4, five were errors in the spec item the orchestrator had
written itself, including a worked contrast that named the wrong matrix
and a sentence about a LAPACK routine that was the opposite of what the
code must pass. The three tasks whose spec items were written before the
plan began had no such finding. The lesson is that a spec item written by
the session that is also running the plan wants the same first reader and
the same spec reviewer as any other, and got neither until the review.
