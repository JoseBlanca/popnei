# Work report: the linear algebra the association study needs

The plan `docs/plans/linalg-gwas.md` is under way, on the branch
`plan/linalg-gwas`, in the worktree `.claude/worktrees/linalg-gwas`,
where it started on 23 September 2026. It builds, in
`crates/popnei-linalg`, the seven operations that a genome wide
association study needs: a Cholesky factorization, which is the
factorization of a symmetric matrix that no vector makes negative, and
the solve, the log of the determinant and the inverse that come off it;
the thin QR of a design, which is what fits a linear model to more
individuals than coefficients; the solve against an upper triangular
matrix; and the rank, how many of a design's columns are independent. It
also gives `product` a typed first operand, so that it computes all four
of `a b`, `a b'`, `a' b` and `a' b'`, and adds one case to the crate's
error enum, `Singular`. The spec behind it is `docs/specs/linalg.md`.

Where the plan stands on 23 September 2026: work package 1 is done,
reviewed and fixed, and work package 2 is done, reviewed and fixed. Work package 3 has not
started.
Work package 3 has not started. Two things
are waiting on the owner and neither stops the plan; the last section of
this report says what they are.

This report is written as the work goes. Each work package gets a section
below when it is done, with the command that checked each deliverable and
what it gave, what was changed in the plan and why, what the review
found, and what the owner should know. When the plan is done, what the
owner reads first goes at the top of this file.

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
| both wasm targets build | `cargo wasm-check` | clean |
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
backend call each. It is the right call for one commit inside a work
package, and nothing outside the crate can reach it: no caller in popnei
names that case.

That subagent used 123480 tokens.

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

That subagent used 144252 tokens.

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

That subagent used 180622 tokens for the fixes.

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

Two things the task did differently, both right. It checked the buffer
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

That subagent used 176640 tokens.

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

That subagent used 155982 tokens.

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

| Deliverable | Filter | What it gave |
| --- | --- | --- |
| 1, the factorization and the row it stops at | `cholesky` | 34 tests |
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

Those subagents used 176640, 155982 and 174752 tokens.

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

That subagent used 229755 tokens for the fixes.

## What is waiting on the owner

Neither of these stops the plan, and work packages 2 and 3 do not depend
on them.

**How the decision about the vector instructions of WebAssembly is
worded.** "Open points" of `docs/specs/linalg.md` records the owner's
decision of 23 September 2026 as "the rustc flag stays off". The flag is
on: `.cargo/config.toml` sets it for both wasm targets because
`sums_of_two` of the `dists` module counts the bits of a pair sixteen
bytes at a time behind `cfg(target_feature = "simd128")`, 1.8 ms against
5.2 ms, and that code went in on 22 September 2026 in `7f3b6cc`, after the
byte comparison the decision rests on was measured and before the spec was
written. Built again on 23 September 2026 on the same machine, the
WebAssembly of `crates/popnei-js` in release is 2249734 bytes with the
flag and 2246645 without it, with different md5 sums. What the owner chose
is untouched, since for the linear algebra the flag still changes no file,
and `68b29fd` corrected every statement of fact around it; the one
sentence that is theirs to write is how the decision itself is put.

**Whether `a' b'` should replace a buffer and a loop in the principal
component analysis.** The `architecture` reviewer found that
`the_components_of_the_product_of_the_rows` of `crates/popnei/src/pca.rs`
writes a matrix of the traits by the components and then copies it entry
by entry into its transpose, and that the fourth combination writes that
transpose directly. It ran both on both backends and got the same numbers.
It costs an allocation and a copy of the traits times the components once
per call, not once per variant, so no result and no time of a whole
analysis is known to change. It changes the principal component analysis,
which work package 1 is not allowed to do, so it is left for the owner to
put in a task or an issue. This repository has no issue open and none has
been filed, so none was filed for this.
