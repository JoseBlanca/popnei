# Plan: the linear algebra the association study needs

23 September 2026. State: under way. Approved by the owner on 23
September 2026 and carried out the same day; work packages 1 to 3 are
done, reviewed and fixed, and work package 4 was added by the owner on the
same day, after the rest was finished, and is what is left. It
builds from `docs/specs/linalg.md`, which went through its first reader
and its review and whose four open points the owner decided on 23
September 2026. Three parts of that spec are what is built: "The seven
operations of the GWAS", a Cholesky factorization and the three things
that come off it, the thin QR of a design, the solve against an upper
triangular matrix and the rank of a matrix; "The product with its first
operand turned", which lets `product` multiply a matrix by another that
has one row for each of the same individuals, as pyNei's association
study writes it fourteen times; and, in "The Rust interface", "The
signatures of the seven of the GWAS" and the one case the error enum
gains, `Singular`.

It is carried out in the worktree `.claude/worktrees/linalg-gwas` on the
branch `plan/linalg-gwas`, which starts from `spec/linalg-gwas` and not
from `main`, because that branch holds the spec; it merges into `main`
with no conflict, checked on 23 September 2026. The report is
`docs/reports/linalg-gwas.md`.

## In and out

Built, all of it in `crates/popnei-linalg` but for four calls in the core
crate:

- The seven operations of the association study, each on both backends,
  the BLAS and LAPACK of the system and faer.
- `product` computes all four of `a b`, `a b'`, `a' b` and `a' b'`, chosen
  by the two types that say how each operand's buffer is laid out. Its
  signature changes, so the three calls of it in `crates/popnei/src/pca.rs`
  and the one in `crates/popnei/src/ld.rs` are rewritten to name
  `ByTheRowsOfTheResult` for their first operand. They compute what they
  compute today.
- One case of the crate's error enum, `Singular`, for a matrix that could
  not be factored, with the row it stopped at.

Not built, with where it goes:

- **The association study**, which is what calls the seven: which of its
  four null models calls what, in what order, and what a `Singular` means
  where it is caught, is `docs/specs/gwas.md`, which is not written.
- **The two functions that turn a statistic into a p-value**, the
  complementary error function and the incomplete beta function, and **the
  kinship**: the same spec and `docs/specs/kinship.md`, neither written.
- **Anything in Python or in TypeScript.** No function this plan builds is
  reached from either: the binding crates and the packages expose the
  calculations of the core crate, and nothing of the core crate calls the
  seven until the association study exists. So no work package of this
  plan ends at a Python function, and the comparison with GMMAT and
  rrBLUP that `docs/objectives.md` names as the reference outside the
  project for an association study is made by the plan that builds
  `docs/specs/gwas.md`. What this plan compares against instead is numpy
  2.5.3, whose numbers are the literals of "How the seven are verified"
  and of "How it is verified" of "The product with its first operand
  turned".
- **The Linux and Windows builds of OpenBLAS**, and **faer against
  OpenBLAS on x86**: this machine is a Mac, its BLAS is Accelerate, and
  both are for the plan that first builds those wheels.

The spec has no open point. The four it had are decided, and "Open points"
of the spec says where each decision now lives.

## What has to be in place

Every layer of popnei exists, so every check of the `coding` skill runs
during this plan and none is reported as not there. On the commit this
starts from, `7ccc709` of `spec/linalg-gwas`, each of them was run in the
worktree on 23 September 2026:

- `cargo fmt --all --check` and `cargo clippy --workspace --all-targets --
  -D warnings`: clean.
- `cargo test --workspace`: 472 tests in the core crate and 42 in the
  linear algebra crate, 2 ignored, all passing.
- `cargo test -p popnei-linalg`: 42. The same with
  `--no-default-features`, which is the faer backend: 37, the five that
  are fewer being tests that only the BLAS backend has.
- `cargo wasm-check`: clean. It checks both crates for both wasm targets
  with the warnings denied, and on a wasm target the cargo feature `blas`
  turns nothing on, so it covers the two commands of "Checks" below as
  well; those are run at the end of the plan anyway.
- `uv run ruff format --check && uv run ruff check`: clean.
- `uv run maturin develop && uv run pytest`: 257 passed, 0 failed.
- `npm run build && npm test` in `js/popnei`: 180 pass.

Two things to know before running them. `uv run maturin develop
--release` makes one pytest test fail,
`test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file`: it
sends the process a ctrl-C 0.1 s into a write that a release build has
already finished. The check of the `coding` skill is the debug build, and
with it the 257 pass. And `npm run build` in `js/popnei` needs
`wasm-bindgen` on the PATH; it was there on 23 September 2026.

The working reference for every task is a trial crate that was written
before the spec and that runs all seven operations on both backends:
`tmp/linalg_gwas_trial/` in this worktree, copied there from the spec's
worktree on 23 September 2026 and not in git, `tmp/` being ignored. It
holds every routine call with the transposes and the layouts worked out,
every literal of the spec asserted, and the timings. `cargo run --release`
inside it prints the checks on the LAPACK backend and `cargo run --release
--no-default-features` on faer; both were run in this worktree on 23
September 2026 and every check of both printed `ok`. It is a reference and
not the deliverable: it has neither the checks of `lib.rs` nor the error
enum, it factors with `unwrap` and `assert_eq!`, and its BLAS thin QR and
rank are the slow route through `dgelqf` and `dgesdd` on the buffer as it
lies, which "What the seven of the GWAS cost" measured at 2.98 ms and 1.54
ms for a design of 10000 x 5 against 0.165 ms and 0.145 ms for the route
through a column major copy that the spec chose. The functions named
`_by_transposing` in its `src/blas_backend.rs` are the ones to copy: the
name is what they do first, which is to write the transpose of the matrix
into a buffer of their own, and that buffer is the matrix held column
after column, which is what those two routines are fast on.

## Work package 1: the product with its first operand turned

### What it gives

`product` computes all four of `a b`, `a b'`, `a' b` and `a' b'`, so
that the association study can multiply a matrix by another laid out
over the same individuals without writing the transpose of one of them
into a buffer of its own, which costs more than the product it feeds:
writing the transpose of a 1000 x 1000 takes 0.361 ms, and the product
of that matrix with a 1000 x 5 takes 0.102 ms on Accelerate and 0.069 ms
on faer, from "What it costs" of the spec. It stops inside the core
crate, as every work package of this plan does. Its caller today is the
principal component analysis and the r² of the linkage disequilibrium,
which keep their numbers.

### Its deliverables

1. `product` takes a `TheFirstOperand` beside its `TheSecondOperand`,
   with the two cases and the doc comments of "The Rust interface" of
   the spec, and the four calls of it in the core crate name
   `ByTheRowsOfTheResult`, the case that says a buffer holds one row for
   each row of the result, which is what all four of them hold today.
   The check: `grep -c ByTheRowsOfTheResult crates/popnei/src/pca.rs
   crates/popnei/src/ld.rs` gives 3 and 1, where today `grep -rc
   TheFirstOperand crates/` gives 0 in every file; and `cargo test
   --workspace` passes with 472 tests in the core crate, the same tests
   as today and none of them changed.
2. The four combinations give the same 2 x 2 matrix, rows (5, 1) and
   (2, 9), from the four pairs of operands that "How it is verified" of
   "The product with its first operand turned" lists, asserted exactly,
   and the test runs on both backends. The check: `cargo test -p
   popnei-linalg --lib four_ways -- --list` names the test, where today it
   prints `0 tests`, and `cargo test -p popnei-linalg` and the same with
   `--no-default-features` both pass.
3. The two combinations that are new refuse the same wrong dimensions and
   the same value that is not finite as the two that exist, each buffer
   checked against the dimensions of the way it is read. The check: `cargo
   test -p popnei-linalg --lib first_operand -- --list` names the tests,
   where today it prints `0 tests`.
4. The comment of `crates/popnei-linalg/Cargo.toml` about the cargo
   feature `wasm-simd128-enable` of `gemm` says what the owner decided on
   23 September 2026, that the feature stays on and the rustc flag `-C
   target-feature=+simd128` stays off, and names no open point. The check:
   `grep -c "Open 1" crates/popnei-linalg/Cargo.toml` is 0, where today it
   is 1.

### What it stands on

Nothing of this plan. `product`, its two backends and its tests are on the
branch this starts from.

### Its tasks

- [x] 1.1 `TheFirstOperand` and the new signature of `product` in
      `crates/popnei-linalg/src/lib.rs`, with the two combinations that
      exist today wired to the backend functions they already call, the
      27 tests that name `product` moved to the new call form, the four
      calls in `crates/popnei/src/pca.rs` and `crates/popnei/src/ld.rs`
      rewritten, and the comment of `crates/popnei-linalg/Cargo.toml`. It
      is one task and one commit because the signature does not compile
      without its callers. Built from "The Rust interface" and "The
      product with its first operand turned" of `docs/specs/linalg.md`.
      Serves deliverables 1 and 4. Needs nothing.
- [x] 1.2 The two combinations that are new: `c = a' b` and `c = a' b'`
      on both backends, one function each in
      `crates/popnei-linalg/src/blas.rs` and in
      `crates/popnei-linalg/src/faer.rs`, reached from the `match` in
      `product`, with the length of each buffer checked against the way
      that operand is read. Its tests are the four combinations of "How it
      is verified" of "The product with its first operand turned" and the
      dimensions each one refuses. Built from that section and from "What
      it costs", which says which `dgemm` and which `matmul` each one is.
      Serves deliverables 2 and 3. Needs 1.1.

### What could go wrong

The principal component analysis and the r² are what say that 1.1 changed
no number, and they say it only if their tests were not touched. A task
that edits an assertion of `pca.rs` or of `ld.rs` has removed the check it
was being checked by; the four calls change how they name their first
operand and nothing else.

The four combinations differ by which buffer each routine reads the other
way round, and each is given the operands that make it compute the one
matrix, rows (5, 1) and (2, 9). That matrix is not symmetric, so a backend
that wrote the transpose of the result fails on it, and a backend that
read one operand the way another combination reads it gives something else
for at least one of the four. `a' b'` has no caller in popnei and is built
and tested like the other three, which a reviewer would otherwise ask
about.

## Work package 2: the Cholesky factorization and the three things off it

### What it gives

The core crate can factor a symmetric positive definite matrix and then
solve against it, take the log of its determinant and invert it, each from
the one factorization, which is half the arithmetic of the LU numpy uses
and which refuses exactly the matrices that are not positive definite. A
matrix that is refused gives `Singular` with the row the factorization
stopped at, counting from 0. It stops inside the core crate: the module
that calls these is the association study, `docs/specs/gwas.md`, which is
not written.

### Its deliverables

1. `cholesky_lower` gives the factorization of the 3 x 3 of "How the
   seven are verified" exactly, gives `Singular` at the row 1 for the 3
   x 3 of "The errors the seven add", and gives the first and the last
   entries of the diagonal of the 1000 x 1000 within 1e-12 relative, on
   both backends. That 1000 x 1000 is the matrix "How the seven are
   verified" builds with the xorshift generator the spec writes out, the
   same one the eigendecomposition is already checked on, and every
   large literal of this work package is of it. The check: `cargo test
   -p popnei-linalg --lib cholesky -- --list` names the tests, where
   today it prints `0 tests`, and `cargo test -p popnei-linalg` and the
   same with `--no-default-features` both pass.
2. `solve_with_cholesky` gives (1, 2, 3) for the right hand side (8, 40,
   27) and gives (1, 2, 3) and (1, 0, 0) for the two right hand sides of
   "How the seven are verified", one row each, within 1e-14; and the first
   three entries and the sum of the solution of the 1000 x 1000 within
   1e-11 relative. The two right hand sides are what catches a backend
   that read the rows of `b` as its columns. Same check command with the
   filter `solve_with`, `0 tests` today.
3. `log_determinant_with_cholesky` gives 3.58351893845611 for the 3 x 3
   exactly and 3963.7986384485084 for the 1000 x 1000 within 1e-13, and
   gives `Singular` and `NotFinite` for a diagonal that holds a value that
   is not above 0 or that is not finite. Filter `determinant`, `0 tests`
   today.
4. `invert_with_cholesky` gives the six entries of the lower half of the
   inverse of the 3 x 3 within 1e-15 relative, against the exact values
   7/18, -5/18, 5/9, 1/3, -2/3 and 1 and not against numpy's, which its LU
   leaves up to 2 units in the last place away; and the first and last
   entries of the diagonal of the inverse of the 1000 x 1000 and its
   trace, within 1e-11. Filter `invert`, `0 tests` today.
5. The error enum has `Singular` with the fields and the doc comment of
   "The Rust interface", and its message is the one "The errors the
   seven add" writes out, "the matrix a is singular: the factorization
   stopped at its row 3, counting from 0", where the 3 is the spec's
   example of the shape and not a case of this plan. A cargo test
   asserts that shape on the row the case it builds gives. Filter
   `singular`, `0 tests` today.

### What it stands on

Nothing of work package 1: the four operations here take no `product`.
They are after it because `product` changes a function the core crate
already runs and this does not, and a regression there is what would stop
the plan.

### Its tasks

- [x] 2.1 `Singular` in the error enum of
      `crates/popnei-linalg/src/lib.rs` and `cholesky_lower` above both
      backends, with `dpotrf` in `blas.rs` and
      `llt::factor::cholesky_in_place` in `faer.rs`, the row that LAPACK
      counts from 1 and faer from 0 turned into the one the spec gives,
      and the checks of the dimensions and of the values of the lower half
      in `lib.rs`. Built from "The Cholesky factorization" of "What the
      seven give", "The errors the seven add" and "How the seven are
      verified". Serves deliverables 1 and 5. Needs nothing.
- [x] 2.2 `solve_with_cholesky` and `log_determinant_with_cholesky`. The
      solve is `dpotrs` and `llt::solve::solve_in_place_with_conj`, with
      `b` of `sides` rows of `n`, one row for each right hand side, which
      is the layout the callers of lines 483 and 693 of pyNei's `gwas.py`
      already have. The log determinant calls neither library: it is twice
      the sum of the logs of the diagonal of the factorization, and the
      diagonal is what it checks for a value that is not finite and for
      the `Singular` a `cholesky_lower` would have given first. Built from
      "The solve of a factorized matrix" and "The log of the determinant"
      of "What the seven give". Serves deliverables 2 and 3. Needs 2.1 for
      the tests, which factor before they solve.
- [x] 2.3 `invert_with_cholesky`: `dpotri` on a copy of the factorization
      in the buffer the caller gave, and faer's
      `llt::inverse::inverse`, whose scratch of n x n is asked for with
      `try_new` of `dyn_stack` so that a machine without the memory gets
      `Memory` and not the end of the process. Built from "The inverse of
      a factorized matrix". Serves deliverable 4. Needs 2.1.

### What could go wrong

The spec names `Memory` for one allocation of the seven, faer's scratch of
the inverse, and for no other. Task 2.2 uses faer's
`llt::solve::solve_in_place_scratch`, which the trial crate takes with
`MemBuffer::new`, and that ends the process if the allocation fails. How
large that scratch is at the size the association study solves at, `n` of
5 with 10000 right hand sides, has not been measured. The task measures it
and reports the number; if it grows with the right hand sides, whether the
solve also gives `Memory` is a point for the spec and for the owner, and
the meanwhile is what the spec says, which is `MemBuffer::new`.

Both backends have to give the same row for a matrix that is not positive
definite, and they count it differently: `dpotrf` gives the leading corner
counting from 1 and faer gives an index from 0. Both were checked on the 3
x 3 of "The errors the seven add" and both gave the row 1 once turned.

## Work package 3: the thin QR, the triangular solve and the rank

### What it gives

The core crate can fit a linear model to more individuals than
coefficients, which is the thin QR of the design and then the solve
against the upper triangular matrix it gives, and can ask how many of a
design's columns are independent, which is what refuses a design whose
covariates repeat each other before any model is fitted. It stops inside
the core crate, as the two work packages before it do.

### Its deliverables

1. `thin_qr` of the 4 x 2 design of an intercept and one covariate gives
   the upper triangular matrix with rows (2, 5) and (0, 2.23606797749979)
   and the two columns of the orthogonal one that "How the seven are
   verified" lists, within 1e-14, with the sign of each column taken so
   that the diagonal is positive, since the sign is the backend's. The
   check: `cargo test -p popnei-linalg --lib qr -- --list` names the
   tests, where today it prints `0 tests`, and `cargo test -p
   popnei-linalg` and the same with `--no-default-features` both pass.
2. `solve_upper_triangular` gives the coefficients (-1, 2) for the trait
   (1, 3, 5, 7) within 1e-14, which is an exact fit, so a backend that
   read the matrix the wrong way round gives something else; and gives
   `Singular` at the row 1 for the matrix with rows (2, 5) and (0, 0),
   which the crate reads the diagonal for above the backends, since faer
   divides by that 0 and gives an infinity while `dtrtrs` gives an `info`.
   Filter `triangular`, `0 tests` today.
3. `rank` gives 2, 2 and 1 for the three designs of "How the seven are
   verified" and 2 and 1 for the two 2 x 2 matrices that sit either side
   of the tolerance, which is what pins the tolerance itself and not only
   the singular values. Filter `rank`, `0 tests` today.
4. The BLAS backend gives `dgeqrf`, `dorgqr` and `dgesdd` a column major
   copy of the matrix, which is the route "What the seven of the GWAS
   cost" measured at 0.165 ms and 0.145 ms for a design of 10000 x 5
   against 2.98 ms and 1.54 ms on the buffer as it lies. The check: `grep
   -c dgeqrf crates/popnei-linalg/src/blas.rs` is 1 at least, where today
   it is 0, and `grep -c dgelqf` of the same file is 0, which is the other
   route; and the task times `thin_qr` and `rank` on a design of 10000 x
   5, once for each backend, and reports the four numbers against the four
   the spec has. The slow route is eighteen times the fast one for the
   QR and eleven times for the rank, so a timing tells them apart.

### What it stands on

`Singular` of task 2.1, which `solve_upper_triangular` also gives. Nothing
else.

### Its tasks

- [x] 3.1 `thin_qr`: `dgeqrf` and then `dorgqr` on a column major copy of
      the matrix in `blas.rs`, which is the reference's
      `thin_qr_by_transposing` and not its `thin_qr`, and faer's `qr` with
      `compute_thin_Q` and `thin_R` in `faer.rs`, both writing the two
      results row after row with the lower half of the triangular one set
      to 0. Built from "The thin QR of the design" of "What the seven
      give" and its part of "How the seven are verified". Serves
      deliverables 1 and 4. Needs nothing of this work package.
- [x] 3.2 `solve_upper_triangular`: `dtrtrs` with the halves and the
      transpose turned as the reference has them, and faer's
      `solve_upper_triangular_in_place`, with `b` laid out as it is for
      the Cholesky solve, and the diagonal read for a 0 in `lib.rs`, above
      both backends. Built from "The solve against an upper triangular
      matrix" and from the paragraph of "The errors the seven add" that
      says why the diagonal is read here. Serves deliverable 2. Needs
      2.1 for `Singular`.
- [x] 3.3 `rank`: `dgesdd` with `jobz` `N` on a column major copy in
      `blas.rs` and faer's `singular_values` in `faer.rs`, and above them
      the count of the values strictly above numpy's tolerance, `s *
      max(rows, cols) * 2.220446049250313e-16`, with `s` the largest
      singular value and the last number the distance from 1 to the next
      `f64` above it. Built from "The rank" of
      "What the seven give". Serves deliverables 3 and 4. Needs nothing of
      this work package.

The three tasks need nothing of each other, and they are run in order
anyway: each adds a function to the same three files, and one tree has one
writer per file.

### What could go wrong

`rank` gives a count and no singular values, so a test at it cannot assert
what the decomposition gave. The five matrices of deliverable 3 are chosen
so that the counts alone pin both the values and the tolerance, and the
two 2 x 2 matrices are the pair that fails if the tolerance is taken at
the wrong threshold. So the five are asserted as the spec gives them, and
one swapped for a matrix of the implementer's own would leave the
tolerance unchecked.

`thin_qr` and `rank` are the two places where the BLAS backend copies, and
the copy is what makes them eighteen and eleven times faster on a tall
design.
A task that gives the routines the buffer as it lies will still pass every
deliverable above, since the trial measured the two routes to agree to
1e-13, and will be slow at the size the association study runs at.

## Work package 4: the solve against a lower triangular matrix

Added on 23 September 2026, after work packages 1 to 3 were done, by the
owner's order, at the request of the session writing `docs/specs/gwas.md`.

### What it gives

The core crate can solve `l x = b` for a lower triangular `l` as well as
for an upper triangular one, which is what the fit of the null model of
the logistic mixed model needs: the trace that each step of it wants comes
from the Cholesky factor of the covariance solved against with one right
hand side for each individual, and a Cholesky factor is lower triangular,
so the upper form cannot serve it. It stops inside the core crate, as
every work package of this plan does.

The session writing `docs/specs/gwas.md` measured what it is worth, with
numpy 2.5.3 on Accelerate on the owner's Apple M5 Pro at 4000 individuals,
against pyNei's fit of 9.959 s: 5.75 s with the lower form, 7.0 s with
`solve_with_cholesky`, which is two triangular solves where one is wanted,
and 5.29 s with an inverse of a triangular matrix, which would be another
operation written twice. Those numbers and the fit they belong to are
`docs/reports/glmm-method/README.md` of the branch `spec/gwas`, and this
plan neither checks nor owns them: what it builds is the operation.

### Its deliverables

1. The spec says what the lower form gives, in "The signatures of the
   seven of the GWAS" and in "The errors the seven add", in a commit that
   comes before the code. The check: `grep -c "lower triangular" docs/specs/linalg.md`
   is above what it is today, and the commit of the spec is an ancestor of
   the commit of the code.
2. The crate solves `l x = b` for a lower triangular `l` on both backends,
   with the numbers of "How the seven are verified", and reads the
   diagonal for a 0 above the backends as the upper form does, since faer
   divides by it and gives an infinity where `dtrtrs` gives an `info`. The
   check: `cargo test -p popnei-linalg --lib triangular -- --list` names
   more tests than the 12 it names today, among them the lower form's, and
   `cargo test -p popnei-linalg` and the same with `--no-default-features`
   both pass.
3. The upper form computes what it computed: the 12 tests that name
   `triangular` today pass unchanged.

### What it stands on

Work package 3, whose `solve_upper_triangular` this extends, and task
2.1's `Singular`, which both forms give.

### Its tasks

- [ ] 4.1 The spec item: what the lower form gives, its signature and its
      errors, and which of the two shapes the interface takes, written
      into `docs/specs/linalg.md` in a commit of its own. Built from the
      request of the session writing `docs/specs/gwas.md` and from the
      existing "The solve against an upper triangular matrix". Serves
      deliverable 1. Needs nothing.
- [ ] 4.2 The code and its tests, in `crates/popnei-linalg/src/lib.rs`,
      `blas.rs` and `faer.rs`: `dtrtrs` takes a `uplo` already and faer
      has `solve_lower_triangular_in_place` beside the upper one. Serves
      deliverables 2 and 3. Needs 4.1.

### What could go wrong

The two forms take the same arguments and give different answers, which is
the shape of mistake "The Rust interface" made `product` one function with
typed operands to stop: a caller that names the wrong one gets a wrong
matrix and no error, since no length can tell them apart. Task 4.1 decides
which shape the interface takes and says why.

## How the whole plan is checked

The sum of the work packages, and these four besides, run at the end on
the last commit of the branch:

- `cargo test -p popnei-linalg` and `cargo test -p popnei-linalg
  --no-default-features`, which are the two backends and which "How the
  seven are verified" asks for by name. Both pass, and the second runs
  every test of the first but the five of the workspace of `dsyevd`.
- `cargo check -p popnei-linalg --target wasm32-unknown-unknown
  --no-default-features` and the same for
  `wasm32-unknown-emscripten`, which are the two targets popnei ships to
  the browser and to pyodide.
- `uv run maturin develop && uv run pytest`, 257 passed, and `npm run
  build && npm test` in `js/popnei`, 180 passing: the principal component
  analysis and the r² reach Python and TypeScript through the `product`
  that work package 1 changes, and these are where a change in what they
  compute would be seen by a user.
- The trial crate is left where it is and is not committed. What replaces
  it is the cargo tests of the three work packages, which assert the same
  literals through the crate's own checks and error enum.
