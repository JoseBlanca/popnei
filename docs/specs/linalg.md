# The linalg crate: the linear algebra and its two backends

September 2026. Every calculation of popnei that works on a block as a
matrix, the PCA, the kinship, the GWAS, needs a few operations of linear
algebra, and this crate is the one place that has them: natively they run
on the BLAS and LAPACK libraries of the system, the same ones numpy
calls, and in WebAssembly, where there is none, on faer, a linear
algebra library written in Rust. There is no code. This spec develops the
row `linalg` of the table in section 9 of `docs/architecture.md` and
decision 5 of `docs/rust_core.md`, which chose the two backends. It
covers the crate, its backends, the builds, and thirteen operations. Two
are the product of a matrix with itself and the eigendecomposition of a
symmetric matrix, which `docs/specs/pca.md` calls. Four are the product
of two matrices, one for each way round the two can be laid out, of which
`docs/specs/pca.md` and the r² of `docs/specs/ld.md` call two and the
genome wide association study calls a third: "The product with its first
operand turned" is the one that adds them. The last seven are what the
association study calls besides, read from pyNei, and they are under "The
seven operations of the GWAS": a Cholesky factorization and the solve,
the log determinant and the inverse that come off it, the thin QR of the
design, the solve against an upper triangular matrix, and the rank of a
matrix.

The routines of BLAS and LAPACK are reached through the crates `blas`
0.23 and `lapack` 0.20, and each is an `unsafe fn` over slices whose
lengths nothing checks against the dimensions it is given, while the
core crate forbids `unsafe`. The owner decided on 22 September 2026 that
the linear algebra is a crate of its own, `crates/popnei-linalg`, where
the `unsafe` blocks of those calls live, each with a comment that says
why the slices are long enough for the dimensions, and that the core
crate keeps its forbid. The options not taken were the `linalg` module
inside the core crate with the forbid lowered to a deny, and a crate that
wraps the calls safely, `ndarray-linalg` or `nalgebra-lapack`, which
brings an array library of its own for five routines. The owner decided
the same day that faer is also a native backend, chosen by a cargo
feature; the option not taken was faer for the wasm targets alone.

Four documents say something else and are changed by the plan that
builds this crate: section 8 of `docs/architecture.md`, whose layout has
three crates, and its section 9, whose row names a module; the line of
`docs/objectives.md` that says "one small linear algebra module"; two
lines of the `coding` skill, that linear algebra goes through the
`linalg` module and that the backend is chosen by the target and not by
a feature; and `docs/glossary.md`, which gets the crate under "The
layers".

## The crate, its backends and the builds

### What it gives

One interface of a few functions over matrices held as `&[f64]`, row
after row, the layout of the blocks and of everything the core crate
holds, and two backends behind it that give the same numbers within the
tolerance of "How it is verified". A caller of the core crate does not
know which backend ran.

The backend is BLAS and LAPACK when the target is native and the cargo
feature `blas` is on, which it is by default, and faer otherwise: on
both wasm targets, where there is no BLAS, and natively when the crate
is built with `--no-default-features`. A cargo feature can only add, so
the feature that a user turns off is the one that links BLAS, and not
one that turns faer on: the crates of BLAS and LAPACK are optional
dependencies of the targets that are not wasm, and `blas` is the
feature that enables them. faer is a dependency in every build. With
the feature off the crate links faer alone, and the interface is the
same. So `cargo test` runs the whole crate on either backend, and a
machine with no BLAS and no Fortran compiler builds popnei. It is also
how the open question of `docs/rust_core.md`, whether faer is level with
OpenBLAS on x86, gets measured: the same binary built twice.

The BLAS backend natively, from decision 5 of `docs/rust_core.md`:
Accelerate on macOS, a framework of the system that `accelerate-src`
links and that ships nothing, and OpenBLAS in the Linux and Windows
wheels through `openblas-src` built statically, which builds OpenBLAS
from source on the machine that builds and needs a Fortran compiler for
its LAPACK. The Linux and Windows builds have not been made; this machine
is a Mac and the Linux development container does not exist yet. What
they need is settled in the plan that first builds them.

Accelerate holds two LAPACKs, the one of 2009 under the plain names,
`dsyevd_`, which `accelerate-src` links, and a newer one under names
with a suffix, which numpy 2.5 links. Both were timed below and the
older one is not slower on the eigendecomposition; nothing else was
compared.

faer, 0.24.4, the version the spike of `docs/rust_core.md` used, with its
default features off and `std` and `linalg` on. Its `rayon` feature does
not build for emscripten, section 3.2 of that document, so it is on only
off the wasm targets.

### The wasm builds and the vector instructions

WebAssembly has a set of 128 bit vector instructions, which work on
sixteen bytes at a time, and popnei's wasm builds use them. Two different
switches ask for them and both are named after them, which is why they are
kept apart here and never called `simd128` alone:

- **The feature of `gemm`**, `wasm-simd128-enable`, a cargo feature of the
  crate that does faer's products. It is **on**, named by
  `crates/popnei-linalg/Cargo.toml` under `cfg(target_family = "wasm")`,
  and it is what the speed comes from.
- **The rustc flag**, `-C target-feature=+simd128`, which tells the
  compiler to emit those instructions for the crates it compiles. It is
  **off**.

The feature is what pays. Measured on 21 September 2026 with the trial
crate of `docs/specs/pca.md`, the product of a block of 5000 variants x
1000 individuals with itself, its lower half, under node 26: 306 ms with
neither switch and 187 ms with both. With the feature off and the flag on,
the analysis of 100000 variants x 1000 individuals under node takes 7.066
s against 4.500 s with the feature on, which task 4.2 of
`docs/plans/pca.md` measured on 22 September 2026.

The flag adds nothing to the products. With the feature on, the module
built with the flag and the module built without it were the same file,
byte for byte, checked on three builds from empty target directories with
rustc 1.98 on 22 September 2026, the flag being on the rustc command line
of every crate of the build that carries it. Nor would the flag cost any
browser: the feature already emits those instructions, so the floor of
Chrome and Edge 91, Firefox 89 and Safari 16.4 that `docs/objectives.md`
sets is what popnei ships today either way, and Open 6 of
`docs/specs/pca.md` details it.

The flag is on all the same, and for another module. It is set for both
wasm targets by `.cargo/config.toml`, and by the `wasm-check` alias of
that file, because `sums_of_two` of the `dists` module counts the bits of
a pair of individuals sixteen bytes at a time behind
`cfg(target_feature = "simd128")`, and without the flag the scalar loop
beside it is what compiles: 5.2 ms against 1.8 ms, which
`docs/specs/pca.md` has. That code went in on 22 September 2026, in the
commit `7f3b6cc`, after the byte comparison above was made and before this
spec was written, so the comparison stopped holding for the module popnei
ships. Measured again on 23 September 2026 on the same machine, the
WebAssembly of `crates/popnei-js` built in release with the flag is
2249734 bytes and without it 2246645, and their md5 sums differ: the 3089
bytes are that bit count. What was measured of the linear algebra is
unchanged, since nothing of this crate reads that feature and faer's
kernels carry their own annotation.

The owner decided on 23 September 2026 that the feature stays on and that
nothing is added for the sake of the linear algebra, since on this rustc
the flag changes no file of it. The options not taken were to set the flag
for the products in both wasm builds and in the wasm package alone; on
rustc 1.98 all three give the same products. What is not known is whether
a later rustc, or a dependency that stops annotating its own functions,
would make the flag matter to them, and the measurement that would say so
is the byte comparison above, run again on a build with the `dists` module
held still.

The decision was recorded on 23 September 2026 as "the flag stays off",
which was read back on the same day against the repository and is not
what the repository does: the paragraph above has where the flag is set
and what it changes. What the owner chose is unaffected, since the choice
was about the products, and how the decision is worded is theirs to
settle.

Where the flag would go if it were ever wanted for the products alone:
into the two commands that build for wasm, `build:wasm` of
`js/popnei/package.json` and `scripts/build_pyodide_wheel.sh`, through
`RUSTFLAGS`, which is how it was tried. Putting it in `.cargo/config.toml`
has a cost the `wasm-check` alias already pays: a `rustflags` for a target
in that file makes cargo ignore the `build.rustflags` that the alias would
otherwise pass to deny warnings, so the alias carries the deny and the
flag together in its own `--config`, and a flag added to a target of that
file has to be added to the alias too. The pyodide wheel built with the
flag in `RUSTFLAGS` loaded under pyodide on node 26 on 22 September 2026
and passed the smoke test of `tests/pyodide/`.

### Threads

Each backend has threads of its own and the core crate never nests them
inside rayon, section 3 of the architecture: a product or a
decomposition is called from outside a rayon loop. No calculation of
popnei takes a `num_threads`, which the owner decided on 22 September
2026 in `docs/specs/dists.md`, so the crate has no function for it
either, and each backend takes the threads its own rule gives: Accelerate
the cores it finds, unless `VECLIB_MAXIMUM_THREADS` is set when the
process starts, and the crate cannot pin it for one call; OpenBLAS the
same with `OPENBLAS_NUM_THREADS`; faer natively the global pool of rayon,
the one popnei's own loops run on, which `RAYON_NUM_THREADS` sizes, and
one thread in wasm. The product of the block above takes Accelerate 10.5
ms with the variable at 1 and 7.4 ms without it.

## The four operations of the PCA and the LD

### What they give

**The product of a matrix with itself.** For A of r rows and c columns,
`G += A'A`, the c x c matrix whose entry i, j is the sum over the rows of
A of the value in column i times the value in column j. Only the lower
half of G is written, j ≤ i, and the upper half is left as it was; a
caller that needs the whole matrix mirrors it. It is `dsyrk` in BLAS and
faer's `matmul` for triangular results, told that the result is
`TriangularLower`, both of which compute the half alone: 95 ms against
174 ms for the block above with faer on one thread. The PCA adds the
product of every block to one G this way, and a block whose rows all had
no variance gives an A of no rows, which adds nothing and is not an
error.

**The eigendecomposition of a symmetric matrix.** For G of n x n, given
by its lower half, its n eigenvalues from the largest and, for each, the
eigenvector of length 1, so that G v = λ v and the vectors are at right
angles to each other. It is `dsyevd` in LAPACK, the routine numpy's
`eigh` calls, and `self_adjoint_eigen` of faer. Both give the eigenvalues
from the smallest and the crate turns them round. The sign of an
eigenvector is whatever the backend gave, and a caller that needs a
fixed sign fixes it, as the PCA does.

**The product of two matrices.** For A of r x k and B of k x c, `C = AB`,
r x c. `dgemm` in BLAS and `matmul` of faer. The PCA uses it for the
weights of the variants and for the projections of a table with more
rows than columns.

**The product of a matrix with the transpose of another.** For A of r x
k and B of c x k, both with one row for each of the r and the c things
and one column for each of the k things they are described by, `C = AB'`,
r x c, whose entry i, j is the sum over the k columns of row i of A times
row j of B. It is not a function of its own: it is the same product as
above, given a second operand that says it holds one row for each column
of the result. The routine is the same `dgemm` and the same `matmul`,
told that the second operand is to be read the other way round, which
both libraries do inside the routine and neither pays a copy for,
measured at 1.174 ms against 1.179 ms for 512 x 1000 on Accelerate on
one thread and within 0.3 % on faer, on 23 September 2026.

Without it a caller whose two matrices are both laid out with one row
for each thing has to write the transpose of one of them into a buffer
of its own, which is a matrix operation in a crate that is not this one
and a copy of the whole matrix, 4.1 MB for 512 x 1000. What that copy
costs was timed four ways on the machine of "Speed" on 23 September 2026
and came out between 0.14 and 0.47 ms, a spread too wide to quote a
figure from; what was measured end to end is that the r² of one pair of
512 variants of 1000 individuals went from 8.25 ms to 7.40 ms when the
three copies it made were dropped for this operation, and from 5.97 ms
to 5.02 ms for a set of variants against itself. The r² of
`docs/specs/ld.md` is the caller: its three matrices hold one row for
each variant and one column for each individual, and the product of one
set of variants with another sums over the individuals.

### Layout, half and the backends

Every matrix crosses the interface row after row, and the Fortran
routines of BLAS and LAPACK read a matrix column after column. The
buffer of an r x c matrix read that way is its transpose, c x r, so the
BLAS backend calls each routine on the transposes: the lower half of G
in popnei's layout is the upper half for the routine, `uplo` is `U`;
`A'A` is `A A'` of the transposed view, `trans` `N`; and `C = AB` is `C'
= B'A'`, so `dgemm` gets the buffer of B as its first operand and that
of A as its second, both with `trans` `N`; "The product with its first
operand turned" has the two calls where the first operand is turned as
well. The faer backend tells faer that the buffers are row major and
calls its functions as written. A test of each product on matrices that
are not square, "How it is verified", is what catches a backend that
mixed the two.

### Errors

Each of these is an error of the crate before any routine runs, and the
core crate turns it into its own error and Python into a `RuntimeError`,
since they are defects of the caller: a dimension that does not match,
a `g` that is not c x c for the product nor n x n for the
eigendecomposition, or an `a`, a `b` or a `c` shorter than its rows times
its columns, a longer one being taken by its first rows times columns
values; a c or an n of 0, and in a product an `inner` of 0 as
well, while the rows of a product may be 0; a dimension, or a number of
values of a matrix, above 2147483647, which is what the routines of
BLAS and LAPACK count in, checked for both backends so that the two
refuse the same calls; and a value that is not finite in a matrix. The
last is checked here because the backends do not agree on it: `dsyevd`
on a matrix with a NaN gives NaN eigenvalues and an `info` of 0, measured
through numpy 2.5 on 22 September 2026, and faer's `self_adjoint_eigen`
gives its error of no convergence. A routine of LAPACK that stops,
`info` other than 0, which for `dsyevd` is an eigendecomposition that
did not converge when the `info` is positive and an argument the routine
refused, a defect of popnei, when it is negative, is an error with the
routine and the `info`, and faer's error for the same case is the same
error of the crate, with the `info` at 0, since faer gives none.

What none of the seven checks is its own result. A matrix that is positive
definite and nearly not factors, and the inverse or the solution that comes
off it can hold an infinity: measured on 23 September 2026, the 2 x 2 with
1e-320 and 1 on its diagonal factors on both backends and its inverse is an
infinity and a 1, with no error. numpy does the same and says nothing. That
is the case "Why a Cholesky where numpy uses an LU" calls a fit that is
running away, which pyNei reaches by another road and gives NaN for, so the
caller is what notices it and `docs/specs/gwas.md` is where that is decided.
The crate refuses what it is given, not what it produced.

One error is not a defect of the caller: the workspace of the
eigendecomposition is 2n² floats and more, 1.6 GB at n = 10000, and a
machine that has not the memory for it would abort the process where the
crate asks for it. So it is asked for with `try_reserve_exact`, and a
machine that refuses gives an error that says what could not be
allocated and how many values it was. The core crate wraps it, and which
exception it becomes is for the spec of the module that calls it, since
it is neither a wrong argument nor a defect of popnei. No test reaches
it: the case is written and read, not run.

### How they run

`dsyevd` overwrites the matrix it is given with the eigenvectors, so
`eigh_lower` takes `g` by value and gives its buffer back as the
eigenvectors, with no copy; a caller that needs G afterwards copies it
first. With eigenvectors the routine needs a workspace of 1 + 6n + 2n²
floats and 3 + 5n integers besides, which the crate asks the routine for,
as LAPACK is asked, with the two lengths at -1, and allocates at each
call: 16 MB at n = 1000 and 1.6 GB at n = 10000, twice the matrix. The
routine writes the number of floats as an `f64`, so a value that is not a
count, a NaN, an infinity, a negative number or one above 2147483647, is
left out and the minimum above stands, and a value that is a count is
taken as its whole part. faer allocates what it needs inside the call.
The crate holds nothing between calls.

Both backends give the eigenvalues from the smallest, LAPACK and faer
alike, and both give each eigenvector where popnei's row major buffer
reads it as a row. So the turning round, the values and the rows of the
vectors together, is done once, above the backends, and neither of them
does it.

### How it is verified

The reference outside the project is numpy 2.5 on Accelerate, whose
`eigh` is `dsyevd` of the newer LAPACK of the framework; the examples
below were run with it on 22 September 2026 and their numbers are the
literals of the first cargo tests. The tests run twice, `cargo test -p
popnei-linalg` and the same with `--no-default-features`, and both have
to pass. Each check is made at the function of "The Rust interface"
that it names.

At `add_self_product_lower`: A of 2 x 3, rows (1, 2, 0) and (0, 1, 3).
A'A is

    1  2  0
    2  5  3
    0  3  9

and the test asserts the lower half exactly, since every entry is a sum
of at most two products of small integers, and that the upper half of G
kept the values it was given; and that an A of no rows leaves G as it
was. At `product`, three cases, exactly, all with the same A. Times B of
3 x 2 with rows (1, 0), (2, 1) and (0, 3) it is the 2 x 2 matrix with
rows (5, 2) and (2, 10). That one is symmetric, so a backend that wrote
the transpose of C would pass it; times B of 3 x 2 with rows (1, 1),
(2, 0) and (0, 3) it is the 2 x 2 matrix with rows (5, 1) and (2, 9),
which is not. And times B of 3 x 1 with rows (1), (0) and (2) it is the
2 x 1 matrix with rows (1) and (6), where the three dimensions are
different, so a backend that swapped two of them would pass neither this
nor the self product of the A above. C holds values other than 0 before
each call, which an operation that added to C instead of overwriting it
would leave in the result.

At `product` with its second operand read by the columns of the result,
the same A of 2 x 3 and three cases, exactly.
Times the transpose of B of 2 x 3 with rows (1, 1, 0) and (0, 2, 1) it
is the 2 x 2 matrix with rows (3, 4) and (1, 5), which is not symmetric,
so a backend that wrote the transpose of C would fail it. Times the
transpose of B of 1 x 3 with the row (2, 0, 1) it is the 2 x 1 matrix
with rows (2) and (3), where the three dimensions differ. And with A
given for both operands it is A A', the 2 x 2 matrix with rows (5, 2)
and (2, 10), which is the same matrix as the first case of `product`
above, since the B of 3 x 2 there, with rows (1, 0), (2, 1) and (0, 3),
is this A written the other way round: the two functions are asserted to
give it alike, which is what catches one of them reading an operand the
way the other does.

At `eigh_lower`: the 3 x 3 matrix with rows (4, 1, 0), (1, 3, 0) and
(0, 0, 1). Its eigenvalues are (7 + √5)/2 = 4.618033988749895,
(7 - √5)/2 = 2.381966011250105 and 1, and the eigenvectors, each with
the sign that makes its entry of largest absolute value positive, are
(0.8506508083520399, 0.5257311121191335, 0), (-0.5257311121191335,
0.8506508083520399, 0) and (0, 0, 1); the first two are (1, λ - 4)
divided by their length. Every number here and below is written with the
digits that name the `f64` numpy computed and no more, which is what
Python's `repr` prints and what Rust reads back as the same value, so
that the tolerance of a test is spent on the code and not on the
rounding of a literal. The test gives each vector that sign and compares
the eigenvalues and the vectors within 1e-12.

At `eigh_lower`, on a matrix too large to write down, so that the two
backends are checked on the size they will run at: G = ZZ' for Z of 1000
rows and 1200 columns, with z(i, c) the (c · 1000 + i)-th number of the
xorshift generator below, started at 7. G then has full rank and its
eigenvalues are apart, between 0.79 and 361.9, no two closer than 0.003.
The literals, from numpy 2.5.3 on 22 September 2026: the trace
99996.3873081677, which is the sum of the eigenvalues, and which numpy
adds up to 99996.38730816769 when it adds the eigenvalues instead; the
three largest eigenvalues 361.9125119011332, 359.66517178659313 and
356.4439329949563, and the smallest 0.7933289215408484; and the first
three entries of the eigenvector of the largest, with the sign of its
largest entry positive, 0.018111301861995926, -0.004576022169100455 and
-0.00518748985091827. The test compares the eigenvalues within 1e-12
relative and the entries within 1e-9, because an eigenvector is less
well determined than its eigenvalue by the gap to its neighbours, and a
matrix of a real dataset has closer ones than this. Measured on 22
September 2026 with the trial of "Speed": the two backends agree to
2.9e-14 relative on the eigenvalues and to 1.3e-12 on the entries of the
eigenvectors after the sign, and their products of the same Z agree to
the bit. A matrix of rank below n is not such a test: with Z of 1000 x
200 the 800 eigenvalues that are 0 come out between -8e-14 and 8e-14 and
their eigenvectors are any base of that space, different in each
backend.

The generator, so that the test and numpy make the same Z: a 64 bit
state s that starts at 7 with its lowest bit set, and for each number
`s ^= s << 13; s ^= s >> 7; s ^= s << 17`, the shifts dropping the bits
that leave the 64, and the number is `(s >> 11) / 2^53 - 0.5`.

The checks at the level above, that a PCA made from these operations
gives the numbers of R, are in `docs/specs/pca.md`.

## The seven operations of the GWAS

The genome wide association study fits one model of the trait to the
covariates alone, the null model, and then tests every variant against
what that model left unexplained. pyNei fits four null models, one for
each pairing of a continuous or a binomial trait with and without a
kinship, and `src/pynei/gwas.py` at commit ef0ca6e calls five functions
of `numpy.linalg` that this crate has not, at fifteen places: `qr` at
376, `solve` at 378, 412, 483, 512, 556, 639, 646, 692 and 693,
`slogdet` at 418, `inv` at 586, 666 and 689, and `matrix_rank` at 841.
This section says which operations those become and what each one gives. It
says nothing about which of the four models calls what, nor about the
order in which a model calls them, which is for `docs/specs/gwas.md`.

What the GWAS works on, before the list. The **design**, `d` below, is
the matrix of the covariates with one row for each individual and one
column for each number the model fits, its **coefficients**: a column of
ones for the intercept and one column for each covariate the user gave.
For the datasets of `docs/objectives.md` it has up to 10000 rows and a
handful of columns. A **right hand side** is one of the vectors `b` of a
system `a x = b`, and one matrix `a` is often solved against many of them
at once. A symmetric matrix `a` is **positive definite** when `v' a v` is
above 0 for every vector `v` that is not 0, which for the matrices here
means that the columns of the design are independent; every such matrix,
and no other, has a **Cholesky factorization**, the lower triangular `l`
with `l l' = a`.

### What the GWAS calls of numpy, and where

The lines are of `src/pynei/gwas.py` at commit ef0ca6e, and n is the
individuals.

| what pyNei calls | at | the matrix, and how large |
|---|---|---|
| `solve` | 556 | the design multiplied by itself with one weight for each individual, `d' w d`, coefficients square; one right hand side |
| `solve` | 412, 512 | the same, of the design turned by the eigenvectors of the kinship, `u' w u` for `u = e' d` with `e` those eigenvectors; one right hand side |
| `solve` | 483 | `d' V⁻¹ d`, where `V` is n x n and dense, the variance of the trait under the mixed model; one right hand side for each individual |
| `solve` | 692, 693 | `d' Σ⁻¹ d`, where `Σ` is the n x n of line 689; one right hand side at 692 and one for each individual at 693 |
| `solve` | 639, 646 | one matrix for each variant of a block, of the coefficients plus one, square |
| `solve` | 378 | the `r` of the QR of the design, coefficients square and upper triangular; one right hand side |
| `slogdet` | 418 | the matrix that line 412 solves against |
| `inv` | 586 | `d' w d`, coefficients square |
| `inv` | 666 | one matrix for each variant, as at 639, of which one entry of each inverse is read |
| `inv` | 689 | n x n, a number not below 0 times the kinship, plus a positive number on each diagonal entry |
| `qr` | 376 | the design, n x coefficients |
| `matrix_rank` | 841 | the design |

Of the square matrices in that table, all but two are symmetric positive
definite by construction. `d' w d` and `u' w u` have a positive weight for
each individual, so they are positive definite whenever the columns of `d`
are independent, which line 841 refuses a design without; and `d' V⁻¹ d`
and `d' Σ⁻¹ d` are positive definite for the same reason whenever `V` and
`Σ` are.

The two that are not. The `r` of line 378 is upper triangular and not
symmetric at all, and it gets a solve of its own below. And `Σ` itself,
the matrix of line 689, is a number that is not below 0 times the kinship
plus a positive number on each diagonal entry, which is positive definite
only while the kinship never makes `v' k v` negative. A kinship of
genotypes none of which is missing never does, being a matrix multiplied
by itself. pyNei's does, because `_KinshipCalc` of `gwas.py` divides that
product entry by entry by the variants called in both individuals of each
pair, at lines 208 and 228, and an entry by entry division does not keep a
matrix positive semidefinite. Simulated on 23 September 2026 on 200
individuals and 2000 variants, with allele frequencies drawn between 0.05
and 0.5 and genotypes missing at random: the smallest eigenvalue of that
kinship is 0 with nothing missing, -0.02 at 2 per cent missing, -0.11 at
10 per cent and -1.00 at 50 per cent, against a largest of 1.7 to 2.4.

`Σ` is that eigenvalue times the variance component the fit is searching
over, plus one over the weight of each individual on the diagonal, and a
weight is at most 0.25, so the diagonal adds at least 4: `Σ` goes
indefinite when the variance component passes about 4 divided by the size
of that negative eigenvalue, which is about 200 at 2 per cent missing and
4 at 50 per cent. pyNei never meets this, because its LU factors an
indefinite matrix without complaining and gives an inverse that is an
inverse. A Cholesky refuses, so popnei stops where pyNei gave a fit. pyNei
clamps the kinship's negative eigenvalues at 0 at line 506, calling them
rounding, but only in the continuous mixed model, which eigendecomposes
the kinship anyway; the logistic one never does.

The owner decided on 23 September 2026 to leave the kinship as pyNei
computes it and let the fit fail with `Singular`, which tells the user
that the kinship is not one, and to reconsider if the implementation of
the GWAS meets it. The options not taken were to clamp the kinship's
eigenvalues at 0 before the logistic mixed model uses it, as line 506 does
for the continuous one, which costs an eigendecomposition of n x n per
fit, 0.035 s at 1000 individuals and 6.3 s at 5000 from "Speed", and
changes the fit for every dataset with missing genotypes in a direction
nothing here has measured; and to give the kinship a per pair denominator
that keeps it positive semidefinite, which is `docs/specs/kinship.md`'s
decision and not this one. What is not known is how often a real dataset
reaches it: no run of pyNei's logistic mixed model on missing genotypes
was made for this spec.

The sizes are of two kinds. Everything but lines 689, 483 and 693 is a
matrix of the size of the coefficients, a handful square, and the cost of
a call is the cost of the call and not of the arithmetic. Lines 689, 483
and 693 grow with the individuals, and they are paid while a null model is
fitted and not once for each variant.

### What the seven give

Each entry says what the operation gives and then names the routine of
LAPACK and the entry point of faer 0.24.4 that runs it, which the
implementer needs and a reader deciding whether these are the right
operations can pass over. Every one of the named routines and entry
points was called on both backends and on both wasm targets before this
was written.

**The Cholesky factorization.** For the symmetric positive definite `a`
of n x n, read by its lower half, the lower triangular `l` with `l l' =
a`, which overwrites the lower half of `a`; the upper half is left as it
was. It is `dpotrf` in LAPACK and `cholesky_in_place` of faer's
`linalg::cholesky::llt::factor`, which is given a scratch of n values
and the regularization faer calls its default. That default is the one
that refuses: told to regularize, faer patches a pivot that is not
positive and factors on, and `Singular` would never be raised. A matrix
that is not positive definite gives the error `Singular` of "The errors
the seven add", and this factorization is the test for it: it stops at
the first row whose diagonal entry, once the rows above it have been
taken out, is not above 0.

**The solve of a factorized matrix.** For the `l` above, the `x` of `a x
= b`. `b` is `sides` x n, one row for each right hand side, and it comes
back holding the solutions the same way. It is `dpotrs` in LAPACK and
`solve_in_place_with_conj` of faer's `linalg::cholesky::llt::solve`.

That layout, one row for each right hand side, is what the caller of
lines 483 and 693 already has, and it was chosen for it. There `a` is the
coefficients square, so each right hand side holds one number for each
coefficient, and there is one of them for each individual: the matrix of
them is the individuals by the coefficients, which is what those two lines
hold, and numpy is handed its transpose because numpy takes the right hand
sides as columns. popnei takes them as rows, so the caller passes the
buffer it has with `sides` the individuals and `n` the coefficients, and
nothing is copied and nothing is moved. Checked against numpy 2.5.3 on 23
September 2026 on 3 individuals and 2 coefficients: `solve` of a 2 x 2
against a 2 x 3 gives a 2 x 3, which is 3 right hand sides of 2 numbers
each. The solutions come back the same way, one row for each column of the
matrix the next product needs, so that product is `product` with its
second operand `ByTheColumnsOfTheResult`.

**The log of the determinant.** For the `l` above, twice the sum of the
logs of its diagonal, which is the log of the determinant of `a`. It is
arithmetic over n numbers and neither library is called for it. Line 418
takes the second value of numpy's `slogdet` and throws the first, the
sign, away; the determinant of a positive definite matrix is above 0, so
that sign is always 1 and nothing is lost.

**The inverse of a factorized matrix.** For the `l` above, the lower
half of the inverse of `a`, written into a buffer of n x n that the
caller gives; its upper half is left as it was. It is `dpotri` in
LAPACK, which the crate calls on a copy of `l` in that buffer and which
needs no workspace, and `inverse` of faer's
`linalg::cholesky::llt::inverse`, which writes into the buffer and asks
for a scratch of n x n of its own, 8 bytes to a value, which faer 0.24.4
says and which it was asked for at n = 1000, where it wanted 8000000
bytes. So an inverse at 10000 individuals needs the two buffers the
caller gives, 800 MB each, and in faer 800 MB more, which is
`docs/rust_core.md`'s open question about the memory in the browser,
where a kinship of 10000 individuals is already 800 MB. The crate asks
for faer's scratch with `try_new` of `dyn_stack`, so that a machine
without the memory gets `Memory` and not the end of the process.

**The thin QR of the design.** For `a` of `rows` x `cols` with `rows` at
least `cols`, the `q` of `rows` x `cols` whose columns are of length 1
and at right angles to each other, and the upper triangular `r` of
`cols` x `cols`, with `a = q r`. It is `dgeqrf` and then `dorgqr` in
LAPACK, the two routines numpy's `qr` calls, and `qr` of faer with
`compute_thin_Q` and `thin_R`. The BLAS backend gives the two routines a
column major copy of `a`, which is its transpose written out, because in
their view the buffer as it lies is the wide matrix and they are much
slower on it than on the tall one; "Speed" has the two numbers. The GWAS
keeps both: `r` for the coefficients of the plain linear model, and `q`
for every variant of every block, whose dosages line 391 takes the
covariates out of with `q q'`.

The sign of a column of `q`, and of the row of `r` that goes with it, is
whatever the backend gave, as it is for an eigenvector. Nothing downstream
sees it: line 391 and the residuals of line 379 use `q` only as `q q'`,
and the coefficients of line 378 are the `c` of `r c = q' y`, and
flipping the sign of column j of `q` flips row j of `r` and entry j of
`q' y` with it, which leaves `c` alone. Measured on 23 September 2026:
LAPACK, faer and numpy 2.5.3 all gave the negative diagonal of `r` for the
design of "How the seven are verified" below, so no case has yet been
found that makes them differ, and the test still fixes the sign so that a
backend which chose the other one would pass.

**The solve against an upper triangular matrix.** For the upper triangular
`r` of n x n, whose lower half is not read, the `x` of `r x = b`, with
`b` laid out as for the Cholesky solve. It is `dtrtrs` in LAPACK and
`solve_upper_triangular_in_place` of faer. It exists for line 378 alone.

**The rank.** For `a` of `rows` x `cols`, how many of its singular
values are strictly above the tolerance `s * max(rows, cols) *
2.220446049250313e-16`, where `s` is the largest singular value and the
last number is the distance from 1 to the next `f64` above it. That is
numpy's tolerance, read from `matrix_rank` of numpy 2.5.3, and popnei
takes it so that a design popnei refuses is a design pyNei refuses. The
singular values of a matrix are the factors by which it stretches space
along as many directions at right angles to each other as it has
columns, from the largest, and one of them is 0 exactly when a column of
the matrix is a combination of the others, so the rank is the columns
that are independent. It is `dgesdd` with `jobz` `N`, which computes the
values and no vectors, in LAPACK, and `singular_values` of faer, and the
BLAS backend copies the transpose for it too.

The options not taken were a QR with column pivoting, which takes the
columns of the design in the order of how much of them is left once the
columns before have been taken out and so brings the nearly dependent
ones last, and whose diagonal of `r` is then a cheaper stand-in for the
singular values but is not them, so that its rank at the same tolerance
can differ from numpy's on a design whose columns are nearly dependent; and a
test the caller makes, trying the Cholesky of the design multiplied by
itself, which is cheaper and refuses a band of designs that pyNei accepts,
measured under "Why a Cholesky where numpy uses an LU". Neither is worth
taking for its speed: the rank of
a design of 10000 x 5 took 0.145 ms on Accelerate and 0.159 ms on faer,
once for a whole study.

### Why a Cholesky where numpy uses an LU

numpy's `solve`, `inv` and `slogdet` all factor a general square matrix
into a lower and an upper triangular one with row exchanges, an LU, which
works for any matrix that can be solved against. Every argument in the
table above but the `r` of line 378 is symmetric positive definite, and
for those a Cholesky is the factorization to make: it is half the
arithmetic of an LU, it gives the log of the determinant of line 418 from
the same factorization the solve of line 412 uses, and, unlike an LU, it
fails on exactly the matrices that are not positive definite instead of
carrying on with a tiny pivot. It is also 2.1 times faster than numpy's
`inv` at 5000 individuals, which "What the seven of the GWAS cost" has
the two numbers of.

The owner decided on 23 September 2026 that the factorization is a
Cholesky. The option not taken was an LU, `dgesv` and faer's
`partial_piv_lu`, which matches numpy input for input at twice the
arithmetic, gives no log determinant off the same factorization, and is
no test for a design whose columns are not independent.

What it changes for a user is which inputs are refused, and there are two
of them. An LU carries on with a pivot that rounding made tiny and gives a
number; a Cholesky stops, so popnei refuses some inputs pyNei answers for.

The first is a design whose covariates are independent on paper and nearly
not in the numbers, which the rank of line 841 lets through and the
Cholesky then refuses. The band was measured on 23 September 2026 with
numpy 2.5.3 on a design of 10000 rows and 4 columns built with its
smallest singular value at a chosen fraction of its largest:
`matrix_rank` gave 4 down to a fraction of 1e-11 and 3 at 1e-12, which is
where its tolerance of `max(rows, cols)` times 2.2e-16 lies, while the
Cholesky of that design multiplied by itself succeeded down to 1e-8 and
refused at 1e-9, the fractions in between going either way as rounding
decided. So the band is about three orders wide, between 1e-9 and 1e-12.
A design that falls in it has its largest singular value between 1e9 and
1e12 times its smallest, and the relative error of a solve is about that
ratio times 2.2e-16, so the coefficients numpy gives for such a design
keep roughly between seven and four of the sixteen digits an `f64` holds.

The second is a variant of the logistic Wald test, the test that fits the
model again with the variant in it and asks how many of its own
uncertainties the variant's effect is away from 0, whose fit is running
away and whose matrix goes flat as the weights of the individuals go to 0.
pyNei reaches those variants by another road and gives them the same
answer: lines 649 to 651 mark a variant whose step is not finite and line
658 one whose coefficient passes 30, and line 670 gives both NaN.

What each refusal becomes where it is caught is for `docs/specs/gwas.md`.

### The systems of one block, and what the fallback becomes

Line 639 solves a stack of systems at once, one for each variant of a
block that is still iterating, each of the size of the coefficients plus
one; line 666 inverts the same stack and reads one entry of each inverse,
the last of the diagonal. numpy solves a stack by looping over it, and no
routine of LAPACK and no entry point of faer takes a stack, so an
operation of this crate for it would be a loop over the two operations
above and nothing more. It is a loop in the caller: the crate gives the
factorization and the solve of one matrix, and the core crate runs them
over the variants of the block. Measured on 23 September 2026: one system
of 7 x 7 factored and solved takes 0.173 µs on Accelerate and 0.159 µs on
faer, so 5000 of them, a block, take 0.9 ms and 0.8 ms, for one iteration
of the fit.

The fallback of lines 642 to 648 goes away. It is there because numpy's
stacked `solve` raises on the whole stack when any one matrix of it is
singular, so pyNei catches that and solves the stack again one matrix at a
time, keeping the answers of the good ones and leaving NaN for the bad.
Checked on 23 September 2026 with numpy 2.5.3 on a stack of two matrices
of which one was singular: the stacked call raised `LinAlgError` and gave
nothing, and the two calls one at a time gave the answer for the good
matrix and raised for the other. A loop over a function that gives a
`Result` has nothing to fall back from: the variants whose system is
refused are the ones the loop gets an error for, and the caller gives them
the step of 0 and the mark of having diverged that lines 649 to 651 give
them, which line 670 turns into the NaN a user sees.

The last entry of the diagonal of the inverse, line 666, needs no
inverse: for `v` the vector that is 1 in its last entry and 0 everywhere
else, that entry is `v'` times the inverse of `a` times `v`, which is the
solution of `a x = v` read at its last entry. So it is one solve, or it is
the inverse and one read of it. Which of the two the GWAS takes is for its
spec; the operations here give either.

### The errors the seven add

Every one of the seven adds its own cases of `Dimension`, a buffer that
does not hold the values of the dimensions it was given or a dimension of
0, and of `NotFinite`, checked over the values the operation reads and
over the lower half alone where only the lower half is read. The rank adds
`NoConvergence`: `dgesdd` gives an `info` other than 0 and faer gives
`SvdError::NoConvergence`, which is the case the crate already has for the
eigendecomposition, with the routine and the `info`. The inverse adds
`Memory`, for faer's scratch of n x n. The log of the determinant reads
the diagonal of `l` alone, so `NotFinite` and the `Singular` below are
what it checks that diagonal for, and it is the same `Singular` a
`cholesky_lower` that gave that `l` would have given first.

The solve and the inverse read that diagonal for the same `Singular`, and
the inverse for the reason the triangular solve below does: the two
backends do not agree on an `l` whose diagonal holds an entry that is not
above 0. Measured on 23
September 2026 on an `l` with a 0 at its row 1, `dpotri` gave an `info` of
2 and faer's `inverse` gave no error at all and wrote infinities and NaN
into the buffer. No `l` that `cholesky_lower` gave is such a matrix, since
that is what it stops at, and a caller holds the two buffers apart and can
pass one that never was a factorization. So the diagonal is read in the
crate, above the backends, where it holds for both, and the caller gets
the error instead of a matrix of NaN.

The solve is refused for a different reason, and it is the one that makes
the rule general: there the two backends agree, and both are wrong. On the
same `l`, measured on 23 September 2026, `dpotrs` gave a solution of NaN,
an infinity and an infinity with the sign turned round, and faer's solve
gave three NaN, each with no error at all. numpy refuses the same system,
`LinAlgError: Singular matrix`. So the three operations that read an `l`
read its diagonal first, and each gives the `Singular` that the
`cholesky_lower` which would have produced that `l` gives at the same row.
What is left to the caller is an `l` that is the factorization of some
other matrix, which no check here can see: "The Rust interface" says so at
each of the three.

One case is new. A matrix that cannot be factored is neither a wrong
dimension nor a value that is not finite nor a routine that ran out of
iterations, and it is not a defect of the caller either: it is what the
data was.

```rust
/// A matrix that could not be factored at the row the value names, ///
counting from 0: the Cholesky reached a diagonal entry that is not ///
above 0 there, or the solve against an upper triangular matrix ///
reached one that is 0. Singular { argument: &'static str, at: usize },
```

Its message is "the matrix a is singular: the factorization stopped at its
row 3, counting from 0". The core crate wraps it and the module that
called it decides what it means, because the same error is two different
things to a user: a design whose covariates are not independent, which the
rank of line 841 is meant to catch before any model is fitted, and a
variant whose fit has run away, which the GWAS gives NaN for.

Both backends give the row, and both count it from 0 once the crate has
turned LAPACK's round: `dpotrf` gives an `info` of k for the leading k x k
corner, counting from 1, and faer gives `LltError::NonPositivePivot` with
an index from 0. Checked on 23 September 2026 on the 3 x 3 with rows (4,
2, 0), (2, 1, 0) and (0, 0, 1), whose leading 2 x 2 has a determinant of
0: both gave the row 1. numpy's `cholesky` raises `LinAlgError` on the
same matrix, and its eigenvalues are 0, 1 and 5: none below 0 and one
exactly 0, which is the case a factorization has to catch and an LU does
not.

faer's `solve_upper_triangular_in_place` gives no error and divides by a
diagonal entry of 0 as it finds it, giving an infinity, while `dtrtrs`
gives an `info`. So the diagonal of `r` is read for a 0 in the crate,
above the backends, where it holds for both, as the check for a value that
is not finite already is.

### How the seven are verified

The reference outside the project is numpy 2.5.3 on Accelerate, the same
one the operations above use. Every number below was taken from it on
23 September 2026 and every one was then run on both backends, `cargo test
-p popnei-linalg` and the same with `--no-default-features`, in a trial
crate that is not in git. Each check is made at the function of "The Rust
interface" that it names. That a GWAS made from these operations gives
the numbers of GMMAT and of rrBLUP, the two programs `docs/objectives.md`
names as the reference outside the project for an association study, is
for `docs/specs/gwas.md`, as the checks of the PCA are for its own spec:
no function of this crate is reached from Python.

The small matrix the Cholesky and what comes off it are checked on is the
3 x 3 with rows (4, 2, 0), (2, 10, 6) and (0, 6, 5), whose factorization is
the lower triangular with rows (2, 0, 0), (1, 3, 0) and (0, 2, 1), every
entry a small integer, and whose determinant is 36.

At `cholesky_lower`, that factorization exactly. At
`log_determinant_with_cholesky`, 3.58351893845611, which is the log of 36
and the number numpy's `slogdet` gives for the same matrix, within 1e-15
relative. It is not asserted to the bit, although the sum of the logs of 2,
3 and 1, doubled, does land on that `f64` here: `ln` is not rounded the
same on every platform, which is why the `coding` skill does not let a test
assert the bits of a value that went through it, and this sum is one bit
wide. Measured on 23 September 2026, five of the nine ways of moving `ln 2`
and `ln 3` by one unit in the last place give another `f64`, and the libm
Rust uses for `wasm32-unknown-unknown` already differs from this machine's
`ln 3` by one of those units. One unit in the last place of the answer is
1.2e-16 relative, so the tolerance leaves about eight of them. At
`solve_with_cholesky`, the right hand side (8, 40, 27) gives (1, 2, 3),
and the two right hand sides (8, 40, 27) and (4, 2, 0), one row each, give
(1, 2, 3) and (1, 0, 0), both within 1e-14, which is what catches a
backend that read the rows of `b` as its columns. At
`invert_with_cholesky`, the lower half is 7/18, -5/18, 5/9, 1/3, -2/3 and
1, within 1e-15 relative. Those are the exact values, and the literals are
they and not numpy's, because numpy's `inv` goes through an LU and gives
0.38888888888888884 where 7/18 is 0.3888888888888889 and
0.9999999999999998 for the entry that is 1, up to 2 units in the last
place away from the number a Cholesky gives.

At `cholesky_lower` again, the 3 x 3 with rows (4, 2, 0), (2, 1, 0) and
(0, 0, 1) is `Singular` at the row 1, as above.

Then a matrix too large to write down, so that the two backends are
checked at a size they will run at, and it is the G of the
eigendecomposition above: `G = ZZ'` for Z of 1000 rows and 1200 columns
built with the xorshift generator of "How it is verified" of "The four
operations of the PCA and the LD", started at 7, whose trace is
99996.3873081677. Its eigenvalues are between 0.79 and 361.9, so it is
positive definite. The literals, from numpy 2.5.3: the first and the
last entries of the diagonal of its factorization, 10.135944716832457
and 4.115426436421405, within 1e-12 relative; its log determinant,
3963.7986384485084, within 1e-13 relative, which numpy's `slogdet` gives as
3963.7986384485057 through its LU, the tolerance being relative and not
absolute since the two backends are 2.4e-12 absolute away from it and 6.1e-16
relative; the solution of `G x = v` for `v` the
vector of 1000 ones, of which the test asserts the first three entries,
-0.3054837936349659, -0.04576083734778211 and -0.21314634692119025, and
their sum over the 1000, 73.9565335781636, within 1e-11 relative; and of
its inverse the first and last entries of the diagonal,
0.06230734831937398 and 0.059043258015696806, and the trace,
59.78707893196584, within 1e-11. What the two backends were measured to
be away from those numbers on 23 September 2026 is 5.7e-16 and 8.0e-16
relative on the log determinant, 1.8e-14 and 1.1e-14 on the four numbers
of the solve, and 5.4e-15 and 4.9e-15 on the three of the inverse,
Accelerate first and faer second, so each tolerance above has at least
two orders of room.

The design the QR, the triangular solve and the rank are checked on is
the 4 x 2 of an intercept and one covariate, with rows (1, 1), (1, 2),
(1, 3) and (1, 4). At `thin_qr`, with the sign of each column of `q`
taken so that the diagonal of `r` is positive, `r` has rows (2, 5) and
(0, 2.23606797749979), which is the square root of 5, and `q` has the
first column (0.5, 0.5, 0.5, 0.5) and the second (-0.6708203932499368,
-0.22360679774997894, 0.223606797749979, 0.6708203932499369), which is
the covariate less its mean, divided by the length of that; both within
1e-14, because Accelerate gave -2.0 for the first entry of `r` and faer
-1.9999999999999998. At `solve_upper_triangular`, the trait (1, 3, 5,
7), which is twice the covariate less 1, gives the coefficients (-1, 2)
from `r c = q' y`, within 1e-14: an exact fit, so a backend that read
`r` the wrong way round gives something else. The `q' y` of that trait
is (8, 4.47213595499958), which is the right hand side the check passes.
Two more traits go with it: (4, 7, 10, 13), three times the covariate
plus 1, whose `q' y` is (17, 6.70820393249937) and whose coefficients
are (1, 3); and the covariate itself, (1, 2, 3, 4), whose `q' y` is
(5, 2.23606797749979) and whose coefficients are (0, 1). The three
passed together, one row each, are what catches a backend that read the
rows of `b` as its columns, as the two right hand sides of the Cholesky
solve above are, and three of them against an `r` of 2 x 2 is what tells
`sides` from `n`. All three fits are exact, and numpy 2.5.3 gives
0.9999999999999991 and 3.0000000000000004 for the second of them, 9e-16
relative away, which the 1e-14 holds. And at
`solve_upper_triangular` again, the `r` with rows (2, 5) and (0, 0) is
`Singular` at the row 1, which is the case the crate reads the diagonal
for, since faer would divide by that 0 and give an infinity.

At `rank`, five matrices and their ranks, which is all `rank` gives: it
returns a count, so its singular values cannot be asserted at it, and the
five are chosen so that the count alone pins both the values and the
tolerance. The 4 x 2 above is rank 2. The 4 x 3 with rows (1, 1, 2),
(1, 2, 3), (1, 3, 4) and (1, 4, 5), whose third column is the sum of the
first two, is rank 2. The 4 x 2 whose covariate is the constant 5 is rank
1. Those last two are the two designs line 841 refuses, the collinear one
and the constant one.

The other two pin the tolerance itself, since a rank that counted at the
wrong threshold would still give 2, 2 and 1 for the three above: the 2 x 2
with rows (1, 0) and (0, 5e-16) is rank 2 and the one with rows (1, 0) and
(0, 4e-16) is rank 1, either side of the tolerance for a 2 x 2 whose
largest singular value is 1, which is 4.440892098500626e-16. numpy 2.5.3
gives 2 and 1 for them, and both backends were run on them on 23 September
2026 and gave the same two counts and the same tolerance.

The singular values themselves, from numpy 2.5.3, are 5.779378813233887
and 0.773809106397227 for the first; 9.344132686098556,
0.8289658283575813 and 3.651382431893325e-17 against a tolerance of
8.299257002607302e-15 for the second; and 10.19803902718557 and 0 for the
third. They are here so that a reader can see where each count comes from,
and no test asserts them. The last two are the
two designs line 841 refuses, the collinear one and the constant one.

## The product with its first operand turned

pyNei's GWAS writes `a.T @ b`, a matrix multiplied by another that has
one row for each of the same individuals, at fourteen places of
`gwas.py`, and the `product` of "The four operations of the PCA and the
LD" can do none of them: it reads its first operand as one row for each
row of the result, and there is no way to ask it for the other. Four of
the fourteen are the design multiplied by itself with one weight for
each individual, at 411, 511, 554 and 586, which
`add_self_product_lower` does give, since the weights are above 0 and
the caller can multiply each row of the design by the square root of its
weight first. The other ten are `q' y` at 378 and 379, the design by the
weighted trait at 412 and 512, the design by the projected design at 482
and 691, the eigenvectors of the kinship by the trait and by the design
at 507 and 508, the design by the residuals at 555, and the projected
design by the working trait at 692.

The owner left this decision to the writer of this spec on 23 September
2026, saying that changing the signature of `product` is not a cost to
weigh against the tool being the best one popnei can build. So the first
operand gets a type, the way the fourth operation above gave the second
one a type for the r² of `docs/specs/ld.md`. The options not taken were
one function for `c = a' b` beside `product`, which brings back what the
typed second operand was for, two functions over the same four numbers
either of which takes the other's call and gives a different matrix with
no error; and leaving it out, so that the caller writes the transpose of
one operand into a buffer of its own, which is a matrix operation in a
crate that is not this one and which costs more than the product it
feeds.

### What it gives

`product` takes a `TheFirstOperand` beside its `TheSecondOperand`, and the
two together say which of four matrices it computes. Each operand says how
its buffer is laid out, never what the routine should do to it, so a
caller states what it has and not what it wants done.

| the first operand | the second | what `product` writes into `c` |
|---|---|---|
| `ByTheRowsOfTheResult` | `ByTheValuesSummedOver` | `c = a b` |
| `ByTheRowsOfTheResult` | `ByTheColumnsOfTheResult` | `c = a b'` |
| `ByTheValuesSummedOver` | `ByTheValuesSummedOver` | `c = a' b` |
| `ByTheValuesSummedOver` | `ByTheColumnsOfTheResult` | `c = a' b'` |

The first two are what `product` computes today and the callers of them do
not change what they compute, only how they say it: the four calls of
`product` in the core crate, three in `pca.rs` and one in `ld.rs`, name
`ByTheRowsOfTheResult` for their first operand and are otherwise as they
were. The third is the ten places above. The fourth has no caller in
popnei and is built anyway, because the two enums make it a call that
compiles, and an interface that takes a call it then refuses at run time is
worse than one that answers it; it costs one `dgemm` with both `trans`
flags set and one `matmul` over two transposed references, and it is
tested like the other three.

`add_self_product_lower` stays and does not become a case of `product`,
although `a' a` is now one: it is `dsyrk`, which writes the lower half
alone and which "What they give" measured at 95 ms against 174 ms for
`dgemm` on a block of 5000 x 1000 with faer on one thread.

### What it costs

Neither backend copies anything for it and neither routine is a different
routine. In the BLAS backend `c = a' b` is the same `dgemm` with `transb`
set to `T`, since in the routine's column major view `c' = b' a`, and
`c = a' b'` is `c' = b a` with both `trans` flags `T`. In faer both are
`matmul` over a transposed `MatRef`, which is another reference over the
same values read the other way round.

Measured on 23 September 2026 on the machine of "Speed", a 1000 x 1000 by
a 1000 x 5: `c = a b` took 0.102 ms and `c = a' b` 0.030 ms on Accelerate,
and 0.069 ms and 0.091 ms on faer, the two giving the same matrix to the
bit when the first operand of one is the transpose of the first operand of
the other. Writing that transpose of the 1000 x 1000 into a buffer
instead, which is what a caller without this would do, took 0.361 ms,
three to twelve times the product it would feed.

### How it is verified

At `product`, on the A of 2 x 3 with rows (1, 2, 0) and (0, 1, 3) and
the B of 3 x 2 with rows (1, 1), (2, 0) and (0, 3), which "How it is
verified" of "The four operations of the PCA and the LD" already uses,
and whose product `a b` is the 2 x 2 with rows (5, 1) and (2, 9). That
matrix is not symmetric, so a backend that wrote the transpose of `c`
fails on it.

The same three numbers come out of all four combinations when each operand
is given as the matrix the combination expects, and the test asserts that,
exactly: `a b` from A and B; `a' b` from the 3 x 2 with rows (1, 0),
(2, 1) and (0, 3), which is A written the other way round, and B; `a b'`
from A and the 2 x 3 with rows (1, 2, 0) and (1, 0, 3), which is B written
the other way round; and `a' b'` from those last two. Run on both backends
on 23 September 2026 and all four gave (5, 1) and (2, 9). A backend that
read one operand the way another combination reads it gives a different
matrix for at least one of the four, which is what the four assertions
together catch.

Those four matrices are all 2 x 2, so the rows and the columns of the
result are the same number and a call that swapped the two would give them
all. So each of the two combinations that turn the first operand is
asserted once more on a result whose three dimensions differ, which is
what "Layout, half and the backends" asks of every product and what the
third case of `product` above is for. Both use the A of 2 x 3 written the
other way round, the 3 x 2 with rows (1, 0), (2, 1) and (0, 3), with
`inner` 3 and one column: `a' b` with the B of 3 x 1 with rows (1), (0)
and (2) is the 2 x 1 with rows (1) and (6), and `a' b'` with the B of 1 x
3 with the row (2, 0, 1) is the 2 x 1 with rows (2) and (3). They are the
two matrices the second and third cases of "How it is verified" of "The
four operations of the PCA and the LD" already give for the same A, since
each of the four combinations computes the same product as one of those
cases, and they are asserted exactly. Run on both backends on 23 September
2026, and a backend whose result had its rows and its columns the other
way round was checked to fail them and to pass every other test of the
crate.

## The Rust interface

The lower half of `g += a'a`; `a` is `rows` x `cols`, and `g` is `cols`
x `cols`. `rows` may be 0.

```rust
pub fn add_self_product_lower(a: &[f64], rows: usize, cols: usize, g:
&mut [f64]) -> Result<()>;
```

The four products are one function, and each operand says how its own
buffer is laid out, which together says which of the four the call is, as
"The product with its first operand turned" lays out. `c`, which is
overwritten, is `rows` x `cols`. `rows` may be 0, and then nothing is
written, as an `a` of no rows adds nothing to `g` above: the second pass
of the PCA multiplies the block it has standardized by the eigenvectors,
and a block whose rows all had no variance leaves an `a` of no rows here
too. `inner` and `cols` are 1 at least.

The cases of an operand carry the same two values and differ in what its
rows are, so a caller cannot reach for the wrong one without writing the
name of the wrong one. They are one function and not several because
functions of the same arguments would take each other's call and give a
different matrix with no error: the check of a length is `rows` times
`inner`, or `rows` times `cols`, whichever way the buffer is read.

Giving the same slice for both operands, the first by the rows of the
result and the second by the columns of it, with `rows` equal to `cols`,
is the product of a matrix with its own transpose, which is what a set of
variants against itself asks for.

```rust
pub enum TheFirstOperand<'a> {
    /// `rows` x `inner`, one row for each row of the result, which is
    /// the first matrix of a product as it is usually written.
    ByTheRowsOfTheResult { values: &'a [f64], rows: usize },
    /// `inner` x `rows`, one row for each of the values the product
    /// sums over, which gives `c = a' b`.
    ByTheValuesSummedOver { values: &'a [f64], rows: usize },
}

pub enum TheSecondOperand<'a> {
    /// `inner` x `cols`, one row for each of the values the product
    /// sums over, which gives `c = a b`.
    ByTheValuesSummedOver { values: &'a [f64], cols: usize },
    /// `cols` x `inner`, one row for each column of the result, which
    /// gives `c = a b'`.
    ByTheColumnsOfTheResult { values: &'a [f64], cols: usize },
}

pub fn product(a: TheFirstOperand<'_>, inner: usize, b:
TheSecondOperand<'_>, c: &mut [f64]) -> Result<()>;
```

The eigendecomposition of the symmetric `g` of `n` x `n`, whose lower
half is read and whose buffer comes back holding the eigenvectors. The
values are from the largest, and row `j` of `vectors` is the eigenvector
of `values[j]`.

```rust
pub struct Eigen {
    pub values: Vec<f64>,
    /// n x n, row after row.
    pub vectors: Vec<f64>,
}

pub fn eigh_lower(g: Vec<f64>, n: usize) -> Result<Eigen>;
```

The error of the crate, which the core crate wraps.

```rust
pub enum Error {
    /// The name of the argument and what was expected of it.
    Dimension { argument: &'static str, expected: String },
    /// The name of the argument that holds a value that is not finite.
    NotFinite { argument: &'static str },
    /// A routine that stopped: the routine of LAPACK with its `info`,
    /// which did not converge when the `info` is positive and refused
    /// that argument of its own, a defect of popnei, when it is
    /// negative; or faer, with the `info` at 0, since faer gives none.
    NoConvergence { routine: &'static str, info: i32 },
    /// What could not be allocated and how many values it holds. The
    /// eigendecomposition asks for a workspace of floats and one of
    /// integers, so the count is of values and not of bytes.
    Memory { what: &'static str, values: usize },
}
```

### The signatures of the seven of the GWAS

The Cholesky factorization of the symmetric positive definite `a` of `n`
x `n`, read by its lower half: the lower triangular `l` with `l l' = a`,
which overwrites that lower half. The upper half is left as it was, and a
caller that needs `a` afterwards copies it first.

```rust
pub fn cholesky_lower(a: &mut [f64], n: usize) -> Result<()>;
```

The `x` of `a x = b`, where `l` of `n` x `n` is the factorization above.
`b` is `sides` x `n`, row after row, one row for each right hand side,
and it comes back holding the solutions the same way. `sides` is 1 at
least.

```rust
pub fn solve_with_cholesky(l: &[f64], n: usize, b: &mut [f64], sides:
usize) -> Result<()>;
```

The log of the determinant of the `a` whose factorization is `l`. There
is no sign to give: the determinant of a positive definite matrix is
above 0.

```rust
pub fn log_determinant_with_cholesky(l: &[f64], n: usize) -> Result<f64>;
```

The lower half of the inverse of that same `a`, into `inverse` of `n` x
`n`, whose upper half is left as it was. `l` and `inverse` are two
buffers and not one.

```rust
pub fn invert_with_cholesky(l: &[f64], n: usize, inverse: &mut [f64]) ->
Result<()>;
```

The thin QR of `a` of `rows` x `cols`, with `rows` at least `cols` and
`cols` 1 at least. It and `solve_upper_triangular` below are the "least
squares" that the `linalg` row of section 9 of `docs/architecture.md`
names: fitting a linear model to more individuals than coefficients is
`thin_qr` of the design and then `solve_upper_triangular` of its `r`.

```rust
pub struct ThinQr {
    /// `rows` x `cols`, row after row: the columns are of length 1 and at
    /// right angles to each other. The sign of a column is the backend's.
    pub q: Vec<f64>,
    /// `cols` x `cols`, row after row: upper triangular, its lower half
    /// 0, with `a = q r`.
    pub r: Vec<f64>,
}

pub fn thin_qr(a: &[f64], rows: usize, cols: usize) -> Result<ThinQr>;
```

The `x` of `r x = b` for the upper triangular `r` of `n` x `n`, whose
lower half is not read. `b` is laid out as it is for the Cholesky solve.

```rust
pub fn solve_upper_triangular(r: &[f64], n: usize, b: &mut [f64], sides:
usize) -> Result<()>;
```

How many singular values of `a` of `rows` x `cols` are above numpy's
tolerance: its rank, which is at most the smaller of the two dimensions.

```rust
pub fn rank(a: &[f64], rows: usize, cols: usize) -> Result<usize>;
```

The error enum gains one case, which "The errors the seven add"
describes.

```rust
    /// A matrix that could not be factored at the row the value names,
    /// counting from 0: the Cholesky reached a diagonal entry that is
    /// not above 0 there, or the solve against an upper triangular
    /// matrix reached one that is 0.
    Singular {
        /// The name of the argument, as this section spells it.
        argument: &'static str,
        /// The row the factorization stopped at, counting from 0.
        at: usize,
    },
```

## Speed

Measured on 22 September 2026 on the owner's Apple M5 Pro, rustc 1.98,
from a trial crate kept in `tmp/linalg_trial/`, not in git, on G = ZZ'
for Z of n x 200 with random values, one run each, which is why the
digits are few; the numbers of faer are of one thread and those of
Accelerate of the threads it takes by itself, which for the
eigendecomposition changes nothing measurable. G has rank 200, which does
not change the cost of the routines, all of which reduce the whole
matrix.

| n | `dsyevd`, Accelerate | faer |
|---|---|---|
| 1000 | 0.035 s | 0.096 s |
| 2000 | 0.27 s | 0.68 s |
| 5000 | 6.3 s | 11.4 s |

`dsyevr`, the routine of LAPACK that computes a chosen range of
eigenvalues and their vectors, asked for the 10 largest took 0.025 s,
0.18 s and 4.9 s. It is not taken now: the PCA gives every component
until Open 1 of `docs/specs/pca.md` is answered. The eigendecomposition
is 0.04 s of the 0.3 s of the PCA at 1000 individuals; at 5000 it is 6.3
s beside products of 5.8 s for 100000 variants, from the 58 ms per block
of `docs/specs/pca.md`, and with faer 11.4 s. Under pyodide on the same
machine the eigendecomposition with faer took 0.31 s at n = 1000, 2.4 s
at 2000 and 13.6 s at 3000, section 3.2 of `docs/rust_core.md`, on 19
September 2026.

The numbers to reach are those of the backends themselves, since the
crate adds an allocation, a call and the reading of the matrices it
checks: the times above for the eigendecomposition, and for the product
of a block of 5000 x 1000 with itself, 10.5 ms with Accelerate on one
thread, 95 ms with faer natively and 187 ms in wasm with `simd128`.

What the PCA measured through this crate on 22 September 2026, in
`docs/reports/pca-measurement.md`: natively a block of 5000 x 1000 takes
12.7 ms in `add_self_product_lower`, the 12.05 ms of `dsyrk` and the
checks; and in wasm a block of the analysis takes 215 ms, which is the
product and the standardizing together, with what is left of the
analysis once the blocks are paid, the eigendecomposition of 1000
individuals and the building of the result, at 0.197 s under node and
0.187 s under pyodide. Those two are a difference of two sizes and not a
call timed by itself, and they are below the 0.31 s that section 3.2 of
`docs/rust_core.md` measured for the eigendecomposition alone under
pyodide on 19 September 2026, with another build.

What the checks cost is the two readings of `add_self_product_lower`,
every value of A and the lower half of G, 44 MB for that block, which
nothing of the product needs. Measured on the same machine on 22
September 2026 with a scratch crate, the best of 10 runs of each, a
release build with `VECLIB_MAXIMUM_THREADS=1`: 1.68 ms of a call of
12.05 ms, an eighth of it. They stay because the alternative is a NaN
that one backend turns into an error and the other into a result, and
because a block the PCA standardized cannot be checked once for the
whole pass: each block is a new matrix.

### What the seven of the GWAS cost

Measured on 23 September 2026 on the owner's Apple M5 Pro, rustc 1.98,
release, from a trial crate kept in `tmp/linalg_gwas_trial/`, not in
git, the best of 3 runs at 5000 and above and of 5 below. The matrix of
each size is `ZZ'` for Z of n rows and as many columns of the xorshift
generator, 200 at n = 1000 and 20 above it, with n added to each of its
diagonal entries, which makes it positive definite; its rank before the
addition does not change the cost, since both routines reduce the whole
matrix. Accelerate takes the threads it finds, and faer is given twice,
once for each of its two rules: the global pool of rayon, which is what
it runs on natively, and one thread, which is what it runs on in the
browser.

| operation | n | Accelerate | faer, the pool | faer, one thread |
|---|---|---|---|---|
| `cholesky_lower` | 1000 | 0.0012 s | 0.0029 s | 0.0073 s |
| | 2000 | 0.0089 s | 0.0114 s | 0.0523 s |
| | 5000 | 0.135 s | 0.118 s | 0.762 s |
| | 10000 | 1.02 s | 0.816 s | 5.93 s |
| `invert_with_cholesky` | 1000 | 0.0025 s | 0.0041 s | 0.0140 s |
| | 2000 | 0.0148 s | 0.0182 s | 0.1020 s |
| | 5000 | 0.205 s | 0.183 s | 1.51 s |
| | 10000 | 2.05 s | 1.29 s | 11.8 s |

numpy 2.5.3 on the same machine and day took 0.723 s for `inv` of the
5000 x 5000, against the 0.340 s of `cholesky_lower` and
`invert_with_cholesky` together on Accelerate. The two numbers are what
each library takes from the matrix to its inverse: numpy's `inv` makes
its own LU inside the call and is one call, not two. numpy's `cholesky`
of the same matrix, which factors and stops, took 0.192 s, against the
0.135 s of `cholesky_lower`.

Those two together are line 689, which the fit of the logistic mixed model
runs a few dozen times, and they are what the open question of
`docs/rust_core.md` about that fit at many individuals is about: at 5000
individuals a few dozen turns of 0.34 s is a few tens of seconds, beside
the 25 s that document measured for pyNei there, and no backend of this
crate changes that.

Everything else is paid once for a whole study or once for a block.
`log_determinant_with_cholesky` reads n numbers. One system of the size of
the coefficients, `cholesky_lower` and then `solve_with_cholesky` with one
right hand side, took 0.060, 0.173 and 0.320 µs for systems of 3, 7 and 11
on Accelerate and 0.060, 0.159 and 0.277 µs on faer, so the 5000 systems
of a block at 7 take 0.9 ms and 0.8 ms, for one iteration of a fit.
`solve_with_cholesky` on a 5 x 5 with one right hand side for each of
10000 individuals took 0.070 ms on Accelerate and 0.140 ms on faer. On a
design of 10000 x 5, `thin_qr` took 0.165 ms on
Accelerate, 0.256 ms on faer and 0.209 ms on faer at one thread, and
`rank` 0.145 ms, 0.159 ms and 0.130 ms; numpy took 0.229 ms and 0.136 ms
for the same two.

Two of the LAPACK routines are much slower on the wide matrix than on the
tall one, so the BLAS backend writes the transpose of the design into a
buffer of its own, `rows` x `cols` values, 400 KB at 10000 x 5, and calls
them on that. Measured on the same day at 10000 x 5: `thin_qr` takes 0.165
ms that way and 2.98 ms called on the buffer as it lies, through `dgelqf`
and `dorglq`, which give the same `q` and `r` to 1e-13; and `rank` takes
0.145 ms against 1.54 ms. faer is given the buffer as it lies and copies
nothing.

## Open points

None. The four this spec had were decided by the owner on 23 September
2026 and each is in the text where its subject is, with the option that
was not taken and what is still unknown about it:

- The feature of `gemm` stays on and nothing is added for the sake of the
  linear algebra: "The wasm builds and the vector instructions", which
  also has where the rustc flag is set, for the `dists` module, and how
  the decision was worded on the day.
- Every square matrix but the `r` of the QR is factored with a Cholesky
  and not an LU: "Why a Cholesky where numpy uses an LU", which has the
  band of designs that changes.
- `product` gets a typed first operand, so that it computes all four of
  `a b`, `a b'`, `a' b` and `a' b'`: "The product with its first operand
  turned".
- The kinship is left as pyNei computes it, and the fit of the logistic
  mixed model fails with `Singular` when missing genotypes have made it
  not positive semidefinite: "What the GWAS calls of numpy, and where".
  This one is reconsidered if the implementation of the GWAS meets it.

Two of them leave work outside this spec. The rustc flag is named in a
comment of `crates/popnei-linalg/Cargo.toml` as the "meanwhile" of an open
point, which the commit that first changes code updates. The typed first
operand changes the signature of `product`, so the four calls of it in the
core crate, three in `pca.rs` and one in `ld.rs`, are rewritten with it,
computing what they computed before.

## Not in this spec

- Which of the four null models of the GWAS calls which of the seven
  operations, in what order, and what a `Singular` means where it is
  caught: `docs/specs/gwas.md`, which is not written. The seven are here
  because they are operations of linear algebra; the models are not.
- Two things this spec assumed about the GWAS so that it could name the
  seven, which the GWAS spec can contradict without changing any of them.
  That the four null models stay as pyNei fits them: another way of
  fitting the logistic mixed model, which `docs/rust_core.md`'s open
  question about that fit at many individuals invites, would want other
  operations, most likely fewer of them. And that the one system for each
  variant of a block is a loop in the caller, which holds while those
  matrices are of the size of the coefficients; a model with many
  coefficients would make a routine that takes a whole stack worth
  looking for, and neither LAPACK nor faer has one.
- The two functions the GWAS turns a statistic into a p-value with: the
  complementary error function, which gives how much of a normal
  distribution lies past a point, and the incomplete beta function, which
  gives the same for the t distribution the tests with an estimated
  variance use. Both are arithmetic over one number and not linear
  algebra: `docs/specs/gwas.md`.
- Whether the logistic mixed model can be fitted without inverting an n
  x n matrix a few dozen times, which is an algorithmic question and not
  one of this crate: the open question of `docs/rust_core.md` about that
  fit at many individuals.
- Which routines the kinship calls: the product above;
  `docs/specs/kinship.md`, which is not written.
- The threads of popnei, which no argument sets: `docs/specs/dists.md`.
- The Linux and Windows builds of OpenBLAS: the plan that first makes
  them.
- faer against OpenBLAS on x86: `docs/rust_core.md`, measured with the
  feature of this crate on a Linux machine.
