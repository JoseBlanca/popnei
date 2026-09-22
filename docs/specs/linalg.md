# The linalg crate: the linear algebra and its two backends

September 2026. Every calculation of popnei that works on a block as a
matrix, the PCA, the kinship, the GWAS, needs a few operations of linear
algebra, and this crate is the one place that has them: natively they run
on the BLAS and LAPACK libraries of the system, the same ones numpy
calls, and in WebAssembly, where there is none, on faer, a linear
algebra library written in Rust. There is no code. This spec develops the
row `linalg` of the table in section 9 of `docs/architecture.md` and
decision 5 of `docs/rust_core.md`, which chose the two backends. It
covers the crate, its backends, the builds, and the three operations that
`docs/specs/pca.md` calls: the product of a matrix with itself, the
eigendecomposition of a symmetric matrix, and the product of two
matrices. The operations of the GWAS are later items.

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

WebAssembly has a set of 128 bit vector instructions, `simd128`, which a
build uses only when it is told to, with `-C target-feature=+simd128`,
and which faer's products use only when, besides, the crate that does
them, `gemm`, has its feature `wasm-simd128-enable` on. Measured on 21
September 2026 with the trial crate of `docs/specs/pca.md`, the product
of a block of 5000 variants x 1000 individuals with itself, its lower
half, under node 26: 306 ms without and 187 ms with both. The pyodide
wheel built with the flag in `RUSTFLAGS` loaded under pyodide on node 26
on 22 September 2026 and passed the smoke test of `tests/pyodide/`.
Whether the flag is turned on is **Open 1**, below.

Where the flag goes when it is: into the two commands that build for
wasm, `build:wasm` of `js/popnei/package.json` and
`scripts/build_pyodide_wheel.sh`, through `RUSTFLAGS`, which is how it
was tried, and not into `.cargo/config.toml`. A `rustflags` for a target
in that file makes cargo ignore the `build.rustflags` that the
`wasm-check` alias of the same file passes to deny warnings, and the
alias would stop denying them without a word. The feature of `gemm` is
a dependency of this crate under `cfg(target_family = "wasm")`.

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

## The three operations

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

### Layout, half and the backends

Every matrix crosses the interface row after row, and the Fortran
routines of BLAS and LAPACK read a matrix column after column. The buffer
of an r x c matrix read that way is its transpose, c x r, so the BLAS
backend calls each routine on the transposes: the lower half of G in
popnei's layout is the upper half for the routine, `uplo` is `U`; `A'A`
is `A A'` of the transposed view, `trans` `N`; and `C = AB` is
`C' = B'A'`, so `dgemm` gets the buffer of B as its first operand and
that of A as its second, both with `trans` `N`. The faer backend tells
faer that the buffers are row major and calls its functions as written.
A test of each product on matrices that are not square, "How it is
verified", is what catches a backend that mixed the two.

### Errors

Each of these is an error of the crate before any routine runs, and the
core crate turns it into its own error and Python into a `RuntimeError`,
since they are defects of the caller: a dimension that does not match,
a `g` that is not c x c for the product nor n x n for the
eigendecomposition, or an `a`, a `b` or a `c` shorter than its rows times
its columns, a longer one being taken by its first rows times columns
values; a c or an n of 0, and in `product` an `inner` of 0 as well, while
the rows of either product may be 0; a dimension, or a number of values
of a matrix, above 2147483647, which is what the routines of BLAS and
LAPACK count in, checked for both backends so that the two refuse the
same calls; and a value that is not finite in a matrix. The last is
checked here because the backends do not agree on it: `dsyevd` on a
matrix with a NaN gives NaN eigenvalues and an `info` of 0, measured
through numpy 2.5 on 22 September 2026, and faer's `self_adjoint_eigen`
gives its error of no convergence. A routine of LAPACK that stops,
`info` other than 0, which for `dsyevd` is an eigendecomposition that
did not converge when the `info` is positive and an argument the routine
refused, a defect of popnei, when it is negative, is an error with the
routine and the `info`, and faer's error for the same case is the same
error of the crate, with the `info` at 0, since faer gives none.

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
backends are checked on the size they will run at: G = ZZ' for Z of
1000 rows and 1200 columns, with z(i, c) the (c · 1000 + i)-th number of
the xorshift generator below, started at 7. G then has full rank and its
eigenvalues are apart, between 0.79 and 361.9, no two closer than 0.003.
The literals, from numpy 2.5.3 on 22 September 2026: the trace
99996.3873081677, which is the sum of the
eigenvalues, and which numpy adds up to 99996.38730816769 when it adds
the eigenvalues instead; the three largest eigenvalues 361.9125119011332,
359.66517178659313 and 356.4439329949563, and the smallest
0.7933289215408484; and the first three entries of the eigenvector of
the largest, with the sign of its largest entry positive,
0.018111301861995926, -0.004576022169100455 and
-0.00518748985091827. The test compares the
eigenvalues within 1e-12 relative and the entries within 1e-9, because
an eigenvector is less well determined than its eigenvalue by the gap to
its neighbours, and a matrix of a real dataset has closer ones than
this. Measured on 22 September 2026 with the trial of "Speed": the two
backends agree to 2.9e-14 relative on the eigenvalues and to 1.3e-12 on
the entries of the eigenvectors after the sign, and their products of
the same Z agree to the bit. A matrix of rank below n is not such a
test: with Z of 1000 x 200 the 800 eigenvalues that are 0 come out
between -8e-14 and 8e-14 and their eigenvectors are any base of that
space, different in each backend.

The generator, so that the test and numpy make the same Z: a 64 bit
state s that starts at 7 with its lowest bit set, and for each number
`s ^= s << 13; s ^= s >> 7; s ^= s << 17`, the shifts dropping the bits
that leave the 64, and the number is `(s >> 11) / 2^53 - 0.5`.

The checks at the level above, that a PCA made from these operations
gives the numbers of R, are in `docs/specs/pca.md`.

## The Rust interface

The lower half of `g += a'a`; `a` is `rows` x `cols`, and `g` is `cols`
x `cols`. `rows` may be 0.

```rust
pub fn add_self_product_lower(a: &[f64], rows: usize, cols: usize, g: &mut [f64]) -> Result<()>;
```

`c = a b`; `a` is `rows` x `inner`, `b` is `inner` x `cols`, and `c`,
which is overwritten, `rows` x `cols`. `rows` may be 0, and then nothing
is written, as an `a` of no rows adds nothing to `g` above: the second
pass of the PCA multiplies the block it has standardized by the
eigenvectors, and a block whose rows all had no variance leaves an `a`
of no rows here too. `inner` and `cols` are 1 at least.

```rust
pub fn product(a: &[f64], rows: usize, inner: usize, b: &[f64], cols: usize, c: &mut [f64]) -> Result<()>;
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
until Open 1 of its spec is answered. The eigendecomposition is 0.04 s
of the 0.3 s of the PCA at 1000 individuals; at 5000 it is 6.3 s beside
products of 5.8 s for 100000 variants, from the 58 ms per block of
`docs/specs/pca.md`, and with faer 11.4 s. Under pyodide on the same
machine the eigendecomposition with faer took 0.31 s at n = 1000, 2.4 s
at 2000 and 13.6 s at 3000, section 3.2 of `docs/rust_core.md`, on 19
September 2026.

The numbers to reach are those of the backends themselves, since the
crate adds an allocation, a call and the reading of the matrices it
checks: the times above for the eigendecomposition, and for the product
of a block of 5000 x 1000 with itself, 10.5 ms with Accelerate on one
thread, 95 ms with faer natively and 187 ms in wasm with `simd128`.

What the checks cost is the two readings of `add_self_product_lower`,
every value of A and the lower half of G, 44 MB for that block, which
nothing of the product needs. Measured on the same machine on 22
September 2026 with a scratch crate, the best of 10 runs of each, a
release build with `VECLIB_MAXIMUM_THREADS=1`: 1.68 ms of a call of
12.05 ms, an eighth of it. They stay because the alternative is a NaN
that one backend turns into an error and the other into a result, and
because a block the PCA standardized cannot be checked once for the
whole pass: each block is a new matrix.

## Open points

The owner decides this one. Until then the implementer follows its
"meanwhile".

**Open 1: `simd128` in the two wasm builds.** With the flag and the
feature of `gemm`, the product of a block takes 187 ms in wasm instead
of 306 ms, and the product is nearly all the time of the PCA in the
browser. A build with the flag does not load in a browser without those
instructions, which every browser of 2023 and later has, as Open 6 of
`docs/specs/pca.md` details; the pyodide wheel built with it loaded. The
options are to turn it on in both builds, in the wasm package alone, or
in neither. Recommendation: on in both. Meanwhile the flag is not set and
the feature of `gemm` is.

## Not in this spec

- What `pynei/gwas.py` calls besides the three operations here:
  `solve`, `inv`, `qr`, `slogdet` and `matrix_rank` of `numpy.linalg`.
  They are later items of this spec, written with `docs/specs/gwas.md`,
  which says which the four null models of the GWAS need and how each is
  computed.
- Which routines the kinship calls: the product above;
  `docs/specs/kinship.md`.
- The threads of popnei, which no argument sets: `docs/specs/dists.md`.
- The Linux and Windows builds of OpenBLAS: the plan that first makes
  them.
- faer against OpenBLAS on x86: `docs/rust_core.md`, measured with the
  feature of this crate on a Linux machine.
