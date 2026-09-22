# Work report: the principal component analysis and the linalg crate

The plan `docs/plans/pca.md` is under way, on the branch `plan/pca`, in
the worktree `.claude/worktrees/pca`, since 22 September 2026. The
orchestrator, in this report, is the session of the assistant that runs
the plan: it sends each task to a subagent on Opus, checks what comes
back and has each work package reviewed.

## Before the first task

The owner approved the plan in chat on 22 September 2026. The branch
stands on `main` at 7d8366f, the merge of the plan filters, with the two
specs and the reference data at 5852540 and the plan at 25562b4.

Run by the orchestrator in the worktree at 25562b4: `git log --oneline
-1 -- docs/specs/pca.md` gives 5852540; `which` finds `plink2`,
`Rscript` and `node` in `/opt/homebrew/bin` and `wasm-bindgen` in
`~/.cargo/bin`, and `~/devel/emsdk` and `~/devel/pyodide-venv` are
there; `pyproject.toml` names pyNei at ef0ca6e. `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `306 passed`, 2 ignored; `cargo
wasm-check` finished; ruff `18 files already formatted` and `All checks
passed!`; `uv run maturin develop && uv run pytest` `174 passed`; `npm
test` in `js/popnei` `tests 126`, `fail 0`; `scripts/build_pyodide_wheel.sh`
built the wheel and `node tests/pyodide/smoke.mjs` exited with 0, after
an `npm install` in `tests/pyodide`, which a new worktree lacks.

## Work package 1: the linalg crate

Task 1.4, the documents, ran side by side with task 1.1 and came first,
at 9ff88e7: the layout of section 8 of `docs/architecture.md` has the
crate and a paragraph on its two backends, section 9 names it in its
heading, its column header and its row, `docs/objectives.md` says "one
small linear algebra crate", and three bullets of the `coding` skill,
on threads, on the linear algebra and on `unsafe`, say what the spec
says. The message of the board is `2026-09-22T1148-plan-pca.md`; the
README of the board asks for UTC in the name, and the messages there
are named in local time, so this one is too, so that `ls` keeps them in
order. Two things for the owner from that task: the rule of the `coding`
skill that a new dependency of the core crate is pure Rust still stands,
and the linalg crate, which natively links BLAS and LAPACK, is one; and
the branch `plan/dists-kosman` announced on the board that it adds a
module to `lib.rs` of the core, cases to its error enum, an export to
`python/popnei/__init__.py` and to `js/popnei/src/index.ts`, and a
function to each binding crate, which are the files work packages 2 and
3 of this plan change too, so the two branches meet there at the merge.
The subagent of task 1.4 used 91183 tokens.

Task 1.1, the crate with the BLAS backend, is at ad170a2: the crate
with its `Error`, the three functions, the checks of the arguments in
`lib.rs` and the four calls of BLAS and LAPACK in `src/blas.rs`, each
`unsafe` block with its `SAFETY` comment and an `#[expect(unsafe_code)]`
that the workspace lint table takes without a change. `cargo test -p
popnei-linalg` `19 passed`, 16 of which failed against the stubs before
the code. One choice the spec did not state, that `product` takes a
`rows` of 0 and writes nothing, as the self product does, because the
second pass of the PCA gives such a matrix for a block whose rows all
had no variance, goes into the spec with task 1.2. The subagent used
144815 tokens.

Task 1.2, the faer backend and the feature `blas`, is at 6bb3685, after
cda2579, which put into the spec what task 1.1 chose for a product of
no rows. The four BLAS crates are optional dependencies of the targets
that are not wasm behind the feature, on by default, and faer 0.24.4 is
in every build; `cargo tree -p popnei-linalg -e normal --depth 1` lists
blas, blas-src, faer, lapack, lapack-src and thiserror, and with
`--no-default-features` faer and thiserror alone, so cargo takes an
optional dependency under a target table as the spec assumed. `cargo
test -p popnei-linalg` `20 passed` and the same with
`--no-default-features` `20 passed`, run by the orchestrator. The 1000 x
1000 test builds G with the crate's own self product, so it checks that
operation at size too, and asserts the trace from the diagonal and from
the sum of the eigenvalues. The subagent, the one of task 1.1, had used
210394 tokens at the end of the two tasks.

Task 1.3, the wasm builds, is at 0e4791d: the `wasm-check` alias names
the crate beside the core, the feature `wasm-simd128-enable` of `gemm`
0.19 is a dependency of the crate under the wasm cfg, and one `gemm`
alone is in the tree of the wasm target, under faer, carrying the
feature; the core crate depends on the crate with a `use popnei_linalg
as _;`, which is what keeps it in the link until work package 2 calls
it. The extension module that maturin links natively carries
`-framework Accelerate` through the build script of `accelerate-src`,
and popnei needs none of its own:
`otool -L` on `python/popnei/_core.cpython-314-darwin.so` lists the
framework. So the first "What could go wrong" of the work package did
not happen. The subagent had used 232493 tokens at the end of its three
tasks.

The deliverables, checked by the orchestrator at 0e4791d:

1. `cargo test -p popnei-linalg` `20 passed`, and with
   `--no-default-features` `20 passed`; `-- --list` names 20 tests, among
   them the self product with an A of no rows, the product of two
   matrices that are not square, the 3 x 3 eigendecomposition and the
   1000 x 1000 one of the generator. 8 were asked.
2. The errors: a dimension that does not match (a `g` that is not cols x
   cols, a `b` without the inner dimension, a `g` shorter than n x n), a
   buffer too short (an `a`, a `c`), a `cols`, an `inner` and an `n` of
   0, a value that is not finite in `a`, in the lower half of `g` and in
   the operands of the product, and the decomposition that did not
   converge, each a test in that list.
3. `cargo wasm-check` finished; the manifest of the crate has `gemm`
   under `[target.'cfg(target_family = "wasm")'.dependencies]`;
   `scripts/build_pyodide_wheel.sh` exit 0 and `node
   tests/pyodide/smoke.mjs` exit 0; `npm run build && npm test` in
   `js/popnei` `tests 126`, `fail 0`.
4. `grep -n "linalg"` on the three documents finds no "module" beside
   it, and section 8 of the architecture lists `crates/popnei-linalg/`.
5. `cargo fmt --all --check` exit 0; `cargo clippy --workspace
   --all-targets -- -D warnings` no warning; `cargo test --workspace`
   `306 passed`, 2 ignored, and `20 passed`; ruff `18 files already
   formatted` and `All checks passed!`; `uv run maturin develop && uv
   run pytest` `174 passed`.

The review, six reviewers at 0e4791d, spec, tests, numbers, errors, api
and architecture, and the fixes at 27f66e1 to dbcecb3. What was found
and is fixed:

- The `gemm` line of the workspace kept gemm's default features, which
  put rayon and ten more crates into both wasm trees, against "Threads"
  of the spec; `default-features = false` takes them out and the
  feature `wasm-simd128-enable` still reaches the one gemm faer uses.
- Nothing above the linalg crate could turn BLAS off, so the spec's
  sentence that a machine with no BLAS builds popnei held for the crate
  alone; the core and the Python binding crate now have a feature `blas`,
  on by default, that forwards to the crate's, and `cargo test -p popnei
  --no-default-features` passes on faer.
- Four tests could not fail: the product test had a symmetric result on
  a 2 x 2, so a transposed `c` passed; `c` was always zeros, so a product
  that accumulated passed; no test put a value that is not finite on the
  diagonal of `g`, so a check cut one entry short passed; and the two
  comparison helpers of the tests returned false on NaN, so an all NaN
  eigendecomposition passed every assertion. Each has its test now.
- Two silent returns in the reversal of the eigenvectors of the BLAS
  backend, unreachable today, would have paired values from the largest
  with vectors from the smallest; the reversal is made once in the
  crate's own function over a length it has checked.
- A `g` longer than its dimensions was truncated in silence by the
  eigendecomposition and written in part by the product; both refuse it,
  and the spec says so.
- The workspace of `dsyevd`, 2n² floats, was allocated with `vec!`,
  which ends the process when the machine has not the memory; it is
  asked for with `try_reserve` and the crate's `Error` has a case
  `Memory`, added to the spec.
- The message of a decomposition that did not converge said "it gave
  the info 0" for faer, which gives none, and "did not converge" for a
  negative info of LAPACK, which is an argument the routine refused, a
  defect of popnei; the three cases have their own message.
- A dimension above `i32::MAX`, the integer of BLAS, was an error on one
  backend only, and the count of values of a matrix was not bounded
  against it; both backends refuse them in the crate, and the spec
  names the limit.
- The 12 digit literals of the spec left the test of the smallest
  eigenvalue, compared within 1e-12 relative, 3.6e-13 of its budget for
  the rounding of the literal; the spec and the test carry the digits
  that name numpy's `f64`.
- `Eigen` derived `PartialEq`, an `==` on floats that nothing used;
  dropped. Doc comments, one function name and a comment corrected.

Not taken, with the reason: the panic inside faer's
`get_global_parallelism` when a program disabled it, which nothing in
popnei does; a test of no convergence coming out of `dsyevd`, for which
no input is known; a test of the threads under wasm, which work package
4 runs; `Error` not `#[non_exhaustive]`, which is the spec's enum.

What the owner should know:

- The two scans for values that are not finite, which "Errors" of the
  spec asks for before any routine runs, take 1.68 ms of the 12.05 ms of
  the product of a 5000 x 1000 block with itself on one thread, best of
  10, and the spec's "Speed" now says so. Dropping the scan of `g`, the
  matrix the crate itself wrote, would save about half; that is a change
  to "Errors" of the spec.
- A finite input can give a result that is not finite with no error:
  `eigh_lower` of a 2 x 2 of 1e308 gives an infinite eigenvalue on both
  backends, and the second one is 0 on BLAS and NaN on faer. The spec
  checks the inputs and says nothing of the outputs.
- faer's smallest eigenvalue of the 1000 x 1000 matrix moves by 3.6e-17
  relative between one thread of rayon and the pool, so a test of the
  PCA compares within the tolerance of the spec and not to the bit.
- The `accelerate` feature of the BLAS crates is unconditional, and
  `accelerate-src` emits `-framework Accelerate` on every OS; the Linux
  build is the next plan's and the workspace manifest says the two lines
  are macOS only.

After the fixes, at dbcecb3: `cargo test -p popnei-linalg` `33 passed`,
with `--no-default-features` `28 passed`, the five of the workspace of
`dsyevd` being of the BLAS backend; `cargo test --workspace` `306
passed`, 2 ignored, and `33 passed`; `cargo test -p popnei
--no-default-features` `306 passed`; clippy, fmt and `cargo wasm-check`
clean; ruff clean; pytest `174 passed`; the wheel built and the smoke
test exited 0; `npm run build && npm test` `tests 126`, `fail 0`; rayon
in neither wasm tree, and no BLAS crate under `-p popnei
--no-default-features`.

How the work went: the four tasks and the fixes went to one subagent,
371966 tokens at the end; the six reviewers used 93145, 123151, 95743,
91363, 95193 and 109045 tokens. No task had to be sent twice.

## Work package 2: the PCA of a table

Task 2.1, the PCA of a table in the core, is at 6ca1199, after 1b0cc72,
which put into the spec two tables with fewer rows than traits and
their numbers from numpy, so that the side of the product that iris does
not take has a test. The module `pca` of the core copies the table once,
transposed when it has fewer rows than traits so that both sides use the
same product, drops the components under λ₁ · max(n, p) · 2.2e-16 and
fixes the sign. `cargo test -p popnei --lib pca` `13 passed`, 6 asked;
broken on purpose, n - 1 in the divisor fails 2 of them, every component
given fails 2, no sign rule fails 4. Two errors the spec did not list
were added to it, a buffer that is not rows x columns and the wrap of an
error of the linalg crate, both a `RuntimeError`. Two things for the
owner. pyNei's 3 x 3 table of `test_pca_refuses_traits_with_no_variance`
has a near tie in its second component, two projections of ±0.7071 one
bit apart, so the sign of that component is decided by rounding and
differs between the two backends: the spec's promise that the three
builds give the same numbers does not hold on an exact tie, and the test
compares magnitudes there, as the spec now says. And the `Memory` error
of the linalg crate arrives in the wrapped case, so task 2.2 makes it a
`RuntimeError` and not a `MemoryError`. The subagent used 213330 tokens.

Tasks 2.2 and 2.3 ran side by side in the one tree, each on its own
files. Task 2.3, `doPca` in TypeScript, is at a07c2ce: the binding in
`crates/popnei-js/src/pca.rs`, `doPca` and `PcaResult` in
`js/popnei/src/pca.ts`, exported from the two entry points of the
package, `node.ts` and `web.ts`, since it has no `index.ts`, and 12
tests in `test/pca.test.ts`, iris and the 3 x 5 table of the spec. The
TypeScript layer checks that `data` holds exactly numRows x numCols
values, because the core reads the first values of a longer buffer; that
goes to the review. The subagent used 142411 tokens. Task 2.2, `do_pca`
in Python, is at b3a2575: `_core.pca` reads the float64 array without a
copy and runs the core with the interpreter released; `do_pca` and
`PCAResult` in `python/popnei/pca.py`; the error of the traits with no
variance crosses as a subclass of `ValueError` that carries the
positions, and the Python layer names the traits as pyNei's message
does; 8 tests in `tests/test_pca.py` against pyNei at ef0ca6e; pandas
3.0.2 became a dependency of the package, which returns frames and
declared only numpy. One trap found: `as_slice()` of the numpy crate
accepts a Fortran contiguous array and hands its values column after
column, where `.claude/skills/coding/pyo3.md` says it fails, so the
binding asks the array for its layout; the skill is corrected with the
review. The subagent used 168021 tokens.

The deliverables, checked by the orchestrator at b3a2575:

1. `cargo test -p popnei --lib pca` `13 passed`, 6 asked, with iris
   standardized and not, the sign rule, the tie, the component with no
   variance on pyNei's 3 x 3 table and the errors.
2. `uv run pytest tests/test_pca.py` `8 passed`: iris with both values
   of `standardize_data` and with `center_data` false against pyNei,
   all 4 components within 1e-9 after the sign rule, the index and the
   columns of the frame, `pass_stats` `None`, and the error that names
   the trait.
3. `npm run build && npm test` `tests 138`, `fail 0`, 12 of them in
   `test/pca.test.ts`.
4. `cargo fmt --all --check` exit 0; clippy no warning; `cargo test
   --workspace` `319 passed`, 2 ignored, and `33 passed`; ruff `20 files
   already formatted` and `All checks passed!`; `uv run maturin develop
   && uv run pytest` `182 passed`; `cargo wasm-check` finished.
