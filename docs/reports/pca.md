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

The review, seven reviewers at b3a2575, spec, tests, numbers, errors,
api, binding and architecture, and the fixes at d512968 to 61b49d3.
What was found and is fixed:

- A trait whose squared deviations overflow, values above about 1e154,
  was zeroed in silence and left the analysis with a weight of 0 and the
  whole variance given to the other traits; one whose squares underflow,
  values below about 2e-162, or whose sum overflows, reached the linalg
  crate as a value that is not finite and gave a `RuntimeError` blaming
  popnei. Each is a `ValueError` naming the trait and which of the three
  happened, in the spec, the core, Python and TypeScript.
- Two products overflowed before the division that would have kept them
  finite: the percentages, 100 times λ before dividing by the sum, were
  infinite for λ above 1.8e306 where pyNei divides first, and the
  threshold of a component with no variance, λ₁ times max(n, p), dropped
  every component for λ₁ above 1.8e308 / max(n, p). Both divide first.
- A table in which no component has variance gave empty frames and no
  word; it is a `ValueError`, "no trait has variance, there is nothing
  to do a PCA with", the wording of the variants half, and pyNei gives
  NaN percentages there. This is a value a user sees that the spec did
  not state, and the owner can reverse it.
- The core read the first values of a buffer longer than rows x columns
  while its error case said the buffer had another size; it refuses it.
- Four tests could not fail: no fixture had more rows than traits and a
  component dropped, so the truncation of the weights on that side was
  unguarded; three of the four `ValueError`s the spec promises in Python
  had no pytest, so moving them to the `RuntimeError` arm failed
  nothing; the ten-name cutoff of the no-variance message was
  unguarded; the branches of the TypeScript arguments for an undefined
  option and a transferred buffer had no test. Each has its test now.
- The exception of the traits with no variance carried the positions
  and no message on the Python side; it carries the core's message and
  the positions, and the new error of a trait out of range does the
  same.
- A missing value of a nullable pandas dtype gave numpy's `TypeError`;
  the frame is converted with NaN for it, so the core names the row and
  the trait in a `ValueError`.
- The defaults `center_data` and `standardize_data` were written in the
  Python and the TypeScript layers and nowhere in the core; they are two
  `pub const` of the core that both bindings export and read, as the
  ploidy of the VCF reader is.
- The JavaScript binding's export was `do_pca` where the Python one was
  `pca`; both are `pca`. A `bool` parameter of a private function became
  an enum. Doc comments, docstrings, comments with a wrong number, and
  the `pca` row of section 9 of the architecture, which did not name
  the PCA of a table, corrected. The spec says which three fields the
  TypeScript result of `doPca` leaves out, what the memory of the
  analysis is as the code holds it, and that the sign of the second
  component of pyNei's 3 x 3 table is decided by rounding.
- The `pyo3.md` of the coding skill said `as_slice()` fails on an array
  that is not C contiguous; it takes a Fortran contiguous one and hands
  its values column after column, and a pandas frame of one dtype is
  Fortran contiguous, so a binding that trusted the skill would read
  every table transposed. The skill says to ask `is_c_contiguous()`
  first.

Not taken, with the reason: naming the row and the trait in the message
of a value that is not finite, which the spec asks of the no-variance
message only; a message naming the argument for a wrong dtype at
`_core.pca`, unreachable from `do_pca`; pinning the threshold constant
with a fixture near it; a `popnei:` prefix on the core's messages in
TypeScript, the package's existing convention.

What the owner should know:

- With more traits than rows the core holds the weights twice, so the
  memory of `do_pca` on a 1000 x 8000 frame is 216 MB in the core plus
  the 64 MB contiguous copy Python makes, where the spec said the table
  plus the square of its smaller side; the spec says now what the code
  uses, and a strided write would save 64 MB there.
- `do_pca` is one call into the core with the interpreter released, and
  a Ctrl-C reaches it only when it returns; the core has no callback
  for it. The PCA of the variants runs block by block and can check
  for signals between blocks.
- The threshold of a component with no variance is pinned by no test:
  a threshold 1e6 times larger passes every test, since the smallest
  real eigenvalue of every fixture is far above it.
- The Python tests of `do_pca` compare with pyNei at run time and hold
  no literal of R; the core's tests hold every row of R's files.

After the fixes, at 61b49d3: `cargo test -p popnei --lib pca` `19
passed`; `cargo test --workspace` `325 passed`, 2 ignored, and `33
passed`; clippy, fmt and `cargo wasm-check` clean; ruff clean; `uv run
maturin develop && uv run pytest` `192 passed`, 18 in
`tests/test_pca.py`; `npm run build && npm test` `tests 142`, `fail 0`;
the wheel built and the smoke test exited 0 with pandas now a
dependency of the package, which micropip fetches under pyodide.

How the work went: the three tasks went to three subagents, 213330,
168021 and 142411 tokens at the end of the tasks and 290523, 224367 and
173885 after the fixes; the seven reviewers used 132110, 119847,
120233, 110675, 115533, 104737 and 94049 tokens. No task had to be sent
twice.

## Work package 3: the PCA of the variants

Task 3.1, the first pass, is at 5a14396, after 292d96b, which put into
the spec three sizes the analysis refuses: a ploidy above 254, since a
genotype's code is one byte; more than 46340 individuals, since the
linalg crate counts a matrix in the 32 bit integer of BLAS and 46341²
is above it; and more variants than a `usize` counts. `cargo test -p
popnei --lib pca` `35 passed`, 16 of them of the variants, the panel's
in 12 ms. Three things to know. `reblock` inside `pca_of_variants`
takes the size popnei chooses, 10000 variants at 5 or at 200
individuals, so the reader's block size never reaches the product
through the public function: the test of the block size is at the
first pass, with blocks of 1, 2 and 3, comparing G within 1e-10, and the
test through `pca_of_variants` with blocks of 1, 2 and 5 of the reader
stays. The standardizing loop is three passes over a row and not two,
because the major allele has to be known before a dosage is: one writes
a byte per genotype, the dosage or 255 for a missing one, one counts
the codes in runs of 255 with counters of one byte, and one looks each
code up in a table of 256 values. The used variants are kept as the
positions the result carries and not as one bit per variant, which
would be a second copy. The subagent used 323610 tokens.

Task 3.2, the second pass, is at bd3b000, after 9fbcd86, which put into
the spec that weights of more values than a `usize` counts are refused:
`num_prin_comps` times the variants used is the one product of this
analysis whose sides are not both at most the individuals. Each block
writes its weights into the columns of its own variants, so nothing
holds them twice. A second pass over another dataset, other individuals
or another ploidy, is refused before its first block, beside the three
the spec named, a variant more, a variant fewer and a variant that lost
its variance; the four are one error, a `RuntimeError` in Python.
`num_prin_comps` of the result is what was asked for capped at the
components with variance, which is what the bindings read for the names
of the components. `cargo test -p popnei --lib pca` `41 passed`, all in
30 ms. The subagent, the one of task 3.1, had used 375363 tokens at the
end of the two tasks.

Tasks 3.3 and 3.4 ran side by side. Task 3.4, `doPcaFromVariants`, is at
7113c7c: the binding opens both chains over the source and the steps and
lends them to the core, as `writeVars` does for one, and each source, a
VCF and a vars file, gets the method; the result extends the one of
`doPca` with the individuals, the used variants and the pass stats; 8
tests. The names of the individuals do not cross wasm a second time, the
package having read them when the file was opened. The subagent used
250839 tokens. Task 3.3, `do_pca_from_variants`, is at 1f7e516: the
binding opens a reader per pass and reads the filters' counts from the
first; the four errors of a dataset the analysis cannot read are a
`ValueError` naming the file and the two of a second pass that differs a
`RuntimeError`; 10 tests, 27 in the file. One decision to know: the
components are named for the components of the projections and the
weights take the first of those names, so the panel's first weight row
is `PC000` and not `PC00`, which is what makes one component have one
name in both frames. A negative `num_prin_comps` is refused in the
Python layer, where the message names the argument. Over the panel's
first 10 components popnei and pyNei differ by 1.0e-12 at most, where
the projections reach 17.10; without the sign rule the same comparison
differs by 34.2. The subagent used 304934 tokens.

The deliverables, checked by the orchestrator at 1f7e516:

1. `cargo test -p popnei --lib pca -- --list` names 41 tests, 22 of them
   of the variants, 12 asked: the worked example with `num_prin_comps` 3
   and 0, `worked3` without and with `transform_to_biallelic`, the panel
   with its five literals, the five cases pyNei asserts, the result that
   does not change with the size of the blocks, each error, and the
   second pass whose variants differ. They pass.
2. `uv run pytest tests/test_pca.py` `27 passed`: the panel and the
   worked example against pyNei at ef0ca6e, 10 and 3 components within
   1e-9 after the sign rule, `pass_stats` with `num_vars` 1200 and the
   count a `filter_by_maf` step kept, and `num_prin_comps` 0 giving
   `princomps` with no rows and the used variants as columns.
3. `npm run build && npm test` `tests 150`, `fail 0`, with the worked
   example at `doPcaFromVariants`, its projections, `usedVars` and
   `passStats`.
4. `cargo fmt --all --check` exit 0; clippy no warning; `cargo test
   --workspace` `347 passed`, 2 ignored, and `33 passed`; ruff clean;
   `uv run maturin develop && uv run pytest` `201 passed`; `cargo
   wasm-check` finished.
