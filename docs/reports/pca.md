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
