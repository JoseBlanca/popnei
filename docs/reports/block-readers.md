# Work report: readers that give blocks

The plan `docs/plans/block-readers.md` is under way, on the branch
`plan/block-readers`, in the worktree `.claude/worktrees/block-readers`,
since 20 September 2026. The orchestrator, in this report, is the session
of the assistant that runs the plan: it sends each task to a subagent,
checks what comes back and has it reviewed. What the owner reads first,
whether the plan is done, what exists now and what is asked of them, is
written here when the plan ends.

## Before the first task

The owner approved the plan in chat on 20 September 2026, and it was
committed on `main` as a11aea9, which the branch starts from. Everything
of "What has to be in place" was checked by running it, in the worktree:

- The three specs, the architecture and the glossary are on `main` since
  8f19650 and cae6a4c, and `git status --short` was empty. None of the
  specs has an open point.
- `cargo fmt --all --check` exit 0; `cargo clippy --workspace
  --all-targets -- -D warnings` no warning; `cargo test --workspace` `92
  passed`; `cargo wasm-check` finished; `uv run ruff format --check` `8
  files already formatted`; `uv run ruff check` `All checks passed!`;
  `uv run maturin develop && uv run pytest` `38 passed`; `npm run build`
  and `npm test` in `js/popnei` `tests 39`, `pass 39`, `fail 0`;
  `bash scripts/build_pyodide_wheel.sh` built the wheel and `node
  tests/pyodide/smoke.mjs` exited with 0.
- `which` finds cargo, rustc, uv, node, npm, wasm-bindgen, bcftools, bgzip
  and plink2; `~/devel/emsdk/emsdk_env.sh` and
  `~/devel/pyodide-venv/bin/pyodide` are there; both wasm targets are
  installed.
- `git ls-remote https://github.com/JoseBlanca/pynei main` gives ef0ca6e.
- The checks that say the work is done fail, as they must: 98 mentions of
  `read_variant`, `VariantReader` or `BlockCollector` in seven files of
  `crates`, none of `BlockReader`, `VariantRef` or `reblock`, no `LICENSE`,
  pyNei at `/Users/jose/devel/pynei` in `pyproject.toml`. The second gzip
  member of `tests/reference/vcf/many.vcf.gz` ends at byte 12336 with 280
  variants before it, recomputed from the file.
- The file of the bench was made again with
  `crates/popnei/benches/make_big_vcf.py`, 403572954 bytes, and its
  `bgzip -k`, 37695707 bytes, in `/Users/jose/devel/popnei-bench/`, a
  directory that the orchestrator made outside the repository and that
  the owner can delete when the plan is merged.

The reader as built, timed again before any change, `cargo bench --bench
read_vcf -- big.vcf --threads N --runs 5`, the owner's Apple M5 Pro,
release, the file in the page cache, the genotypes asked for: a median of
1.307 s on one thread and of 0.161 s on 18, where the report of the plan
before this one has 1.24 s and 0.160 s. The one thread number of today is
the one that the new reader is compared with.

What the orchestrator changed in the plan before it was approved, none of
which changes what is built: the opening names a third change that a user
sees, a column that was not asked for is no longer checked, which "How it
runs" of `docs/specs/io_vcf.md` decides; task 1.1 corrects section 8 of
`docs/architecture.md` too, which has the path dependency on pyNei as
`docs/objectives.md` has; tasks 2.2 and 2.3 run one after the other,
because `memchr` is not a dependency today and 2.3 writes the manifest of
the core crate and `Cargo.lock`, and 2.6 runs after 2.5; tasks 2.4 and
2.5 name "The Rust interface" of the VCF reader spec.

## Work package 1: the five decisions

Tasks 1.1 and 1.2 went to one subagent in one prompt, since both are
small: commits 8ca8935 and a16ea06, 113 thousand tokens and 7 minutes.

### The deliverables, run by the orchestrator at a16ea06

1. The license. `LICENSE` at the root starts with `MIT License` and
   `Copyright (c) 2026 Jose Blanca`; `cargo metadata --no-deps
   --format-version 1` gives `MIT` for `popnei`, `popnei-python` and
   `popnei-js`, which take it from the workspace manifest;
   `pyproject.toml` has `license = "MIT"` and `license-files`, and
   `js/popnei/package.json` `"license": "MIT"`; the wheel of pyodide holds
   `popnei-0.1.0.dist-info/licenses/LICENSE`, and `npm pack --dry-run` in
   `js/popnei` lists `LICENSE` among its 21 files. npm packs the files of
   the directory of the package and leaves a symbolic link out, so a
   `prepack` script copies the `LICENSE` of the root beside
   `package.json` before every pack, and `.gitignore` holds that copy:
   git has one `LICENSE`.
2. pyNei by its repository. `pyproject.toml` has pyNei as a git source of
   uv at ef0ca6e177be5a18c2dba1cc78901940ad847cf9; `grep -n "Users/jose"
   pyproject.toml uv.lock` finds nothing; `uv sync`, `uv run maturin
   develop` and `uv run pytest -k pynei` give `8 passed, 30 deselected`,
   and the whole of pytest `38 passed`; `grep -rn "path dependency"
   docs/objectives.md docs/architecture.md` finds nothing.
3. The `coding` skill. "Errors, and no panics" says that the core crate
   has one error enum, `non_exhaustive`, to which each module adds its
   cases, with the owner's two reasons, and how each binding crate turns
   it into what its language throws, in `errors.rs` of each. The
   orchestrator read the paragraph against the code: the two files, the
   enums `PyPopneiError` and `JsPopneiError`, and the cases it names are
   there.

### What was changed in the plan, and what the tasks did beyond it

`.claude/skills/code-review/categories.md` asked a reviewer of errors for
one error type per operation, so a reviewer would have reported the single
enum as a defect; the orchestrator added its correction to task 1.2.
The subagent found two more places that said what the tasks correct, and
corrected them: "Tests" of the `coding` skill, which had pyNei as a path
dependency, and `pyo3.md` beside it, which described a newtype with a
`From` for each error type of the core where the binding crate has an
enum with one case for the one error of the core.
