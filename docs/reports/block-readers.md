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
