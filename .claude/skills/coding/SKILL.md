---
name: coding
description: How code is written in popnei, in the Rust core crate, in the pyo3 binding crate and in the Python package. Use it before writing or changing any code or test of popnei. It covers the order of the work, the arithmetic that keeps the numbers right, errors without panics, the types and the names, what the architecture asks of the code, the binding and Python layers, the tests and the checks to run before the work is called done.
---

# Coding

popnei has three layers, and `docs/architecture.md` describes them: the
core crate `crates/popnei`, pure Rust with no pyo3, where every
calculation lives; the binding crate `crates/popnei-python`, which
translates between Python and the core and does nothing else; and the
Python package `python/popnei`, which is the API, the result objects and
the tests. The goals are in `docs/objectives.md`, in order: right, usable
from Python as pyNei is, in the browser, fast where it matters, streaming.
When two of them pull against each other the earlier one wins.

Most of what follows is enforced by the lint table in `lints.toml`, beside
this file, which goes into the `Cargo.toml` of the workspace. A rule that a
lint enforces is given here with its reason and not repeated in detail.
The prose is for what no lint catches.

## Before the code

Code is written from a spec, `docs/specs/<module>.md`, and usually from a
step of an implementation plan. Read the item of the spec, the pyNei
function it names, and the part of `docs/architecture.md` it stands on.

When the spec can be read in two ways that give different numbers, run
pyNei, or the reference program, on a case where the two readings differ,
before choosing. It takes a minute, and the choice made without it is a
guess that the tests will not catch, because they were written from the
same reading.

When the spec does not say what should happen in a case, the choice is not
made in silence. A choice that changes a value a user sees or the public
API goes to the owner as an open point of the spec. A smaller one is made
and written in the spec, in a commit of its own that comes before the
commit of the code and the tests, as the `writing-specs` skill says, and
it is named in the commit message of the code.

## The order of the work

1. The test first. It asserts the numbers the spec gives, as literals, and
   it fails before the change. To make it compile, the new function gets a
   body that returns a wrong value, `None`, `Ok(0)`, and never `todo!()`,
   which the lints deny. Run the test and see it fail on the assertion.
   When a change cannot have such a test, a rename, a move, say so in the
   commit message.
2. The code, as small as the step asks for.
3. The checks at the end of this file, all of them, with their real
   output. A check that could not be run is reported as not run.
4. The commit, with the message the `writing` skill describes.

## Integers: no operator that can be silently wrong

In a release build `+`, `-`, `*`, `<<` on integers wrap when they overflow,
and the program goes on with a wrong number. In a debug build the same
line panics, and `/` and `%` by zero panic in every build. `as` between
integer types truncates in both. Neither is
acceptable in a library whose first goal is to be right, and the sizes of
popnei reach the limits: the largest dataset of the objectives, a million
variants of 10000 samples, has 1e10 genotypes, which does not fit in a
`u32`, nor in a `usize` under wasm, where it is 32 bits.

- A count that runs over the whole dataset is a `u64`. A `usize` is for
  indexing memory and for nothing that has to be the same in wasm and
  natively.
- Integer arithmetic uses the methods that say what happens on overflow:
  `checked_add` and its family, with the `None` turned into an error, as
  the default; `saturating_` and `wrapping_` only when that is the meaning
  wanted. The lint `arithmetic_side_effects` is denied, so a plain
  operator on integers does not compile.
- Where a bound makes overflow impossible, the plain operator is allowed
  with the bound written down: `#[expect(clippy::arithmetic_side_effects,
  reason = "at most num_samples * ploidy alleles in one variant, checked
  to fit in u32 when the reader is built")]`. This is how the loop over the
  genotypes of one variant stays free of checks: the per variant count is
  bounded, and the one checked addition is the one that adds it to the
  total of the dataset. clippy follows the range of some expressions,
  `2 * u64::from(x)`, and does not fire on them; an `#[expect]` there is
  itself an error, so the bound goes in a plain comment.
- The lint does not look at floats. Nothing checks the arithmetic that
  gives the result except the tests.
- Anything that comes from a file, a length, an offset, a position, is
  checked before it is used to index or to allocate.
- To widen, `u64::from(x)`. To narrow, `u32::try_from(x)?`. No `as`
  between integer types. From a float to an integer, check for NaN and for
  the range first, because `f64::NAN as usize` is 0.
- A count to `f64` is `as f64` and is exact below 2^53, which every count
  of popnei is. That cast is the one `as` that stays.

## Floats

- Results are `f64`. An `f32` appears only where a file format has one.
- A missing value is `Option<f64>` or an explicit enum inside the core,
  and becomes NaN only at the boundary with Python, where pandas expects
  it. NaN that travels through Rust arithmetic hides where it was born.
- A total of floats must not depend on the number of threads. rayon's
  `sum` and `reduce` join the parts in an order it chooses at run time, so
  the last bits change with the pool. Reduce over chunks of a fixed size,
  collect the partial results in index order and add them serially. The
  same for the order of a `HashMap`: it never feeds a float total. Sort
  the keys or use a `BTreeMap`. Integer totals are exact and reduce in any
  order.
- popnei does not promise the same bits on every platform. The owner
  decided that in September 2026: results agree with the reference
  programs and between platforms within the tolerance of the spec. `exp`,
  `ln`, `powf` and the trigonometric functions are not rounded the same on
  every platform, and macOS, Linux and wasm will differ in the last
  place. `sqrt` and the four operations are exact everywhere. So a
  test never asserts the bits of a result that went through the first
  group, and a test that compares two platforms uses the tolerance of the
  spec. `a.mul_add(b, c)` rounds differently from `a * b + c`; it is not
  swapped in when the number has to match R or plink2.
- Floats are compared within a tolerance that has a reason, the digits the
  reference program prints, and ordered with `total_cmp`. A comparison
  that decides a result, which variant passes, which allele is the major
  one, breaks its ties by a rule the spec gives, the lower index, and not
  by the last bit.

## Errors, and no panics

A panic in the core becomes a `PanicException` in Python, which derives
from `BaseException` so that it is not caught, and in a notebook or a
browser tab that ends the session. So library code does not panic:
`unwrap`, `expect`, `panic!`, `todo!`, `unimplemented!` and indexing with
`[]` are denied outside the tests. Slices are walked with iterators,
`chunks_exact(ploidy)`, `zip`, `get`, which are also what lets the
compiler drop the bounds checks.

- `Result` everywhere, fail fast, as section 7 of the architecture says. A
  malformed line of a VCF is an error with the line number and the field,
  not a warning and a skipped record.
- Errors are typed, with `thiserror`. One error type for each operation
  that fails in its own way, not one for the crate. The variants name what
  was being done, `ReadHeader`, `ParseGenotype`, and carry what is needed
  to find the cause: the path, the line, the field, the value.
- A public error type is `#[non_exhaustive]`, and it does not hold the
  error type of a dependency in a public variant.
- A `Result` is never dropped. `let _ =` on one needs a reason in a
  comment.
- The conversion to Python exceptions lives in the binding crate, in one
  newtype, as `pyo3.md` describes, and chooses the exception a pyNei user
  would expect: `ValueError` for a bad argument, `OSError` for a file.

## Types, names and defaults

- A name says what the value is: `called_alleles`, `max_missing_rate`,
  never `n`, `data`, `tmp`, `val`. The symbols of a formula, `p`, `k`,
  stay in the doc comment that gives the formula. Where pyNei has a name
  for the thing, that is the name, and `docs/glossary.md` has the names
  of the things of the domain, `pop`, `gts`, a block and not a chunk.
- The same thing has the same name everywhere in the three layers.
- Where `docs/architecture.md` or the spec gives a signature or a field,
  that is the signature, a bare `usize` index and a `bool` field included.
  When it looks wrong, that is a point for the owner and not a silent
  change.
- In what they leave free: two scalars of the same primitive that meet in
  one signature, an index of a sample and an index of a variant, get
  newtypes. No `bool` parameters, an enum with two named variants. No value of a
  finite set passed as a string inside Rust.
- A default that changes a result is a named `pub const` with a doc
  comment that says where it comes from: "20, inherited from pyNei, nobody
  has measured whether it is the right threshold". No `Default` derived on
  a struct whose fields change results, and no `Option` parameter that
  quietly becomes a value.
- A `match` on an enum of popnei names every variant, so that a new one
  does not fall into a `_` arm.
- Private by default, `pub(crate)` between modules, `pub` for what the
  binding crate or a user of the core calls. `MISSING_ALLELE`, never a
  literal -1.
- Every `pub` item has a doc comment as the `writing` skill describes it,
  with `# Errors` when it returns a `Result`.
- No `unsafe` in the core crate, `#![forbid(unsafe_code)]`. In the binding
  crate an `unsafe` block carries a `// SAFETY:` comment that names each
  condition and why it holds there.
- A lint is silenced with `#[expect(lint, reason = "...")]` on the
  smallest item, never with a bare `#[allow]`.
- A new dependency of the core crate is pure Rust, builds for
  `wasm32-unknown-emscripten`, and is justified in the commit that adds
  it: what it gives and why it is not a few lines of our own.

## What the architecture asks of the code

- The variants flow in blocks, and one variant is a view into a block.
  Nothing is allocated per variant: a `Vec` created inside the loop over
  the rows of a block is a defect, not a style point. A filter compacts
  the block it was given and does not build another.
- A calculation asks for the fields it reads with `Needs`, and checks that
  the column is there when it depends on one.
- rayon runs over the rows of a block. The matrix
  products run on BLAS natively. The two are never nested: a rayon worker
  that calls BLAS has it pinned to one thread, and a big product is called
  from outside rayon. The library never builds the global pool of rayon.
- Everything builds and runs with one thread, because wasm has no
  threads. Code that needs threads is behind
  `#[cfg(not(target_family = "wasm"))]` with a serial version beside it,
  and the choice of the linear algebra backend is a `cfg` on the target,
  not a pair of cargo features that exclude each other.
- Linear algebra goes through the `linalg` module and nowhere else.
- Readers take `impl Read` or `impl BufRead`, so that a test feeds them
  bytes from memory.

## The binding crate and the Python package

The binding crate translates and holds no logic. If a function there has
an `if` about genetics, it is in the wrong crate. Before touching it read
`pyo3.md`, beside this file: the current names of pyo3, which are not the
ones of a year ago, how arrays cross without a copy, releasing the
interpreter around long work with `py.detach`, classes that are `frozen`,
the one newtype that turns the errors of the core into Python exceptions,
and the two builds, native and wasm.

The Python package mirrors the signatures of pyNei, holds the defaults'
docstrings and the type hints of every public function, makes the arrays
contiguous and of the right dtype before they go to `_core`, builds the
frozen result dataclasses with their pandas frames, and does no
calculation. `ruff format` and `ruff check` clean.

## Tests

- The numbers of a test are literals, from the spec, which has them from a
  reference program or from pyNei. A test never computes its expected
  value with the code under test, nor with a second copy of the same
  formula. When the spec names a case and gives no number for it, the
  number is got by running pyNei or the reference program and is added to
  the spec with how it was got. Until then the test asserts only what the
  spec states, that a value exists or does not.
- The cargo tests cover the core on their own. The pytest tests run pyNei
  and popnei on the same input where they overlap, pyNei being a
  development dependency at the commit that `pyproject.toml` names, and
  compare as the spec says: exactly for counts and sets of
  variants, within the spec's tolerance for floats.
- Every field and every parameter takes, in some test, a value that
  differs from the others and from its default. A suite in which the
  exponent is always the ploidy, and the ploidy is always 2, cannot tell
  which of the two the code reads, nor whether it reads the ploidy at all.
- A test has to be able to fail. Check the fixture against the usual ways
  it cannot: a condition no fixture reaches; a fixture in which several
  wrong implementations give the same number, all weights equal, all
  populations of the same size, no missing genotype; a fixture in a regime
  where the thing tested does not happen. When in doubt, break the code on
  purpose and see the test fail.
- A parser gets the malformed inputs as cases of their own: a truncated
  line, a missing field, a genotype with the wrong ploidy, a bad gzip.
- A calculation that runs in parallel is tested with one thread and with
  several, over data that spans several chunks, and the results are equal.
- The name of a test says the behaviour and the outcome:
  `missing_filter_keeps_a_variant_exactly_at_the_threshold`.
- The tests are numerical work, so the test profile is built with
  `opt-level = 2`. Overflow checks and debug assertions stay on in it.

## Before the work is called done

```
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
uv run ruff format --check && uv run ruff check
uv run maturin develop && uv run pytest
```

The three cargo commands run for every change. The two Python ones run
from the moment the binding crate and the package exist, also for a change
in the core alone, because the pytest tests are the ones that compare with
pyNei. A layer that does not exist yet is reported as not there, not as
passed.

When the change touches what wasm builds differently, threads, the linear
algebra backend, a dependency, the wasm wheel is built as well, with the
steps the walking skeleton leaves in the repository.

Report what each command printed when it failed and that it passed when it
passed. Speed is not claimed without a measurement, with the dataset and
the machine, before and after, as the objectives ask.

## When this skill is wrong

No code existed when this was written, in September 2026, and neither did
the development container the project is going to have. When it exists,
this file says which of the commands above run inside it and which on the
host, the timings and the native build against Accelerate among them. The lint table
and the commands above are to be tried on the walking skeleton. A lint
that proves too noisy is changed in `lints.toml` with the count that
showed it, and a command that is not the right one is corrected here.
