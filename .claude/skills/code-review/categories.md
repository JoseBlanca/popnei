# The categories of a code review

A reviewer reads the section for its category and the parts of
`.claude/skills/coding/SKILL.md` it names. The rules are there with their
reasons and are not repeated here. What is here is what to look for and
how to look.

## spec

Does the code do what the spec item says, and what pyNei does?

- Go through the spec item part by part: what it gives, the cases a
  reader would not guess, what pyNei does that is odd and whether popnei
  reproduces it, the open points and their "meanwhile". For each
  statement find the code that makes it true and the test that would fail
  if it stopped being true. A statement with no code, or with code and no
  test, is a finding.
- Run the cases. The worked example of the spec, the edges it names, and
  two or three of your own that it does not name: nothing called, one
  sample, one allele, a half called genotype, a ploidy other than 2. Where
  pyNei covers the case, run pyNei from `/Users/jose/devel/pynei` on the
  same input and compare as the spec says to compare.
- Look for what the code does that the spec does not say: a default, a
  clamp, a skipped record, an early return. Each is either a gap of the
  spec or a defect of the code, and the finding says which you think it
  is.
- When the code and the spec disagree and the code looks right, the
  finding is about the spec.

## tests

Can each test fail, and are the numbers the change claims true?

- For every test of the change, break the code it guards and run it: flip
  a comparison, drop a term, return early, swap two arguments. A test that
  still passes guards nothing, and that is a finding. Before reporting
  one, show that your change did alter the behaviour on some input,
  because a change that alters nothing proves nothing about the test.
- The ways a fixture hides a defect: a condition no fixture reaches;
  values under which several wrong implementations give the same number,
  all populations of one size, all weights equal, nothing missing, the
  exponent equal to the ploidy; a regime where the thing tested does not
  happen.
- The expected values are literals from the spec. A test that computes its
  expectation with the code under test, or with a second copy of the
  formula, is a finding.
- Parallel code: is it run with one thread and with several, over data
  that spans several chunks, and are the results compared?
- Every number that a doc comment, a test comment or the commit message
  gives about this change, how many cases, how large a difference, how
  much faster, is computed again, not read again. Report each as right or
  as wrong with the right value. This is where wrong numbers gather: the
  ones copied from a source are usually right, the ones the writer states
  about their own work often are not. An explanation of why something
  happens is a claim too, and is checked by reproducing it.
- Restore the code after each experiment and end with `git status` clean.

## numbers

The sections "Integers" and "Floats" of the coding skill.

- Every `#[expect(clippy::arithmetic_side_effects)]` and every comment
  that gives a bound: is the bound true, and is it established where the
  reason says, for every caller? This is the first thing to check,
  because a false bound is an overflow with a signature on it.
- Counts that run over the dataset in anything narrower than `u64`, and
  `usize` in anything that must be the same under wasm.
- `as` between integer types, and from float to integer without the check
  for NaN and range.
- Values read from a file that reach an index, a length or an allocation
  unchecked.
- The lint does not look at floats, so read the float arithmetic as
  arithmetic: a division whose divisor can be 0, a subtraction of nearly
  equal numbers, a sum of many small terms, the order of operations
  against the formula of the spec.
- Float totals that depend on the threads or on the order of a `HashMap`.
  `mul_add` where the number has to match a reference. NaN going through
  arithmetic inside the core. Floats compared with `==` or ordered with
  `partial_cmp().unwrap()`. A tie that the last bit decides.

## errors

The section "Errors, and no panics" of the coding skill.

- Every path to a panic outside the tests: `unwrap`, `expect`, `[]`,
  `panic!`, a division by zero, a `RefCell` or a lock, an `#[expect]` of
  one of those lints whose reason does not hold.
- For each way the input can be wrong, what does the user of Python see?
  Is it the right exception, and does the message say which file, which
  line, which field, which value?
- Errors: the one `#[non_exhaustive]` enum of the core crate, to which
  each module adds its cases, with the cases named for what was being
  done and no dependency's error in a public one.
- A `Result` that is dropped, an error turned into a default, a record
  skipped without a word.

## api

The section "Types, names and defaults" of the coding skill, and the doc
comments as `.claude/skills/writing/SKILL.md` asks for them.

- Read every new name as someone who has not seen the code: does it say
  what the value is? Is the same thing called the same in Rust, in the
  binding and in Python, as pyNei calls it and as `docs/glossary.md`
  names it?
- Signatures and fields against `docs/architecture.md` and the spec.
- Defaults: a named constant with its source, or hidden in a constructor,
  a `Default`, an `Option` that becomes a value?
- `bool` parameters, strings for a finite set, a `_` arm on an enum of
  popnei, a `pub` that could be `pub(crate)`.
- Doc comments: what the item is, in the words of the domain, the units
  and the shape, `# Errors`. A comment that says what the next line does
  and not why is noise.
- A function longer than a screen or with many branches, the same logic
  in three places, dead code, a `TODO` with no issue.

## architecture

The section "What the architecture asks of the code" of the coding skill,
against `docs/architecture.md`.

- An allocation inside the loop over the variants: a `Vec`, a `String`, a
  `collect`, a `clone`, a `to_vec`.
- `Needs` asked for and `filled` checked.
- rayon and BLAS nested. The global pool of rayon built by the library.
- Code that needs threads and is not under
  `cfg(not(target_family = "wasm"))`, or has no serial version.
- Linear algebra outside `linalg`. pyo3 in the core crate.
- A new dependency: is it pure Rust, does it build for
  `wasm32-unknown-emscripten`, is it justified?
- A reader that takes a path where it could take `impl BufRead`.

## binding

`.claude/skills/coding/pyo3.md`, for the binding crate and the Python
package.

- Logic in the binding crate. A default or a docstring in `_core`.
- Names of pyo3 from an older version.
- Long work outside `py.detach`, and rayon used with the interpreter
  held. Something that is not `Send` moved into the closure by a trick.
  A loop over a dataset with no `check_signals`.
- Arrays: `unsafe` access to a `Bound<PyArray>`, a copy where the
  allocation could be handed over, a list where an array belongs, the
  dtype and the layout not made right on the Python side.
- A `#[pyclass]` that is not `frozen`, or is `unsendable`.
- `map_err` at the call sites in place of the one newtype. An `unwrap`.
- The Python package: calculation done in Python, a signature that
  differs from pyNei's with the difference not in the spec, missing type
  hints, a result that is not a frozen dataclass.
- The tests the binding needs: a wrong dtype, shape and layout give a
  `ValueError` that names the argument, a core error arrives as the right
  exception.
