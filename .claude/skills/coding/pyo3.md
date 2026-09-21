# The binding crate: pyo3, numpy and maturin

Read this before touching `crates/popnei-python`. It was written in
September 2026 from the guide of pyo3 0.29.2, the docs of the numpy crate
0.29.0, the maturin and pyodide-build guides, and the spike that built for
pyodide, `spike/pynei_spike` in the pyNei repository. No binding code of
popnei existed yet, so what the walking skeleton finds to be different is
corrected here.

The names of pyo3 change from version to version, and most of what a model
remembers is from older ones. These are the current names: `Bound<'py, T>`
and `Py<T>` for references to Python objects, `Python::attach` where it
was `with_gil`, `py.detach` where it was `allow_threads`, `cast` where it
was `downcast`, `IntoPyObject` where it was `IntoPy` and `ToPyObject`,
`PyOnceLock` where it was `GILOnceCell`, `Py<PyAny>` where it was
`PyObject`. When the version of pyo3 moves, read its migration guide,
https://pyo3.rs/latest/migration.html, before writing anything.

## What the crate is

It translates and holds no logic. A function here checks its arguments,
gets slices or views out of the arrays, releases the interpreter, calls
one function of the core crate, and turns the result into arrays. An `if`
about genetics is in the wrong crate, and so is a loop over variants.

The native module is private, `popnei._core`: in `pyproject.toml`,
`module-name = "popnei._core"` and `python-source = "python"`; in the
`Cargo.toml` of the crate, `lib.name = "_core"`; and the `#[pymodule]` has
the same name. The Python package imports from `_core` and is the only
thing a user imports. The module is written in the declarative form,
`#[pymodule] mod _core { ... }` with `#[pymodule_export]`.

The core crate never depends on pyo3, and `cargo test -p popnei` runs with
no Python. The `extension-module` feature of pyo3 is deprecated and is not
used; maturin tells pyo3 what it needs. This is from the docs and has to
be confirmed on the walking skeleton, because the spike still used it.

## Signatures

The Python package gives the API its pyNei signatures, its defaults and
its docstrings. The functions of `_core` are plain: required positional
arguments, no defaults, so that a default exists in one place only, the
`pub const` of the core that the Python package reads.

A path comes in as `PathBuf`, which takes `str` and `pathlib.Path`. Sample
names come in as `Vec<String>`. Anything with one value per genotype, per
variant or per sample comes in as an array, never as a list: a list is
converted element by element.

## Arrays

- In: `PyReadonlyArray3<'py, i8>` for genotypes, and the like. Get the
  data with `as_slice()`, which fails when the array is not C contiguous,
  or with `as_array()`. Both are safe on a readonly array. The same
  methods on a `Bound<PyArray>` are `unsafe`, because Python can change
  the data meanwhile, and are not used.
- The Python package makes the array right before the call,
  `numpy.ascontiguousarray(gts, dtype=numpy.int8)`, and checks the number
  of dimensions, because numpy 0.29.0 reports a wrong dtype with a message
  that says "'ndarray' is not an instance of 'ndarray'". The binding crate
  still turns a failed `as_slice` into a `ValueError` that names the
  argument.
- Out: `vec.into_pyarray(py)` and `array.into_pyarray(py)` hand the Rust
  allocation to numpy without a copy. `to_pyarray` and `from_slice` copy,
  and are for small things only. The return type is
  `Bound<'py, PyArray2<f64>>`.
- An array that borrows memory owned by a Rust object,
  `PyArray::borrow_from_array`, is `unsafe` and needs the owner to outlive
  it and never reallocate. Not used unless a measurement shows the copy
  matters, and then with its `// SAFETY:` comment.
- `ndarray` is taken from `numpy::ndarray`, so that the binding crate and
  the numpy crate agree on its version.

## Releasing the interpreter

Work in Rust that lasts more than a few milliseconds runs inside
`py.detach(|| ...)`. It is required whenever the work uses rayon: a worker
thread that needs the interpreter while the caller holds it deadlocks.

What goes into the closure has to be `Send`. `Python`, `Bound` and
`PyReadonlyArray` are not, so take the slice or the view first and move
only that in:

```rust
let gts = gts.as_slice().map_err(|_| not_contiguous("gts"))?;
let result = py.detach(|| popnei::stats::exp_het(gts, ...))?;
```

Ctrl-C works only if the code asks: `py.check_signals()?` between blocks,
outside `detach`, in anything that runs over a whole dataset.

A module of pyo3 0.28 or later declares itself safe without the GIL by
default. That is a promise about our code. It holds while every
`#[pyclass]` is `frozen` with its state behind a `Mutex`, and there are no
`static mut` and no std `OnceLock` that calls into Python; use
`PyOnceLock` for that.

## Classes

A `#[pyclass]` cannot have lifetime or type parameters and has to be
`Send + Sync`. popnei needs few: the object that holds a reader or a
source of blocks for the Python `Variants`, and little else. Results are
not classes; they go out as arrays and tuples, and the Python package
builds the frozen dataclasses.

- `#[pyclass(frozen)]`, with what changes inside a `Mutex`:
  `reader: Mutex<Box<dyn VariantReader + Send>>`. A class that is not
  frozen gets a borrow check at run time, and two overlapping uses raise
  `RuntimeError: Already borrowed` in the user's session.
- Never `unsendable`, which panics when another thread touches the object.
- Python does not read one variant at a time, as section 5 of the
  architecture says. `__next__` gives a chunk: the genotypes as an array
  of variants x samples x ploidy and the columns that were asked for.

## Errors and panics

`impl From<popnei::Error> for PyErr` cannot be written in this crate,
because neither type is ours. So the crate has an error type of its own,
`enum PyPopneiError`, with a `From<popnei::Error>`, a `From<PyErr>` and one
`From<PyPopneiError> for PyErr`. Its eight cases are the error of the core;
that same error with the file it happened in, which the core was not given;
an argument that says how many of something there are and counts nothing,
which this crate refuses before the core sees it; the threshold of a filter
that is not a number from 0 to 1, under the name of the argument the user
wrote it in, which the core does not know, since it names a filter by its
kind, `maf`, and a user wrote `max_allowed_maf`; a path that a file is
already at, given to a call that writes one, which this crate also refuses
before the core sees it; another of those errors with the file that the
call was writing and could not take away afterwards, which becomes a note
on the exception, the text that Python keeps in `__notes__` and prints under
the message; a defect of this crate, a lock that a panic left broken or a
chromosome that is not in the table it came from; and an exception the
interpreter itself raised, the `KeyboardInterrupt` that `py.check_signals`
finds between two blocks, which travels back as it is.

Every function of the crate returns `Result<T, PyPopneiError>`, the
`#[pyfunction]` and the `#[pymethods]` that pyo3 exports among them, so
that `?` carries an error of the core across. A call site maps one by hand
only to add what the core does not have: `PyPopneiError::of_the_file` puts
the file into an error that came from reading one.

The one conversion chooses the exception by the owner's convention of 21
September 2026, which "Errors, and no panics" of `SKILL.md` gives: a
`ValueError` for a wrong input of a function, which a file whose content is
not what a VCF holds is, and which a case nobody has written yet gets; a
`RuntimeError` for a defect of popnei, which is a defect of this crate, one
of the three cases with which `docs/specs/block.md` says that a reader has
one, the number of values a filter gave `retain_vars`, or the parse of a
batch of lines that did not come back; and an `OSError` for a file that
cannot be read, that was cut short or that is corrupted. The `OSError`
is built with the number the system gave, so that it is the
`FileNotFoundError`, the `IsADirectoryError` or the `PermissionError` of
that number, and with `None` in its place when nothing of the system
refused anything, which leaves an `OSError` whose `errno` is `None`. Either
way it carries the file in `filename`.

What the message is depends on the exception, and a new case follows the
one it is:

- A `ValueError` and a `RuntimeError` of a file: the `Display` of the core
  error, which has the line, the column and the value, with the path before
  it, `<path>: <what the core says>`. The path is `to_string_lossy` there,
  since a message is text and a name of a file is bytes.
- An `OSError`: the `Display` of the core error with no path, because
  Python prints `filename` after the message and would say it twice. The
  file goes in as an `OsString`, which arrives as the text the standard
  library would give, so `error.filename` is the path the caller wrote.
- The two `OSError`s that wrap a `std::io::Error`, a file that could not be
  opened and a read that failed: this crate writes the message itself, "the
  file could not be opened: " and what the system said, without the number
  that Rust puts at the end of it, `(os error 2)`, which Python prints of
  its own in `[Errno 2]`.
- An argument that this crate or the core refuses, `fields`, `ploidy`,
  `num_vars_per_block`, names no file: what a user wrote is wrong whatever
  file is read.

A panic in Rust reaches Python as `PanicException`, which derives from
`BaseException`, is not caught by `except Exception`, and usually ends the
session. That is the reason the core does not panic, and it holds here
too: no `unwrap` on `as_slice_mut`, on a lock or on a conversion. A
poisoned `Mutex` is an error, not an `unwrap`.

## Two builds

Everything here builds twice: natively, and for
`wasm32-unknown-emscripten` as the wheel for pyodide, where there are no
threads. rayon is a dependency only off wasm,
`[target.'cfg(not(target_family = "wasm"))'.dependencies]`, and code that
uses it is under the same `cfg` with a serial version beside it. Nothing
is compiled or linked with `-pthread`, or the module will not load.

The wasm wheel is tied to one version of pyodide and its emscripten, and
the versions are read from `pyodide config get`, not written into a
script. The host Python of pyodide-build must not be the free threaded
one. The steps that worked are in `spike/README.md` of pyNei until the
walking skeleton puts them in this repository.

## Tests

The logic is tested in the core with `cargo test`. What this crate adds is
tested from Python with pytest, through the package: that an array of the
wrong dtype, shape or layout gives a `ValueError` that names the argument,
that an error of the core arrives as the right exception with its message,
that a big result comes back without a copy where that was the intent,
and that a long call can be interrupted.

`maturin develop` builds without optimization. Any timing is taken after
`maturin develop --release`.
