# The Python boundary

**Purpose.** A user of popnei calls Python. What they wait for is the Rust work plus everything that happens on the way in and out: conversions, copies, the interpreter held while Rust works, pandas objects built from the results. This checklist is for that path. The rules for writing the binding are in `.claude/skills/coding/pyo3.md`.

**Triggers.** The scope includes the binding crate `crates/popnei-python`, the Python package, or a timing taken from Python.

**Skip when.** The scope is the core crate alone and every timing is a cargo benchmark.

## Rules

- **Time the release build.** `maturin develop` builds without optimization. A timing from Python is taken after `maturin develop --release`, and the report says so. A comparison of a debug popnei against pyNei's numpy is not a finding, it is a mistake.
- **Measure at the boundary the user sees, and then inside it.** Time the Python call, and time the same work as a cargo benchmark of the core function. The difference is the boundary. When it is small against the whole, stop looking at it.
- **Copies on the way in.** `numpy.ascontiguousarray(x, dtype=numpy.int8)` copies when `x` is not already contiguous `int8`: a slice with a step, a transposed array, an `int64` array of genotypes. One copy of the genotypes of a chunk is about 5 MB and cheap; a copy of every chunk of a million variants, made because an upstream step produced the wrong dtype, is a finding. The count of bytes copied is deterministic, so gate on it.
- **Copies on the way out.** `into_pyarray` hands the Rust allocation to numpy. `to_pyarray`, `from_slice` and building a Python list from a `Vec` copy, the last one element by element. A result with one value per variant or per genotype leaves as an array.
- **Lists are converted element by element.** A `Vec<String>` of sample names is fine, once. A `Vec<f64>` or a list of lists with one entry per variant is not: each element is a Python object and a conversion.
- **The interpreter is released around the work.** Long Rust work outside `py.detach` holds every other Python thread, including pyNei-style read ahead threads, for its whole length. With rayon inside, it can deadlock. A sampling profile of the Python process shows this as one busy thread and the rest waiting.
- **Calls per chunk, not per variant.** `docs/architecture.md` section 5: Python never reads one variant at a time. A loop in the Python package that crosses into Rust once per variant, or once per sample, pays the call overhead a million times. The number of crossings is countable: gate on it.
- **Building the results.** A pandas `DataFrame` built from a dict of numpy arrays is cheap. One built row by row, or from Python lists, or re-indexed and sorted after the fact, can cost more than the calculation. Time the construction of the result object separately when the Rust part is already fast.
- **Compare like with like against pyNei.** Same input file, same chunk size where it applies, same number of threads, both after a warm run so that the file is in the page cache, and the same thing timed: pyNei's numbers in `docs/rust_core.md` say what was included.
