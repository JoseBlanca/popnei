# Task B: a reply in chat

You have spent the session measuring linear algebra libraries for popnei.
The owner's message, which you now have to answer, was:

> Do we use faer everywhere, so that there is one linear algebra
> implementation, or BLAS natively and faer only in wasm? I'd rather have
> one code path if it doesn't cost much.

Below are your session notes. They are all you know. Do not read
`docs/rust_core.md` or `docs/architecture.md`, which already hold the
answer. Write the reply you would send in chat, as markdown, to the output
path you were given.

## Session notes

```
linalg bakeoff -- notes

setup: the Mac (M series, 18 cores). numpy here links Accelerate.
arm A = numpy/Accelerate, arm B = faer (pure rust), via the pyo3 spike crate.
all float64. z = dosage matrix, variants x samples.

LA-1  z'z 5000x1000
   A 0.01 s | B 1 thread 0.17 s | B 6 threads 0.03 s

LA-2  eigh sweep (symmetric eigendecomposition, n = samples)
   n=1000  A 0.05  | B1 0.13 | B6 0.11
   n=2000  A 0.33  | B1 0.90 | B6 0.44
   n=5000  A 6.35  | B1 14.3 | B6 4.57
   -> correctness: z'z bit exact vs numpy, eigenvalues agree to 1e-14

note: Accelerate runs f64 matmul on the AMX units. that's why LA-1 is 17x.
apple only. faer's OWN published benchmarks say level w/ OpenBLAS and MKL on
x86 avx2/avx512. NOT measured by us. no x86 box at hand today.

pyo-run (pyodide 314.0.7, emscripten 5.0.3, 1 thread, same data)
   z'z 5000x1000   native B 0.169 | pyodide B 0.58 | pyodide numpy 5.09
   eigh 1000       0.124 | 0.31 | 0.69
   eigh 2000       0.90  | 2.39 | 5.57
   eigh 3000       2.97  | 13.6 | 19.2   (32 bit heap starts to hurt)
   pyodide numpy has no BLAS, reference LAPACK.
   there is no BLAS to link to in wasm anyway + objective says no C dependency there.

lost 40 min: uv venv --python 3.14 gives the free threaded python, pyodide-build
then puts sysconfigdata under python3.14t. fixed by asking for the normal build.
also faer "rayon" feature doesn't build on emscripten (spindle -> atomic-wait),
gated it with cfg(not(target_family = "wasm")).

what the core needs from linalg, as far as I can see from pyNei: z'z style
products (kinship, PCA), symmetric eigh (kinship, PCoA, mixed models),
a cholesky or solve for the mixed model null fits, maybe svd for PCA.
4-5 functions. if two backends: one small module, same signatures, picked
with cfg on the target. each function written twice, tests run on both natively
(faer builds natively too, so the wasm path can be tested with cargo test
without a browser).

where it matters: objectives say up to 10000 samples. eigh n=5000 one thread
6 s vs 14 s. kinship z'z over 1e6 variants x 1000 samples = 200 blocks of 5000:
A 2 s total, B1 34 s, B6 6 s.
mixed models in pyNei were level with GMMAT because of BLAS.

threads: rule in the design is rayon for records, BLAS for products, never nested.
with faer native multi thread we'd have faer's own rayon pool inside ours. needs care.

leaning: BLAS native + faer wasm. cost = the 4-5 functions twice. revisit if x86
numbers show faer level w/ OpenBLAS, then faer everywhere on linux might be fine,
but mac users still lose 2-17x.
```
