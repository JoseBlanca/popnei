# Linear algebra and wasm

**Purpose.** The kinship, the PCA, the PCoA, the LD and the mixed models of the GWAS are bound by matrix products and decompositions, and there the library underneath decides the speed, not popnei's loops. And everything has to run in a browser tab, single threaded, in a 32 bit address space. This checklist is for those two.

**Triggers.** The scope touches the `linalg` module or a module that calls it, `kinship`, `pca`, `ld`, `gwas`, `dists`, the size or the layout of a block, or the wasm build.

**Skip when.** The scope is per variant work with no matrix and the question is not about wasm.

## What is already measured

`docs/rust_core.md`, on an M series Mac, 100000 variants x 1000 samples unless said: numpy on Accelerate against faer, z'z of 5000 x 1000 in 0.01 s against 0.17 s on one thread and 0.03 s on six; the symmetric eigendecomposition at n = 5000 in 6.35 s against 14.3 s and 4.57 s. Accelerate runs float64 products on the AMX units, so that gap is Apple's; faer against OpenBLAS or MKL on x86 has not been measured here. Under pyodide faer is the fast one: z'z in 0.58 s against 5.09 s for pyodide's numpy, which has no BLAS, and the parser loses 18% against native. Read those sections before proposing anything here, so that the review starts from them.

## Rules

- **Below the `linalg` module the cost is the library's.** A profile that puts the time inside `dgemm`, `dsyrk`, `dsyevd` or faer's kernels is not improved by changing the loop around the call. What can be changed is how often it is called, on how much data, and what is computed at all.
- **Call the product that does less.** z'z is symmetric: `syrk` computes half of what `gemm` does. A kinship accumulated block by block needs the products and never the whole z. A decomposition of which only the first components are used does not need all of them. These are algorithm findings and they outrank any tuning.
- **The block is the unit of the product.** A product over a block of a few thousand variants runs 5x to 10x faster than the same work variant by variant (`docs/architecture.md` section 2), and a block that does not fit in the cache of the BLAS kernel loses part of that. The block size is a parameter to measure, with the dataset stated, not to reason about.
- **The conversion to float is part of the cost.** BLAS needs `f64` and the genotypes are `i8`: the dosage matrix of a block is 8 bytes per genotype written, 40 MB for 5 million genotypes. Count how many times a block is converted. The same block converted once for the kinship and again for the PCA is a finding.
- **One pool at a time.** rayon for the records, the BLAS pool for the products, never nested: a rayon worker that calls BLAS pins it to one thread, and a big product is called from outside rayon. Accelerate, OpenBLAS and faer each bring a pool of their own, and the oversubscription playbook of `concurrency.md` applies to them as it does to two rayon pools. faer's `rayon` feature starts faer's own use of the pool inside ours.
- **Packed genotypes are an experiment, not a finding.** `docs/rust_core.md` section 2.2: in numpy, 2 bit packing lost to BLAS on floats by 4x to 20x, and in compiled code the packed path is the one plink uses. `docs/architecture.md` leaves it as an option behind the same row views. A proposal to pack comes with the fused pass that uses it and its benchmark against the `i8` layout.
- **Two backends, two measurements.** A change in `linalg` is measured on BLAS natively and on faer, natively too, since faer builds there and the wasm path can be timed without a browser. A change that helps one and hurts the other says which one popnei's users are on.

## wasm

- **Measure in pyodide when the claim is about the browser.** Native faer on one thread is a stand-in for the algorithm, not for the target: wasm has no SIMD wider than 128 bits, a 32 bit heap, and its own allocator. The steps to build and run the wasm wheel are in the repository once the walking skeleton leaves them there, and until then in `spike/README.md` of pyNei.
- **One thread.** Everything parallel has a serial version under `cfg(target_family = "wasm")`, and that is the one that runs. A finding about rayon says nothing about the browser, and a serial path that allocates or copies to feed a parallel one that is not there is a finding.
- **Memory is the limit before time is.** A tab has at most 4 GB of address space, in practice less, and growing the heap copies it. `docs/rust_core.md` saw it at n = 3000: the eigendecomposition takes 2.7x its native time at n = 2000 and 4.6x at n = 3000. Peak memory is a first class metric here: state it for the largest dataset the objectives name, 10000 samples, where a samples x samples `f64` matrix is 800 MB and two of them do not fit.
- **`usize` is 32 bits.** A count or an offset that overflows there is a correctness finding; route it to the code review.
