# How the Kosman distances of a block are computed: the trial and the reference

21 September 2026. The programs behind two things that
`docs/specs/dists.md` states: the times of three ways of computing the
Kosman distances of every pair of individuals over one block, on which the
owner decided that popnei computes them with sets of bits, and the
comparison of pyNei, and of the paper's formula for other ploidies, with
an R function, which is where the literals of that spec come from. The Kosman distance of two diploid individuals at a
variant is 0 for the same genotype, 1 for no allele in common and 0.5
otherwise, and a block is 5 million genotypes, the unit the variants flow
in. The spec has the tables and what was concluded; this page says what
each file is and how to run it. Nothing here is part of popnei's build or
of its tests.

The machine was the owner's Apple M5 Pro, 18 cores, macOS 27.0, with
rustc 1.98.0, node 26.8.2, R 4.6.1, and the Python environment of popnei,
numpy 2.5.3 on Accelerate and pyNei at commit ef0ca6e. The load average
was 2 to 3 during the runs. A first set of runs was taken while R
packages compiled on the same machine and came out about a quarter
slower; the spec has the second set.

## The trial of the three ways

`kspike/` is a crate of its own, outside the workspace of popnei. Its
`src/lib.rs` makes a block at random, each variant with its own allele
frequency and 3 in 100 genotypes missing, and computes, for every pair,
twice the sum of the distances and the number of variants called in both,
in four ways: `bits`, sets of bits and counts of ones, and `bits_par`, the
same with the pairs on the threads of rayon, the Rust library popnei runs
its parallel loops on; `matmul`, the products of matrices of 0 and 1 in
float32 that pyNei makes, with faer 0.22.6, the linear algebra library in
pure Rust that popnei uses in wasm; and `naive`, a loop over the pairs and
the variants on the int8 genotypes. The binary takes the variants, the
individuals and the alleles of the block, times the first three ways, and
the fourth too when `naive` is written after them, which takes about a
second for the first block and much longer for many individuals. It
asserts that all the ways it ran give the same integers.

    cd kspike
    cargo build --release
    ./target/release/native 5000 1000 2 naive
    ./target/release/native 5000 1000 4
    ./target/release/native 500 10000 2

In wasm, under node, with the function `run` that the library exports:

    cargo build --release --lib --target wasm32-unknown-unknown
    node ../run_wasm.mjs target/wasm32-unknown-unknown/release/kspike.wasm

`run_wasm.mjs` takes the time of making the block out of each time. The
wasm was built without the SIMD instructions of wasm.

pyNei on blocks of the same sizes, and its whole calculation on 100000
variants x 1000 individuals from memory, from the root of popnei. Here
numpy makes the products with Accelerate, the BLAS of macOS:

    uv run python docs/reports/kosman-method/time_pynei.py
    uv run python docs/reports/kosman-method/time_full.py

numpy makes its own random blocks, with the same frequencies and the same
missing rate, and not the ones of the crate.

## The reference in R

`gd.kosman` of the R package PopGenReport 3.1.3 computes the Kosman
distance of every pair. The package did not install, because its
dependency terra needs the GDAL library, so the function is read from the
source package with `source()`. It needs adegenet; `ref.R` also loads
poppr and mmod, to show that their `diss.dist` and `dist.codom` treat
missing genotypes in another way and are not references. From the root of
popnei, with `S` a working directory:

    Rscript -e 'install.packages(c("adegenet", "poppr", "mmod"), lib="'$S'/rlib")'
    Rscript -e 'download.packages("PopGenReport", destdir="'$S'", type="source")'
    mkdir $S/pgr && tar -xzf $S/PopGenReport_3.1.3.tar.gz -C $S/pgr
    uv run python docs/reports/kosman-method/ref_export.py $S
    S=$S Rscript docs/reports/kosman-method/ref.R
    uv run python docs/reports/kosman-method/literals.py $S

`ref_export.py` writes three sets of genotypes as csv files for R, and
pyNei's distances for the last two: the worked example of the spec; the
panel of 200 individuals and 1200 variants, which it reads from
`sim_missing.vars` of pyNei's `test/gwas_reference/` at
`/Users/jose/devel/pynei`, a path written in the script; and 300 variants
of 40 individuals with four alleles from `numpy.random.default_rng(3)`.
`ref.R` runs `gd.kosman` on the three and prints the largest difference
from pyNei, 2.8e-17 on the panel and 0 on the other. `literals.py` prints,
for five pairs of the panel, twice the sum of the distances and the number
of variants from a loop in Python, beside the value of R.

`poly_export.py` and `poly.R` are the same for the other ploidies, added
on 22 September 2026 when the owner asked for them: the first writes a
tetraploid and a haploid dataset of 200 variants and 12 individuals, and
the two worked examples of the spec for those ploidies, with the
distances of a loop in Python that pairs the alleles as formula 2 of the
paper says; the second runs `gd.kosman` on them and prints the largest
difference, 5.6e-17 and 1.1e-16, and the first pairs.

    uv run python docs/reports/kosman-method/poly_export.py $S
    S=$S Rscript docs/reports/kosman-method/poly.R

`odd.py` runs pyNei on the cases of the spec where the result is not what
the formula would make a reader expect: a pair with no variant in common,
`min_num_snps`, the fewest variants a pair needs to get a distance, the
ploidies, no variants, one individual, and chunks of another size. It
takes no argument: `uv run python docs/reports/kosman-method/odd.py`.
