# I/O & syscalls checklist

**Purpose.** Surface I/O patterns that pay per-record syscall overhead, miss readahead, or hold large buffers in memory longer than necessary.

**Triggers.** `File::open`, `std::io::Read`, `std::io::Write`, `mmap`, a gzip decoder, an arrow IPC reader or writer, anything iterating over records in a file or stream.

**popnei.** The two readers are the VCF parser, text, often gzipped, parsed in parallel by batches, and the vars file, arrow IPC with zstd, one record batch per chunk. `docs/rust_core.md` has the baselines: 400 MB of VCF, 100000 variants x 1000 samples, in 0.55 s on one thread and 0.11 s with rayon, where the inflate is the ceiling for a gzipped one, and 0.17 s per 100 million genotypes from a vars file in pyNei.

**Skip when.** Code does no I/O. If dispatched anyway, write `No findings.`

## Rules

- **`File::read` and `File::write` without buffering are a per-call syscall.** Wrap in `BufReader::with_capacity(N, file)` / `BufWriter::with_capacity(N, file)`. Default capacity is 8 KiB; for large sequential reads (genomic files, log streams), bump to 64 KiB or 1 MiB and benchmark. Filed for any `File` used directly with line / record reads.
- **`BufWriter::drop` may swallow flush errors.** For correctness-critical writes, call `flush()` explicitly and propagate the error. Performance-relevant when the implicit flush at drop is the I/O cost.
- **Reading a whole file into memory is the fastest path — when the file fits.** `std::fs::read` / `read_to_string` performs one syscall and one allocation, sized to the file, and skips the per-call buffering machinery. Trade-off: peak RSS equals file size; not applicable to streaming workloads.
- **`mmap` wins for random access in large files.** A memory-mapped file lets the kernel page in only what is touched and avoids copies into user-space buffers. Caveats: SIGBUS on truncation, page-cache pressure, no safe cross-process semantics. Random-access patterns over multi-GB files using `seek` / `read` are filed Likely with the `mmap` (`memmap2` crate) experiment.
- **Compressed I/O: decompress once into a reusable buffer.** Reading a `.gz` line-by-line through a fresh decompressor per record is a finding; reuse a single `flate2::read::GzDecoder` wrapping a `BufReader`. Per-record decompression context creation is a hidden allocation. Cross-reference: `allocations.md`.
- **`fsync` is not free.** A `flush` on `BufWriter` does not call `fsync`. If durability is required, the right place is at end-of-stage, not per-record.
- **Endianness conversion is autovectorizable.** Loops doing `u32::from_le_bytes` per record vectorize when fed a `chunks_exact(4)` iterator. Manual byte-by-byte shifts do not. Cross-reference: `hot_loops.md` on `chunks_exact`.
- **Stdin / stdout: lock once.** `io::stdin()` returns a new lock per call; in a tight read loop, hoist `let stdin = io::stdin(); let mut h = stdin.lock();` and wrap in `BufReader`. Same for `stdout` in high-throughput writers.
