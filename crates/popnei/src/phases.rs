//! The two clocks of one turn of a pass over the blocks: how long the pass
//! was inside its reader and how long it was working on the block the
//! reader gave.
//!
//! A pass asks its reader for a block and then works on it, so the read of
//! the next block and the work on the one in hand never overlap unless
//! something puts them on two threads, which is what
//! [`with_one_block_ahead`](crate::block::with_one_block_ahead) does. What
//! that reader can save a pass is the smaller of the two clocks below, per
//! block, and a sampling profile cannot give either: the reader
//! decompresses on the thread that then computes, so the samples of the
//! read and the samples of the work are of one thread and one stack. This
//! is what says whether a pass is worth the reader's thread and its second
//! block of memory, and `docs/reports/perf-read-ahead-2026-09-25.md` has
//! what it said for each of them.
//!
//! Nothing here is compiled into a build popnei ships: the clocks are
//! behind the cargo feature `bench-phases`, which `crates/popnei/Cargo.toml`
//! describes, and with the feature off [`timed`] is the work it is given
//! and nothing else. The benchmarks under `crates/popnei/benches/` read
//! them with [`taken`].
//!
//! The `gwas` module has three clocks of its own, `gwas::phases`, which cut
//! the work of its pass into the dosages of a block and the test of its
//! variants. They are read by the benchmark of that module alone and never
//! together with these two.

/// Which of the two phases of one turn of a pass a clock is of.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Phase {
    /// Inside [`BlockReader::next_block`](crate::block::BlockReader::next_block)
    /// of the chain of readers, which is the read of the file and the
    /// decompression of a block.
    NextBlock,
    /// Everything the pass does with the block that call gave it.
    Work,
}

/// What `work` gives, with how long it took added to the clock of `phase`.
///
/// Without the cargo feature `bench-phases` nothing is timed and this is
/// `work()` and nothing else. With it, one `Instant::now` and one
/// `fetch_add` per call, which is twice per block of the pass.
pub(crate) fn timed<T>(phase: Phase, work: impl FnOnce() -> T) -> T {
    #[cfg(not(feature = "bench-phases"))]
    {
        let _ = phase;
        work()
    }
    #[cfg(feature = "bench-phases")]
    {
        let started = std::time::Instant::now();
        let value = work();
        add(phase, started.elapsed());
        value
    }
}

/// The nanoseconds spent inside the reader since the clocks were last
/// taken. A pass runs on one thread, so this and the one below are counters
/// and not a point of contention; they are atomics because a static that is
/// written needs to be one.
#[cfg(feature = "bench-phases")]
static NEXT_BLOCK: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// The nanoseconds spent working on the blocks since the clocks were last
/// taken.
#[cfg(feature = "bench-phases")]
static WORK: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// How long each of the two phases of a pass took.
#[cfg(feature = "bench-phases")]
#[derive(Debug, Clone, Copy)]
pub struct Phases {
    /// Inside `next_block` of the chain of readers.
    pub next_block: std::time::Duration,
    /// Working on the blocks it gave.
    pub work: std::time::Duration,
}

/// The two clocks since they were last taken, and both of them back to
/// zero, so that the next pass is timed on its own.
///
/// A calculation that makes two passes over the blocks, which the principal
/// components of the variants of a reader do, adds the clocks of both into
/// these two: what they say is how long that calculation spent reading and
/// how long working, over every pass it made.
#[cfg(feature = "bench-phases")]
#[must_use]
pub fn taken() -> Phases {
    use std::sync::atomic::Ordering;
    Phases {
        next_block: std::time::Duration::from_nanos(NEXT_BLOCK.swap(0, Ordering::Relaxed)),
        work: std::time::Duration::from_nanos(WORK.swap(0, Ordering::Relaxed)),
    }
}

/// `took` added to the clock of `phase`. A time longer than 584 years
/// saturates, which no phase of a pass reaches.
#[cfg(feature = "bench-phases")]
fn add(phase: Phase, took: std::time::Duration) {
    use std::sync::atomic::Ordering;
    let nanos = u64::try_from(took.as_nanos()).unwrap_or(u64::MAX);
    let clock = match phase {
        Phase::NextBlock => &NEXT_BLOCK,
        Phase::Work => &WORK,
    };
    clock.fetch_add(nanos, Ordering::Relaxed);
}
