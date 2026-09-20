# The format of a finding and its severity

Every reviewer of a performance review reads this once before writing. It is ported from pop_var_caller.

## Severity

A performance finding is not a defect. Its severity says how confident we are that the change gives a gain that can be measured and is worth its complexity.

- **Hot-path**: strong evidence that the site is on a hot path, it is named in a profile, at the top of a flamegraph, or inside a function with a reproducible benchmark, and the change has a clear mechanism for a gain. The site accounts for a meaningful share of the time, the allocations, the lock waits or the syscalls in the measurement that is quoted.
- **Likely**: an improvement matched by pattern, with a plausible call frequency and no profile yet. Worth running its measurement before anything is changed. Most findings start here.
- **Speculative**: the pattern is there and it is unclear that the site is hot or the gain real. Filed so that it is known. Not acted on without an experiment that contradicts "this will not matter".
- **Note**: a pointer for future code, or a cold path mentioned for completeness. Not actionable as written.

Hot-path needs high confidence and the quoted output of a profile or a benchmark that names the site. The benchmark is one that can be run again: it is in `benches/`, or its command is in the report. A timing whose harness was not kept is an anecdote, and a finding that rests on it is Likely. A finding from a pattern alone is at most Likely. Code that runs once, at start up, in the handling of an error, in the building of a result object for a small result, is at most a Note unless the orchestrator asked for it.

When the same pattern appears more than about five times in the same form, report it once with a recommendation for a sweep, not five times.

## Each finding

```
- `path/to/file.rs:LINE` — [Severity] Title
- Confidence: high / medium / low
- Hot-path evidence: the quoted profile or benchmark, or "pattern only"
- Pattern: which rule of the checklist it matches
- Mechanism: one to three sentences that name the cost, an allocation, a cache miss, a lock wait, a syscall, a mispredicted branch, a copy across the Python boundary, a conversion to f64, and why the change removes it
- Measurement plan: the benchmark or the profile that would confirm or refute the gain: the command, the dataset, the metric, and the threshold that makes the change worth keeping. When the mechanism is something that can be counted, allocations, syscalls, bytes copied, crossings into Rust, gate on the count, which is the same on every run, and use the wall time as the second check. Name only tools that `profiling_environment.md` lists as present.
- Effect on the numbers: none, or what changes in the result and how it stays within the tolerance of the spec and the same at every number of threads (`checklists/numbers.md`)
- Complexity cost: what the change adds, a type, a lifetime, a dependency, an `#[expect]` of a lint, a build flag, an invariant somebody has to keep
- Suggested experiment: a diff, a replacement, or numbered steps, complete enough to run
```

Group the findings by severity, Hot-path first. With nothing to report, write "No findings." and stop.

## What belongs to another category

Do not file it. Put one line at the end, under "For another category": `file:line`, what it looks like, and which category. A lock taken inside a tight loop is `concurrency` even when the call shape suggests a clone. False sharing is `data_layout` although concurrent state triggers it. A `format!` in a hot serializer is `allocations` when the cost is the allocation and `hot_loops` when it is the formatting.

## What is never done

- A path or a line that was not read is not cited.
- The output of a tool is quoted as it was printed. A benchmark that was not run is not described as if it had been: say "pattern only".
- No speed-up is estimated as a percentage or a multiple unless the experiment was run. Say what is removed, "one allocation per variant", "the lock out of the inner loop".
- Cold code is not a hot path.
- A finding made by reading a file outside the scope is probably out of scope: note it and stop.
