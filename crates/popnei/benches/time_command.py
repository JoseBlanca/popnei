"""It runs a command a number of times and prints each wall time and the
median.

It is the clock that task 4.1 of docs/plans/filters.md put on bcftools, so
that what the missing data filter costs bcftools is measured as it is
measured for popnei: the median of five runs of the whole command with the
filter, less the median of five runs of it without, the two run back to
back. `docs/reports/filters-measurement.md` has the numbers and the load
averages they were taken at. The two commands of that report, from the
directory that holds `big.vcf`:

    uv run --no-project python time_command.py 5 -- \\
        bcftools view -H -o /dev/null big.vcf
    uv run --no-project python time_command.py 5 -- \\
        bcftools view -H -i "F_MISSING<=0.1" -o /dev/null big.vcf

The records go to `/dev/null` in both, so that the destination costs the
same in each; `bcftools view` writes the records it keeps as text, so the
pass with the filter writes fewer of them and its difference is not the
filter alone. The report says what that does to the numbers.

One run that is not timed comes first, as the benchmark `filter_vars.rs`
has one, so that the timed runs read the file from the page cache. What the
command writes is the command's own business: this script redirects
nothing, and a command that writes to a terminal is timed writing to a
terminal. A run that exits with anything but 0 stops the script, which
exits with what the command exited with: a command that fails every time
would otherwise be timed and reported as if it had worked.

    uv run --no-project python time_command.py <runs> -- <command> [argument ...]
"""

import statistics
import subprocess
import sys
import time


def main() -> int:
    runs = int(sys.argv[1])
    command = sys.argv[sys.argv.index("--") + 1 :]
    print(" ".join(command), f", {runs} runs", sep="")
    started = time.perf_counter()
    done = subprocess.run(command, check=False)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"exit {done.returncode}"
    )
    if done.returncode != 0:
        return done.returncode
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        done = subprocess.run(command, check=False)
        took = time.perf_counter() - started
        if done.returncode != 0:
            print(f"run {run}: exit {done.returncode}")
            return done.returncode
        times.append(took)
        print(f"run {run}: {took:.3f} s")
    print(
        f"best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
