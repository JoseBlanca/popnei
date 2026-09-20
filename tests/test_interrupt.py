"""A Ctrl-C while a block is read: what the session of the user gets.

`.claude/skills/coding/pyo3.md` asks, under "Tests", that a long call can be
interrupted. A block of a big VCF takes seconds, so the SIGINT of a Ctrl-C
usually arrives while the core is reading one, with the interpreter released
and no bytecode running to raise it. What the user must get then is the
`KeyboardInterrupt` of any other interrupted call, and not the
`PanicException` of a panic in Rust, which derives from `BaseException` and
ends the session.

The case is run in a process of its own because it depends on being the
first numpy array of the process: the pending interrupt is raised inside the
import of the numpy C API, and the array that follows it in the same process
finds that API already loaded.
"""

import subprocess
import sys
from pathlib import Path

# How many variants the VCF of the test holds, and how long the core takes
# to read them as one block: 400000 lines of three individuals are 17 MB and
# 0.3 s with the core that `maturin develop` builds, which is unoptimised.
# The interrupt is sent 0.05 s after the read starts, so it lands inside it.
_NUM_VARS = 400_000
_SECONDS_BEFORE_THE_INTERRUPT = 0.05

# What the process prints when it got the exception a user expects.
_INTERRUPTED = "KeyboardInterrupt"

# The process that reads one block and sends itself the SIGINT of a Ctrl-C
# while it is being read. A signal that arrives after the block was given is
# raised by the sleep, and is the same outcome for a user; the failure this
# guards against is a `PanicException`, which no `except KeyboardInterrupt`
# catches and which leaves the process with another exit code.
_READ_A_BLOCK_AND_INTERRUPT_IT = f"""
import os
import signal
import sys
import threading
import time

import popnei

variants = popnei.open_vcf(sys.argv[1])
blocks = variants.iter_blocks(num_vars_per_block={_NUM_VARS})
threading.Timer(
    {_SECONDS_BEFORE_THE_INTERRUPT}, lambda: os.kill(os.getpid(), signal.SIGINT)
).start()
try:
    next(blocks)
    time.sleep(2)
except KeyboardInterrupt:
    print("{_INTERRUPTED}")
    sys.exit(0)
print("the interrupt was never raised")
sys.exit(2)
"""


def _vcf_of_many_variants(path: Path) -> Path:
    """A VCF of three individuals and `_NUM_VARS` variants, at `path`."""
    header = (
        "##fileformat=VCFv4.4",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3",
    )
    lines = (
        f"chr1\t{variant + 1}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
        for variant in range(_NUM_VARS)
    )
    path.write_text("\n".join((*header, *lines)) + "\n")
    return path


def test_a_ctrl_c_while_a_block_is_read_raises_keyboard_interrupt(tmp_path: Path):
    path = _vcf_of_many_variants(tmp_path / "many_variants.vcf")
    read = subprocess.run(
        [sys.executable, "-c", _READ_A_BLOCK_AND_INTERRUPT_IT, str(path)],
        capture_output=True,
        text=True,
        timeout=120,
        # What the process ended with is what the test reads, so a code
        # other than 0 is an assertion of this test and not an exception of
        # `subprocess`, whose message would not say what was expected.
        check=False,
    )
    assert read.returncode == 0, (
        f"the process ended with {read.returncode} and not with the "
        f"KeyboardInterrupt of a Ctrl-C\nstdout: {read.stdout}\n"
        f"stderr: {read.stderr}"
    )
    assert read.stdout.strip() == _INTERRUPTED, read.stdout
