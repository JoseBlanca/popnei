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

import pytest

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


# What the process of the second test prints in each of the three cases: the
# pass was over after the interrupt, it gave the block that follows the one
# that was lost, and the interrupt arrived after the first block was given,
# where no block was lost and the pass goes on.
_THE_PASS_IS_OVER = "the pass is over"
_ANOTHER_BLOCK = "the pass gave the block after the one that was lost"
_AFTER_THE_BLOCK = "the interrupt arrived after the block"

# The process that is interrupted while the first of two blocks is read and
# then asks the pass for another block. The block the interrupt happened in
# is lost, so a pass that gave the next one would hand out the variants that
# follow the lost ones as if nothing had happened.
_ASK_FOR_A_BLOCK_AFTER_AN_INTERRUPT = f"""
import os
import signal
import sys
import threading
import time

import popnei

variants = popnei.open_vcf(sys.argv[1])
blocks = variants.iter_blocks(num_vars_per_block={_NUM_VARS // 2})
threading.Timer(
    {_SECONDS_BEFORE_THE_INTERRUPT}, lambda: os.kill(os.getpid(), signal.SIGINT)
).start()
the_block_was_given = False
try:
    next(blocks)
    the_block_was_given = True
    time.sleep(2)
except KeyboardInterrupt:
    pass
if the_block_was_given:
    print("{_AFTER_THE_BLOCK}")
    sys.exit(0)
try:
    next(blocks)
except StopIteration:
    print("{_THE_PASS_IS_OVER}")
    sys.exit(0)
print("{_ANOTHER_BLOCK}")
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


def test_a_pass_gives_no_block_after_the_interrupt_that_lost_one(tmp_path: Path):
    path = _vcf_of_many_variants(tmp_path / "many_variants.vcf")
    read = subprocess.run(
        [sys.executable, "-c", _ASK_FOR_A_BLOCK_AFTER_AN_INTERRUPT, str(path)],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert read.returncode == 0, (
        f"the process ended with {read.returncode}\nstdout: {read.stdout}\n"
        f"stderr: {read.stderr}"
    )
    what_happened = read.stdout.strip()
    if what_happened == _AFTER_THE_BLOCK:
        # The interrupt lost no block, so the pass was right to go on and
        # this run says nothing about the case. It has not happened in 20
        # runs: the read of half of the file takes 0.15 s and the interrupt
        # is sent 0.05 s after it starts.
        pytest.skip(_AFTER_THE_BLOCK)
    assert what_happened == _THE_PASS_IS_OVER, read.stdout
