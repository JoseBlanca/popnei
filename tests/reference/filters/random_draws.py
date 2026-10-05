"""It prints the variants that the filter that keeps variants at random keeps of many.vcf.

Run from the root of the repository:

    uv run --no-project python tests/reference/filters/random_draws.py

There is no reference program for a sample drawn with the rule of "The filter
that keeps variants at random" of docs/specs/filters.md, so this script is a
second version of that rule, written in Python for the spec: SplitMix64
started at the seed, one draw for each variant in the order of the file, the
top 53 bits of the draw divided by 2^53, and the variant kept when that number
is below the keep rate. It first checks the generator against the five draws
from a seed of 1234567 that Java's java.util.SplittableRandom gives, printed by
SplitMix.java beside this file, and stops when they differ.

It prints the ten numbers of the worked example of the spec, from a seed of 42,
and, for each row of the table of the spec, how many of the 500 variants of
tests/reference/vcf/many.vcf are kept and the positions of the first five. Every
line of the file is a variant, those that failed their FILTER among them, as
popnei reads the file with only_passed false.
"""

from pathlib import Path

MASK = (1 << 64) - 1
VCF = Path("tests/reference/vcf/many.vcf")
PUBLISHED = [
    6457827717110365317,
    3203168211198807973,
    9817491932198370423,
    4593380528125082431,
    16408922859458223821,
]


class SplitMix64:
    """The generator of Steele, Lea and Flood (2014), as splitmix64.c has it."""

    def __init__(self, seed: int) -> None:
        self.state = seed & MASK

    def next_draw(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def next_number(self) -> float:
        """A number from 0 to 1, below 1, from the top 53 bits of a draw."""
        return (self.next_draw() >> 11) * 2.0**-53


def main() -> None:
    generator = SplitMix64(1234567)
    draws = [generator.next_draw() for _ in PUBLISHED]
    if draws != PUBLISHED:
        raise SystemExit(f"the draws from 1234567 are {draws}, not {PUBLISHED}")

    generator = SplitMix64(42)
    numbers = [generator.next_number() for _ in range(10)]
    print("seed 42:", ", ".join(f"{number:.6f}" for number in numbers))
    print(
        "kept at 0.5:",
        [variant for variant, number in enumerate(numbers, start=1) if number < 0.5],
    )

    positions = [
        (line.split("\t", 2)[0], int(line.split("\t", 2)[1]))
        for line in VCF.read_text().splitlines()
        if not line.startswith("#")
    ]
    for keep_rate, seed in [(0.1, 42), (0.5, 42), (0.1, 7)]:
        generator = SplitMix64(seed)
        kept = [pos for pos in positions if generator.next_number() < keep_rate]
        print(f"keep rate {keep_rate}, seed {seed}: {len(kept)} of {len(positions)}")
        print("  the first five:", kept[:5])


if __name__ == "__main__":
    main()
