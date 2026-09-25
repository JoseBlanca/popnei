"""Every number that every pass over the blocks gives, written at
seventeen digits, so that a change made for speed can be shown to have moved
none of them.

A change that is kept for its wall time has to leave the results as they
were, and "as they were" means byte for byte here and not within a
tolerance: the passes below add floats, and a change that reorders an
addition moves the last bits without moving anything a tolerance would
catch. So this writes the numbers of a commit to a file, the same command is
run on the commit before it, and the two files are compared with `cmp`. It
is run at one thread and at the threads of the machine, since a pass whose
total depends on how many threads read a block would differ between the two
and not between the commits.

    RAYON_NUM_THREADS=1 uv run python the_numbers_of_every_pass.py \
        <path> <populations file> <output>

The path is a vars file or a VCF, plain or gzipped, as every other timing
script of this directory reads them: a path that ends in `.vars` is opened
with `open_vars` and anything else with `open_vcf`. The populations file
holds one line for each individual, its name, a tab and the name of its
population, which is what `make_pops.py` writes.

The passes it runs: the kinship, the principal components of the variants
with the weights of ten of them, which is the one calculation that reads the
source twice, the Kosman distance of every pair of individuals, the
distances between populations, the five statistics of a variant with the
populations and without them, the two counts of every individual, and the
diversity of every population with every statistic. Each of them reads the
whole file.

Two passes over blocks are not here. The association study writes its own
three columns, which `crates/popnei/benches/gwas.rs` does with
`--write-answers`. The r² matrix holds one value for each pair of the
variants it was given and reads whole blocks, so over this panel, whose
blocks are 5000 variants, the smallest matrix it can be asked for is 25
million values: a file of them would be larger than the panel itself, and
the pass is checked by the cargo and pytest tests instead.

It writes one line for each number, `pass field index value`, the value with
the seventeen digits that tell two `f64` apart. Nothing of it is a timing:
it opens the file as many times as there are passes and the passes run one
after another.
"""

import dataclasses
import sys

import numpy
import pandas
import popnei


def the_pops(path: str) -> dict[str, list[str]]:
    """The populations of the file at `path`, as the dict of population name
    to the names of its individuals that the calculations take.

    It is `the_pops` of `time_diversity.py`, beside this file, and reads the
    same files: one line for each individual, its name, a tab and the name of
    its population. The populations come out in the order in which the file
    first names each of them.
    """
    pops: dict[str, list[str]] = {}
    with open(path) as fhand:
        for number, line in enumerate(fhand, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise ValueError(
                    f"{path}, line {number}: a line of a populations file is "
                    f"the name of an individual, a tab and the name of its "
                    f"population, and this one is {line.rstrip(chr(10))!r}"
                )
            pops.setdefault(fields[1], []).append(fields[0])
    if not pops:
        raise ValueError(f"{path}: the populations file names no individual")
    return pops


def open_the_file(path: str) -> popnei.Variants:
    """The variants of the file at `path`, opened as the other timing scripts
    of this directory open theirs."""
    if path.lower().endswith(".vars"):
        return popnei.open_vars(path)
    return popnei.open_vcf(path)


def the_lines_of(name: str, value: object, field: str = "") -> list[str]:
    """One line for each number inside `value`, `name field index value`.

    It walks what the results of popnei are made of: a frame, a series, a
    numpy array, a dataclass, a dict, a list or a tuple, and a number. A
    value that is none of them is written as its `repr`, so that a field
    which stops being a number shows itself as a difference instead of
    passing unwritten.
    """
    if isinstance(value, pandas.DataFrame):
        return the_lines_of(name, value.to_numpy(), field)
    if isinstance(value, pandas.Series):
        return the_lines_of(name, value.to_numpy(), field)
    if isinstance(value, numpy.ndarray):
        flat = value.reshape(-1)
        return [
            f"{name} {field} {index} {one!r}"
            if not isinstance(one, (float, numpy.floating))
            else f"{name} {field} {index} {float(one):.17g}"
            for index, one in enumerate(flat.tolist())
        ]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        lines = []
        for one in dataclasses.fields(value):
            inside = getattr(value, one.name)
            lines.extend(the_lines_of(name, inside, f"{field}.{one.name}"))
        return lines
    if isinstance(value, dict):
        lines = []
        for key in value:
            lines.extend(the_lines_of(name, value[key], f"{field}[{key}]"))
        return lines
    if isinstance(value, (list, tuple)):
        lines = []
        for index, one in enumerate(value):
            lines.extend(the_lines_of(name, one, f"{field}[{index}]"))
        return lines
    if isinstance(value, float):
        return [f"{name} {field} . {value:.17g}"]
    return [f"{name} {field} . {value!r}"]


def the_numbers(path: str, pops_path: str) -> list[str]:
    """The numbers of every pass over the file at `path`, in the order the
    docstring of this file lists them."""
    pops = the_pops(pops_path)
    lines: list[str] = []

    lines.extend(the_lines_of("kinship", popnei.calc_kinship(open_the_file(path))))
    lines.extend(
        the_lines_of(
            "pca",
            popnei.do_pca_from_variants(open_the_file(path), num_prin_comps=10),
        )
    )
    lines.extend(
        the_lines_of("kosman", popnei.calc_pairwise_kosman_dists(open_the_file(path)))
    )
    lines.extend(
        the_lines_of(
            "pop_dists",
            popnei.calc_pop_dists(open_the_file(path), pops, jackknife_group=1000),
        )
    )
    lines.extend(
        the_lines_of("per_var", popnei.calc_per_var_distribs(open_the_file(path)))
    )
    lines.extend(
        the_lines_of(
            "per_var_pops",
            popnei.calc_per_var_distribs(open_the_file(path), pops=pops),
        )
    )
    lines.extend(
        the_lines_of(
            "per_individual", popnei.calc_per_individual_stats(open_the_file(path))
        )
    )
    lines.extend(
        the_lines_of(
            "diversity",
            popnei.calc_pop_diversity(
                open_the_file(path),
                pops,
                stats=tuple(popnei.PopDiversityStat),
                num_called_alleles=200,
            ),
        )
    )
    return lines


def main() -> None:
    """The three arguments, or what is wrong with them."""
    if len(sys.argv) != 4:
        raise SystemExit(
            "the arguments are the path of a vars file or a VCF, the path of "
            "a populations file and the path to write the numbers to"
        )
    path, pops_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
    lines = the_numbers(path, pops_path)
    with open(out, "w") as fhand:
        fhand.writelines(f"{line}\n" for line in lines)
    print(f"{out}: {len(lines)} numbers of the passes over {path}")


if __name__ == "__main__":
    main()
