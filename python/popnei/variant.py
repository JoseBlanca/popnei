"""The handle a user holds: a source of variants, its individuals and the
steps that were put on it, and the counts of a pass over it."""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from popnei import _core
from popnei.block import Block, Field, _block_of
from popnei.filters import FilteringStats, Step


@dataclass(frozen=True)
class PassStats:
    """The counts of one pass over a source of variants.

    Every consumer of a :class:`Variants` gives them back with its result,
    and so does the iterator of :meth:`Variants.iter_blocks`. A pass is one
    reading of the source from its start, through the steps the ``Variants``
    had when it started, so these counts are of that reading alone and of
    what the ``Variants`` holds afterwards nothing reaches them.
    """

    num_vars: int
    """How many variants the consumer took, after the steps."""

    filtering: dict[str, FilteringStats]
    """How many variants each filter of the pass was given and kept, under
    the kind of the filter, ``"missing_data"``, ``"maf"`` or ``"obs_het"``,
    in the order of the steps. It is empty for a pass with no filter."""


def _pass_stats_of(counts) -> PassStats:
    """The counts that `popnei._core` gives for one pass.

    The core gives the filters of the chain of readers, the outermost first,
    and a user reads them in the order of the steps, which is the one the
    filters were put on the ``Variants`` in and the reverse of the chain's.
    """
    num_vars, filtering = counts
    return PassStats(
        num_vars=num_vars,
        filtering={
            kind: FilteringStats(vars_processed=vars_processed, vars_kept=vars_kept)
            for kind, vars_processed, vars_kept in reversed(filtering)
        },
    )


class Blocks:
    """The blocks of one pass, one after another, and the counts of it.

    It is what :meth:`Variants.iter_blocks` gives: an iterator of
    :class:`popnei.Block`, and a `pass_stats` that says how many variants
    have come out of it and what each filter of the pass has been given and
    kept.
    """

    def __init__(self, blocks: _core.Blocks):
        """The pass that `popnei._core` started, as the blocks of it."""
        self._blocks = blocks

    def __iter__(self) -> Blocks:
        """The iterator itself: a pass is read once, from its start."""
        return self

    def __next__(self) -> Block:
        """The next block of the pass.

        When a variant cannot be read, the error comes in the place of the
        block that would have held it, and every call after that one says
        that the pass is over.
        """
        return _block_of(next(self._blocks))

    @property
    def pass_stats(self) -> PassStats:
        """The counts of the pass as it stands.

        Read when the pass is over, they are of everything it gave. Read
        between two blocks, they are of the blocks that came out so far,
        which can be fewer variants than the filters of the pass have kept:
        the reader that cuts the blocks to the size the user asked for keeps
        the variants of the next block.
        """
        return _pass_stats_of(self._blocks.pass_stats())


class Variants:
    """A source of variants: a VCF with the options it is read with, or a
    vars file, and the steps that were put on it.

    It holds no genotypes. A user gets one from :func:`popnei.open_vcf` or
    from :func:`popnei.open_vars` and gives it to as many calculations as
    they want: each one opens the source again and runs its loop over the
    variants inside the Rust core, so the dataset is never in memory as a
    whole.

    What is done with it is of two kinds, and what a call gives back says
    which. A step, a filter of :mod:`popnei.filters`, is a method that adds
    itself to the list of steps, reads nothing and returns nothing, and
    :attr:`steps` is that list. A consumer, :meth:`iter_blocks`,
    :func:`popnei.write_vars` or the function of a calculation, gives
    something back, and it runs the steps: it makes as many passes over the
    source as it needs, each one built from the steps the ``Variants`` has
    when that pass starts. So a step added between two consumers holds for
    the second, and one added while a pass runs holds from the next pass.

    It is pyNei's ``Variants`` under the word of ``docs/glossary.md``: what
    pyNei calls a sample is here an individual, one organism that was
    genotyped. The genotypes come out of it through :meth:`iter_blocks` and
    through nothing else.

    It cannot be pickled, and ``copy.deepcopy`` of one is a ``TypeError``
    for the same reason: what it holds is an object of Rust with the path
    and the options of the source, which no pickle carries. What travels
    between processes is the path and the arguments of
    :func:`popnei.open_vcf` or of :func:`popnei.open_vars`, and a
    ``Variants`` is opened again at the other end. ``copy.copy`` gives a
    second handle over the same source, which reads the same variants and
    shares its steps: a filter put on either of them is on both, and a
    handle holds nothing of a pass, so the two are used as one is. A user
    who wants a second set of thresholds over one file opens it again,
    which reads the header and nothing else.
    """

    def __init__(self, source: _core.VcfSource | _core.VarsSource):
        """The handle over `source`, which :func:`popnei.open_vcf` and
        :func:`popnei.open_vars` build."""
        self._source = source
        # The names of the individuals of the source come from the header,
        # which was read once. The steps are given them: they resolve the
        # names of a filter of individuals against them, and they are what
        # says which individuals the next pass gives.
        self._steps = _core.Steps(list(source.individuals()))

    def __repr__(self) -> str:
        """The source the variants are read from, the options it is read
        with and the steps on it.

        A second filter of one kind is refused, so a user has to be able to
        see which are set, and a notebook whose cells were run out of order
        is where they most need it. The options are part of the source: two
        handles over one VCF that differ in them give different variants. A
        vars file is read with none, since what its genotypes hold is
        written in the file.
        """
        what = [f"of {self._source.path()}"]
        if isinstance(self._source, _core.VcfSource):
            what.append(f"ploidy={self._source.ploidy()}")
            what.append(f"only_passed={self._source.only_passed()}")
        steps = ", ".join(f"{step.kind}({_arguments_of(step)})" for step in self.steps)
        what.append(steps or "no steps")
        return f"<Variants {', '.join(what)}>"

    @property
    def individuals(self) -> tuple[str, ...]:
        """The names of the individuals the next pass gives, in its order.

        They are those of the source, in the order the source has them,
        until :meth:`filter_individuals` is put on the ``Variants``: from
        then on they are the ones that filter keeps, in the order they were
        named, which is the order of the genotypes of every block. A pass
        changes nothing of them, so they are the same read before one and
        after one.
        """
        return tuple(self._steps.individuals())

    @property
    def num_individuals(self) -> int:
        """How many individuals the next pass gives the genotypes of."""
        return len(self.individuals)

    @property
    def ploidy(self) -> int:
        """How many alleles the genotype of one individual holds."""
        return self._source.ploidy()

    @property
    def steps(self) -> tuple[Step, ...]:
        """The steps that were put on this ``Variants``, in their order.

        Each one is a :class:`popnei.Step` with the kind of the step and the
        arguments it was given. A pass takes the steps that are there when
        it starts, so this is what the next consumer will run.

        The tuple and the ``args`` dict of every step in it are built at
        each read, out of what the ``Variants`` holds, so writing into one
        of those dicts changes nothing of the steps: a step is added by one
        of the four filter methods and by nothing else.
        """
        return tuple(
            Step(kind=kind, args=dict(args)) for kind, args in self._steps.steps()
        )

    def filter_by_missing_data(self, max_allowed_missing_rate: float) -> None:
        """Keep the variants whose missing rate is at most
        `max_allowed_missing_rate`.

        The missing rate of a variant is its missing genotypes divided by
        all the individuals of the dataset, and not by the ones that were
        called at it. A genotype is missing when one of its alleles at least
        was not called, so ``0/.`` in a VCF is a missing genotype, as it is
        in pyNei and in bcftools.

        The call adds a step and gives nothing back. What runs it is the
        next pass over the source, which every consumer makes: a filter
        added between two of them holds for the second, and one added while
        a pass runs holds from the pass after it.

        `max_allowed_missing_rate` has no default, where pyNei's is 0.0,
        which keeps only the variants with every genotype called. What is
        no number, a string, ``None`` and a truth value among them, is a
        ``TypeError`` that names the argument and what was given, and so is
        a call with no threshold; a number that is not from 0 to 1, both
        included, NaN and a whole number too large for a float among them,
        is a ``ValueError`` that names the argument and the value. A second
        filter of this kind on the same ``Variants`` is a ``ValueError``
        too, with the threshold that is set: two thresholds of one kind
        keep what the stricter of them keeps alone, so the second says that
        the steps are not what their user thinks, which running the cell of
        a notebook twice gives. :attr:`steps` is what they hold. After any
        of them the steps are as they were.
        """
        self._steps.filter_by_missing_data(max_allowed_missing_rate)

    def filter_by_maf(self, max_allowed_maf: float) -> None:
        """Keep the variants whose major allele frequency is at most
        `max_allowed_maf`.

        The major allele frequency of a variant, "maf" in pyNei and in
        popnei, is the count of its commonest allele divided by its called
        alleles, where most of the literature and plink2 give those letters
        to the minor allele. Every allele of a multiallelic variant has its
        own count, and an allele is counted wherever it was called, in a
        half called genotype too. A filter at 0.95 takes out the variants
        that hardly vary among these individuals. A variant with no called
        allele has no major allele frequency and is not kept, whatever the
        threshold.

        It asks for no minimum of called data, as pyNei does not: a variant
        with one called genotype has the frequency of the alleles of that
        genotype. A user who does not want the variants that have little
        called data puts :meth:`filter_by_missing_data` before this one.

        The call adds a step and gives nothing back, and it refuses what
        :meth:`filter_by_missing_data` refuses: what is no number and a call
        with no threshold are a ``TypeError``, and a number that is not from
        0 to 1 and a second filter of this kind are a ``ValueError``.
        """
        self._steps.filter_by_maf(max_allowed_maf)

    def filter_by_obs_het(self, max_allowed_obs_het: float) -> None:
        """Keep the variants whose observed heterozygosity is at most
        `max_allowed_obs_het`.

        The observed heterozygosity of a variant is its heterozygous
        genotypes divided by its called ones, where a genotype is
        heterozygous when it is called and its alleles are not all the same,
        at any ploidy. It takes out the variants in which too many
        individuals are heterozygous, which in most datasets are paralogous
        regions read as one site. A variant with no called genotype has no
        observed heterozygosity and is not kept, whatever the threshold.

        It asks for no minimum of called data, as pyNei does not: a variant
        with one called genotype, heterozygous, has an observed
        heterozygosity of 1.

        The call adds a step and gives nothing back, and it refuses what
        :meth:`filter_by_missing_data` refuses: what is no number and a call
        with no threshold are a ``TypeError``, and a number that is not from
        0 to 1 and a second filter of this kind are a ``ValueError``.
        """
        self._steps.filter_by_obs_het(max_allowed_obs_het)

    def filter_individuals(self, individuals: Sequence[str]) -> None:
        """Keep the genotypes of `individuals` at every variant and drop
        those of the rest.

        Every variant stays: the step takes columns of the genotypes away
        and no row, so it has no entry in the counts of a pass. The
        individuals are kept in the order they are named here, which is the
        order of the genotypes of every block and of the rows of every
        result over individuals, so it is also the way to put a dataset's
        individuals in the order a user wants. pyNei's ``filter_samples``
        keeps them in the order of the source instead.

        A step of it is what every step that comes after it sees:
        :meth:`filter_by_missing_data` before the call divides by all the
        individuals of the source, and after it by the kept ones alone.
        :attr:`individuals` and :attr:`num_individuals` are the kept ones
        from the call on, since they are what the next pass gives.

        A name that is not an individual of the source is a ``ValueError``
        that names it, where pyNei drops it in silence and gives the
        individuals it did find; a name that is there twice is a
        ``ValueError`` too, since one individual is kept once; and so is a
        call with no name, because variants of nobody are no dataset. A
        second filter of individuals on the same ``Variants`` is a
        ``ValueError`` as well: two lists keep the individuals that are in
        both, which is one list, so the second says that the steps are not
        what their user thinks. A user who wants two sets of individuals
        over one file opens it twice. After any of them the steps are as
        they were.

        One name written as a string, ``filter_individuals("ind05")``, is a
        ``TypeError`` that says to write ``("ind05",)``: a string is a
        sequence of its letters and the call would ask for the individuals
        ``i``, ``n``, ``d``, ``0`` and ``5``. What is no sequence at all,
        and an element of it that is no name, are a ``TypeError`` that
        names `individuals` and what was given.
        """
        if isinstance(individuals, str):
            # A string is a sequence of its letters, and one name written
            # without its comma would ask for the individuals `i`, `n`,
            # `d`...
            raise TypeError(
                f"`individuals` is a sequence of names and not one name: write "
                f'individuals=("{individuals}",) for that one individual'
            )
        try:
            names = list(individuals)
        except TypeError:
            # What Python says of its own here, `'int' object is not
            # iterable`, names neither the argument nor the call.
            raise TypeError(
                f"`individuals` is a sequence of the names of the individuals "
                f"to keep, and {individuals!r}, a "
                f"{type(individuals).__name__}, was given"
            ) from None
        for name in names:
            if not isinstance(name, str):
                # pyo3 refuses it with `'int' object is not an instance of
                # 'str'`, which names neither the argument nor which of the
                # names it is.
                raise TypeError(
                    f"`individuals` is a sequence of the names of the "
                    f"individuals to keep, and {name!r}, a "
                    f"{type(name).__name__}, is not one of them"
                )
        self._steps.filter_individuals(names)

    def iter_blocks(
        self,
        fields: Iterable[Field] = ("chrom", "pos"),
        num_vars_per_block: int | None = None,
    ) -> Blocks:
        """The variants of the source, block by block, from its start.

        `fields` names what each block carries besides the genotypes, among
        ``"chrom"``, ``"pos"``, ``"id"``, ``"alleles"`` and ``"qual"``, and
        any other name is a ``ValueError``. The chromosome and the position
        are one field of the reader, so asking for one fills both. A field
        that is not asked for is never parsed, and the block has ``None``
        where it would be.

        `num_vars_per_block` is how many variants a block holds, and
        ``None`` asks for the size the core works out from the number of
        individuals, so that a block is a few million genotypes. The size
        changes nothing but where the cuts fall: the blocks of a source,
        joined, are the same for any size, and only the last one can be
        shorter than the rest.

        What it gives is a :class:`Blocks`, the blocks one after another and
        the counts of the pass in its ``pass_stats``: how many variants have
        come out of it, and what each filter of the ``Variants`` was given
        and kept.

        Every call reads the source from its start. When a variant cannot be
        read, the error comes in the place of the block that would have held
        it, and the variants of that block that were read are lost with it.
        A variant popnei cannot read is a ``ValueError`` whose message names
        the file, the line of a VCF and what is wrong, and so is a batch of
        a vars file whose buffers are compressed with zstd, which popnei
        writes in no build; a file that was cut short, a file that bgzip
        wrote and whose bytes were damaged, and a batch of a vars file that
        cannot be decoded are an ``OSError`` with the path in ``filename``
        and no ``errno``, because nothing of the file system refused
        anything.

        What a user has received when that error comes depends on
        `num_vars_per_block`: the variants that were read and had not filled
        a block are lost with it, so fewer variants come out than the file
        holds before the cut, and with the default size, which is thousands
        of variants, a file that was cut may give none at all. A user who
        has to know how far a damaged file was read asks for small blocks.
        """
        if isinstance(fields, str):
            # A string is a sequence of its letters, and asking for one
            # field would be read as asking for `a`, `l`, `l`, `e`...
            raise TypeError(
                f"`fields` is a sequence of names and not one name: write "
                f'fields=("{fields}",) for that one field'
            )
        return Blocks(
            self._source.blocks(list(fields), num_vars_per_block, self._steps)
        )


def _arguments_of(step: Step) -> str:
    """The arguments of one step, as the user wrote them in the call."""
    return ", ".join(f"{name}={value!r}" for name, value in step.args.items())
