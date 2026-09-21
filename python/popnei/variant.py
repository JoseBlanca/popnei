"""The handle a user holds: a source of variants and its individuals."""

from collections.abc import Iterable, Iterator

from popnei import _core
from popnei.block import Block, Field, _block_of


class Variants:
    """A source of variants: a VCF with the options it is read with.

    It holds no genotypes. A user gets one from :func:`popnei.open_vcf` and
    gives it to as many calculations as they want: each one opens the source
    again and runs its loop over the variants inside the Rust core, so the
    dataset is never in memory as a whole.

    It is pyNei's ``Variants`` under the word of ``docs/glossary.md``: what
    pyNei calls a sample is here an individual, one organism that was
    genotyped. The genotypes come out of it through :meth:`iter_blocks` and
    through nothing else.

    It cannot be pickled or copied: what it holds is an object of Rust with
    the path and the options of the source. What travels between processes
    is the path and the arguments of :func:`popnei.open_vcf`, and a
    ``Variants`` is opened again at the other end.
    """

    def __init__(self, source: _core.VcfSource):
        """The handle over `source`, which :func:`popnei.open_vcf` builds."""
        self._source = source
        # The names come from the header, which was read once, so they are
        # taken out of the core here and not at every use.
        self._individuals = tuple(source.individuals())

    @property
    def individuals(self) -> tuple[str, ...]:
        """The names of the individuals, in the order the source has them."""
        return self._individuals

    @property
    def num_individuals(self) -> int:
        """How many individuals the source holds."""
        return len(self._individuals)

    @property
    def ploidy(self) -> int:
        """How many alleles the genotype of one individual holds."""
        return self._source.ploidy()

    def iter_blocks(
        self,
        fields: Iterable[Field] = ("chrom", "pos"),
        num_vars_per_block: int | None = None,
    ) -> Iterator[Block]:
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

        Every call reads the source from its start. When a variant cannot be
        read, the error comes in the place of the block that would have held
        it, and the variants of that block that were read are lost with it.
        A variant popnei cannot read is a ``ValueError`` whose message names
        the file, the line and what is wrong; a file that was cut short, or
        a file that bgzip wrote and whose bytes were damaged, is an
        ``OSError`` with the path in ``filename`` and no ``errno``, because
        nothing of the file system refused anything.

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
        return map(_block_of, self._source.blocks(list(fields), num_vars_per_block))
