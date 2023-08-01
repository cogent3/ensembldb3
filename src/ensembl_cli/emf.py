# parser for Ensembl Multi Format (EMF) FLATFILE DUMPS
# we limit this to the Compara flavoured version
import os
import typing

from cogent3 import open_

from src.ensembl_cli.name import EmfName


# TODO spaces are optional between columns representing SEQ and SCORE lines
# gah discuss with Ensembl
def _get_block_seqnames(data) -> dict[str, str]:
    names = []
    for i, line in enumerate(data):
        if line.startswith("SEQ"):
            names.append(EmfName(*line.strip().split()[1:]))
        elif line.startswith("DATA"):
            break

    seq_data = [aln_col.split()[0].strip() for aln_col in data[i + 1 :]]
    # we ignore the ancestral sequences
    return {
        n: "".join(s)
        for n, *s in zip(names, *seq_data)
        if n.species != "ancestral_sequences"
    }


def _iter_blocks(data: typing.Iterable[str]) -> list[tuple[int, int]]:
    # find block boundaries
    start = 0
    blocks = []
    for i, line in enumerate(data):
        if line.startswith("//"):
            blocks.append((start, i))
            start = i + 1

    return blocks


# we need a raw parser
def parse_emf(path: typing.Union[str, os.PathLike]) -> dict[EmfName, str]:
    """yield data for alignment from EMF files

    Parameters
    ----------
    path
        location of emf file

    Returns
    -------
    {EmfName(): <seq string>, ...}

    Notes
    -----
    The key (EmfName) has useful attributes, including the python
    coordinates for the sequence, coord name, species, etc...

    Raises
    ------
    NotImplementedError if not compara emf format
    """
    with open_(path) as infile:
        data = infile.readlines()
        if not data[0].startswith("##FORMAT (compara)"):
            raise NotImplementedError(
                f"only compara format supported, not {data[0].strip()!r}"
            )

    blocks = _iter_blocks(data)
    for start, end in blocks:
        yield _get_block_seqnames(data[start:end])
