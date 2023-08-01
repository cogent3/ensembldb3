import pytest

from cogent3 import make_seq


@pytest.mark.parametrize(
    "seq", ("----", "---AGC--TGC--", "AGC--TGC--", "---AGC--TGC", "AGCTGC")
)
def test_roundtrip_gapped_seqs(seq):
    seq = make_seq(seq, moltype="dna")
    ug = seq.degap()
    c = seq_to_gap_coords(seq)
    aligned = gap_coords_to_seq(c, ug)
    assert str(aligned) == str(seq)
