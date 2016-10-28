"""Test sequence module."""

from nose2.tools import params
from nose2.tools.such import helper
from pyrates.sequence import SequenceWithQuality

def test_swq_new():
    """Create sequence objects"""
    read = "ACTGGGTGTAT"
    qual = "IIIIIIIIIII"
    name = 'test'
    seq = SequenceWithQuality(read, qual)
    seq = SequenceWithQuality(read, qual, name)
    assert len(seq) == len(read), "Lengths of sequence and read don't match"
    assert seq.sequence == read, "Sequence doesn't match input read"
    assert seq.quality == qual, "Quality doesn't match input"
    assert seq.name == name, "Name doesn't match input"

def test_swq_access():
    """Access to sequence object members"""
    read = "ACTGGGTGTAT"
    read2 = "TATGTGGGTCA"
    qual = "IIIIIIIIIII"
    qual2 = "DDDDDDDDDDD"
    name = 'test'
    seq = SequenceWithQuality(read, qual, name)
    assert seq[0] == ('A', 'I'), "Failed to extract sequence/quality pair"
    seq.sequence = read2
    assert seq.sequence == read2, "Failed to change read sequence"
    seq.quality = qual2
    assert seq.quality == qual2, "Failed to change qualities"

    ## check that attempts to set members to new values that
    ## would result in an inconsistent state raise an exception
    with helper.assertRaises(ValueError):
        SequenceWithQuality("ACTG", qual)
    with helper.assertRaises(ValueError):
        seq.sequence = "ACTG"
    with helper.assertRaises(ValueError):
        SequenceWithQuality(read, "!!!!")
    with helper.assertRaises(ValueError):
        seq.quality = "!!!!!"

def test_swq_str():
    """String representation of sequence objects"""
    read = "ACTGGGTGTAT"
    qual = "IIIIIIIIIII"
    name = 'test'
    seq = SequenceWithQuality(read, qual, name)
    expect_str = "@test\nACTGGGTGTAT\n+\nIIIIIIIIIII"
    expect_repr = "SequenceWithQuality(sequence='ACTGGGTGTAT', quality='IIIIIIIIIII', name='test')"
    assert str(seq) == expect_str, "%r != %r" % (str(seq), expect_str)
    assert repr(seq) == expect_repr, "%r != %r" % (repr(seq), expect_repr)

def test_grosslydifferent_equal():
    """Recognise identical sequences as similar"""
    seq1 = SequenceWithQuality('AACTGTGAGTGTAGATGTTC', 'I'*20)
    assert not seq1.grosslydifferent(seq1), \
        "Sequence %r considered too different from itself" % seq1.sequence

@params(('AACTGTGAGTGTAGATGTTC', 'AACTTTGAGTGTAGATGTTC'), \
        ('AACTGTGAGTGTAGATGTTC', 'TTTTTTTTTTGTAGATGTTC'))
def test_grosslydifferent_similar(seq1, seq2):
    """Tolerate small differences between sequences."""
    seq1 = SequenceWithQuality(seq1, 'I'*len(seq1))
    seq2 = SequenceWithQuality(seq2, 'I'*len(seq2))
    assert not seq1.grosslydifferent(seq2), \
        "Sequences %r and %r considered too different" % (seq1.sequence, seq2.sequence)

def test_grosslydifferent_diff():
    """Tolerate small differences between sequences."""
    seq1 = SequenceWithQuality('AACTGTGAGTGTAGATGTTC', 'I'*20)
    seq2 = SequenceWithQuality('GTAGATGTTCGTAGATGTTC', 'I'*20)
    assert seq1.grosslydifferent(seq2), \
        "Sequences %s and %s considered similar" % (seq1.sequence, seq2.sequence)
