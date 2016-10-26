"""
Test consensus module.
"""
from nose2.tools import params
import pyrates.consensus as cons

def test_grosslydifferent_equal():
    """Recognise identical sequences as similar"""
    seq1 = 'AACTGTGAGTGTAGATGTTC'
    assert not cons.grosslydifferent(seq1, seq1), \
        "Sequence %s considered too different from itself" % seq1

@params(('AACTGTGAGTGTAGATGTTC', 'AACTTTGAGTGTAGATGTTC'), \
        ('AACTGTGAGTGTAGATGTTC', 'TTTTTTTTTTGTAGATGTTC'))
def test_grosslydifferent_similar(seq1, seq2):
    """Tolerate small differences between sequences."""
    assert not cons.grosslydifferent(seq1, seq2), \
        "Sequences %s and %s considered too different" % (seq1, seq2)

def test_grosslydifferent_diff():
    """Tolerate small differences between sequences."""
    seq1 = 'AACTGTGAGTGTAGATGTTC'
    seq2 = 'GTAGATGTTCGTAGATGTTC'
    assert cons.grosslydifferent(seq1, seq2), \
        "Sequences %s and %s considered similar" % (seq1, seq2)

def test_consensus_idlen():
    """Skip sequences with incompatible IDs"""
    qual_id1 = "IIII"
    qual_id2 = "IIIII"
    (seq1, qual_seq1) = ("AACTGTGAGTGTAGATGTTC", "I"*20)
    consensus = cons.consensus(qual_id1, qual_id2, seq1, seq1, qual_seq1, qual_seq1, 1, {})
    expect = {'qid':qual_id1, 'seq':seq1, 'qseq':qual_seq1}
    assert consensus == expect, "%r != %r" % (consensus, expect)

def test_consensus_seqlen():
    """Skip sequences with different lengths"""
    qual_id1 = "IIII"
    (seq1, qual_seq1) = ("AACTGTGAGTGTAGATGTTC", "I"*20)
    (seq2, qual_seq2) = ("AACTGTGAGTGTAGATGTTCTGTA", "I"*24)
    consensus = cons.consensus(qual_id1, qual_id1, seq1, seq2, qual_seq1, qual_seq2, 1, {})
    expect = {'qid':qual_id1, 'seq':seq1, 'qseq':qual_seq1}
    assert consensus == expect, "%r != %r" % (consensus, expect)

@params(('qqqqq', 'IIIII', 'qqqqq'), \
        ('IIIII', 'qqqqq', 'qqqqq'), \
        ('abcde', 'edcba', 'edcde'))
def test_update_qual(qual1, qual2, expect):
    """Retain highest quality"""
    updated = cons.update_qual(qual1, qual2)
    assert updated == expect, \
           "Failed to retain high quality sequence (%s != %s)" % (updated, expect)

def test_consensus_seq():
    """Compute consensus sequence"""
    seq1 = "ACTGTTTGTCTAAGC"
    seq2 = "ACTTTTTGTCTTAGC"
    qual1 = "IIIDIIIIIIIIIII"
    qual2 = "IIIIIIIIIDIDIII"
    diffs = {}
    consensus = cons.consensus('IIII', 'IIII', seq1, seq2, qual1, qual2, 1, diffs)
    seq_expect = "ACTTTTTGTCTAAGC"
    qual_expect = "I"*len(seq_expect)
    diff_expect = {3:{'A':0, 'C':0, 'T':1, 'G':1, 'N':0}, 11:{'A':1, 'C':0, 'T':1, 'G':0, 'N':0}}
    assert consensus['seq'] == seq_expect, \
           "Failed to update consensus (%s != %s)" % (consensus['seq'], seq_expect)
    assert consensus['qseq'] == qual_expect, \
           "Failed to update qualities (%s != %s)" % (consensus['qseq'], qual_expect)
    assert diffs == diff_expect, \
           "Incorrect sequence diff (%r != %r)" % (diffs, diff_expect)

def test_consensus_diff():
    """Update sequence diff"""
    seq1 = "ACTTTTTGTCTAAGC"
    seq2 = "ACTTTTTGTGTTAGC"
    qual1 = "IIIIIIIIIIIIIII"
    qual2 = "IIIIIIIIIqIDIII"
    diffs = {3:{'A':0, 'C':0, 'T':1, 'G':1, 'N':0}, 11:{'A':1, 'C':0, 'T':1, 'G':0, 'N':0}}
    consensus = cons.consensus('IIII', 'IIII', seq1, seq2, qual1, qual2, 2, diffs)
    seq_expect = "ACTTTTTGTGTAAGC"
    qual_expect = "IIIIIIIIIqIIIII"
    diff_expect = {3:{'A':0, 'C':0, 'T':2, 'G':1, 'N':0},
                   11:{'A':1, 'C':0, 'T':2, 'G':0, 'N':0},
                   9:{'A':0, 'C':2, 'T':0, 'G':1, 'N':0}}
    assert consensus['seq'] == seq_expect, \
           "Failed to update consensus (%s != %s)" % (consensus['seq'], seq_expect)
    assert consensus['qseq'] == qual_expect, \
           "Failed to update qualities (%s != %s)" % (consensus['qseq'], qual_expect)
    assert diffs == diff_expect, \
           "Incorrect sequence diff (%r != %r)" % (diffs, diff_expect)
