"""
Test consensus module.
"""
from nose2.tools import params

import pyrates.consensus as cons
import pyrates.sequence as sequence

def test_consensus_new():
    """Create objects of class Consensus"""
    seq = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    id_seq = sequence.SequenceWithQuality("AAA", "III")
    consensus = cons.Consensus(id_seq, seq)
    assert consensus.sequence == seq
    assert consensus.uid == id_seq
    assert consensus.size == 1

def test_consensus_idlen():
    """Skip sequences with incompatible IDs"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    id2 = sequence.SequenceWithQuality("AAAAA", "IIIII")
    seq = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")

    consensus = cons.Consensus(id1, seq)
    success = consensus.update(id2, seq)
    assert not success
    assert consensus.uid == id1, "%r != %r" % (consensus.uid, id1)

def test_consensus_seqlen():
    """Skip shorter sequences"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    seq1 = sequence.SequenceWithQuality("AACTGTGAGTGTAGATGTTCTGTA", "I"*24)
    seq2 = sequence.SequenceWithQuality("AACTGTGAGTGTAGATGTTC", "I"*20)
    consensus = cons.Consensus(id1, seq1)
    success = consensus.update(id1, seq2)
    assert not success
    assert consensus.sequence == seq1, "%r != %r" % (consensus.sequence, seq1)
    assert consensus.shorter == 1, "Skipped sequence not recorded"

    consensus = cons.Consensus(id1, seq2)
    success = consensus.update(id1, seq1)
    assert not success
    assert consensus.sequence == seq1, "%r != %r" % (consensus.sequence, seq1)
    assert consensus.shorter == 1, "Skipped sequence not recorded"

    consensus = cons.Consensus(id1, seq2)
    success = consensus.update(id1, seq2)
    assert success
    success = consensus.update(id1, seq1)
    assert not success
    assert consensus.sequence == seq2, "%r != %r" % (consensus.sequence, seq2)
    assert consensus.longer == 1, "Skipped sequence not recorded"

def test_consensus_skip():
    """Reject sequences that are too different"""
    uid = sequence.SequenceWithQuality("AAA", "III")
    seq1 = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    seq2 = sequence.SequenceWithQuality("TTCTCCCTGGTAAGC", "IIIDIIIIIIIIIII")
    consensus = cons.Consensus(uid, seq1)
    success = consensus.update(uid, seq2)
    assert not success
    assert consensus.sequence == seq1, "%r != %r" % (consensus.sequence, seq1)
    assert consensus.different == 1, "Skipped sequence not counted"

@params(('qqqqq', 'IIIII', 'qqqqq'), \
        ('IIIII', 'qqqqq', 'qqqqq'), \
        ('abcde', 'edcba', 'edcde'))
def test_update_uid(qual1, qual2, expect):
    """Retain highest quality"""
    id1 = sequence.SequenceWithQuality("A"*len(qual1), qual1)
    id2 = sequence.SequenceWithQuality("A"*len(qual2), qual2)
    seq = sequence.SequenceWithQuality("A"*20, "I"*20)
    consensus = cons.Consensus(id1, seq)
    consensus._update_uid(id2)
    assert consensus.uid.quality == expect, \
           "Failed to retain high quality sequence (%r != %r)" % (consensus.uid.quality, expect)

def test_consensus_seq():
    """Compute consensus sequence"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    seq1 = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    seq2 = sequence.SequenceWithQuality("ACTTTTTGTCTTAGC", "IIIIIIIIIDIDIII")
    consensus = cons.Consensus(id1, seq2)
    success = consensus.update(id1, seq1)

    seq_expect = "ACTTTTTGTCTAAGC"
    qual_expect = "I"*len(seq_expect)
    diff_expect = {3:{'T':1, 'G':1}, 11:{'A':1, 'T':1}}
    assert success, "Sequence %r was rejected" % seq1
    assert consensus.sequence.sequence == seq_expect, \
           "Failed to update consensus (%s != %s)" % (consensus.sequence.sequence, seq_expect)
    assert consensus.sequence.quality == qual_expect, \
           "Failed to update qualities (%s != %s)" % (consensus.sequence.quality, qual_expect)
    assert consensus.diffs == diff_expect, \
           "Incorrect sequence diff (%r != %r)" % (consensus.diffs, diff_expect)

def test_consensus_diff():
    """Update sequence diff"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    seq1 = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    seq2 = sequence.SequenceWithQuality("ACTTTTTGTCTTAGC", "IIIIIIIIIDIDIII")
    seq3 = sequence.SequenceWithQuality("ACTTTTTGTGTTAGC", "IIIIIIIIIqIDIII")
    consensus = cons.Consensus(id1, seq2)
    success = consensus.update(id1, seq1)

    assert success, "Sequence %r was rejected" % seq1
    success = consensus.update(id1, seq3)

    seq_expect = "ACTTTTTGTGTAAGC"
    qual_expect = "IIIIIIIIIqIIIII"
    diff_expect = {3:{'T':2, 'G':1},
                   11:{'A':1, 'T':2},
                   9:{'C':2, 'G':1}}
    assert success, "Sequence %r was rejected" % seq3
    assert consensus.sequence.sequence == seq_expect, \
           "Failed to update consensus (%s != %s)" % (consensus.sequence.sequence, seq_expect)
    assert consensus.sequence.quality == qual_expect, \
           "Failed to update qualities (%s != %s)" % (consensus.sequence.quality, qual_expect)
    assert consensus.diffs == diff_expect, \
           "Incorrect sequence diff (%r != %r)" % (consensus.diffs, diff_expect)

def test_consensus_str():
    """String representation of consensus sequences"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    seq1 = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    seq2 = sequence.SequenceWithQuality("ACTTTTTGTCTTAGC", "IIIIIIIIIDIDIII")
    consensus = cons.Consensus(id1, seq1)
    expect_str1 = "@1\nAAAAACTGTTTGTCTAAGC\n+\nIIIIIIIDIIIIIIIIIII"
    expect_repr1 = "Consensus(uid=SequenceWithQuality(sequence='AAAA', " + \
                                                     "quality='IIII', name=''), " + \
                   "sequence=SequenceWithQuality(sequence='ACTGTTTGTCTAAGC', " + \
                                                "quality='IIIDIIIIIIIIIII', name=''), " + \
                   "diffs={}, size=1)"
    expect_str2 = "@2\nAAAAACTTTTTGTCTAAGC\n+4G1T1 12A1T1\nIIIIIIIIIIIIIIIIIII"

    assert str(consensus) == expect_str1, "\n%s\n!=\n%s" % (consensus, expect_str1)
    assert repr(consensus) == expect_repr1, "\n%r\n!=\n%r" % (consensus, expect_repr1)
    consensus.update(id1, seq2)
    assert str(consensus) == expect_str2, "\n%s\n!=\n%s" % (str(consensus), expect_str2)

def test_merge_simple():
    """Combine two consensus sequences"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    id2 = sequence.SequenceWithQuality("AACA", "IIII")
    seq = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    cons1 = cons.Consensus(id1, seq)
    cons2 = cons.Consensus(id2, seq)
    merged = cons1.merge(cons2, 1)
    assert merged, "Merging failed unexpectedly"
    assert cons1.size == 2, "Incorrect size for merged cluster (%d != %d)" % (cons1.size, 2)
    assert cons1.sequence.sequence == seq.sequence, "Incorrect merged sequence (%r != %r)" % \
                                           (cons1.sequence.sequence, seq.sequence)

def test_merge_fail_uid():
    """Don't merge sequences with very different UIDs'"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    id2 = sequence.SequenceWithQuality("CCAA", "IIII")
    seq = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    cons1 = cons.Consensus(id1, seq)
    cons2 = cons.Consensus(id2, seq)
    merged = cons1.merge(cons2, 1)
    assert not merged, "Merging succeeded unecpectedly"

def test_merge_size():
    """Update size of merged clusters"""
    id1 = sequence.SequenceWithQuality("AAAA", "IIII")
    id2 = sequence.SequenceWithQuality("AACA", "IIII")
    seq = sequence.SequenceWithQuality("ACTGTTTGTCTAAGC", "IIIDIIIIIIIIIII")
    cons1 = cons.Consensus(id1, seq)
    cons1.update(id1, seq)
    cons2 = cons.Consensus(id2, seq)
    cons2.update(id2, seq)
    merged = cons1.merge(cons2, 1)
    assert merged, "Merging failed unexpectedly"
    assert cons1.size == 4, "Incorrect size for merged cluster (%d != %d)" % (cons1.size, 4)
