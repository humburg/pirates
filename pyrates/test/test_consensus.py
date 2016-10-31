"""
Test consensus module.
"""
import os

from nose2.tools import params
from nose2.tools.decorators import with_setup, with_teardown

import pyrates.consensus as cons
import pyrates.sequence as sequence
from pyrates.test import TMP
from pyrates.test.fixtures import (setup_fastq_mismatch, setup_fastq_simple,
                                   teardown_fastq_mismatch, teardown_fastq_simple,
                                   create_consensus)


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
    expect_str2 = "@2 3G1T1 11A1T1\nAAAAACTTTTTGTCTAAGC\n+\nIIIIIIIIIIIIIIIIIII"

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

@params(([0, 1, 2, 3]), ([1, 0, 2, 3]), ([1, 2, 0, 3]), ([2, 1, 0, 3]), ([2, 0, 1, 3]),
        ([0, 1, 3, 2]), ([1, 0, 3, 2]), ([1, 2, 3, 0]), ([2, 1, 3, 0]), ([2, 0, 3, 1]),
        ([0, 3, 1, 2]), ([1, 3, 0, 2]), ([1, 3, 2, 0]), ([2, 3, 1, 0]), ([2, 3, 0, 1]),
        ([3, 0, 1, 2]), ([3, 1, 0, 2]), ([3, 1, 2, 0]), ([3, 2, 1, 0]), ([3, 2, 0, 1]))
def test_merge_diff(idx):
    """Propagate diffs when merging clusters"""
    uid1 = "ACCT"
    uid2 = "ACTT"

    seq1 = ["ACTGTTTGTCTAAGC"]*2
    qual1 = ['I'*len(seq1[0])]*len(seq1)
    seq2 = ["ACTGTTTTTCTAAGC"]*5
    qual2 = ['I'*len(seq2[0])]*len(seq2)
    seq3 = ["ACTGTTTTTCTAAGC"]*2
    qual3 = ['I'*len(seq3[0])]*len(seq3)
    seq4 = ["ACTGTTTGTGTAAGC", "ACTGTTTGTGTAAGC", "ACTGTTTGTATAAGC"]
    qual4 = ['I'*len(seq4[0])]*len(seq4)

    consensus = create_consensus([uid1 + uid2]*len(seq1) + \
                                   [uid2 + uid1]*len(seq2) + \
                                   [uid2 + uid2]*len(seq3) + \
                                   [uid1 + uid1]*len(seq4),
                                 ['I'*(len(uid1)*2)]*(len(seq1) + len(seq2) + \
                                   len(seq3) + len(seq4)),
                                 seq1 + seq2 + seq3 + seq4,
                                 qual1 + qual2 + qual3 + qual4)
    ids = [uid1+ uid2, uid2 + uid1, uid2 + uid2, uid1 + uid1]
    clusters = [consensus[ids[i]] for i in idx]
    merged = clusters[0]
    for i in range(1, len(clusters)):
        success = merged.merge(clusters[i], 2)
        assert success
    expect = "@12 7G5T7 9A1C9G2"
    obs = str(merged).splitlines()
    assert merged.size == 12, "%r != %r" % (merged.size, 12)
    assert obs[0] == expect, "%r != %r" % (obs[0], expect)

@with_setup(setup_fastq_simple)
@with_teardown(teardown_fastq_simple)
def test_fastq_simple():
    """Create consensus from fastq file."""
    cluster = cons.from_fastq(TMP + 'simple.fastq', 4, 'ACGT')
    uid1_expect = 'AAAACCCC'
    uid2_expect = 'CCCCAAAA'
    seq1_expect = 'ACCTCTCCCTGTGGGTCATGTGACT'
    seq2_expect = 'TTGTTTGAAAAACCTCGAAAGTAAC'

    assert uid1_expect in cluster, "%r not in %r" % (uid1_expect, list(cluster.keys()))
    assert uid2_expect in cluster, "%r not in %r" % (uid2_expect, list(cluster.keys()))
    assert cluster[uid1_expect].sequence.sequence == seq1_expect, \
           "%r != %r" % (cluster[uid1_expect].sequence.sequence, seq1_expect)
    assert cluster[uid2_expect].sequence.sequence == seq2_expect, \
           "%r != %r" % (cluster[uid2_expect].sequence.sequence, seq2_expect)

@with_setup(setup_fastq_mismatch)
@with_teardown(teardown_fastq_mismatch)
def test_fastq_mismatch():
    """Create consensus from reads with mismatches in cluster."""
    cluster = cons.from_fastq(TMP + 'mismatch.fastq', 4, 'ACGT')
    uid1_expect = 'AAAACCCC'
    uid2_expect = 'CCCCAAAA'
    uid3_expect = 'AAAAAAAA'
    seq1_expect = 'ACCTCTCCCTGTGGGTCATGTGACT'
    seq2_expect = 'TTGTTTGAAAAACCTCGAAAGTAAC'
    seq3_expect = 'CATTTTTGTGTCCAATGCCTAAATT'

    assert uid1_expect in cluster, "%r not in %r" % (uid1_expect, list(cluster.keys()))
    assert uid2_expect in cluster, "%r not in %r" % (uid2_expect, list(cluster.keys()))
    assert uid3_expect in cluster, "%r not in %r" % (uid3_expect, list(cluster.keys()))
    assert cluster[uid1_expect].sequence.sequence == seq1_expect, \
           "%r != %r" % (cluster[uid1_expect].sequence.sequence, seq1_expect)
    assert cluster[uid2_expect].sequence.sequence == seq2_expect, \
           "%r != %r" % (cluster[uid2_expect].sequence.sequence, seq2_expect)
    assert cluster[uid3_expect].sequence.sequence == seq3_expect, \
           "%r != %r" % (cluster[uid3_expect].sequence.sequence, seq3_expect)

@params((0, 0, 0), (1, 3, 5))
def test_output_simple(tol, size, target):
    """Write output without merging of clusters"""
    uid1 = "ACCT"
    uid2 = "ACTT"
    seq1 = ["ACTGTTTGTCTAAGC"]*3
    qual1 = ['I'*len(seq1[0])]*len(seq1)
    seq2 = ["ACTGTTTTTCTAAGC"]*5
    qual2 = ['I'*len(seq2[0])]*len(seq2)
    consensus = create_consensus([uid1 + uid2]*len(seq1) + [uid2 + uid1]*len(seq2),
                                 ['I'*(len(uid1) + len(uid2))]*(len(seq1) + len(seq2)),
                                 seq1 + seq2, qual1 + qual2)
    cons.to_fastq(consensus, TMP + 'simple_out.fastq', tol, size, target)
    with open(TMP + 'simple_out.fastq') as fastq:
        lines = fastq.readlines()
        lines = [line.rstrip() for line in lines]

    assert len(lines) == 8
    seqs = [lines[0:4], lines[4:8]]
    seqs.sort(key=lambda x: x[0])
    assert seqs[0][0] == '@3', "%r != %r" % (seqs[0][0], '@3')
    assert seqs[0][1] == uid1 + uid2 + seq1[0], "%r != %r" % (seqs[0][1], uid1 + uid2 + seq1[0])
    assert seqs[0][2] == '+'
    assert seqs[0][3] == 'I'*(len(uid1) + len(uid2)) + qual1[0], \
           "%r != %r" % (seqs[0][3] == 'I'*(len(uid1) + len(uid2)) + qual1[0])
    assert seqs[1][0] == '@5', "%r != %r" % (seqs[1][0], '@5')
    assert seqs[1][1] == uid2 + uid1 + seq2[0], "%r != %r" % (seqs[1][1], uid2 + uid1 + seq2[0])
    assert seqs[1][2] == '+'
    assert seqs[1][3] == 'I'*(len(uid1) + len(uid2)) + qual2[0], \
           "%r != %r" % (seqs[1][3] == 'I'*(len(uid2) + len(uid1)) + qual2[0])
    os.remove(TMP + 'simple_out.fastq')

@with_teardown(lambda: os.remove(TMP + 'merge_out.fastq'))
def test_output_merge():
    """Write output, allow merging of clusters."""
    uid1 = "ACCT"
    uid2 = "ACTT"
    uid3 = "GGGG"

    seq1 = ["ACTGTTTGTCTAAGC"]*2
    qual1 = ['I'*len(seq1[0])]*len(seq1)
    seq2 = ["ACTGTTTTTCTAAGC"]*5
    qual2 = ['I'*len(seq2[0])]*len(seq2)
    seq5 = ["ACTGTTTTTCTAAGC"]*2
    qual5 = ['I'*len(seq5[0])]*len(seq5)

    seq3 = ["GGACGGGGCAATTTA"]
    qual3 = ['I'*len(seq3[0])]

    seq4 = ["ACTGTTTTTCTAAGC"]*10
    qual4 = ['I'*len(seq4[0])]*len(seq4)

    consensus = create_consensus([uid3 + uid3]*len(seq4) + [uid1 + uid2]*len(seq1) + \
                                   [uid1 + uid1] + [uid2 + uid1]*len(seq2) + \
                                   [uid2 + uid2]*len(seq5),
                                 ['I'*(len(uid1) + len(uid2))]*(len(seq1) + len(seq2) + \
                                   len(seq3) + len(seq4) + len(seq5)),
                                 seq4 + seq1 + seq3 + seq2 + seq5,
                                 qual4 + qual1 + qual3 + qual2 + qual5)
    cons.to_fastq(consensus, TMP + 'merge_out.fastq', 2, 4, 7)
    with open(TMP + 'merge_out.fastq') as fastq:
        lines = fastq.readlines()
        lines = [line.rstrip() for line in lines]

    expect1 = ['@9 7G2T7', uid2 + uid1 + seq2[0], '+', 'I'*(len(uid1) + len(uid2)) + qual1[0]]
    expect2 = ['@1', uid1 + uid1 + seq3[0], '+', 'I'*(2*len(uid1)) + qual3[0]]
    expect3 = ['@10', uid3 + uid3 + seq4[0], '+', 'I'*(len(uid3) + len(uid3)) + qual4[0]]
    assert len(lines) == 12, "%r != %r" % (len(lines), 12)
    seqs = [lines[0:4], lines[4:8], lines[8:12]]
    seqs.sort(key=lambda x: int(x[0][1:].split()[0]))
    for obs, exp in zip(seqs[1], expect1):
        assert obs == exp, "%r != %r" % (obs, exp)
    for obs, exp in zip(seqs[0], expect2):
        assert obs == exp, "%r != %r" % (obs, exp)
    for obs, exp in zip(seqs[2], expect3):
        assert obs == exp, "%r != %r" % (obs, exp)
