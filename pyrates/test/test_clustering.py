"""Tests for sequence clustering"""

import os
from nose2.tools import params
from nose2.tools.decorators import with_setup, with_teardown

import pyrates.clustering as clust
from pyrates.test import TMP
from pyrates.test.fixtures import (setup_fastq_mismatch, setup_fastq_simple,
                                   teardown_fastq_mismatch, teardown_fastq_simple,
                                   create_consensus)


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

    clusters = create_consensus([uid1 + uid2]*len(seq1) + \
                                   [uid2 + uid1]*len(seq2) + \
                                   [uid2 + uid2]*len(seq3) + \
                                   [uid1 + uid1]*len(seq4),
                                ['I'*(len(uid1)*2)]*(len(seq1) + len(seq2) + \
                                   len(seq3) + len(seq4)),
                                seq1 + seq2 + seq3 + seq4,
                                qual1 + qual2 + qual3 + qual4)
    ids = [uid1+ uid2, uid2 + uid1, uid2 + uid2, uid1 + uid1]
    centres = [clusters[ids[i]] for i in idx]
    merged = centres[0]
    for i in range(1, len(clusters)):
        success = merged.merge(centres[i], 2)
        assert success
    expect = "@12 7G5T7 9A1C9G2"
    obs = str(merged).splitlines()
    assert merged.size == 12, "%r != %r" % (merged.size, 12)
    assert obs[0] == expect, "%r != %r" % (obs[0], expect)

@with_setup(setup_fastq_simple)
@with_teardown(teardown_fastq_simple)
def test_fastq_simple():
    """Create consensus from fastq file."""
    cluster = clust.Clustering.from_fastq(TMP + 'simple.fastq', 4, 'ACGT')
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
    cluster = clust.Clustering.from_fastq(TMP + 'mismatch.fastq', 4, 'ACGT')
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

@params((0, 0), (1, 3))
def test_output_simple(tol, size):
    """Write output without merging of clusters"""
    uid1 = "ACCT"
    uid2 = "ACTT"
    seq1 = ["ACTGTTTGTCTAAGC"]*3
    qual1 = ['I'*len(seq1[0])]*len(seq1)
    seq2 = ["ACTGTTTTTCTAAGC"]*5
    qual2 = ['I'*len(seq2[0])]*len(seq2)
    clustering = create_consensus([uid1 + uid2]*len(seq1) + [uid2 + uid1]*len(seq2),
                                  ['I'*(len(uid1) + len(uid2))]*(len(seq1) + len(seq2)),
                                  seq1 + seq2, qual1 + qual2)
    clustering.merge(tol, size)
    clustering.write(TMP + 'simple_out.fastq')
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

    clustering = create_consensus([uid3 + uid3]*len(seq4) + [uid1 + uid2]*len(seq1) + \
                                   [uid1 + uid1] + [uid2 + uid1]*len(seq2) + \
                                   [uid2 + uid2]*len(seq5),
                                  ['I'*(len(uid1) + len(uid2))]*(len(seq1) + len(seq2) + \
                                   len(seq3) + len(seq4) + len(seq5)),
                                  seq4 + seq1 + seq3 + seq2 + seq5,
                                  qual4 + qual1 + qual3 + qual2 + qual5)
    clustering.merge(2, 4, 7)
    clustering.write(TMP + 'merge_out.fastq')
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
