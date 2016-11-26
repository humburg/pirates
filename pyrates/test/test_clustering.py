"""Tests for sequence clustering"""

import os
from nose2.tools import params
from nose2.tools.decorators import with_setup, with_teardown

import pyrates.clustering as clust
import pyrates.sequence as pseq
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

def test_merge_targets():
    """Identify cluster for merging"""
    uid1 = "ACCT"
    uid2 = "GGGG"
    uid3 = "AAGG"

    seq1 = ["ACTGTTTGTCTAAGC"]*2
    qual1 = ['I'*len(seq1[0])]*len(seq1)
    seq2 = ["ACTGTTTTTCTAAGC"]*5
    qual2 = ['I'*len(seq2[0])]*len(seq2)
    seq3 = ["ACTGTTTTTCTAAGC"]*2
    qual3 = ['I'*len(seq3[0])]*len(seq3)

    clusters = create_consensus([uid1 + uid1]*len(seq1) + \
                                   [uid2 + uid2]*len(seq2),
                                ['I'*(len(uid1)*2)]*(len(seq1) + len(seq2)),
                                seq1 + seq2, qual1 + qual2)
    seq3 = [pseq.SequenceWithQuality(seq, qual) for seq, qual in zip(seq3, qual3)]
    uid = pseq.SequenceWithQuality(uid2 + uid3, 'I'*(len(uid2) + len(uid3)))
    cand = clusters.merge_target(uid, seq3[0], {}, 2, None)
    assert cand == uid2 + uid2, "%r != %r" % (cand, uid2 + uid2)
    cand = clusters.merge_target(uid, seq3[0], {}, 1, None)
    assert cand is None, "%r != %r" % (cand, None)


@with_setup(setup_fastq_simple)
@with_teardown(teardown_fastq_simple)
def test_fastq_simple():
    """Create consensus from fastq file."""
    cluster = clust.Clustering.from_fastq(TMP + 'simple.fastq', 4, 'ACGT', threshold=0)
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
    cluster = clust.Clustering.from_fastq(TMP + 'mismatch.fastq',
                                          id_length=4,
                                          adapter='ACGT',
                                          threshold=0)
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
