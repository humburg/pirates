"""Fixtures for unit tests"""

import os
import os.path
import pyrates.sequence as sequence
import pyrates.consensus as cons
from pyrates.test import TMP

def create_fastq(seq, qual, file_name, name=None):
    """Create a fastq file.

    Args:
        seq (:obj:`list`): Sequences for reads to be written to file.
        qual (:obj:`list`): ASCII encoded qualities for reads to be written to file.
        file_name (:obj:`str`): Name of output file.
        name (:obj:`list`, optional): Names to use for reads in output. If this is `None`
            reads will be named *test_1*, *test_2*, ...
    """
    ## Don't overwrite existing file
    file_name = TMP + file_name
    if os.path.isfile(file_name):
        return
    if not os.path.isdir(TMP):
        os.makedirs(TMP)

    if name is None:
        name = ['test_%d' % i for i in range(len(seq))]
    with open(file_name, 'w') as fastq:
        for data in zip(name, seq, qual):
            fastq.write("@%s\n%s\n+\n%s\n" % data)

def setup_fastq_simple():
    """Create fastq file with reads from two clusters without UID errors."""
    uid1 = 'AAAA'
    uid2 = 'CCCC'
    adapter = 'ACGT'
    adapter_rev = 'ACGT'
    read1 = ['ACCTCTCCCTGTGGGTCATGTGACT']*3
    read1 = [uid1 + adapter + r + adapter_rev + uid2 for r in read1]
    read2 = ['TTGTTTGAAAAACCTCGAAAGTAAC']*5
    read2 = [uid2 + adapter + r + adapter_rev + uid1 for r in read2]
    qual = ['I'*len(read1[0])]*(len(read1) + len(read2))
    create_fastq(read1 + read2, qual, 'simple.fastq')

def teardown_fastq_simple():
    """Remove files created for simple fastq test"""
    os.remove(TMP + 'simple.fastq')

def setup_fastq_mismatch():
    """Create fastq file with reads from three clusters with errors in sequence."""
    uid1 = 'AAAA'
    uid2 = 'CCCC'
    adapter = 'ACGT'
    adapter_rev = 'ACGT'
    read1 = ['ACCTCTCCCTGTGGGTCATGTGACT', 'ACCTCTCCCTGTGTGTCATGTGACT', 'ACCTCTCCCTGTGGGTCATGTGACT']
    qual1 = ['IIIIIIIIIIIIIIIIIIIIIIIII', 'IIIIIIIIIIIIIDIIIIIIIIIII', 'IIIIIIIIIIIIIIIIIIIIIIIII']
    qual1 = ['I'*len(uid1) + 'I'*len(adapter) + q + 'I'*len(adapter) + 'I'*len(uid2) for q in qual1]
    read1 = [uid1 + adapter + r + adapter_rev + uid2 for r in read1]
    read2 = ['TTGTTTGAAAAACCTCGAAAGTAAC', 'TTGTTTGAATAACCTCGAAAGTAAC', 'TTGTTTGAAAAACCTCGAAAGTAAC',
             'ACCTCTCCCTGTGGGTCATGTGACT']
    read2 = [uid2 + adapter + r + adapter_rev + uid1 for r in read2]
    qual2 = ['I'*len(read2[0])]*len(read2)
    read3 = ['CCTTTTTGTGTCCAATGCCTAAATT', 'CATTTTTGTGTCCAATGCCTAAATT', 'CCTTTTTGTGTCCAATGCCTAAATT']
    qual3 = ['I!IIIIIIIIIIIIIIIIIIIIIII', 'IIIIIIIIIIIIIIIIIIIIIIIII', 'I!IIIIIIIIIIIIIIIIIIIIIII']
    read3 = [uid1 + adapter + r + adapter_rev + uid1 for r in read3]
    qual3 = ['I'*len(uid1) + 'I'*len(adapter) + q + 'I'*len(uid1) + 'I'*len(adapter) for q in qual3]
    create_fastq(read1 + read2 + read3, qual1 + qual2 + qual3, 'mismatch.fastq')

def teardown_fastq_mismatch():
    """Remove files created for simple fastq test"""
    os.remove(TMP + 'mismatch.fastq')

def create_consensus(uids, uid_qual, seqs, seq_qual):
    """Create consensus dictionary from raw sequences.import

    Args:
        uids (:obj:`list`): UID sequences.
        seqs (:obj:`list`): Read sequences.

    Returns:
        :obj:`dict`: Consensus sequences.
    """
    uid_with_qual = [sequence.SequenceWithQuality(seq, qual) for seq, qual in zip(uids, uid_qual)]
    seq_with_qual = [sequence.SequenceWithQuality(seq, qual) for seq, qual in zip(seqs, seq_qual)]
    cluster = {}
    for uid, seq in zip(uid_with_qual, seq_with_qual):
        if uid.sequence not in cluster:
            cluster[uid.sequence] = cons.Consensus(uid, seq)
        else:
            cluster[uid.sequence].update(uid, seq)
    return cluster
