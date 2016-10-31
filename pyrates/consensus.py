"""Functions to facilitate consensus calling.
"""

import logging
from collections import defaultdict

import pyrates.sequence as pseq
import pyrates.utils as utils

class Consensus(object):
    """Consensus sequence inferred from observed read sequences.

    Args:
        uid (:obj:`pyrates.sequence.SequenceWithQuality`): The
            unique molecular identifier associated with this
            sequence.
        sequence (:obj:`pyrates.sequence.SequenceWithQuality`): The
            sequence of the first read that serves as the basis for the
            subsequent consensus computations.

    Attributes:
        uid (:obj:`pyrates.sequence.SequenceWithQuality`): The
            unique molecular identifier associated with the
            consensus sequence.
        sequence (:obj:`pyrates.sequence.SequenceWithQuality`): The
            consensus sequence.
        diffs (:obj:`dict`): A collection of all sequence
            differences observed between the consensus and the
            underlying reads.
        size (:obj:`int`): Number of reads used to compute the consensus.
        different (:obj:`int`): Number of times consensus computation failed
            because the other sequence was too different.
        shorter (:obj:`int`): Number of times consensus computation failed
            because the other read was too short.
        longer (:obj:`int`): Number of times consensus computation failed
            because the other read was too long.
    """
    __slots__ = 'uid', 'sequence', 'diffs', 'size', 'different', 'shorter', 'longer'
    _logger = utils.get_logger(__name__)
    def __init__(self, uid, sequence):
        self.uid = uid
        self.sequence = sequence
        self.diffs = defaultdict(lambda: defaultdict(int))
        self.size = 1
        self.different = 0
        self.shorter = 0
        self.longer = 0

    def _update_uid(self, uid_other):
        """Update uid sequence and qualities.

        The provided :obj:`pyrates.sequence.SequenceWithQuality` is used
        to update the current UID of this consensus sequence.

        Note:
            The current implementation assumes that the sequence is identical
            and simply chooses the highest quality for each position.

        Args:
            uid_other (:obj:`pyrates.sequence.SequenceWithQuality`): UID
                sequence with associated qualities.
        """
        qual_update = list(self.uid.quality)
        for (i, qual_other) in enumerate(uid_other.quality):
            if qual_other > qual_update[i]:
                qual_update[i] = qual_other
        self.uid.quality = ''.join(qual_update)

    def update(self, uid_other, seq_other, size_other=1, discard=True):
        """Update consensus sequence.

        The read represented by `seq_other` is added to the consensus.

        Args:
            uid_other (:obj:`pyrates.sequence.SequenceWithQuality`): UID
                sequence with associated qualities.
            seq_other (:obj:`pyrates.sequence.SequenceWithQuality`): Read
                sequence with associated qualities for read that is to be
                added to the consensus.
            size_other (:obj:`int`, optional): Treat the sequence provided for
                updating as a representative of this many sequences.
            discard (:obj:`bool`, optional): If this is `True` sequences that
                are rejected for consensus computations are counted and are
                assumed to be excluded from further concideration. They are
                essentially assigned to this cluster but don't affect the
                consensus sequence.

        Returns:
            :obj:`bool`: `True` if the sequence was successfully added to the
                consensus, `False` otherwise.
        """
        # better do some sanity checking
        if len(self.uid) != len(uid_other):
            self._logger.error("Mismatch in id quality length, this should not happen. " +
                               "Check your input.")
            self._logger.debug("Mismatching quality strings were '%s' and '%s'",
                               self.uid.sequence, uid_other.sequence)
            return False

        # if grossly different then just count this and move on
        if self.sequence.grosslydifferent(seq_other):
            if discard:
                self.different += size_other
            return False

        # if sequence length is shorter, count this occurance, abandon this
        # sequence and move on
        if len(self.sequence) > len(seq_other):
            if discard:
                self.shorter += size_other
            return False

        # if new sequence is longer, count this occurance
        # replace consensus sequence if built from only one other sequence
        if len(self.sequence) < len(seq_other):
            if discard:
                if self.size == 1:
                    self.sequence = seq_other
                    self.shorter += size_other
                else: self.longer += size_other
            return False

        self._update_uid(uid_other)

        # Step through sequence, record highest quality at each step, want to save diffs
        # for changes to the sequence but not the quality.
        # If we encounter a mismatch between consensus and newly observed read,
        # keep the current sequence
        seq_update = self.sequence.sequence
        qual_update = self.sequence.quality
        qual_other = seq_other.quality
        seq_other = seq_other.sequence
        max_qual = list(map(max, zip(zip(qual_other, [0]*len(qual_other), seq_other),
                                     zip(qual_update, [1]*len(qual_update), seq_update))))
        qual_update = [q[0] for q in max_qual]
        diff = [s != o for s, o in zip(seq_update, seq_other)]
        if any(diff):
            for (i, is_diff) in enumerate(diff):
                # check if new sequence has different nucleotide at this position
                if is_diff:
                    nuc = seq_update[i]
                    nuc_other = seq_other[i]
                    # update diff to record reading discrepancy at this position
                    if self.diffs[i][nuc] == 0:
                        # update for count seen so far
                        self.diffs[i][nuc] = self.size
                    self.diffs[i][nuc_other] += size_other
                elif i in self.diffs:
                    self.diffs[i][seq_update[i]] += size_other
            seq_update = [q[2] for q in max_qual]
            self.sequence.sequence = ''.join(seq_update)
        self.size += size_other

        # regardless of sequence values we will remember the highest quality value
        self.sequence.quality = ''.join(map(max, zip(qual_update, qual_other)))
        return True

    def merge(self, other, tolerance):
        """Merge two consensus sequences.

        Merging will only be attempted if the UIDs don't differ at more than `tolerance`
        positions.

        Args:
            other (:obj:`pyrates.consensus.Consensus`): The consensus sequence to merge into
                this one.
            tolerance (:obj:`int`): The maximum number of differences between UIDs.

        Returns:
            :obj:`bool`: `True` if the sequences were successfully merged, `False` otherwise.
        """
        if self.uid.grosslydifferent(other.uid, len(self.uid), tolerance):
            return False
        return self.update(other.uid, other.sequence, other.size, discard=False)

    def __str__(self):
        diff_str = ''
        diff_pos = sorted(self.diffs)
        for pos in diff_pos:
            diff_str += " " + str(pos)
            for nuc in sorted(self.diffs[pos]):
                diff_str += nuc + str(self.diffs[pos][nuc])
        return "@%d%s\n%s%s\n+\n%s%s" % (self.size, diff_str, self.uid.sequence,
                                         self.sequence.sequence, self.uid.quality,
                                         self.sequence.quality)

    def __repr__(self):
        return "Consensus(uid=%r, sequence=%r, diffs=%r, size=%r)" % \
                         (self.uid, self.sequence, dict(self.diffs), self.size)

def from_fastq(input_file, id_length, adapter):
    """Read FASTQ file to generate consensus sequences.

    Args:
        input_file (:obj:`str`): Name of input file.
        id_length (:obj:`int`): Length of UID sequence at beginning/end of read.
        adapter (:obj:`str`): Adapter sequence.

    Returns:
        :obj:`dict`: Computed consensus sequences.
    """
    logger = utils.get_logger(__name__)
    line_count = 0
    total_skipped = 0
    seq = {}
    id_length = id_length
    adapt_length = id_length + len(adapter)

    open_fun = utils.smart_open(input_file)
    with open_fun(input_file) as fastq:
        for line in fastq:
            # print out some stats as we go
            if logger.isEnabledFor(logging.DEBUG) and (line_count % 100000) == 0:
                logger.debug("reads: %d clusters: %d skipped: %d",
                             line_count/4, len(seq), total_skipped)
            elif (line_count % 4) == 1:
                line = line.rstrip("\n")
                nameid = line[0:id_length] + line[-id_length:]
                sequence = line[adapt_length:-adapt_length]
            elif (line_count % 4) == 3:
                line = line.rstrip("\n")
                qnameid = line[0:id_length] + line[-id_length:]
                qsequence = line[adapt_length:-adapt_length]

                uid = pseq.SequenceWithQuality(nameid, qnameid)
                read_seq = pseq.SequenceWithQuality(sequence, qsequence)
                if nameid not in seq:
                    seq[nameid] = Consensus(uid, read_seq)
                else:
                    success = seq[nameid].update(uid, read_seq)
                    if not success:
                        total_skipped += 1
            line_count += 1
    return seq

def to_fastq(seq, output_file, id_tolerance, merge_size, merge_target):
    """Write consensus sequences to fastq file, optionally merging clusters.

    If merging of small clusters can be enabled by setting all relevant arguments
    (`id_tolerance`, `merge_size`, `merge_target`) to positive values. If this is
    the case clusters smaller than `merge_size` will be compared to all clusters
    with no more than `merge_target` members. If their UIDs don't differ at more
    than `id_tolerance` positions an attempt will be made to merge the two clusters,
    but the attempt may fail if the sequences are too different (using the same
    logic employed when computing the initial consensus).

    Args:
        seq (:obj:`dict`): Read clusters and their consensus sequences, identified
            by the sequence of their UIDs.
        output_file (:obj:`str`): File name for output. Will be replaced if it exists.
        id_tolerance (:obj:`int`): Maximum number of mismatches between UIDs allowed
            for merging of clusters.
        merge_size (:obj:`int`): Maximum size for *small* clusters that will be considered
            for merging.
        merge_target (:obj:`int`): Only clusters that are not larger than this are
            considered when trying to find matches for small clusters.
    """
    logger = utils.get_logger(__name__)
    output_fun = utils.smart_open(output_file)
    with output_fun(output_file, 'w') as output:
        if merge_size and merge_target:
            logger.info('Merging small clusters')
            candidates = []
            targets = []
            merge_count = 0
            for (seq_count, uid) in enumerate(seq):
                if logger.isEnabledFor(logging.DEBUG) and seq_count % 10000 == 0:
                    logger.debug("clusters: %d, merged: %d, small: %d, targets: %d",
                                 seq_count, merge_count, len(candidates), len(targets))
                if seq[uid].size <= merge_size:
                    merged = False
                    for consensus in targets:
                        merged = consensus.merge(seq[uid], id_tolerance)
                        if merged:
                            merge_count += 1
                            break
                    if not merged:
                        ## attempt merging with other candidates
                        remove = None
                        for (i, cand) in enumerate(candidates):
                            merged = seq[uid].merge(cand, id_tolerance)
                            if merged:
                                merge_count += 1
                                remove = i
                                break
                        if merged:
                            del candidates[remove]
                            if seq[uid].size > merge_size:
                                targets.append(seq[uid])
                            else:
                                candidates.append(seq[uid])
                        else:
                            candidates.append(seq[uid])
                elif seq[uid].size <= merge_target:
                    targets.append(seq[uid])
                    processed = []
                    for (i, cand) in enumerate(candidates):
                        if seq[uid].merge(cand, id_tolerance):
                            processed.insert(0, i)
                    for i in processed:
                        del candidates[i]
                        merge_count += 1
                else:
                    output.write(str(seq[uid]) + "\n")
            logger.info('Clusters merged: %d', merge_count)
            logger.info('Small clusters remaining: %d', len(candidates))
            for consensus in targets:
                output.write(str(consensus) + "\n")
            for consensus in candidates:
                output.write(str(consensus) + "\n")
        else:
            for uid in seq:
                output.write(str(seq[uid]) + "\n")
