"""Functions to facilitate consensus calling.
"""

from collections import defaultdict
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

    def update(self, uid_other, seq_other, size_other=1, diffs_other=None, discard=True):
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
            diffs_other (:obj:`dict`, optional): Differences already recorded
                for other sequence.
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
        if diffs_other is None:
            diffs_other = {}
        if diffs_other:
            for i in diffs_other:
                if i not in self.diffs:
                    self.diffs[i][seq_update[i]] += self.size
                for nuc in diffs_other[i]:
                    self.diffs[i][nuc] += diffs_other[i][nuc]
        if any(diff):
            for (i, is_diff) in enumerate(diff):
                # check if new sequence has different nucleotide at this position
                if is_diff:
                    if i not in diffs_other:
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
        elif self.diffs and not diffs_other:
            for i in self.diffs:
                self.diffs[i][seq_update[i]] += size_other
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
        return self.update(other.uid, other.sequence, other.size, diffs_other=other.diffs,
                           discard=False)

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
