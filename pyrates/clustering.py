"""CLustering of reads to create consensus sequences.
"""

import logging
import pyrates.utils as utils
import pyrates.sequence as pseq
import pyrates.consensus as cons

class Clustering(object):
    """Clustering of reads with UIDs.

    Args:
        centres (:obj:`dict`): Cluster centres represented by consensus
            sequences and identified by the associated UID.

    Attributes:
        clusters (:obj:`dict`): Cluster centres represented by consensus
            sequences and identified by the associated UID.
    """
    __slots__ = 'clusters'
    _logger = utils.get_logger(__name__)

    def __init__(self, centres):
        self.clusters = centres

    @classmethod
    def from_fastq(cls, input_file, id_length, adapter):
        """Read FASTQ file to generate consensus sequences.

        Args:
            input_file (:obj:`str`): Name of input file.
            id_length (:obj:`int`): Length of UID sequence at beginning/end of read.
            adapter (:obj:`str`): Adapter sequence.

        Returns:
            :obj:`dict`: Computed consensus sequences.
        """
        line_count = 0
        total_skipped = 0
        seq = {}
        id_length = id_length
        adapt_length = id_length + len(adapter)

        open_fun = utils.smart_open(input_file)
        with open_fun(input_file) as fastq:
            for line in fastq:
                # print out some stats as we go
                if cls._logger.isEnabledFor(logging.DEBUG) and (line_count % 100000) == 0:
                    cls._logger.debug("reads: %d clusters: %d skipped: %d",
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
                        seq[nameid] = cons.Consensus(uid, read_seq)
                    else:
                        success = seq[nameid].update(uid, read_seq)
                        if not success:
                            total_skipped += 1
                line_count += 1
        return cls(seq)

    def write(self, output_file, id_tolerance, merge_size, merge_target):
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
        output_fun = utils.smart_open(output_file)
        with output_fun(output_file, 'w') as output:
            if merge_size and merge_target:
                self._logger.info('Merging small clusters')
                candidates = []
                targets = []
                merge_count = 0
                for (seq_count, uid) in enumerate(self):
                    if self._logger.isEnabledFor(logging.DEBUG) and seq_count % 10000 == 0:
                        self._logger.debug("clusters: %d, merged: %d, small: %d, targets: %d",
                                           seq_count, merge_count, len(candidates), len(targets))
                    if self[uid].size <= merge_size:
                        merged = False
                        for consensus in targets:
                            merged = consensus.merge(self[uid], id_tolerance)
                            if merged:
                                merge_count += 1
                                break
                        if not merged:
                            ## attempt merging with other candidates
                            remove = None
                            for (i, cand) in enumerate(candidates):
                                merged = self[uid].merge(cand, id_tolerance)
                                if merged:
                                    merge_count += 1
                                    remove = i
                                    break
                            if merged:
                                del candidates[remove]
                                if self[uid].size > merge_size:
                                    targets.append(self[uid])
                                else:
                                    candidates.append(self[uid])
                            else:
                                candidates.append(self[uid])
                    elif self[uid].size <= merge_target:
                        targets.append(self[uid])
                        processed = []
                        for (i, cand) in enumerate(candidates):
                            if self[uid].merge(cand, id_tolerance):
                                processed.insert(0, i)
                        for i in processed:
                            del candidates[i]
                            merge_count += 1
                    else:
                        output.write(str(self[uid]) + "\n")
                self._logger.info('Clusters merged: %d', merge_count)
                self._logger.info('Small clusters remaining: %d', len(candidates))
                for consensus in targets:
                    output.write(str(consensus) + "\n")
                for consensus in candidates:
                    output.write(str(consensus) + "\n")
            else:
                for uid in self:
                    output.write(str(self[uid]) + "\n")

    def keys(self):
        """UIDs used to identify clusters."""
        return self.clusters.keys()

    def values(self):
        """Consensus sequences corresponding to clusters."""
        return self.clusters.values()

    def items(self):
        """UIDs / consensus sequence pairs."""
        return self.clusters.items()

    def iterkeys(self):
        """UIDs used to identify clusters."""
        return self.clusters.iterkeys()

    def itervalues(self):
        """Consensus sequences corresponding to clusters."""
        return self.clusters.itervalues()

    def iteritems(self):
        """UIDs / consensus sequence pairs."""
        return self.clusters.iteritems()

    def has_key(self, key):
        """Test for presence of UID"""
        return key in self.clusters

    def __contains__(self, item):
        return item in self.clusters

    def __getitem__(self, key):
        return self.clusters[key]

    def __iter__(self):
        return iter(self.clusters)

    def __len__(self):
        return len(self.clusters)

    def __str__(self):
        str_lst = [str(c) for c in self]
        return "\n".join(str_lst)

    def __repr__(self):
        return "Clustering(centres=%r)" % self.clusters
