"""CLustering of reads to create consensus sequences.
"""

import logging
from collections import defaultdict
import pyrates.utils as utils
import pyrates.sequence as pseq
import pyrates.consensus as cons

class Clustering(object):
    """Clustering of reads with UIDs.

    Args:
        centres (:obj:`dict`): Cluster centres represented by consensus
            sequences and identified by the associated UID.
        store (:obj:`pyrates.sequence.SequenceStore`, optional): Precomputed set of UID sequences.
            If this is missing it will be computed from the cluster centres.

    Attributes:
        clusters (:obj:`dict`): Cluster centres represented by consensus
            sequences and identified by the associated UID.
    """
    __slots__ = 'clusters', '_store'
    _logger = utils.get_logger(__name__)

    def __init__(self, centres, store=None):
        self.clusters = centres
        if store is not None:
            self._store = store
        else:
            self._store = pseq.SequenceStore.from_list(list(centres.keys()))

    def _filter(self, pattern, candidates, read_seq):
        candidates = [cand for cand in candidates if
                      len(self[cand].sequence) == len(read_seq)]
        candidates = [cand for cand in candidates if
                      not self[cand].sequence.grosslydifferent(read_seq)]
        candidates = [(cand, pseq.SequenceStore.diff(cand, pattern)) for cand in candidates]
        return candidates

    def merge_target(self, uid, read_seq, id_map, threshold, wildcard):
        """Compute set of candidate clusters for a given read.

        Args:
            uid (:obj:`pyrates.sequence.SequenceWithQuality`): UID sequence.
            read_seq (:obj:`pyrates.sequence.SequenceWithQuality`): Read sequence.
            id_map (:obj:`dictionary`): A mapping of known approximate matches for UIDs.
            threshold (:obj:`int`): Maximum number of differences allowed between UIDs.
            wildcard (:obj:`string`): Single character that should be treated as wildcard or
                `None` to disable use of wildcard matching.

        Returns:
            :obj:`string`: Either the best approximate match for the UID or `None`
                if no valid match was found.
        """
        nameid = uid.sequence
        id_cands = self._store.search(nameid, raw=True, wildcard=wildcard)
        id_cands = self._filter(nameid, id_cands, read_seq)
        if id_cands:
            similar_id = min(id_cands, key=lambda x: x[1])
            if similar_id[1] > threshold:
                similar_id = None
            else:
                similar_id = similar_id[0]
        else:
            similar_id = None
        ## Create new cluster or merge with existing consensus
        if similar_id is None:
            self.clusters[nameid] = cons.Consensus(uid, read_seq)
            self._store.add(nameid, wildcard=wildcard)
        else:
            id_map[nameid] = similar_id
        return similar_id

    @classmethod
    def from_fastq(cls, input_file, id_length, adapter, threshold=5):
        """Read FASTQ file to generate consensus sequences.

        Args:
            input_file (:obj:`str`): Name of input file.
            id_length (:obj:`int`): Length of UID sequence at beginning/end of read.
            adapter (:obj:`str`): Adapter sequence.
            threshold (:obj:`int`): Maximum number of differences allowed between UIDs.
        Returns:
            :obj:`dict`: Computed consensus sequences.
        """
        total_skipped = 0
        total_merged = 0

        id_set = pseq.SequenceStore(id_length*2)
        id_wild = defaultdict(list)
        id_map = {}
        seq = cls({}, id_set)

        adapt_length = id_length + len(adapter)
        open_fun = utils.smart_open(input_file)
        with open_fun(input_file) as fastq:
            for (line_count, line) in enumerate(fastq):
                # print out some stats as we go
                if cls._logger.isEnabledFor(logging.DEBUG) and (line_count % 10000) == 0:
                    cls._logger.debug("reads: %d clusters: %d merged: %d skipped: %d",
                                      line_count/4, len(seq), total_merged, total_skipped)
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
                    if nameid.count('N'):
                        id_wild[nameid].append((uid, read_seq))
                        continue
                    ## Look for similar IDs that may be candidates for merging
                    similar_id = None
                    if nameid in id_map:
                        similar_id = id_map[nameid]
                    elif nameid not in seq:
                        similar_id = seq.merge_target(uid, read_seq, id_map, threshold, None)
                    else:
                        similar_id = nameid
                    if similar_id is not None:
                        success = seq[similar_id].update(uid, read_seq)
                        if success:
                            total_merged += 1
                        else:
                            total_skipped += 1
                    else:
                        seq.add(uid, read_seq, None)
            for nameid in id_wild:
                similar_id = seq.merge_target(uid, read_seq, id_map, threshold, 'N')
                if similar_id is not None:
                    success = seq[similar_id].update(uid, read_seq)
                    if success:
                        total_merged += 1
                    else:
                        total_skipped += 1
                else:
                    seq.add(uid, read_seq, 'N')
        return seq

    def find(self, uid, tolerance=3):
        """Find best approximate match for UIDs.

        Args:
            uid (:obj:`str`): UID for which the best match should be returned.
            tolerance (:obj:`int`): Number of mismatcheds that should be tolerated.

        Returns:
            :obj:`str`: If a match was found the best matching UID in the Clustering
                is returned, otherwise :obj:`None`.
        """
        if tolerance is None or tolerance == 0:
            if uid in self:
                return uid
            else:
                return None
        else:
            return self._store.find(uid, tolerance)

    def write(self, output_file):
        """Write consensus sequences to fastq file.

        Args:
            output_file (:obj:`str`): File name for output. Will be replaced if it exists.
        """
        output_fun = utils.smart_open(output_file)
        with output_fun(output_file, 'w') as output:
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

    def add(self, uid, sequence, wildcard=None):
        """Add a new cluster centre.

        Args:
            uid (:obj:`pyrates.sequence.SequenceWithQuality`): UID for the new cluster.
            sequence (:obj:`pyrates.sequence.SequenceWithQuality`): Sequence to represent cluster.
            wildcard (:obj:`str`): Wildcard character to use, `None` to disable wildcard matching.
        """
        nameid = uid.sequence
        self.clusters[nameid] = cons.Consensus(uid, sequence)
        self._store.add(nameid, wildcard)

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
