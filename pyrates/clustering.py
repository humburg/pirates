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
        seq = {}
        id_length = id_length
        adapt_length = id_length + len(adapter)

        id_set = pseq.SequenceStore(id_length)
        id_map = {}

        open_fun = utils.smart_open(input_file)
        with open_fun(input_file) as fastq:
            for (line_count, line) in enumerate(fastq):
                # print out some stats as we go
                if cls._logger.isEnabledFor(logging.DEBUG) and (line_count % 4) == 0:
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
                    if nameid in id_map:
                        nameid = id_map[nameid]
                        total_merged += 1
                    elif nameid not in seq:
                        ## Look for similar IDs that may be candidates for merging
                        id_cands = id_set.search(nameid, raw=True)
                        id_cands = [cand for cand in id_cands if
                                    len(seq[cand].sequence) == len(read_seq)]
                        id_cands = [cand for cand in id_cands if
                                    not seq[cand].sequence.grosslydifferent(read_seq)]
                        id_cands = [(cand, id_set.diff(cand, nameid)) for cand in id_cands]
                        if id_cands:
                            similar_id = min(id_cands, key=lambda x: x[1])
                            if similar_id[1] > threshold:
                                similar_id = None
                            else:
                                similar_id = similar_id[0]
                        else:
                            similar_id = None
                        ## Check that there are no obvious differences between sequences
                        if similar_id is None:
                            seq[nameid] = cons.Consensus(uid, read_seq)
                            id_set.add(nameid)
                            continue
                        else:
                            id_map[nameid] = similar_id
                            nameid = similar_id
                            total_merged += 1
                    success = seq[nameid].update(uid, read_seq)
                    if not success:
                        total_skipped += 1
        return cls(seq, id_set)

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
