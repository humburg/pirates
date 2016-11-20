"""CLustering of reads to create consensus sequences.
"""

import logging
import ngram
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
    __slots__ = 'clusters', '_store'
    _logger = utils.get_logger(__name__)

    def __init__(self, centres):
        self.clusters = centres
        self._store = ngram.NGram(centres.keys(), N=3)

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

    def find(self, uid, tolerance=3):
        """Find best approximate match for UIDs.

        Args:
            uid (:obj:`str`): UID for which the best match should be returned.
                threshold (:obj:`float`): Minimum similarity required for match.
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
            ngram_count = len(uid) - 2
            threshold = (ngram_count - 3*tolerance)/float(ngram_count)
            return self._store.find(uid, threshold)

    def search(self, uid, max_hits=1, tolerance=3):
        """Find all approximate matches.

        Args:
            uid (:obj:`str`): UID for which the best match should be returned.
            threshold (:obj:`float`): Minimum similarity required for match.
            tolerance (:obj:`int`): Number of mismatcheds that should be tolerated.
            max_hits (:obj:`int`): Maximum number of hits to return.

        Returns:
            :obj:`list`: All identified matches.
        """
        ngram_count = len(uid) - 2
        threshold = (ngram_count - 3*tolerance)/float(ngram_count)
        hits = self._store.search(uid, threshold)
        if hits:
            hits = hits[1:]
        if len(hits) > max_hits:
            hits = hits[:max_hits]
        return hits

    def merge(self, id_tolerance, merge_size):
        """Merge small clusters into larger ones.

        For each cluster not larger than `merge_size` an attempt is made to
        identify a larger cluster with similar UID and sequence. If UIDs
        are sufficiently similar the merging may still fail if the actual
        sequences are too different (using the same logic employed when computing
        the initial consensus).

        Args:
            id_tolerance (:obj:`int`): Maximum number of mismatches between UIDs allowed
                for merging of clusters.
            merge_size (:obj:`int`): Maximum size for *small* clusters that will be considered
                for merging.
        """
        remove = set()
        merge_count = 0
        small_count = 0
        for (i, uid) in enumerate(self.clusters, start=1):
            cur_cluster = self.clusters[uid]
            if cur_cluster.size <= merge_size:
                candidate = self.search(uid, 1, id_tolerance)
                candidate = [c[0] for c in candidate if c not in remove]
                if candidate:
                    candidate = candidate[0]
                    cand_cluster = self[candidate]
                    success = cand_cluster.merge(self[uid], id_tolerance)
                    if success and candidate:
                        remove.add(uid)
                        merge_count += 1
                    else:
                        small_count += 1
            if self._logger.isEnabledFor(logging.DEBUG) and i % 10000 == 0:
                self._logger.debug("clusters: %d, small: %d, merged: %d",
                                   i, small_count + merge_count, merge_count)
        for uid in remove:
            del self.clusters[uid]
            self._store.remove(uid)

        self._logger.info('Clusters merged: %d', merge_count)
        self._logger.info('Small clusters remaining: %d', small_count)

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
