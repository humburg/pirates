"""CLustering of reads to create consensus sequences.
"""

import logging
import time
import datetime
import os.path
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
        wildcard (:obj:`string`): Single character that should be treated as wildcard or
                `None` to disable use of wildcard matching.

    Attributes:
        clusters (:obj:`dict`): Cluster centres represented by consensus
            sequences and identified by the associated UID.
    """
    __slots__ = 'clusters', '_store', 'stats'
    _logger = utils.get_logger(__name__)

    def __init__(self, centres, store=None, wildcard=None,
                 alphabet=('A', 'C', 'G', 'T'), tag_size=5, max_diff=3,
                 read_length=None):
        self.clusters = centres
        ## keep track of UID handling for fragments that are shorter/longer than read length
        created_at = time.time()
        self.stats = {
            'read_length':read_length,
            'total_skipped':[0, 0],
            'total_merged':[0, 0],
            'total_fixed': [0, 0],
            'single_count':[0, 0],
            'start_time':created_at,
            'batch_start':created_at
        }
        if store is not None:
            self._store = store
        elif tag_size > 0:
            self._store = pseq.GroupedSequenceStore.from_list(list(centres.keys()),
                                                              alphabet=alphabet,
                                                              tag_size=tag_size,
                                                              max_diff=max_diff,
                                                              wildcard=wildcard)
        else:
            self._store = pseq.SequenceStore.from_list(list(centres.keys()),
                                                       alphabet=alphabet)

    def _filter(self, pattern, candidates, read_seq, threshold):
        candidates = [cand for cand in candidates if
                      len(self[cand].sequence) == len(read_seq)]
        candidates = [cand for cand in candidates if
                      not self[cand].sequence.grosslydifferent(read_seq)]
        candidates = [(cand, pseq.SequenceStore.diff(cand, pattern)) for
                      cand in candidates]
        candidates = [cand for cand in candidates if cand[1] <= threshold]
        return candidates

    def merge_target(self, uid, read_seq, id_map, threshold):
        """Compute set of candidate clusters for a given read.

        Args:
            uid (:obj:`pyrates.sequence.SequenceWithQuality`): UID sequence.
            read_seq (:obj:`pyrates.sequence.SequenceWithQuality`): Read sequence.
            id_map (:obj:`dictionary`): A mapping of known approximate matches for UIDs.
            threshold (:obj:`int`): Maximum number of differences allowed between UIDs.

        Returns:
            :obj:`string`: Either the best approximate match for the UID or `None`
                if no valid match was found.
        """
        nameid = uid.sequence
        id_cands = self._store.search(nameid, max_hits=100, raw=True)
        id_cands = self._filter(nameid, id_cands, read_seq, threshold)
        if id_cands:
            similar_id = min(id_cands, key=lambda x: x[1])
            similar_id = similar_id[0]
        else:
            similar_id = None
        ## Create new cluster or merge with existing consensus
        if similar_id is None:
            self.clusters[nameid] = cons.Consensus(uid, read_seq)
            self._store.add(nameid)
        else:
            id_map[nameid] = similar_id
        return similar_id

    @classmethod
    def from_fastq(cls, input_file, id_length, adapter, threshold=5, prefix=5, read_length=None):
        """Read FASTQ file to generate consensus sequences.

        Args:
            input_file (:obj:`str`): Name of input file.
            id_length (:obj:`int`): Length of UID sequence at beginning/end of read.
            adapter (:obj:`str`): Adapter sequence.
            threshold (:obj:`int`, optional): Maximum number of differences allowed between UIDs.
            prefix (:obj:`int`, optional): Length of UID prefix to use in clustering algorithm.
            read_length (:obj:`int`, optional): Original read length used. If this is set and
                and the logging level is sufficiently high, additional log entries are generated
                to track the number of short and long fragments processed.
        Returns:
            :obj:`dict`: Computed consensus sequences.
        """
        adapt_length = id_length + len(adapter)
        if read_length is not None:
            max_short = read_length - id_length - len(adapter)
        else:
            max_short = 0
        name = os.path.basename(input_file).split('.')[0]

        id_set = pseq.GroupedSequenceStore(id_length*2, tag_size=prefix, max_diff=threshold,
                                           wildcard='N')
        id_map = {}
        seq = cls({}, id_set, read_length=read_length)

        open_fun = utils.smart_open(input_file)
        with open_fun(input_file) as fastq:
            for (line_count, line) in enumerate(fastq):
                # print out some stats as we go
                if cls._logger.isEnabledFor(logging.DEBUG) and line_count > 0 and \
                        (line_count % 10000) == 0:
                    seq.log_progress(line_count)
                elif (line_count % 4) == 1:
                    line = line.rstrip("\n")
                    nameid = line[0:id_length] + line[-id_length:]
                    sequence = line[adapt_length:-adapt_length]
                    is_long = len(sequence) > max_short
                elif (line_count % 4) == 3:
                    line = line.rstrip("\n")
                    qnameid = line[0:id_length] + line[-id_length:]
                    qsequence = line[adapt_length:-adapt_length]

                    uid = pseq.SequenceWithQuality(nameid, qnameid)
                    read_seq = pseq.SequenceWithQuality(sequence, qsequence, name=name)
                    ## Look for similar IDs that may be candidates for merging
                    similar_id = None
                    if nameid in id_map:
                        similar_id = id_map[nameid]
                        seq.stats['total_fixed'][is_long] += 1
                        id_matched = False
                    elif nameid not in seq:
                        similar_id = seq.merge_target(uid, read_seq, id_map, threshold)
                        id_matched = False
                        if similar_id is not None:
                            seq.stats['total_fixed'][is_long] += 1
                    else:
                        similar_id = nameid
                        id_matched = True
                    if similar_id is not None:
                        success = seq[similar_id].update(uid, read_seq)
                        if success:
                            if not id_matched:
                                seq.stats['total_merged'][is_long] += 1
                            if seq[similar_id].size == 2:
                                seq.stats['single_count'][is_long] -= 1
                        else:
                            seq.stats['total_skipped'][is_long] += 1
                    else:
                        seq.add(uid, read_seq)
                        seq.stats['single_count'][is_long] += 1
        return seq

    def log_progress(self, line_count):
        checkpoint = time.time()
        total_time = checkpoint - self.stats['start_time']
        batch_time = checkpoint - self.stats['batch_start']
        self.stats['batch_start'] = checkpoint
        self._logger.debug("reads: %d, clusters: %d, singletons: %d (%.1f%%), " +
                           "corrupted UIDs: %d (%.2f%%)",
                           line_count/4, len(self), sum(self.stats['single_count']),
                           sum(self.stats['single_count'])/(len(self)/100.0),
                           self.fail_count,
                           self.fail_count/(line_count/4.0)*100)
        if self.stats['read_length'] is not None:
            self._logger.debug("singletons (short/long): %d %d",
                               self.stats['single_count'][0], self.stats['single_count'][1])
        self._logger.debug("similar UIDs: %d (%.1f%%), UIDs merged: %d (%.1f%%), " +
                           "merge failures: %d (%.1f%%)",
                           sum(self.stats['total_fixed']),
                           sum(self.stats['total_fixed'])/(line_count/4.0)*100,
                           sum(self.stats['total_merged']),
                           sum(self.stats['total_merged'])/(line_count/4.0)*100,
                           sum(self.stats['total_skipped']),
                           sum(self.stats['total_skipped'])/(line_count/4.0)*100)
        if self.stats['read_length'] is not None:
            self._logger.debug("similar UIDs (short/long): %d %d",
                               self.stats['total_fixed'][0],
                               self.stats['total_fixed'][1])
            self._logger.debug("merged UIDs (short/long): %d %d",
                               self.stats['total_merged'][0],
                               self.stats['total_merged'][1])
            self._logger.debug("merge failures (short/long): %d %d",
                               self.stats['total_skipped'][0],
                               self.stats['total_skipped'][1])
        self._logger.debug("total time: %s, increment: %s, rate: %.1f reads/s",
                           datetime.timedelta(seconds=total_time),
                           datetime.timedelta(seconds=batch_time),
                           line_count/4.0/total_time)

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

    def add(self, uid, sequence):
        """Add a new cluster centre.

        Args:
            uid (:obj:`pyrates.sequence.SequenceWithQuality`): UID for the new cluster.
            sequence (:obj:`pyrates.sequence.SequenceWithQuality`): Sequence to represent cluster.
        """
        nameid = uid.sequence
        self.clusters[nameid] = cons.Consensus(uid, sequence)
        self._store.add(nameid)

    def cluster_fails(self):
        """Sequences that were not assigned to any cluster.

        Returns:
            :obj:`dict`: UID / sequence pairs
        """
        return {key:self.clusters[key] for key in self._store.wild_tags}

    @property
    def fail_count(self):
        """Number of sequences that were not assigned to any cluster.
        """
        return len(self._store.wild_tags)

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
