"""Functions to facilitate consensus calling.
"""

import logging

LOG = logging.getLogger(__name__)

def grosslydifferent(seq1, seq2):
    """Partial comparison of sequences to determine whether
    they are substantially different.

    Args:
        seq1 (string): First sequence.
        seq2 (string): Second sequence.

    Returns:
        bool: True if the sequences are considered too different
        to warrant further comparison.
    """
    diff = 0
    for i in range(10):
        if seq1[i] != seq2[i]:
            diff = diff + 1
    if diff >= 8:
        return True
    return False

def consensus(qidE, qidN, seqE, seqN, qseqE, qseqN, count, diffs):
    """Determine consensus value for this sequence.

    Args:
        qidE (:obj:`str`): Qualities for stored ID sequence.
        qidN (:obj:`str`): Qualities for ID of new sequence.
        seqE (:obj:`str`): Existing consensus sequence.
        seqN (:obj:`str`): New sequence to be merged into consensus.
        qseqE (:obj:`str`): Quality value for existing consensus.
        qseqN (:obj:`str`): Quality values for the new sequence.
        count (int): How many times has this sequence label been seen.
        diffs (:obj:`dict`): Sequence differences from the consensus
            observed so far.

    Returns:
        :obj:`dict`:
            {
                'qid': Cluster ID quality,
                'seq': Consensus sequence,
                'qseq': Consensus quality
            }
    """
    # better do some sanity checking
    if len(qidE) != len(qidN):
        LOG.error("Mismatch in id quality length, this should not happen. Check your input.")
        LOG.debug("Mismatching quality strings were '" + qidE + "' and '" + qidN)
        return {'qid':qidE, 'seq':seqE, 'qseq':qseqE}

    if len(qseqE) != len(qseqN):
        LOG.error("Mismatch in sequence quality length, this should not happen." + \
                     " Check your input.")
        LOG.debug("Mismatching sequence qualities were '" + qseqE + "' and '" + qseqN)
        return {'qid':qidE, 'seq':seqE, 'qseq':qseqE}

    # step through quality and record highest at each step
    id_qual_update = ""
    for (i, qual) in enumerate(qidE):
        if qual > qidN[i]:
            id_qual_update += qual
        else:
            id_qual_update += qidN[i]

    # step through sequence, record highest quality at each step, want to save diffs
    # for changes to the sequence but not the quality
    seq_update = ""
    qual_update = ""
    for (i, nuc) in enumerate(seqE):
        # regardless of sequence values we will remember the highest quality value
        if qseqE[i] > qseqN[i]:
            qual_update += nuc
        else:
            qual_update += qseqN[i]
        # check if new sequence has different nucliotide value at this position
        if seqN[i] != nuc:
            # if this position in the new sequence has higher quality reading than
            # the consensus then swap it in
            if qseqN[i] > qseqE[i]:
                seq_update += seqN[i]
            else:
                seq_update += nuc
            # update diff to record reading discrepancy at this position
            if i not in diffs:
                # create diff entry for this position
                diffs[i] = {}
                diffs[i]['A'] = 0
                diffs[i]['G'] = 0
                diffs[i]['C'] = 0
                diffs[i]['T'] = 0
                diffs[i]['N'] = 0
                # update for count seen so far
                diffs[i][nuc] = count
            diffs[i][seqN[i]] += 1
            continue

        # sequences agree at this psoition, if diff exists then update
        seq_update += nuc
        if i in diffs:
            diffs[i][nuc] += 1

    return {"qid": id_qual_update, "seq": seq_update, "qseq": qual_update}
    