"""Functions to facilitate consensus calling.
"""

import pyrates.utils as utils

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

def update_qual(qual_cur, qual_new):
    """Update quality sequence.

    Args:
        qual_cur (:obj:`str`): Qualities for stored sequence.
        qual_new (:obj:`str`): Qualities for newly observed sequence.

    Returns:
        :obj:`str`: Updated quality scores.
    """
    qual_update = ""
    for (qual1, qual2) in zip(qual_cur, qual_new):
        if qual1 > qual2:
            qual_update += qual1
        else:
            qual_update += qual2
    return qual_update

def consensus(id_qual_cur, id_qual_new, seq_cur, seq_new, qual_cur, qual_new, count, diffs):
    """Determine consensus value for this sequence.

    Args:
        id_qual_cur (:obj:`str`): Qualities for stored ID sequence.
        id_qual_new (:obj:`str`): Qualities for ID of new sequence.
        seq_cur (:obj:`str`): Existing consensus sequence.
        seq_new (:obj:`str`): New sequence to be merged into consensus.
        qual_cur (:obj:`str`): Quality value for existing consensus.
        qual_new (:obj:`str`): Quality values for the new sequence.
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

    logger = utils.get_logger(__name__)
    # better do some sanity checking
    if len(id_qual_cur) != len(id_qual_new):
        logger.error("Mismatch in id quality length, this should not happen. Check your input.")
        logger.debug("Mismatching quality strings were '%s' and '%s'", id_qual_cur, id_qual_new)
        return {'qid':id_qual_cur, 'seq':seq_cur, 'qseq':qual_cur}

    if len(qual_cur) != len(qual_new):
        logger.error("Mismatch in sequence quality length, this should not happen." + \
                     " Check your input.")
        logger.debug("Mismatching sequence qualities were '%s' and '%s'", qual_cur, qual_new)
        return {'qid':id_qual_cur, 'seq':seq_cur, 'qseq':qual_cur}

    # step through quality and record highest at each step
    id_qual_update = update_qual(id_qual_cur, id_qual_new)

    # step through sequence, record highest quality at each step, want to save diffs
    # for changes to the sequence but not the quality
    seq_update = ""
    qual_update = ""
    for (i, (qual, nuc)) in enumerate(zip(qual_cur, seq_cur)):
        # regardless of sequence values we will remember the highest quality value
        if qual_cur[i] > qual_new[i]:
            qual_update += qual
        else:
            qual_update += qual_new[i]
        # check if new sequence has different nucliotide value at this position
        if seq_new[i] != nuc:
            # if this position in the new sequence has higher quality reading than
            # the consensus then swap it in
            if qual_new[i] > qual:
                seq_update += seq_new[i]
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
            diffs[i][seq_new[i]] += 1
            continue

        # sequences agree at this psoition, if diff exists then update
        seq_update += nuc
        if i in diffs:
            diffs[i][nuc] += 1

    return {"qid": id_qual_update, "seq": seq_update, "qseq": qual_update}
    