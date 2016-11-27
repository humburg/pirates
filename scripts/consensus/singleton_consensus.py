# read every line, match singletons to groups

import gzip

fname = "data/fastq/Jurkat_only.assembled.fastq.gz"
f = gzip.open(fname, "r")

c = 0

group_dict = {} 
singleton_dict = {}
diff = {}

def grosslydifferent(seq1, seq2, max_diff=8):
    d = 0
    for i in range(10):
       if seq1[i] != seq2[i]:
            d = d + 1
    if d >= max_diff:
        return True
    return False

for line in f:
    if (c % 100000) == 0:
        print c

    line = line.rstrip("\n")

    # count line
    if (c % 4) == 0: 
        # line[0] is @
        count = line[1]
    # id + sequence line
    elif (c % 4) == 1:
        nameid = line[0:16]
        sequence = line[16:]
    # quality of id + quality of sequence line
    elif (c % 4) == 3:
        qnameid = line[0:16]
        qsequence = line[16:]
        
    # if a group already, add to group_dict
    if count > 1:
        # group_dict field is list with count, sequence,
        # qnameid, qsequence, 
        # count of how many singletons merged into group and consensus updated,
        # count of how many singletons merged into group but with different length
        group_dict[nameid] = [
            count,sequence,
            qnameid,qsequence,
            0,0
        ]
        continue
    # if a singleton
    for ref_nameid in group_dict:
        maxdiff = 5
        # compare if nameid's are grossly different
        if grosslydifferent(ref_nameid,nameid,maxdiff):
            # put into pool of singletons
            # will compare these singletons after group comparisons
            singleton_dict[nameid] = [
                count,sequence,
                qnameid,qsequence,
                0
            ]
            continue
        ## check if sequences are different lengths
        ref_seq = group_dict[ref_nameid][1]
        ref_qual = group_dict[ref_nameid][3]
        # if sequences are the same lengths, create consensus
        if len(ref_seq) == len(sequence):
           create_consensus(ref_nameid,ref_seq,group_qual,sequence,qsequence)
           group_dict[ref_nameid][4] += 1 
           continue
        # if the sequences are different lengths...
        elif abs(len(ref_seq) -  len(sequence)) < 10:
            # and not too different, update count in group_dict
            if not grosslydifferent(ref_seq,sequence):
                group_dict[ref_nameid][5] += 1 
                continue
            else:
                # put into singleton pool
                singleton_dict[nameid] = [
                    count,sequence,
                    qnameid,qsequence,
                    0
                ]

def create_consensus(ref_nameid,ref_seq,ref_qual,new_seq,new_qual):
    """
    Create consensus sequence out of two sequences.
    - ref_seq is our initial consensus sequence.
    - new_seq is the new sequence, which we update the consensus with.
    """
    # might need to do another check on length
    # what happens if they are within 10 chars in length
    # but new_seq is shorter than ref_seq?
    # error will occur here
    # error will also happen if we do range(len(new_seq))
    # and ref_seq is shorter

    # additionally should be comparing over full length of sequence?
    for i in range(len(ref_seq)):
        # compare sequence chars, if same go to next char
        if ref_seq[i] == new_seq[i]:
            continue
        # compare sequence quals, if 
        if new_qual > ref_qual:
            # update sequence
            group_dict[ref_nameid][1] = ref_seq[0:i] + new_seq[i] + ref_seq[i+1:]
            # update quality
            group_dict[ref_nameid][3] = ref_qual[0:i] + new_qual[i] + ref_qual[i+1:]
        if nameid in diff:
            diff[nameid] += 1
        else:
            diff[nameid] = 1
            continue 
        continue
    c =c + 1

print "Length consensus ", len(seqdict)
print "Number grossly different ", len(differentdict)
print "Number of shorter ", len(shorterdict)
print "Number of diffs ", len(diff)

raw_input("How big am I")


