# read every line and create consensus

import gzip
import sys

f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")

c = 0
gd = 0
seq = {}
seqcount = {}
different = {}
shorter = {}
longer = {}
diff = {}

# 
# returns true if the two given sequences are very different
#
def grosslydifferent(seq1, seq2):
    d = 0
    if len(seq2) < 10:
        print seq2
    for i in range(10):
       if seq1[i] != seq2[i]:
            d = d + 1
    if d >= 8:
        return True
    return False
    

# determine consensus value for this sequence
# seqE - existing consensus sequence
# seqN - new sequence to be merged into consensus
# qseqE, qseqN - corresponding quality values for the sequences
# count - how many times has this sequence label been seen
# function returns a dictionary
# {"qid": string, "seq": string, "qseq": string}
#
def consensus(qidE, qidN, seqE, seqN, qseqE, qseqN, count, diffs):
    # better do some sanity checking
    if len(qidE) != len (qidN):
        print "Error, mismatch in id quality length, should not happen"
        exit(1)
    
    if len(qseqE) != len(qseqN):
        print "Error, mismatch in sequence quality length, should not happen"
        exit(1)
        
    # step through quality and record highest at each step
    q = ""
    for i in range(len(qidE)):
        if qidE[i] > qidN[i]:
            q += qidE[i]
        else:
            q += qidN[i]
        
    # step through sequence, record highest quality at each step, want to save diffs
    # for changes to the sequence but not the quality
    s = ""
    sq = ""
    for i in range(len(seqE)):
        # regardless of sequence values we will remember the highest quality value
        if qseqE[i] > qseqN[i]:
            sq += qseqE[i]
        else:
            sq += qseqN[i]
            
        # if new sequence character differs then change if it has a higher quality
        # if we change then record a diff
        if seqN[i] != seqE[i]:
            if qseqN[i] > qseqE[i]:
                # change values
                s += seqN[i]
                # update diff
                if i not in diffs:
                    # create diff entry for this position
                    diffs[i] = {}
                    diffs[i]['A'] = 0
                    diffs[i]['G'] = 0
                    diffs[i]['C'] = 0
                    diffs[i]['T'] = 0
                    # update for count seen so far
                    diffs[i][seqE[i]] = count
                    
                diffs[i][seqN[i]] += 1
                continue
        
        # else keeping value, udpate diffs if they exist for this position
        s += seqE[i]
        if i in diffs:
            diffs[i][seqE[i]] += 1
        
    return {"qid": q, "seq": s, "qseq": sq}
    
    
#
# Main loop for reading file
#
for line in f:
    # print out some stats as we go
    if (c % 101) == 100:
        print >> sys.stderr, c
        print >> sys.stderr, len(seq)
        break

    if (c % 4) == 1:
        line = line.rstrip("\n")
        nameid = line[0:8] + line[-8:]
        sequence = line[12:-12]
        c = c + 1
        continue

    if (c % 4) == 3:
        line = line.rstrip("\n")
        qnameid = line[0:8] + line[-8:]
        qsequence = line[12:-12]
        c = c + 1

        # if not already seen then record and move on, create zero diffs dict
        if nameid not in seq:
            seq[nameid] = [qnameid, sequence, qsequence]
            seqcount[nameid] = 1
            diff[nameid] = {}
            continue

        # we have this id already recorded so must process for possible consensus
        #
        
        # if grossly different then just count this and move on
        if grosslydifferent(sequence, seq[nameid][1]):
            if nameid not in different:
                different[nameid] = 0
                
            different[nameid] += 1
            continue 

        # if sequence length is shorter, count this occurance, abandon this sequence and move on
        if len(sequence) < len(seq[nameid][1]):
            if nameid not in shorter:
                shorter[nameid] = 0
            
            shorter[nameid] += 1
            continue 

        # if new sequence is longer, count this occurance 
        # replace consensus sequence if built from only one other sequence
        if len(sequence) > len(seq[nameid][1]):
            if nameid not in longer:
                longer[nameid] = 0
            longer[nameid] += 1
            
            # if consensus built from just one sequence then replace
            if seqcount[nameid] == 1:
                seq[nameid] = [qnameid, sequence, qsequence]
                continue 
                
            # else abandon this sequence
            continue

        # do consensus
        #

        result = consensus(seq[nameid][0], qnameid, seq[nameid][1], sequence, seq[nameid][2], qsequence, seqcount[nameid], diff[nameid])
        # some prints for consensus debug
        #print "Calling consensus"
        #print "Existing, new and result"
        #print seq[nameid][0], seq[nameid][1], seq[nameid][2]
        #print qnameid, sequence, qsequence        
        #print result['qid'], result['seq'], result['qseq']
        #if nameid in diff:
        #    print diff[nameid]
        #raw_input("Press key to continue")
        
        # {"qid": string, "seq": string, "qseq": string}
        seq[nameid][0] = result['qid']
        seq[nameid][1] = result['seq']
        seq[nameid][2] = result['qseq']
        
        # update count for this sequence label
        seqcount[nameid] += 1
        if len(diff[nameid]) > 1:
            break
        continue
        
    c = c + 1

print >> sys.stderr,  "Number of consensus sequence with unique labels", len(seq)
print >> sys.stderr,  "Number sequences grossly difference from consensus with same label", len(different)
print >> sys.stderr,  "Number of sequences that were shorter than consensus sequence ", len(shorter)
print >> sys.stderr,  "Number of sequences that were longer then consensus sequence ", len(longer)

#
# Print everything out
# '@'int [int,A,int,C,int,T,int,G,int ...]
#
for label in sorted(seqcount, key=lamba x: d[x]):
    name = "@" + str(seqcount[label])
    # loop over any diffs and add to name string
    for pos in diff[label]:
        name += " " + str(pos)
        for c in ['A', 'C', 'T', 'G']:
            name += c + str(diff[label][pos][c])
            
    # print out sequence data
    print name
    print label + seq[label][1]
    print "+"
    print seq[label][0] + seq[label][2]

