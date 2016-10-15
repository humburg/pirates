# read every line and create consensus

import gzip
f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")

c = 0
gd = 0
seqdict = {}
seqcount = {}
differentdict = {}
shorterdict = {}
diff = {}

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

        

for line in f:
    if (c % 100000) == 0:
        print c

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

        # check if in dictionary
        if nameid not in seqdict:
            seqdict[nameid] = [qnameid, sequence, qsequence]
            seqcount[nameid] = 1
            continue

        # if grossly different then record and move on
        if grosslydifferent(sequence, seqdict[nameid][1]):
            if nameid in differentdict:
                differentdict[nameid] += 1
            else:
                differentdict[nameid] = 1
            continue 

        # if sequence length is shorter, record and move on
        if len(sequence) < len(seqdict[nameid][1]):
            if nameid in shorterdict:
                shorterdict[nameid] += 1
            else:
                shorterdict[nameid] = 1
            continue 

        # if sequence is longer, record, replace existing and move on
        if len(sequence) > len(seqdict[nameid][1]):
            seqdict[nameid] = [qnameid, sequence, qsequence]
            if nameid in shorterdict:
                shorterdict[nameid] += 1
            else:
                shorterdict[nameid] = 1
            continue 

        # if duplicate id then do consensus
        seqcount[nameid] += 1
        for i in range(len(sequence)):
            if sequence[i] == seqdict[nameid][1][i]:
                continue
            if qsequence[i] > seqdict[nameid][2][i]:
                seqdict[nameid][1] = seqdict[nameid][1][0:i] + sequence[i] + seqdict[nameid][1][i+1:]
                seqdict[nameid][2] = seqdict[nameid][2][0:i] + qsequence[i] + seqdict[nameid][2][i+1:]

            if nameid in diff:
                diff[nameid] += 1
            else:
                diff[nameid] = 1
            continue 
        continue
        
    c = c + 1

print "Length consensus ", len(seqdict)
print "Number grossly different ", len(differentdict)
print "Number of shorter ", len(shorterdict)
print "Number of diffs ", len(diff)

raw_input("How big am I")


