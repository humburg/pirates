# extract labels and associated quality readings for each entry

import gzip
f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")
c = 0
for line in f:
    if (c % 4) == 1:
        line = line.rstrip("\n")
        nameid = line[0:8] + line[-8:]
        clabelstart = line[8:12]
        clabelend = line[-12:-8]
    if (c % 4) == 3:
        line = line.rstrip("\n")
        qualityid = line[0:8] + line[-8:]
        print nameid, clabelstart, clabelend, qualityid, line[8:12], line[-12:-8]
    c = c + 1

