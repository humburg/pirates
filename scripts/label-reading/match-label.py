import gzip
import sys

usage = "Must have 16 character label to match to on command line"

# get label we will match to from command line
if len(sys.argv) != 2:
    print usage
    exit(0)

match = sys.argv[1]
if len(match) != 16:
    print usage
    exit(0)

f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")
c = 0
for line in f:
    if (c % 4) == 1:
        line = line.rstrip("\n")
        nameid = line[0:8] + line[-8:]
        if nameid == match:
            print line
    c = c + 1

