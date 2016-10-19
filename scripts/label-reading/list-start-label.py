import gzip
f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")
c = 0
for line in f:
    if (c % 4) == 1:
        line = line.rstrip("\n")
        print line[0:8], line[8:12]
    c = c + 1

