import gzip
f = gzip.open("data/fastq/Jurkat_only.assembled.fastq.gz", "r")
c = 0
for line in f:
    if (c % 4) == 1:
        line = line.rstrip("\n")
        print line
        nameid = line[0:8] + line[-8:]
        print nameid, line[8:12], line[-12:-8]
    c = c + 1

