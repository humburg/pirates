# extract the sequence data by stripping the label from the start and end
# reads from either standard input or file name given on command line
# assumes every line is sequence data

import fileinput

for line in fileinput.input():
    line = line.rstrip("\n")
    # nameid = line[0:8] + line[-8:]
    # print nameid, line[8:12], line[-12:-8]
    # print line
    print line[12:-12]
