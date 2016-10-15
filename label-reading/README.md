Python utilities to print out label data from start and end of sequence. These assume that you have a local directory named `data` and will try and read the file `data/fastq/Jurkat_only.assembled.fastq.gz`

`list-labels.py` - lists out the labels for each sequence in the format variable-label, start-constant-label, end-constant-label. The code just prints out characters from the sequence string, no effort is made to see if there might be characters missing.

