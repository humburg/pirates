Python utilities to print out label data from start and end of sequence. These assume that you have a local directory named `data` and will try and read the file `data/fastq/Jurkat_only.assembled.fastq.gz`

`list-labels.py` - lists out the labels for each sequence in the format variable-label, start-constant-label, end-constant-label. The code just prints out characters from the sequence string, no effort is made to see if there might be characters missing.

Sample use

To list all the labels into the one 425MB file:

```
python list-labels.py > all-labels.txt
```

To get a sorted listogram of labels
```
cut -f1 -d' ' all-labels.txt | sort | uniq -c > histogram-labels.txt
sort -k1,1rn histogram-labels.txt > sorted-histogram-labels.txt
```
To get a histogram of that histogram
```
awk '{print $1}' histogram-labels.txt | sort -n | uniq -c > count-of-counts.txt
```


