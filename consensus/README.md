Parses fastq file in one pass and builds consensus about sequences. 

## What it Does

Takes a fastq file of sequence readings and sorts them according to the observed 16 character label. For each sorted group the sequence is examined and at each nucleotide position the reading with the highest quality is kept, along with that quality, to build what we have called a consensus sequence. A list of differences is recorded and printed out with the sequence data. On the data tested we filtered around 16.5 million sequence readings into around 800k unique sequences.

## File Formats

Reads from fastq file and outputs sequences to stdout. Some stats reported to stderr.

Input format is usual fastq format:

```
@name string
[8 char label part 1][4 char constant][variable length sequence][4 char constant][8 char label part 2]
+
[quality readings matching sequence data above]
```

Output format is

```
@[int sequence count] [int char position]A[int count]T[int count']C[int count]G[int count] ...
[16 char label][variable length sequence]
+
[quailty readings match label and sequence above]
```

where:
* 1st line: name - line starting with '@' character, sub fields
  * sequence count - number of sequences that were processed in making this consensus sequence
  * diff - one or more difference counts indicating where nucleotide has been substituted in building the consensus, first value is character position in the sequence and then the count of 'A', 'C', 'G' and 'T' observed at that position
* 2nd line: label id and sequence
* 3rd line: '+' character
* 4th line: quality data associated with label and sequence reading

### Future Work
Current script only records difference where substituation have been made into the consensus sequence which happen when a observed nucleotide has a higher quality reading then the existing value at that position in the consensus sequence. May be better to reord all difference observed. 
