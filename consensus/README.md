Parses fastq file in one pass and builds consensus about sequences. On the data tested we filtered around 16.5 million sequence readings into around 800k unique sequences.

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
* name - line starting with '@' character, sub fields
  * sequence count - number of sequences that were processed in making this consensus sequence
  * diff - one or more difference counts indicating where nucleotide has been substituted in building the consensus, first value is character position in the sequence and then the count of 'A', 'C', 'G' and 'T' observed at that position
