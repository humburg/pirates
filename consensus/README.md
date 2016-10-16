Parses fastq file in one pass and builds consensus about sequences. On the data tested we filtered around 16.5 million sequence readings into around 800k unique sequences.

Input format is usual fastq format:

```
@name string
[8 char label part 1][4 char constant][variable length sequence][4 char constant][8 char label part 2]
+
[quality readings matching sequence data above]
```

