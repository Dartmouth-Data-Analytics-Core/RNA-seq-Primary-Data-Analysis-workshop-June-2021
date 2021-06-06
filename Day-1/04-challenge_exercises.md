## File manipulation

Part of getting comfortable with bash is flexing your muscles and practicing using the skills you learned today. Here are a couple of exercises to encourage you to work with the bash tools that you learned to query the fastq files in `/scratch/rnaseq1/data/raw-fastq/subset`.

1. Since this data was downloaded from the SRA write a code to determine if the adapter sequence has already been trimmed from these reads. You've forgotten if the kit was Nextera or Truseq, so to be on the safe side you should check for both.<br> <br>*Hint: you will need to count the sequences that contain this sequence, remember to use regular expressions to differentiate matches at the beginning of the read*.<br>
<br>**Nextera adapter sequence:** CTGTCTCTTATA<br>
**Truseq adapter sequence:** AGATCGGAAGAG

2. How many of the reads in each file start with a string of A's that could indicate a poly-A tail that needs trimming? 
