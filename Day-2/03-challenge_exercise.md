# Prokaryotic Alignment

Earlier today we used STAR to align human genomes to the hg38 reference genome using a splice aware aligner. In some cases a splice aware aligner is not approprate for the dataset. One example of this are prokaryotic datasets. In these cases you will want to use an aligner like *Bowtie* or *BWA*. 

To practice with non-splice aware aligners we will use transcriptome data from *Staphylococcus aureus* MRSA 1369 exposed to human urine. 

**Experimental Design**
Total RNA was extracted from bacteria grown in TS broth using RibopureTM bacteria RNA extraction kit (Invitrogen). 
Library construction included DNase treatment (TURBO DNase, ThermoFisher Scientific) and rDNA depletion (QIAseq FastSelect, Qiagen) followed by RNA fragmentation and random priming. 
cDNA synthesis (NEBNext Ultra II, New England Biolabs) was followed by end repair, 5 phosphorylation and dA-tailing before. 
Libraries were sequenced on a partial lane of Illumina HiSeq 400 with 150bp PE sequencing.

The SRA accesssion numbers for this data are:

-SRR14057225
-SRR14061515
-SRR14069283
-SRR14069797
-SRR14072892
-SRR14078369

Though the data is paired end the pairs of reads have been joined and these files should be treated as unpaired reads. 

## Bowtie

Before you use any tool it is a good idea to look at the [manual](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner) to get a feel for the options available in the tool that you're using and figure out which of the options are the most appropriate for your dataset. 
There are two main ways that bowtie can align a transcriptome to a reference genome, an end-to-end alignment (which is the default) or a local alignment. The local alignment enables soft clipping of reads to improve the alignment score, while the end-to-end alignment uses all of the bases in the read but this may result in a poorer alignment score. The type of alignment you chose depends on the experiment that you run. For example if the reference that you were aligning to is very distantly realted than a local alignment that accomodates genomic differences expected between distantly related genomes. End-to-end alignment will be more appropriate for aligning reads to reference genomes that are closely related. 

### Bowtie Index 

A reference sequence was downloaded for *S. aureus* from the NCBI refseq data base using the ftp site [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/GCF_000013425.1_ASM1342v1/). If the strain of interest is not available in the refseq database another good resource for prokaryotic genomes is the IMG database hosted by the [Joint Genome Institute](https://img.jgi.doe.gov). 

The reference genome must be indexed in the same way that we indexed the reference before running *STAR*. This is done with the *bowtie2-build* command. This command takes two arguments the reference genome file to be indexed and the prefix to be used for the indexed file. The output of this command is a set of indexed genome files with the suffix `.bt2`.

``` bash

bowtie2-build <reference_genome.fasta> <bt2-base>

```

The reference genome in this case has already been indexed for you and can be accessed at `/scratch/data/prok_data/indexed_ref/S_aureus`.

The basic flags to be used with `bowtie2` command for aligning reads to a reference genome.

**-x** prefix for indexed reference genome
**-1** a fastq file with forward (R1) reads
**-2** a fastq file with reverse (R2) reads
**-U** a fastq file with unpaired reads
**-S** SAM file to write aligned reads to

Optional flags that may enhance your analysis:

**--interleaved** a fastq file with interleaved paired reads (R1, R2, R1, R2, etc.)
**--sra-acc** SRA accession number, bowtie can access the SRA directly using only the accession number
**-b** unaligned reads to be aligned from a BAM file (rather than a fastq file)
**-N** number of mismatches allowed per alignment
**-L** length of seed substring, smaller values make the alignment more sensitive but are slower
**--end-to-end** *default setting* use end to end alignment method
**--local** use local alignment
**--no-unal** don't print unaligned reads to SAM output file
**-p** number of threads to use during the alignment 

These are just a few of the enhanced options available there are many more available that you can (and should) read about [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner).

### Challenge exercises

1. Write your own code to align the reads in SRR14057225 to the reference genome (remember that the paired ends here have been merged and these files should be treated as unpaired reads).
2. Write a loop that aligns each of the samples to the reference genome and writes each sample to its own SAM output file. 
3. Write a script that runs bowtie alignments on samples that are handed to the script.
