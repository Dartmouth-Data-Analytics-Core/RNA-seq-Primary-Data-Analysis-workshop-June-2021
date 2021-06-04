

# change script to use complete raw counts matrix




# Part 4 - Data normalization in RNA-seq

### Learning objectives:
- Understand why read counts must be normalized in RNA-seq data
- Learn the principles behind the major normalization strategies in RNA-seq and when to apply them
- Learn how to perform these normalization strategies

For this lesson we will work in the `quant/` directory we made previously:
```bash
# go to our home worksp dir
rnaw

# move into folder for quantification analysis
cd results/quant
```

If you get lost and need to catch up quickly you can copy the files needed for the next step with this command:

```bash
cp /scratch/rnaseq1/data/raw-fastq/htseq-count/* results/quant
```

## Count normalization in RNA-seq

In order to compare expression levels between genes within a sample, or genes across multiple samples, it is critical the data is normalized to allow appropriate interpretation of the results. Which normalization strategy depends on several factors such as the library type, as well as the comparison you wish to make (e.g. within- vs between-sample).

Below we will discuss the major sources of variation that need to be accounted for during normalization in order to reach appropriate conclusions about your results. Subsequently, we will discuss the normalization approaches that account for these sources of variation and when to use them.

### Gene length

Genes length typically varies a great deal across the genome of virtually all organisms.
In the below example, we see two genes, X and Y. If we were to simply compare the number of reads successfully mapped to gene X and gene Y, we would conclude gene X is expressed ~2x that of gene Y.

However, since gene X is ~2x longer than gene Y, gene X contributed ~2x as many RNA fragments to the original sequencing library. Gene X therefore only has more reads mapped to it because it is longer, **NOT** because it is truly expressed at greater level than gene Y.  

<p align="center">
<img src="../figures/gene_length.png" alt="glength"
	title="" width="80%" height="80%" />
</p>

To address gene length bias, we must normalize raw read counts in a way that accounts for the size of each gene. If we did this for gene X & gene Y, we would conclude their gene expression levels are similar.

Normalization for gene length is critical when comparing between genes **within the same sample**, however when comparing expression of the same gene **across different samples**, correction for gene length is not as important since we assume the gene is of the same length in all samples.

### Library size/sequencing depth  

Although samples are pooled together (each sample is tagged with a barcode and samples are combined) at similar concentrations in a sequencing run, some samples will end up being sequenced more than others, leading to slight differences in how many reads are produced for that sample, and therefore sequencing depth and size. Furthermore, if samples are sequenced on separate runs, their sequencing depths may be very different. If we don't account for this variation in sequencing depth, we might conclude some genes are expressed at greater levels in a sample that has simply been sequenced to a higher depth.  

<p align="center">
<img src="../figures/library_composition.png" alt="lib-composition"
	title="" width="85%" height="85%" />
</p>







### Library composition

The presence of truly differentially expressed genes (in particular, DEGs with very large fold changes) between samples will cause the number of reads for other genes in those samples to be skewed. For example, in the below example, gene C is differentially expressed between the two samples, with much higher expression in sample 1. This high number of reads causes fewer reads to be detected for other genes in this sample, making it appear that these other genes are expressed at lower levels than in sample 2, however this is simply an artifact of library composition differences between the samples.

<p align="center">
<img src="../figures/library_size.png" alt="lib-size"
	title="" width="85%" height="85%" />
</p>




Note on why comparisons between group and within-/between-samples

## Normalization methods

Several normalization methods exist for RNA-seq data. Which method you use depends on the comparison you are trying to make (e.g. between or within samples), therefore it is important to understand how each is calculated and when to use it.

### Counts per million (CPM)

CPM is a simple normalization method that involves scaling the number of reads mapped to a feature by the total number of reads in a sample. This fraction is multiplied by 1 million in order to provide the number of reads per million mapped in the sample.

<p align="center">
<img src="../figures/cpm.png" alt="lib-composition"
	title="" width="65%" height="65%" />
</p>

We will briefly use R to calculate CPM values for our dataset. If you are not familiar with R don't worry, this is not complex R code and many software packages will calculate normalized counts for you.
```r
# read in raw counts matrix
all_counts <- read.table("all_counts_full.txt", sep="\t", stringsAsFactors=F, header=T)

# look at the counts object
head(all_counts)

# write a function that will calculate TPM
cpm <- function(counts) {
	cpm <- c()
	for(i in 1:length(counts)){
		cpm[i] <- counts[i] / sum(counts) * 1e6
	}
	cpm
}

# apply function to the columns of raw counts data
all_counts_cpm <- apply(all_counts[2:ncol(all_counts)], 2, cpm)

# write to file
write.csv(all_counts_cpm, file="all_counts_CPM.csv")
```

**NOTE:** CPM does **NOT** normalize for gene length, therefore cannot be used to compare expression between different genes in the same sample.

### Transcripts per million

<p align="center">
<img src="../figures/tpm.png" alt="lib-composition"
	title="" width="50%" height="50%" />
</p>

Calculate TPM from our raw read counts:
```r
# read in raw counts matrix
all_counts <- read.table("all_counts.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
gene_lengths <- read.csv("lengths-grch38.csv", header=FALSE)

# look at the lengths object
head(lengths)

# write a function that will calculate TPM
cpm <- function(counts, lengths) {
	rate <- counts / lengths
	tpm <- c()
	for(i in 1:length(counts)){
		tpm[i] <- rate[i] / sum(rate) * 1e6
	}
	tpm
}

# apply function to the columns of raw counts data
all_counts_tpm <- apply(all_counts[2:ncol(all_counts)], 2, tpm, gene_lengths)

# write to file
write.csv(all_counts_tpm, file="all_counts_TPM.csv")
```

### Plot expression of example gene



### Reads/fragments per kilobase of exon per million

<p align="center">
<img src="../figures/rpkm-fpkm.png" alt="lib-composition"
	title="" width="90%" height="85%" />
</p>


Since our dataset is paired-end and we counted the number of fragments in the quantification step, we are calculating FPKM. Calculate FPKM from our raw read counts:
```r
# write a function that will calculate TPM
cpm <- function(counts, lengths) {
	rate <- counts / lengths
	fpkm <- c()
	for(i in 1:length(counts)){
		fpkm[i] <- rate[i] / sum(counts) * 1e9
	}
	fpkm
}

# apply function to the columns of raw counts data
all_counts_fpkm <- apply(all_counts[2:ncol(all_counts)], 2, fpkm, gene_lengths)

# write to file
write.csv(all_counts_fpkm, file="all_counts_FPKM.csv")
```




**Normalization method comparison**

**Method** | **Name** | **Accounts for** | **Appropriate comparisons**
-------|-------|-------|-------
CPM | Counts per million | Depth	 | - Between-sample<br>- Within experimental group
TPM | Transcripts per million | Depth & feature length | - Between- and within-sample<br>- Within experimental group
RPKM/FPKM | Reads/fragments per kilobase<br>of exon per million | Depth & feature length | - Within-sample<br>




## add in example of why RPKM is difficult to compare between samples



## Normalization methods for differential expression

size factor based normalization
