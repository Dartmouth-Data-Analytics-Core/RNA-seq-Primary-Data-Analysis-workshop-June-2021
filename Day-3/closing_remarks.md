# Closing remarks

### Workshop goals:
- Develop a working understanding of the analytical workflow for a modern RNA-seq experiments
- Build a working knowledge of sample preparation considerations for RNA-seq experiments
- Learn how to process raw NGS data in FASTQ format to generate a gene expression matrix
- Learn how to perform a detailed quality control analysis

### Workshop overview
<p align="left">
<img src="../figures/analysis_overview.png" alt="lib-composition"
	title="" width="40%" height="40%" />
</p>

### Some final take-aways from the workshop:
- Spend the time to plan, consult, practice, (and money) to generate high quality data sets
- If you are going to do a lot of Bioinformatics, you should get **really** good at the command-line (Bash), otherwise, pre-processing will be slow & painful
- Identify, understand, and check key QC metrics to ensure the quality of your results
- Understand when and how to apply different normalization approaches.


### How to consolidate your learning:
- Re-run the code a week or two after the workshop, as this is a great way to consolidate what you have learned at the command-line
- Edit the code, run sub-sections, read the `man` pages for commands, etc. to build a solid understanding of how everything works
- Practice with the complete dataset (all chromosomes), that is available to you for approx. 1 month on discovery in `/scratch/rnaseq1/data/`. This will give you experience running data from the entire genome, and an appreciation for the computational resources and required time to complete these tasks.
- Read the methods sections of published papers that perform RNA-seq, to gain an appreciation for the range of approaches used in practice and how they are implemented
- Ask us questions! (Bioinformatics office hours: https://dartmouth.zoom.us/s/96998379866, fridays at 1-2 pm, password: bioinfo)

## Suggested reading:

Reading manuscripts that use RNA-seq, or reviews specifically focused on RNA-seq are excellent ways to further consolidate your learning. Additionally, reading the original manuscripts behind common tools will improve your understanding of how that tool works, and allow you to leverage more complicated options and implementations of that tool when required.

Below we provide some suggested reading to help get you on your way:

#### Review articles
- (Stark *et al*, 2019, *Nat. Rev. Genetics*). `RNA Sequencing: The Teenage Years`[https://pubmed.ncbi.nlm.nih.gov/31341269/] from Stark *et al*, 2019, *Nat. Rev. Genetics*, `RNA Sequencing: The Teenage Years`.


### What next?

While multiple downstream applications of RNA-seq exist, differential expression analysis is the most common. As discussed in the last lesson, standard normalization does not provide a statistical assessment of the data which is needed to make inferences (e.g. gene X is differentially expressed between treatment groups).

If you are interested in learning how to perform a DE analysis, sign up for Part 2 in our series of RNA-seq data analysis workshops. You should all have received an email containing the [sign up link](https://sites.dartmouth.edu/cqb/workshops/rnaseq-differential-gene-expression-analysis-workshop/).

#### Differential expression analysis  workshop - July 2020
##### Workshop goals:

- Develop a working understanding of fundamental bioinformatics and statistical concepts for a typical bulk RNA-seq DE analysis
- Learn how to leverage the R/Bioconductor framework to perform DE analysis  
- Learn how to use unsupervised data analysis methods (e.g. principal components analysis) to explore RNA-seq datasets  
- Perform a complete DE analysis on a real RNA-seq dataset  

Monday, July 12, 2021
Wednesday, July 14, 2021
Friday, July 16, 2021

Registration Close Date: July 5, 2021
Registration Limit: 40

### Feedback:

We ask that you all complete the survey that will be sent out over email so that we can gauge what worked well and what we need to improve for our next workshop. If you have additional thoughts that were not addressed in the survey, please feel free to contact any one of us, or reach out to the DAC email directly (*DataAnalyticsCore@groups.dartmouth.edu*).

<img src="../figures/logo.jpg" width="250" height="140" >

We plan to offer this workshop again, as well as workshops covering other types of genomic data analysis and bioinformatics. If you have suggestions for workshops you would like to see, please let us know!

Please feel free to reach out to us with questions about concepts discussed in the workshop, or for a analysis consultations. Our **bioinformatics office hours** on **Fridays 1-2pm** are a great place to do this! (currently on zoom: https://dartmouth.zoom.us/s/96998379866, pword: *bioinfo*)

### Now.. Final questions?
