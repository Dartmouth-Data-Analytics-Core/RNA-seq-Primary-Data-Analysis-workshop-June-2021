# Joining the workflow together  #
Just like experiments in the lab, it is important that you keep careful track of how your work was performed. To do this for computational analysis, we need to keep track of the commands used to perform an analysis.

One way to do this is by storing all your commands for a particular analysis in one file, and writing the commands in such a way that they are completed in succession, with each command using the output from the previous command.

One of the benefits of writing your code this way is the same code can be used over and over by simply providing the code with a different set of sample names. Code written this way is generally written, edited, and revised using a text editor (e.g. [Sublime text](https://www.sublimetext.com/) or [BBedit](https://www.barebones.com/products/bbedit/)).

![](../figures/terminal_commands3.png)

## If/Else statements

A useful way to control operations in a script is with if/else statements.  These allow commands to be run only if certain conditions are true.
if statement example:
```bash
#Just an if statement, without an else statement
a=5
if [ $a == 5 ]
    then
    echo "a is equal to 5"
fi
#Try setting a to 6, and observe there is no output.
```

if/else statement example
```bash
a=5
if [ $a == 5 ]
    then
    echo "a is equal to 5"
    else
    echo "a is not equal to 5"
fi
```

There are many possible conditions to check for variables and files in Bash.  [Here is a link](https://tldp.org/LDP/Bash-Beginners-Guide/html/sect_07_01.html) to a reference containing the syntax for some available conditional tests.


## Exercise: complete the workflow

After you have completed the previous lessons & exercises (or as 'homework' after the workshop), we have provided the foundations of a simple script that links together commands from each part of the RNA-seq analysis we performed. Try filling in the rest of the commands to complete the workflow, and run it on the full dataset.

If you run into trouble leave a slack comment explaining your error(s) and we will do our best to get back to you!

To start out we are going to create an array with our sample names so that we can use the sample names to control the input and output of each of the commands.
Then we will use the array to write several for loops that iterate over the elements of the array to do something with them.

```bash
#!/bin/bash

echo -n "RNA-Seq Pipeline beginning at: "; date

###################################
### Data Gathering ###

echo "Symlinking Raw Data"
mkdir data
cd data
ln -s /dartfs-hpc/scratch/rnaseq1/data/raw-fastq/subset/*.gz ./
echo -n "Symlinks created in: "; pwd
cd ..

sample_list="SRR1039508 SRR1039509 SRR1039512 SRR1039513"
echo "The pipeline will be run for the following samples:"
for i in $sample_list; do echo $i; done

echo "Checking file existence..."
for i in $sample_list
do
if [ -f data/${i}_1.chr20.fastq.gz ]
then
    echo $i "exists."
else
    echo $i "does not exist!!!"
    echo "Exiting pipeline."
    exit
fi
done
echo "File checking complete."



###################################
### FastQC ###

echo "Running FastQC..."
for i in $sample_list; do fastqc data/${i}_1.chr20.fastq.gz; fastqc data/${i}_2.chr20.fastq.gz;done
echo "FastQC complete."

echo "Moving FastQC reports..."
mkdir fastqc_reports
mv data/*fastqc.html fastqc_reports/
mv data/*fastqc.zip fastqc_reports/
echo "Moving reports complete."

###################################
### Cutadapt ###

echo "Running Cutadapt..."
mkdir trim
cd trim
for i in $sample_list; do cutadapt -o ${i}_1.trim.chr20.fastq.gz -p ${i}_2.trim.chr20.fastq.gz ../data/${i}_1.chr20.fastq.gz ../data/${i}_2.chr20.fastq.gz -m1 -q 20 -j4 > ${i}_cutadapt.report; done
echo "Cutadapt complete."
cd ..


###################################
### STAR Alignment ###

mkdir alignment
cd alignment

echo "Running STAR alignment..."
for i in $sample_list; do STAR --genomeDir /dartfs-hpc/scratch/rnaseq1/refs/hg38_chr20_index --readFilesIn ../trim/${i}_1.trim.chr20.fastq.gz ../trim/${i}_2.trim.chr20.fastq.gz --readFilesCommand zcat --runThreadN 4 --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFileNamePrefix ${i}_; done
echo "STAR complete."
cd ..


###################################
### Run MarkDuplicates ###

###################################
### Run CollectRNASeqMetrics ###

###################################
### Move alignment and metrics into a single directory and run multiqc ###

###################################
### Run htseq-count ###

echo -n "RNA-Seq Pipeline finished at: "; date

```
