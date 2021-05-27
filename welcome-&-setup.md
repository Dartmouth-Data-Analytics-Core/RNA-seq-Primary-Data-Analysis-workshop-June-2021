# Welcome to the DAC RNAseq workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup so you are able to participate in all sections of the workshop.  

For those of you that indicated you did not have an account on *discovery* you should have received an email from me explaining how to set that up, please make sure this is done and you are able to log into your account **BEFORE** the workshop begins. YOU WILL NEED A DISCOVERY ACCOUNT!

## Downloading the data ##

For this workshop we will be using a dataset downloaded from the short read archive (SRA), a public repository of genomic data. This dataset comes from [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625), and was collected from human airway smooth muscle cells to test gene pathways effected by exposure to Glucocorticoid drugs, which have been historically used for their anti-inflammatory effects to treat asthma. Four cell lines were treated with either a control vehicle (untreated), dexamethasone (dex), albuterol (alb), or both dexamethasone and albuterol (co-treated) for 18 hours before transcriptomes were extracted.

The commands that you will be following can be found in markdown `(.md)` files where there is a brief description of the command and how it is applied to the data and what it does followed by an example command that you can copy and paste into the terminal window. The majority of analysis will be performed using the terminal window on your local machine, with an open `ssh` connection to discovery7.

In your terminal window navigate to where you want to download the files needed for this workshop onto your local machine. Then execute the following command:

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/RNA-seq_workshop_July2020/
```

## Setting up a Conda Environment ##

Conda is an open source package and environment manager that runs on Windows, MacOS and Linux. Conda allows you to install and update software packages as well as organize them efficiently into environments that you can switch between to manage software collections and versions.

We will be using Conda to make sure everyone has the required software to perform the analyses included in the workshop. To start using conda on Discovery, open an ssh connection using your username & password, then run the following command:


```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add the above line of code to your `.bashrc` file in your home directory, otherwise you will need to run this command each time you start a new session on discovery. You can follow the [tutorial](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=72888) on the *Research Computing* teams website.

Now run the following command to create a .conda/ directory in your home drive to store all of your personal conda environments. You only have to run this command once to make this directory, so it does not need to be added to your .bashrc file.

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

Now create the conda environment that we will be using for the workshop. This takes about 15 minutes to complete. As you will see, many packages are being installed or updated, all managed for you by conda.

```bash
conda env create -f /scratch/rnaseq1/environment.yml
```

When you are ready activate the conda environment, which you will need for the work we are doing for day 1 of the workshop you can use the following command.

```bash
conda activate rnaseq1
```
You will see that the activate command has worked when it reads (rnaseq1) rather than (base) to the left of the prompt.

There is one more conda environment you will need to create in order to run some of the QC metrics for our alignment. We will install the picard program (yes from Star Trek) using the bioconda channel.

When you are finished using a conda environment, it is good practice to deactivate your session with the following command.

```bash
conda deactivate
```

Thats it! This conda environment contains all the software you will need during the workshop. If you run into issues with the setup, please reach out to us at *DataAnalyticsCore@groups.dartmouth.edu* and someone will be in touch to assist you. 
