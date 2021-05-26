# Unix/Linux Shell basics

The Unix/Linux *'Shell'* describes a program that takes commands from a some input (essentially your keyboard) and passes them to an operating system that will execute them. In contrast to a *Graphical User Interface (GUI)* the Shell is both simulatenously a *command line interface (CLI)* and a programming language that allows you to perform tasks on your system.

<p align="center">
  <img src="../figures/terminal.png" height="300" width="500"/>
</p>

Interacting with a system through the Shell has many advantages over a GUI. The Shell allows you to quickly and easily navigate through directories on your computer, make, copy and search files in a systematic way, and construct pipelines that will execute complex tasks on big datasets.

Importantly, the Shell allows us to do each of these in the context of Bioinformatics, and Bioinformatics softwares.

## Why learn Shell?  
Shell can be challenging to learn, however is an absolutely key skill in bioinformatics, as it is used to primary way in which we interface with a lot of bioinformatics software and file types.

Some bioinformatics softwares provide GUIs that enable execute tasks with programs that you would otherwise execute using the Shell. Whuile such softwares can be powerful in the right context, they can also make it very easy to perform tasks in bioinformatics incorrectly, therefore they should be treated with caution.


## The Bash shell

### The absolute basics

There are different different types of Unix shells, however the most popular is Bash (the *Bourne Again Shell*) and is also now the most common. Since the majority of participants will be using the Bash shell, and this is the default shell used on Dartmouth's high performance computing system (which we will be using), this lesson will be introduce the Shell through using the Bash shell, however most, if not all, content should be transferable to other Unix shells.

> Use the Cheat Sheet in the GitHub repo to help you learn commands and available options.

Accessing the (bash) shell:  
- On a mac or linux system, the *Terminal* application provides access to the shell. There are also applications that you can download that provide customizations not present in the Terminal application, such as [iTerm2](https://iterm2.com/).
- On a Windows system, you can use an applicatrion such as [MobaXterm](https://mobaxterm.mobatek.net/).

<p align="center">
  <img src="../figures/shell.png" height="80%" width="80%"/>
</p>

When you open your terminal application you will be presented with the command prompt `$` when you are able to input commands. If the terminal is busy and cannot currently accept new commands, you will not be presented with the prompt.

When the prompt is shown, you can enter commands by typing them in after the prompt. Commands are typically composed of three components:  
- the command itself  
- any flags or options you wish to run the command with (not always required)
- and an argument

In the above example, we are asking the Shell to pass the `mkdir` command to the operating system (for making directories) with the `-p` option (which just lets us make parent and sub directroies at the same time) and the argument detailing what directories we want the command to make.

Manual pages for specific commands can be accessed using the `man` command.
```bash
man mkdir
```

The shell provides us with a number of commands that allow us to list files in our current working directory, as well as change the current working directory to another location. For example:
```bash
# 'ls' command lists files in our current wopkring directory
ls

# run ls with the '-a' option to include hidden files
ls -a

# pwd show you your current working directory
pwd

# cd allows you to change your current working directory ('.' means current directory)
cd ./

# '..' tells the shell to move your current directory up one directory
cd ..

# check you directory again
pwd

# now go back down the tree
cd OwenW/
pwd
```

To go back down the directory structure, we specified a directory that was in our current working directory (cd). This is called a **relative path**, since it is relative to our cd and will only work provided our cd is relative to the directory we are trying to reach in the way written in the command.  

Relative paths are contrasted to **absolute paths** which always starts with a '/' and will start at the root (highest level) of the directory tree, and work from wherever you are in the directory substructure. For example:
```bash
ls /Users/OwenW/
```

By default, your terminal application will start your current directory as your *home directory* (more on that later). No matter where you are, you can always get back to your home directory using the tilde `~` with the `cd` command.
```bash
cd ~/
```

Another useful command is `echo` which will evaluate and print characters provided to it.
```bash
echo 'bla bla bla'
```

We can use the redirect command (>) to redirect the output of commands like echo into a file. As an example, lets save the important note we made above to a text file.
```bash
echo 'bla bla bla' > mynotes.txt
```
## Log on to discovery cluster

Many of the higher level commands for working with NGS data will require a lot of memory and computing power, more than most laptops can handle efficiently.
The discovery cluster is a resource hosted by Dartmouth's Research Computing team. This cluster enables you to execute high level commands wihtout using the memory and computing power on your local machine (more about this tomorrow). Let's log onto the discovery cluster now. We will use a secure shell command `ssh` to log onto the discovery cluster.

```bash

# Establish the secure shell connection
ssh netID@discovery7.dartmouth.edu

# Enter your password at the prompt (when you type no characters will show up to preserve privacy)
netID@discovery7.dartmouth.edu's password:

# You're in!
(base) [netID@discovery7 ~]$

```
All of the commands that you just executed locally in your terminal window work the same way when you are logged into discovery. It is always useful to orient yourself when you're working on an HPC so that you know where the output of all of the commands you run will end up. Lets run our first command to get your location.

```bash

# Check your location on the cluster
pwd

```

You should see something like `/dartfs-hpc/rc/home/h/netID` displayed in response to your command. Initially when you log on you will always be directed to your home directory (the address or path listed above). Your home directory by default will have 50GB of storage space to begin with, if you are running something that requires more storage space it is possible to extend that limit temporarily with the `/scratch/ drive`. This is where we have stored all of the files you will be working with today. Directories and files hosted on the `/scratch/` drive will only be maintained for 45 days, you will receive a notification from research computing before the data is deleted, but it will be maintained longer than 45 days.

It is a good idea when working on projects on an HPC to stay organized so lets start by making a folder, or directory to store all of the work you do today we will call it `fundamentals_of_bioinformatics`. You will notice that I chose a title that has no spaces in it, this is because the space is a special character, special characters need to be *escaped* with the `\` and so `funadmentals_of_bioinformatics` would look like `fundamentals\ of\ bioinformatics` with the escape characters. You can see that file names with spaces become unweildy to type out so most programmers will replace spaces with `_`, `.`, or `-` in their filenames to keep everything neat.

```bash
# navigate to scratch so you can make your own directory there 
cd /scratch/

# make the directory 
mkdir -p omw/fundamentals_of_bioinformatics

# go into it
cd omw/fundamentals_of_bioinformatics

# set an alias so we can get here quicly 
alias biow='cd /scratch/omw/fundamentals_of_bioinformatics'
# NOTE: you can add this line to your .bashrc so it get run everytime you log in, we will cover this below 

# check your location on the cluster
pwd

# list the contents of your directory
ls

```
As expected the new directory that you created is empty there are no files. Lets copy a file from the `/scratch/` directory we created for this workshop to the directory you just created. This file (`all_counts.txt`) provides raw read counts for an RNA-seq experiment, with genes in rows and samples in columns.

```bash

# copy the file from the scratch drive to the fundamentals_of_bioinformatics directory you just created
cp /scratch/fund_of_bioinfo/counts/all_counts.txt ./

```


### Viewing the contents of files

The shell provides us with a number of commands to view the contents of files in define ways. The `cat` command for example (standing for concatenate) will print the entire contents of a file to the terminal. This can be useful for smaller files, but as you will see with larger files can be an unweildy way to look at the contents of a file.
```bash
cat all_counts.txt
```

When working with larger files, which we are usually doing in bioinformatics, you may not wish to print the whole file as it would overrun your terminal. Other commands exist that allow you to explore file contents with more control.
- `more` shows you as much of the file as can be shown in the size of the terminal screen you have open, and you can continue to "scroll" through the rest of the file by using the space bar  
- `head` will print the first 10 lines by default, but this number can be controlled with the `-n` option
- `tail` will print the final lines of a file, and can also be controlled with the `-n` option

We will use a larger text file to show the utility of these commands, as well as other commands in the subsequent parts of this lesson.
```bash
# Show the first 20 lines of the all_counts.txt file
head -n 20 all_counts.txt

# Show the last 50 lines of the all_counts.txt file
tail -n 50 all_counts.txt

# use word count (wc) command with the lines option (-l) to show how many lines (rows) are in the dataset
wc -l all_counts.txt
```

### Renaming and removing files

Sometimes you will need to reorganize your directories or rename a file, which can be achieved wuith the `mv` command. Let's start by copying the all_counts.txt file from the fundamentals_of_bioinformatics directory to your home directory.

```bash
# Copy the all_counts.txt file to your home directory
cp all_counts.txt ~/all_counts.txt
```
Now let's rename the copy of the all_counts.txt file that we just created.
```bash
# Rename the copied all_counts.txt file
mv ~/all_counts.txt ~/all_counts.copy.txt
```
You can also use the `mv` command to move a file to a new location. Let's move the all_counts.copy.txt from your home directory into your fundamentals_of_bioinformatics directory.
```bash
# Move the all_counts.copy.txt into your fundamentals_of_bioinformatics directory
mv ~/all_counts.copy.txt fundamentals_of_bioinformatics

#check the contents of your fundamentals_of_bioinformatics directory
ls
```

Copying the all_counts.copy.txt file was just an exercise to show you how the tools works, in practice you will want to keep your directories as neat as possible as you accumulate a lot of files. Let's remove the all_counts.copy.txt file with the `rm` command.

```bash
# Remove the all_counts.copy.txt file
rm all_counts.copy.txt
```

You will notice that before the file was deleted you were asked if you were sure you wanted this file deleted. You want to be careful not to remove files that you did not create if you are working in shared directories. If you want to bypass this checkpoint, you can use the `-f` flag with `rm -f` to force the removal of a file, but be careful with this, as there is no *Trash* equivalent in the shell.

### Manipulating file contents

Some commands enable you to maniuplate and subset files based on specific paramters. One useful example is the `cut` command, which allows you to 'cut' a file based on the options you select, such as the `f` option, which corresponds to fields (columns). We could use `cut` to obtain read counts for only the first 5 samples in `all_counts.txt`.
```bash
# Look at only the counts from the first four samples
cut -f 1,2,3,4,5 all_counts.txt
```

To prevent all rows being printed to our console, we could combine the `cut` command with the `head` command using a *'pipe'*, specified by a '|'. Pipes send the output from one command an inital command to a subsequent command, all in the same line, such that you do not need to include an argument for the last command.
```bash
# List only the first 100 lines of only samples SRR1039508 (col 2) and SRR1039523 (col 17)
cut -f 1,2,17 all_counts.txt | head -n 100
```

Similarly to how we used the redirect command (>) above, we could redirect the output of the cut command to create a new counts file, that only contains the columns 1 (gene IDs), and samples in columns 2 and 17.
```bash
# Print the counts from SRR1039508 and SRR1039523 to a new file
cut -f 1,2,17 all_counts.txt > all_counts_sub.txt

# look at head of this new file
head all_counts_sub.txt
```

### Pattern matching with *Grep*

Often we will want to pull a specific piece of information from a large file, lets say that we were interested in the read counts for a specific gene, ALDH3B1 (Ensembl ID: ENSG00000006534). We can use the `grep` command to search for this ID, or any other character string we are interested in, in our counts matrix.
```bash
# Get the count data for ENSG00000006534 (ALDH3B1) from all_counts.txt
grep "ENSG00000006534" all_counts.txt
```

`grep` is a pattern recognition tool that searches in files for a character string we can define. We can define the entire character string, as we did above, or combine regular characters with special characters (or 'wildcards') to search for specific types of matches. the most commonly used special characters are included in the table below.

Operator | Effect
---|---
\* | wildcard stands for any number of anything
^ | start of the line
$ | end of the line
[0-9] or \d| any number (0123456789)
[a-z]| any lowercase letter
[A-Z]| any uppercase letter
\t | a tab
\n | a newline
\s | any white space (tab, newline, space)
\S | non-white space (the opposite of \s)

These regular expressions can be used with any of the tools that you have learned thus far, so if we wanted to list all of the files in our directory that end in .txt we could use the following command.

```bash
# List all files that end in .txt
ls *.txt
```

We can even enhance the power of these regular expressions by specifying how many times we expect to see the regular expression with quantifiers.

Quantifier| Operation
---|---
X* | 0 or more repetitions of X
X+ | 1 or more repetitions of X
X? | 0 or 1 instances of X
X{*m*} | exactly *m* instances of X
X{*m*,} | at least *m* instances of X
X{*m*,*n*} | between *m* and *n* instances of X

Now lets use some of these regular expressions in a `grep` command  to see their utility. Let's use regular expressions to see how many genes have no reads expressed for the first four samples.  

```bash
# Count the number of genes with no reads in the first four samples
grep "^ENSG[0-9]*\t0\t0\t0\t0\t" all_counts.txt| wc -l

# Count the number of genes with no reads expressed in any of the samples
grep "^ENSG[0-9]*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0$" all_counts.txt| wc -l
```

### Shell environment variables

The command line *environment* essentially describes a collection of variables that have been set to provide context for the commands that you run. These variable are referred to as *environment variables*. A number of environment variables are set automatically everytime you log into the bash shell. The `env` command will show all environment variables available in the current shell. Try that now:
```bash
env
```

One important environment variable is `$HOME`, which describes your home directory. Variable such as home can be evaluated by placing the `$` in front of them. For example:
```bash
echo $HOME
```

Environment variables can also be set on the fly and then called as needed. These can be virtually anything. For example, perhaps you want to save the name of the genome version you are working with in your current session, so it can be easily called multiple times in some bash code you are writing.
```bash
# set the variable
genv='hg38.patch13'

# call it with echo and the $
echo $genv
```

You can also use environment variables to store commands that you want to run without having to type the entire command out each time. For example, we might run the `ls` command often with the flags `-lah` to show files in a list format, including all hidden files, and with file sizes in human readable format. The entire command would be `ls -lah`, however if we save the command to a variable, and then call the varaible directly, the command will be evaluated by the shell.

```bash
# save to variable
ll="ls -lah"

# call variable to execute command
$ll
```

It is possible to make variables you add to your environment persistent, meaning those changes will define your environment each time you start a new bash session. This can be achieved by adding the variable assignment to one of the *environment files*, which are a set of files that are executed everytime you start a new bash session. These files are typically hidden, so we need to use `ls` wuth the `-a` flag to see them.

List all files in your homoe directory and locate the `.bash_profile` environment file, and view its contents with the `cat` command.
```bash
# navigate to your home directory
cd ../

# view files in current working directory and include hidden files
ls -a

# view contents of bash profile
cat .bash_profile
```

The `.bash_profile` is run everytime you start a bash session and contains variables used to configure the bash environment in a way that is specific to the contents of the `.bash_profile` file. You can add lines to the `.bash_profile` to set environment variables that will be established each time you start a new session. Lets add the command we created above to our `.bash_profile`.
```bash
# use the nano text editor to add the line ' ll="ls -lah" ' to your bash_profile
nano `.bash_profile`

# run the new bash_profile to set the environment variables (or start a new bash session)
source ~/.bash_profile

# now run the command as we did above
$ll
```

Now `$ll` will be set as an environment variable everytime we start a new bash terminal. It is also possible to avoid using the `$` to evaluate this variable by using the `alias` command in bash. `alias` allows you to set command that can be called directly using whatever characters you define, and can be added to your `.bash_profile` in the same way as we did above.
```bash
# make an alias for the ls -lah command
alias ll="ls -lah"

# call command directly with ll
ll
```

Another effective use of an alias is for accessing specific directories quickly. For example, if we had a project sub directory that we regularly want to access, such as `~/project/with/many/directories/`, we would need to write this out everytime to get there from our $HOME directory, using `cd /project/with/many/directories/`. Using an alias, we can save this command so that it is more easily callable.
```bash
# make a long directory that you may want to get to quickly in the future
mkdir -p ~/project/with/many/directories/

# make the alias for it
alias pd="cd ~/project/with/many/directories/"

# now call the alias
pd

# check your current dir
pwd
```

Again, just like above, we could add this line defining the alias command to our `.bash_profile` to make this alias available every time we start a new bash session, without even having to set it (after we have put it in our `.bash_profile`). Do this again with nano:
```bash
nano .bash_profile
```

### The $PATH environment variable

Another very important environment variable is `$PATH`, which stores a list of directories that tells bash where specific programs that we want to be available to us are stored. Programs are all essentially just files, and bash needs to know where these files are in order to run the commands as we call them.

The list is stored as strings separated by colons, so that many directories can be defined. Use `echo` to print `$PATH` variable.
```shell
echo $PATH

# make more readable using 'tr' to swap the colons for newlines
echo $PATH| tr ":" "\n"
```

As you can see, many of the directory names end in `bin` which standards for *binary*, which is a common directory name to store executables (programs).

Importantly, you can add directories to your `$PATH` as you either create or install programs, making them available to you as executables. Since the `$PATH` variable is set each time your `.bash_profile` is run at the start of a new session, the executables you add to `$PATH` will be available for you in a new bash session, without having to add them to your `$PATH` again.

We will create an executable file and add it to our $PATH in another lesson, however below is a toy example of how you would add a new executables directory to your `$PATH` variable:
```
export PATH="~/location/of/new/executables:$PATH"
```

A command for finding where a program lives in the $PATH is the `which` command. This can be useful for debugging environment issues as they arise when trying to use or instal new software. Check where the executable for the `echo` command is located:
```r
which echo
```

Many commands like `ls` will also accept wildcards, which are special character instances that allow you to do things like operate on multiple files at one time, or search for specific patterns (either in files or file names). We don't have time to review all the wildcard characters, however the most commonly used one is the asterisk, which can be used to represent any number of characters.
```bash
# list all files in my current directory with the file extension .txt
ls *.txt
```

### Customizing your environment

You will notice the prompt in your terminal when you are logged onto discovery starts with the term `(base)` what this is indicating is that the environments loaded in your .bash_profile are the tools that are available for you to use. For this workshop (and for most NGS data processing) you will need to extend the software packages that are available to you.

We will do this now by loading a new environment with the tool `conda`. We have pre-built this `conda` environment for you such that all of the tools you will need have been loaded into this environment, you should have created this environment with the commands included in the welcome and setup email. Tomorrow we will talk more about how to create your own custom `conda` environment.

```bash

# Load conda environment
conda activate bioinfo
```
This should change the word at the beginning of your prompt from `(base)` to the name of the conda environment that you just loaded `(bioinfo)`.

> As we move through the subsequent lessons, we will introduce more complex bash commands in order to manipulate common bioinformatics file types. If you are ever confused about what a command does, remember you can always use `man` to check out the manual page (or google it). It you are confused about how commands are used in conjunction with each other, it can also be helpful to break them down and run parts individually, in order to understand what the constituent parts do.

### Breakout room activities

- PRACTICE the bash commands - getting muscle memory for these commands and how to combine them and how they work are going free up your brain power to think about the analysis you want to perform rather than the commands you need to use. 
- Check out the cheat sheet links
