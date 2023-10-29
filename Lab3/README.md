# Laboratory 3
## Lab 3.1: Using the Command Line

## Introduction

Learning how to use the command line has become a crucial skill in phylogenetics. Many phylogenetics programs can only be run through the command line, so it is important to have at least a basic understanding of the syntax. Many of our labs this semester will rely on using the command line. 

You can think of the command line as a more direct way of interacting with your computer compared to using applications in a Graphical User Interface (GUI). Instead of pointing and clicking on an interactive interface, you will give the computer commands through a type of program called a shell. There are several types of languages (specific syntax) for using a shell, but the most common one for phylogenetics (and bioinformatics more generally) is called bash (for Bourne Again Shell). Bash was designed for a type of shell called Unix. For context, Linux and MacOS are two operating systems that use a Unix-like system. Windows uses a system called DOS. 

First, let’s introduce you to some basics. To begin, log on to Jetstream and access the Web Shell from your instance (VM). 

The bottom line on your shell window should end with an ```$``` or ```#``` followed by a white rectangle:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/dollar1.png" width=10% height=10%>

![Shell prompt image](../Images/dollar1.png?raw=true width=50% height=50%)

This means your shell is ready to take commands. If you do not see these symbols, then your computer is likely running something or is otherwise not ready to receive commands. To submit a command, type in your command and hit Enter (Return).

### _Files, directories, and paths_

Like in a GUI interface, you use the command line to access files from folders, subfolders, etc. Folders are often referred to as “directories.” However, in the shell window, you can only be in one location (directory) at a time. Your current location is referred to as a “working directory.” To move around or access files in other directories, you need to specify the path of that file or directory, which is basically a list of directories and subdirectories separated by slashes (/). For example, File X might be located in directory A/A.1/A.1.1. 

**Note: Never label your directories or files with spaces in the names. This can make it difficult to type the path in your command line. Instead, label files and directories with underscores (_) or dots (.) to separate words, etc.**

You will often see/hear paths referred to in directions like “up” or “down” (e.g., move up one directory). This is because using directions is a good way to visualize the path structure:






Thinking about paths and directories in terms of “up” and “down” also allows for using relative paths (relative to your current working directory). In bash, you can utilize relative paths with the following syntax:

```.``` Current working directory

```..``` Up one directory

```~``` Home directory 

Using the above example, if your current working directory is A.1.1, you can access A.1 by referring to ```..```.

In contrast, absolute paths require you to list the entire path starting at the root (the “uppermost” directory). So, to access A.1 from A.1.1 using an absolute path, you would need to reference /A/A.1.

### _Common Bash Commands_

Next, let’s cover some basic bash commands. You will likely use these commands the most frequently throughout the course. Read through the different commands and complete this week’s lab by going through the exercise.
```
pwd
```
Print working directory; writes your working directory to the screen

```
ls
```
List; lists the contents of your current directory

```
cd <path>
```
Change directory; moves you to a different working directory

```
cp <file> <destination>
```
Copy-paste; copies a file from one location to another

```
rm -i <file>
```
Remove; deletes a file; THIS IS PERMANENT (NO TRASH)

**Note: the “-“ (called a “dash” or “flag”) usually indicates specific options for different commands. We will see plenty of “flag” options throughout the semester with various programs. In this instance, the -i flag tells the computer to make you confirm your deletion before proceeding. Without this option, rm will simply remove your file, no questions asked.**

```
mkdir <directory>
```
Creates a directory

```
touch <filename>
```
Creates a blank file

For any bash command, you can type ```man <command>``` to get the manual page for that command. This will give you the basic usage instructions, along with all the possible options (likely more than you would ever use).

























### _Excercise 3.1.1_

Go through the follow exercise to practice working with the command line on Jetstream. Answer the questions on the separate worksheet document and submit to Blackboard to receive credit for the lab.

First, navigate to your mounted volume directory (probably ```/vol_b```).

Create a new directory called ```Lab2```.

Throughout the semester, you will be using git commands to access data for the labs. Git commands allow you to download data from my data folder for the course, which is located on GitHub (https://github.com/; a free data depository usually used to make computer code publicly available). 

Once you are in your volume directory, use the following command:
```
git clone https://github.com/adsweet/A-State_Phylogenetics_Lab.git
````
You should see a new directory in your volume directory. 

:white_check_mark: __1.	What is the name of this new directory?__

Change into the new directory. You should see a directory called ```Lab2```. Now change into this directory. You should see a directory called ```Data```. Change into the ```Data``` directory. You should see a file. 

:white_check_mark: __2.	What is the name of the file?__

Now let’s copy this file into your ```Lab2``` directory. After copying the file, change into your ```Lab2``` directory. This file is a compressed TAR file (short for “tape archive”), which is basically a type of file that can have multiple files and/or directories compressed inside of it (often referred to as a “tar ball”). To extract the contents, we need to use a “tar” command. Use the following command:

```
tar -xvzf <tar file>
```
(x = extract, v = verbose, z = unzip, f = file)

:white_check_mark: __3.	What is the name(s) of the new file(s) in the “Lab2” directory?__

Next, we want to see what is in the file. There are several commands for this. To take a peak into the file, use the following command:
```
head <file>
```

This will print the first 10 lines of a file. If you want to see more or fewer lines, you can include that number after the command (e.g., ```head -5 <file>```).

:white_check_mark: __4.	What do you see in the file? Paste the first two lines below.__

For your reference, you can also see the last ten (or last X #) of lines using “tail” instead of “head.” 

:white_check_mark: __5.	Paste the last two lines below.__

You can also use one of several built-in text editors. We will be using nano (because I think it is the easiest to use), but know that there are other options (vim, emacs, etc.).

Open the file by typing:

```
nano <file>
```

This should print the contents to your screen, along with a cursor (white rectangle) and some commands listed at the bottom of the screen. You can navigate around the file using the arrow keys. Nano also allows you to edit files by typing in characters behind the cursor. Add some text to the first line of the file. Now, exit the file (CTRL + X). You will be asked whether you want to save the file (Y for Yes, N for No). Choose Yes and hit Enter to save. Now check your edited file with “head.”

:white_check_mark: __6.	What does your first line now say?__

Finally, let’s practice creating and deleting files. Create a file called “practice.txt”. Open the file with nano and add some text. Save the file and view the first few lines with “head.”

:white_check_mark: __7.	What is the first line of your new file? Paste below.__

List the content of your “Lab2” file, but this time add the “lh” flag to your command:
```
ls -lh
```
This gives you more information about the contents, including the creator (username), size (number of bytes (K = KB, M = MB, etc.)), creation/edit date and time, and name. 

:white_check_mark: __8.	How large are the three files in the Lab2 directory?__

Download your “practice.txt” file to the Desktop using the Jetstream interface. 

Sadly, your new file was not meant to last long. As your final act in this exercise, delete your “practice.txt” file. List the contents of ```Lab2``` to confirm you removed the file.

## Lab 3.2: NCBI

## Introduction

Many of you in this course do not yet have your own sequence data, so you will need to access the vast stores of information on the Internet in order to complete your project. In the same way that systematists keep type specimens in museums, molecular sequences are also archived in publicly accessible repositories, such as the National Center for Biotechnology Information (NCBI) (http://ncbi.nlm.nih.gov). NCBI’s website houses a set of extensive databases and web-based applications for systematics, bioinformatics, and genetics. Together, the NCBI databases contain an immense wealth of information, including extensive reviews of gene functions, taxonomic authority data, and genome annotations, as well as records of every published molecular sequence (including genomic DNA, cDNA, RNA, and protein). GenBank has collaborative associations with other international databases: DNA DataBank of Japan (DDBJ) and the European Molecular Biology Laboratory (EMBL).

In this exercise, you will learn how to use several of the tools available from NCBI. Knowing how to navigate this resource will be an incredibly useful skill for working with molecular sequence data to estimate phylogenies.

To begin, access the NCBI home page by going to this link: https://www.ncbi.nlm.nih.gov/. You should see a page that looks like this:
















The top of the page has some useful links under the “Resources” and “How To” menus. The left side of the page lists the resources available. The right sidebar lists some popular resources that are frequently accessed, including PubMed (a literature database) and BLAST (a sequence similarity and alignment tool). The search bar, located in the center top of the page, features a pulldown menu for choosing a specific search database and a form for queries.

### _Exercise 3.2.1_

For the first exercise, let’s practice searching for specific sequences on GenBank, one of the databases on NCBI. We are going to search for a nucleotide sequence, so select “Nucleotide” in the pulldown menu. 














The nucleotide database is extensive, so we want to be very specific about our search, both in terms of the gene of interest and the taxon we are studying. Let’s say we want to find out if a sequence exists for the gene Cytochrome Oxidase Subunit 1 (COI) in the species Zenaida macroura. Type the name of the organism (Zenaida macroura) and the name of the gene (COI). The query will pull up a list of results. You will see an accession number, a description of the sequence, and a variety of other data in this format:






The first information given is the sequence name.  This has several parts, including the species, voucher number, gene region, gene origin, and whether or not the gene is complete (partial cds stands for partial coding domain sequence). The accession number (HM033982) is listed in the last line. The GI number is a series of digits assigned consecutively to each sequence record processed by NCBI.

Click on the sequence name to open the record for the sequence, which contains information about this particular sequence, the taxonomy of the organism from which the sequence was obtained, any publications using this sequence, and the sequence itself.

:white_check_mark: __1.	What is the reference for the paper in which this sequence was published?__

:white_check_mark: __2.	How many base pairs exist for this sequence?__

You will also be to see taxonomic information for the organism under “Organism.”

:white_check_mark: __3.	What is the taxonomy of this particular organism (copy and paste from the record)?__

:white_check_mark: __4.	What type of animal is it?__

You can also easily download data from GenBank in several different file formats. Let’s download a sequence in FASTA format. FASTA files are a type of text file for storing genetic information. The basic syntax for a FASTA file is a carrot (“>”) followed by a string of characters (the “header,” usually the sample information), and finally the sequence data on the next line. When in doubt about a file format for sequence data, if you see the “>” and a header, you know it’s a FASTA file.

You can download individual sequences by selecting the “Send” pulldown menu in the upper right-hand corner of the window. Select “Complete Record”, with a destination of “File.” You will now be able to select “FASTA” format. This will download the sequence and label data directly to your computer.












Open the file in a text editor. 

:white_check_mark: __5.	Copy and paste the contents of the downloaded file in your assignment file.__

You can also save multiple records to a single FASTA file. Go back to your search results. You can do this by choosing your search in the Recent Activity field on the right side of the page.     







Choose the last search you performed. This will return you to your search results. You will notice that there are checkboxes next to each item. Click the check boxes of a few items to select them.













Now, at the top of the page, there are a series of drop-down menus. On the drop-down menu labeled “Send to,” select “File” and on the pulldown menu labeled Format, select “FASTA.” Sort the sequences by “Default order.” This will save the selected sequences to a single FASTA file. Open that file in a text editor to verify it is correct.

If you’re interested in finding all available data for a particular species (or higher-level taxon), it’s most efficient to search GenBank with the “Taxonomy” tools. In the pull-down menu left of the search bar, select “Taxonomy” and search for “Zenaida macroura.”













Click the link for that species. On the right-hand side of the screen, you will see a table called “Entrez records.” This table lists the different NCBI databases and the number of links that exist for “Zenaida macroura.” Clicking on the numbers will take you to a page with links to those data.

:white_check_mark: __6.	How many genomes are available for Zenaida macroura?__

:white_check_mark: __7.	What type of genome?__
