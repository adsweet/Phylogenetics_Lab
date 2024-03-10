# Laboratory 8
## Phylogenomics

### Introduction

Until this point, we have been working with relatively small data sets in our labs. By “small data set,” we mean using one or a few genes for estimating phylogenetic relationships. However, recent years has seen the “genomic revolution,” where large data sets have become more cost effective to generate in non-model organisms. Because of this, most phylogenetic analysis seeking to reconstruct evolutionary relationships among different taxa have become phylogenomic studies. “Phylogenomics” can broadly refer to using several different genes for a phylogenetic analysis, but more often it refers to using hundreds or thousands of genes/loci. Note that phylogenetic studies focused on the evolution of genes or gene families can still use genomic approaches to obtain data, but the questions do not involve using hundreds of genes.

Large data sets can be very useful for estimating a well-resolved phylogeny. However, these data present their own challenges. First, it can be challenging to determine if genes are orthologous. Smaller data sets are often generated with primers through Sanger sequencing or are comprised of specific targeted genes, so you are more confident in your assignment of orthology. With larger data sets, you are often assembling genes through “de novo” approaches or mapping them to a reference, rather than targeting with primers or selecting specific genes. Second, when you’re working with hundreds or thousands of genes, it is unfeasible to look through each alignment. You put a certain level of trust in your initial alignment parameters and use summary statistics or additional algorithms to assess your alignments. Essentially, with phylogenomics you rely very heavily on computational algorithms even before you get to the tree-building step.

Today, we will briefly introduce you to the data involved with phylogenomic analyses and cover a few approaches for analyzing and filtering large datasets.

### Exercise 8.1 

### Genomic data

There are several technologies for generating genomic-level data sets, but the most frequently used is the “short read” technology from Illumina. With Illumina sequencing, the genome is fragmented into many little pieces (usually 400-500 base pairs), specific adapters are added to the ends of these fragments, and ~150 base pairs are sequenced from each end of the fragment. This can be done for a specific region of the genome (e.g., mitochondrial genome, 16S gene, conserved nuclear genes, restriction-site associated regions) or for the entire genome (whole genome shotgun sequencing). The result is a file or files with millions of ~100-160 base pair fragments that together makeup the whole genome or targeted region(s). You would then assemble these short “reads” by themselves (de novo assembly) or map them to a closely related reference to obtain your data for phylogenetic analysis.  

Let’s look at some data that was generated from genomic sequencing. First, log on to your Jetstream instance and launch the Web Shell. Navigate to your volume or home directory and create a directory called `Lab8_phylogenomics`.

Next, navigate to the `Phylogenetics_Lab` directory and get the most update data with `git pull`.

Copy the contents of the `Lab8/Data` directory to your `Lab8_phylogenomics` directory and extract it using:

```
tar -xvzf Lab8_phylogenomics_data.tar.gz
```

You should see several files, including FASTA (.fas, .fa) and FASTQ files.

Download one of the FASTQ files and open it in a text editor. These are Illumina whole-genome shot guns data from a louse in the genus _Ibidoecus_, which parasitize ibises.

:white_check_mark: __Have you seen an ibis before (you can use a Google search)?__

:white_check_mark: __What are some main differences that you notice between FASTA and FASTQ files?__

The basic format of a FASTQ file consists of sets of four lines: 1) a header line beginning with “@” , 2) the sequence, 3) a “+”, and 4) quality information. The quality data consists of various ASCII characters which correspond to different PHRED score values, which basically give a confidence level for a particular base.  The sequence and quality lines should be the same length. Each set of four lines is for a single read generated from the Illumina short read sequencing.

:white_check_mark: __How long are the reads?__

We can also use some basic Bash scripting to determine the number of reads in a particular file or files. Use the following command to count the number of lines in your FASTQ files:

```
ls -lh | wc -l *.fastq
```

This command lists the files ending in .fastq and “pipes” (`|`) the output of that command to a second command for counting the number of lines.

:white_check_mark: __How many reads do you have in each FASTQ file (hint: remember how many lines represent a single read)?__

### Trimming genomic data

After getting your raw FASTQ files, you generally want to trim the reads to remove adapters added for the sequencing process (and thus do not represent actual sequence from your organism) and remove bases that might have particularly low quality. There are several programs for trimming FASTQ files, but we will use Trimmomatic. 

Run the following command to use Trimmomatic (the command should be on a single line when you copy it to the shell):

```
fastp --in1 Ibidoecus_1.fastq --in2 Ibidoecus_2.fastq --out1 Ibidoecus_1P.fastq --out2 Ibidoecus_2P.fastq --adapter_fasta TruSeq2-PE.fa --cut_front --cut_tail --length_required 75 --dedup
```

This command will trim paired-end data (PE) by removing the adapters indicated in the file `TruSeq2-PE.fa`, remove bases below a PHRED score of 3 from the 5’ (LEADING) and 3’ (TRAILING) ends of each read, remove chunks of 4 base pair-long regions that have an average PHRED score below 15, and remove any reads less than 75 base pairs after the trimming steps. Download and open the `TruSeq2-PE.fa` file in a text editor to view the adapters and other sequences that wil be trimmed from the reads.

After running Trimmomatic, you should see four new files ending in `P.fastq` and `U.fastq`. These are the trimmed paired (P) and unpaired (U) reads. Download one of the trimmed paired read files and open it in a text editor.

:white_check_mark: __Do you notice any differences between the trimmed and untrimmed file? Are all the reads the same length?__

:white_check_mark: __How many reads are there? What happened to any missing reads?__

### Aligning genomic data

After trimming your genomic reads, you would be all set to assemble or map the reads to obtain genes or other markers for phylogenetic analysis. We’ll skip that step for this lab, and instead work with several genes that have already been assembled. These are in the five files ending in .fas (one file per gene). Download one of these files and open it in SeaView.

:white_check_mark: __How many taxa are in this alignment?__

These are sequences from different species of “hemipteroid” insects.  

:white_check_mark: __What is an example of a “hemipteroid” insect (again, Google searching is fine)?__

Hopefully, you should notice that this gene is not aligned. The same is true for the other four gene files. To align them, we could run each file through your alignment program of choice, but there is a simpler way to align each file with a single command (imagine if you had hundreds or thousands of gene files!). To do this, we will use a for loop.

For loops are a way to “loop” through a series of files and run the same command or set of commands on each file. The basic syntax of a for loop in Bash scripting is: 

```
for <object> in <files>;
  do <command>;
done.
```

Let’s use this handy scripting tool to align each of our ten files:

```
for file in *.fas;do
  mafft --auto $file >$file.aligned.fasta;
done
```

You can type each of these lines separately or paste them all in at once. In this for loop, we assign each file ending in .fas to the variable “file” (“file” can be any word you want; I usually use “file” for simplicity), use each of these files as input into MAFFT, and output the resulting aligned files. In Bash scripting, anytime you set a variable, you recall that variable by putting a “$” before the word or string of characters you chose. In our case, anytime we used the command “$file”, it is actually the name of a file ending in .fas.

This loop will take ~5 minutes, but when it’s finished you should have five new files ending in *aligned.fasta. Download one of these aligned files and open it in SeaView.

:white_check_mark: __Is the new file aligned?__

We will be using for loops in future labs as well; they are very useful for working with large data sets!

### Filtering alignments

It’s still feasible to check five alignment files by-eye, but you can imagine it gets very difficult if you’re working with hundreds or thousands of files. One way to run a quality check on your alignments is to filter them using some criteria. There are several programs you can use to trim based on alignment scores or gaps, but today we will use the program Trimal.  We will keep it simple today and trim alignments based on the proportion of gaps at a particular site. 

To trim your alignments, run the following command:

```
for file in *.aligned.fasta;do
  trimal -in $file -gt 0.5 -out $file.trim.fasta;
done
```

This will use the program Trimal to remove sites that have gaps (gt = gap threshold) in 50% or more of the taxa. Download a trimmed alignment and open it in SeaView.

:white_check_mark: __How much shorter is the trimmed alignment than the untrimmed alignment?__

:white_check_mark: __Were any sequences removed because they only contained gaps?__

Trimal also has several automated methods for selecting filtering parameters. Let’s try filtering again with the “automated1” option, which determines and applies the most appropriate automatic method:

```
for file in *.aligned.fasta;do
  trimal -in $file -automated1 -out $file.trim.auto.fasta;
done
```

Ignore any error messages that print to your screen.

:white_check_mark: __How do the 50% threshold and automated filter alignments compare to each other?__

### Phylogenomic analysis

Once you have trimmed your alignments, you are finally ready to run a phylogenetic estimation with your data. We’ll run a partitioned concatenated alignment, so we’ll need to concatenate all of our trimmed alignment files and generate a partition file.

First, make a new directory called `concat` and copy all of your filtered alignment files to this new directory. You can choose either the 50% or automatically trimmed alignments. __Change into the concat directory__ and run the following command to concatenate your alignments and obtain a partition file:

```
amas concat -f fasta -i *.fasta -d dna -u fasta --part-format raxml
```

Change some of the resulting file name:

```
mv concatenated.out hemipteroid_concat.fasta
mv partitions.txt hemipteroid_concat_partitions.txt
```

The partitioned analysis will take too long to finish in lab, but we’ll want a concatenated tree for next week’s lab. Open new a screen:

```
screen -S iqtree
```

and run the following command:

```
iqtree -s hemipteroid_concat.fasta -p hemipteroid_concat_partitions.txt -m MFP+MERGE -B 1000 -pre hemipteroid_concat
```

:white_check_mark: __Briefly describe what you see on the screen.__

Detach from your screen (`CTRL+a d`). We’ll see how the tree looks next week!
