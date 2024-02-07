# Laboratory 4
## Sequence data and multiple sequence alignment

### Introduction

To this point in the course, you have been primarily estimating phylogenetic trees using morphological data. Although morphology is still used in phylogenetics, the majority of studies now rely on molecular sequence data (either nucleotides or amino acids) to estimate trees. This has become particularly true in the “genomic era,” because researchers can use millions of nucleotide or amino acid characters for their analyses. However, regardless if you are using a single gene with a few hundred characters or thousands of genes with millions of characters, there are several steps that need to happen before you can use proceed with a phylogenetic analysis. First, it is important to established that you are comparing homologous sequences, which we partly covered in last week’s lab on GenBank. A second important step is the multiple sequence alignment (MSA). This step ensures that you are not only comparing homologous genes, but that you are comparing homologous sites within a gene. Today, you will become more familiar with sequence data, including different file types and software for viewing sequences. You will then be introduced to multiple sequence alignment. 

**Objectives**

By the end of the lab, you should be comfortable with:

- Different types of files associated with genetic data
- Opening and viewing sequence data in SeaView
- Converting between different types of files
- Different programs and parameters for multiple sequence alignments
- How alignment effects phylogenetic inference

### _Exercise 4.1: Sequence data_

In this exercise, you will become familiar with different types of files used for storing sequence data. You have already worked with a few of these file types (FASTA and NEXUS), but this lab will introduce you to a few more.

To begin, log on to your Jetstream2 account and access your course instance through the terminal window. Change directories to your mounted volume (```/media/volume/sdb```) or your directory in the root directory, and create a new directory called ```Lab4```. Next, change into the ```Phylogenetics_Lab``` GitHub directory you downloaded last week and enter the following command:

```
git pull
```

This command will update the Phylogenetics_Lab directory from GitHub. For any remaining labs, you will use this command to download data for that lab. This will allow me to make changes to the GitHub content without you having to download the entire directory each time. If you list the contents of your updated GitHub directory, you should see a new directory called ```Lab4```. 

Copy the contents of ```Phylogenetics_Lab/Lab4/Data``` to your ```Lab4``` directory and extract the contents using the following command:

```
tar -xvzf Lab4_Sequences_MSA.tar.gz
```

Once you have extracted your files, let’s start by becoming familiar with a few different types of sequence files.

First, download the file “claravinae_cox1.fasta” and open it in a text editor. This is a FASTA file, which you should be familiar with from previous labs.

:white_check_mark: __Describe the basic structure of this file.__  

Next, download and open the file “apoid_28S_interleaved.fasta.” This is also a FASTA file, but is in the “interleaved” format (the previous file was a “sequential” FASTA file). Note that FASTA files usually have the extensions .fasta, .fa, or .fst.

:white_check_mark: __How does the interleaved FASTA format compare to the sequential format?__

Most phylogenetics programs that use FASTA files are fine with either type of format, but it is useful to know that there can be differences. This can be particularly useful if you ever need to write code to edit a FASTA file.

Next, let’s introduce you to viewing sequence files in a more informative interface. We will use the GUI program SeaView. 

***_SeaView_***

[SeaView](https://doua.prabi.fr/software/seaview) (Gouy et al. 2010) is a free program for viewing and editing molecular sequence (DNA or protein) data. It also allows you to do multiple sequence alignments (Part 2 of this lab) and phylogenetic tree building, but in this course we will primarily use it to view/edit sequence data. In today’s lab, you will get familiar with some of the basic features.

First, let’s view our sequential file in SeaView. Open up the program, select File -> Open Fasta, and select the “claravinae_cox1.fasta” file. You should see something like this:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/alignment.png" width=50% height=50%>

You are seeing the same information as when you opened the file in a text editor, but now there is more information and editing options for the sequences. Your taxon names are on the left, followed by their DNA sequence. Each base has a unique color. Viewing sequences in this way allows you to quickly assess the content of your data file, which will be very useful for viewing sequence alignments. Here are a few useful tools you can use in SeaView:

1. Select taxa by clicking on the name, or select multiple taxa by clicking and dragging across multiple names.
2. Rearrange the order of names by selecting a name or names, holding down CTRL, and clicking on the new location in the order.
3. Edit your sequences by selecting Props -> Allow seq. editing, clicking on a particular nucleotide, and editing with backspace, space, or replacement (A, T, C, G, or an ambiguity code). Use caution when selection this option, because you can easily accidentally delete data. At the very least, always de-select the option when you are finished editing.
4. View your sequence as amino acids by selecting Props -> View as proteins. This will translate your nucleotides to amino acids (although in this instance there will be stop codons, since the data is not in reading frame)

Take a few minutes to get familiar with these options. You can save the file (File -> Save) if you want, but note that SeaView only saves FASTA files in interleaved format.

After becoming familiar with SeaView, let’s use it to introduce you to a few other types of sequence files.

First, let’s look at NEXUS files, which you should be familiar with from using PAUP*. Save your FASTA file as a NEXUS by selecting File -> Save as… and then select “NEXUS (*.nxs)” under the “Save As:” options. 

:white_check_mark: __Save the file as LastName_cox1_nexus.nex.__

NEXUS files usually have the extension .nex or .nxs. Open the new NEXUS file in a text editor. The first line should say “#NEXUS” followed by some information about your data.

:white_check_mark: __How many taxa (NTAX) are in your file? Characters (NCHAR)? What is the data type?__

One of the key characteristics of a NEXUS file are data blocks in the “BEGIN…END” format. In your current file, you should see a “data block;” anything between “BEGIN DATA” and “END” is in the data block. Also notice the semicolons (;), which demarcate different sections of the file.

:white_check_mark: __Other than the block format, what are some differences between FASTA and NEXUS files?__

Finally, let’s introduce you to PHYLIP files. The PHYLIP format was originally introduced by Joe Felsenstein for his software program PHYLIP (surprise!), which was one of the first computer programs for inferring phylogenetic trees. PHYLIP (the program) is not used as frequently anymore, but the associated file format is still used by many different software programs. 

Go back to your sequence file (either FASTA or NEXUS) in SeaView and convert it to PHYLIP by selecting File -> Save as… and choosing “Phylip (*.phy).” 

:white_check_mark: __Save the file as LastName_cox1_phylip.phy.__

Phylip files usually have the extension .phy. Open the new PHYLIP file in a text editor. You should see something that looks like this:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/phylip.png" width=50% height=50%>

The top of a PHYLIP file should have two numbers: the number of taxa (19 in this case) and the number of characters (420). The data starts in the second line with the taxon name followed by the sequence data. Note that PHYLIP files can also have several different formats. Like with FASTA files, there are “interleaved” and “sequential” PHYLIP files. There are also “strict” and “relaxed” PHYLIP files. The strict format only allows a maximum of 10 characters in your taxon names , whereas there is no limitation in the relaxed format. SeaView saves PHYLIP in the interleaved, relaxed format.

### Exercise 4.2: Multiple Sequence Alignment

Now that you are more familiar with viewing, editing, and exporting sequence data, let’s move on to multiple sequence alignments (MSA). MSA is a very important first step for inferring phylogenetic trees from sequence data, as it helps ensure you are comparing homologous characters. However, MSA is often overlooked or quickly breezed over in phylogenetic analyses, which can lead to erroneous trees. In this lab exercise, you will be become familiar with several different programs for performing MSA, along with some best practices for using these programs.

**Alignment Methods**

_Manual Methods_

With morphological data, it can be relatively straight-forward to determine which states in one taxon are homologous with which states in another taxon because, generally, the available states are defined by the context of the character. But with DNA the matter is more complex because each possible character has the same four possible states (A, C, G and T) , and protein sequences have 20 (number of possible amino acids in eukaryotes)  . DNA alignment is the process of assigning a priori hypotheses of homology to your sequence data, which will then be used to generate a phylogeny.

Consider the following example:

```
TaxonA	ACTTCCGAATTTGGCT
TaxonB	ACTCGATTGCCT
```

These data do not immediately indicate how the states in these two taxa should be homologized (or aligned). Three examples of aligned data are provided:

```
ACTTCCGAATTTGG-CT
|||  |||  |||  ||
ACT—-CGA--TTG-CCT
```
```
ACTTCCGAATTTGG-CT
|||*    *||| | ||
ACTC----GATT-GCCT
```
```
ACTTCCGAATTTGGCT
|||*    **|||*||
ACTC----GATTGCCT
```

Which one is correct? In the first example, the alignment was done in a manner that would minimize base substitutions. In the second example, the alignment was done in a manner that would balance the number of insertions and deletions (INDELS) with the number of base substitutions (indicated by asterisks). In the third example, the alignment was done to minimize INDELS.

The problems with by-eye (or manual) alignment are:

1. Different alignments can lead to different phylogenetic trees – however, these differences are merely a matter of opinion, not empirical science.

2. It is hard, if not impossible, to align by eye without allowing one’s pre-conceived notions of relationships to come into play.

Thus, what is needed is some objective method of aligning sequences that is repeatable and logical. By manually aligning sequences there is a trade-off between the cost of a base substitution and the cost of an INDEL. 

_Algorithmic Methods_

A cost can be assigned to base substitutions and gaps. Given the cost of a substitution is 1 and the cost of opening a gap (no matter the length) is 1, the final cost of the alignment given below is 5 [0(1) + 5(1) = 5].

```
ACTTCCGAATTTG-GCT
|||  ||| || |  ||
ACT—-CGA-TT-GC-CT
```

In the next example, the final cost is still 5 [3(1) + 2(1) = 5].  (The three substitutions are noted by asterisks.)

```
ACTTCCGAATTTGGCT
|||*   *||| |*|| 
ACTC---GATT-GCCT
```

If opening a gap costs the same as a base substitution, all of the disagreements between the two sequences can be explained by insertions and deletions (the first alignment immediately above). However, this comes at the expense of hypothesizing no base substitutions and, thus, eliminates any phylogenetic information that could have been present in those parts of the data.  

Usually, you will want to set the cost of opening a gap to be higher than the cost of a substitution. For example, if the cost of a substitution is 1 and the cost of opening a gap is 2, the final cost of the alignment immediately below will be 10 [0(1) + 5(2) = 10].

```
ACTTCCGAATTTG-GCT
|||  ||| || |  ||
ACT—-CGA-TT-GC-CT
```
And the cost of the next alignment will be 7 [3(1) + 2(2) =7].

```
ACTTCCGAATTTGGCT
|||*   *||| |*||
ACTC---GATT-GCCT
```

The second of these two alignments is preferred because it has the lower cost.

Now, let’s consider the alignment that minimizes gaps. A base substitution costs 1 and opening a gap costs 2. The final cost for each of the three alignments below is given after the alignment.

```
ACTTCCGAATTTG-GCT  	
|||  ||| || |  ||
ACT—-CGA-TT-GC-CT
```
Final cost = 0(1) + 5(2) = 10
```
ACTTCCGAATTTGGCT  	
|||*   *||| |*||
ACTC---GATT-GCCT
```
Final cost = 3(1) + 2(2) = 7
```
ACTTCCGAATTTGGCT  	
|||*    **|||*||
ACTC----GATTGCCT
```
Final cost = 4(1) + 1(2) = 6

The last of these three alignments would be preferred based on the assigned costs.

_Computational Alignment_

In the above examples, we considered only short sequences and two taxa. In reality, however, the sequences being compared are many and quite long; thus, one needs an efficient way to discover where the gaps are. Fortunately, automated methods exist which can help you align long sequences for many taxa and determine final cost estimates. Today, there are many programs that can perform multiple sequence alignment (MSA). We will investigate only a small fraction of the available programs today. For a full discussion of the many approaches to MSA, see [this wiki page](https://en.wikipedia.org/wiki/Multiple_sequence_alignment).

Based on many studies, the programs [MAFFT](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html) (Multiple Alignment using Fast 
Fourier Transform has both progressive and iterative alignment capabilities) and [MUSCLE](https://www.drive5.com/muscle/) (an iterative alignment) are relatively accurate. Another widely used program is Clustal. Clustal is updated frequently and they have a new MSA program, [Clustal Omega](http://clustal.org/omega/). Omega uses HMM (hidden Markov Model) analysis. We will discuss Markov processes when we discuss Bayesian methods of phylogenetic analysis later in the semester.

We will explore all three programs today and compare the resulting alignments. Although MUSCLE, MAFFT, and Clustal-Omega can also align amino acid data, other programs exist specifically for the alignment of amino acid data (and thus usually perform better), including [PROBCONS](http://probcons.stanford.edu/), but we will not talk about these programs today.

In Clustal alignments, all sequences are compared to each other using pairwise alignments, which are used to construct a neighbor-joining distance phylogeny. The phylogeny is then used as a guide for subsequent multiple sequence alignment, with the most closely related taxa aligned first.

This is very different from the iterative approach used by MUSCLE, which compares two unaligned sequences and generates a measure of similarity by finding the number of k-mers or “words” in common – that is, how many small fragments of the sequence that the two sequences have in common. This produces a good guide tree (using UPGMA, a distance method) very quickly. 

MAFFT uses a similar iterative approach; however, its refinement process is different. On some datasets MAFFT may perform better than MUSCLE, and vice versa, so exploring both options is a good idea for your project data. It offers a range of multiple alignment methods, L-INS-i (accurate; recommended for <200 sequences), FFT-NS-2 (fast; recommended for >2,000 sequences), etc. (see http://mafft.cbrc.jp/alignment/software/about.html).

Some factors to take into consideration when aligning DNA sequences are:

1. Most alignment software packages will give you a single multiple alignment. However, there may be equally optimal alignments for the same set of data no matter what your optimality criterion is.  

2. Protein-coding sequences may be relatively easy to align using the amino acids themselves as references. In contrast, non-coding sequences (such as ribosomal DNA, introns, etc.) are usually more difficult to align, as are more distantly diverged sequences. In this laboratory, you will obtain experience aligning and working with all of these types of data. 

**Aligning Data**

First, let’s compare different alignment software by aligning the same sequences with each method. Navigate back to the Jetstream terminal window and make sure you are in your ```Lab4``` directory. We will start by aligning the file "apoid_ef1a.fasta". The file contains sequences of the nuclear gene _EF-1α_ (elongation factor 1, subunit alpha) from different species of bees and wasps (superfamily Apoidea). The _EF-1α_ gene is a frequently used locus in phylogenetic studies of eukaryotes.

Download the file "apoid_ef1a.fasta" and open it in SeaView. Notice that the sequences are not aligned.

We will first align these sequences using the default parameters for Clustal-Omega, MUSCLE, and MAFFT.

Run the following command for Clustal-Omega:

```
clustalo -i apoid_ef1a.fasta -o apoid_ef1a_clustalo.fasta 
```

For MUSCLE:

```
muscle -in apoid_ef1a.fasta -out apoid_ef1a_muscle.fasta
```

And finally for MAFFT:

```
mafft --auto apoid_ef1a.fasta >apoid_ef1a_mafft.fasta
```

MAFFT has several different algorithms designed for different sizes and types of data sets. When in doubt about which algorithm to use, run MAFFT with the ```--auto``` parameter and the program will choose for you based on your data. 

:white_check_mark: __Which strategy did MAFFT use (look below “Strategy” near the bottom of the terminal screen)?__

Download all three output files and open them in SeaView.

:white_check_mark: __What differences, if any, are noticeable among these alignments? Note that the name order might be different in different files.__

Now let’s look at a dataset that creates more alignment ambiguities: apoid_28S.fasta. EF-1α, the gene used in the previous analysis, is a nuclear protein-coding gene and the indels occurred exclusively in introns. This new file includes sequences for these taxa based on _28S_, a nuclear ribosomal RNA (rRNA) gene, which has multiple repeats across the genome of eukaryotes and evolves more quickly than protein-coding genes for most organisms. The highly variable areas in this gene correspond to [loop regions](https://rnacentral.org/rna/URS0000ABD82A/9606?tab=2d) of the _28S_ secondary structure. Stem regions tend to be more conserved (i.e., less variable)

Align these _28S_ sequences using the default parameters in MUSCLE:

```
muscle -in apoid_28S.fasta -out apoid_28S_muscle.fasta
```

And with the ```--auto``` setting in MAFFT:

```
mafft --auto apoid_28S.fasta >apoid_28S_mafft.fasta
```

Download the output files and open them in SeaView.

:white_check_mark: __What is a likely example of a loop region (provide position ranges, which can be seen by clicking on a particular base in SeaView)? A stem region?__

Next, let’s see how changing alignment parameters effects an alignment. We will do this by changing the “gap open” parameter in MUSCLE. Higher values of this parameter will increase the penalty of introducing gaps into the alignment, so lower values should result in more gaps. 

Re-run the MUSCLE alignment with the following command, which sets the gap open penalty to -5:

```
muscle -in apoid_28S.fasta -gapopen -5 -out apoid_28S_muscle_gap5.fasta
```

Download the output file and open it in SeaView. 

:white_check_mark: __How does this alignment compare to the alignment from default parameters?__

Try running MUSCLE again with parameters of your choice.

:white_check_mark: __What were your parameters?__  


**Impact of Alignment on Phylogeny**  

We’ve seen that different alignment programs or parameters can produce different alignments from the same sequence data, but do these differences really impact phylogenetic inference? Let’s run some tests to find out.

First, save your MAFFT and MUSCLE 28S alignment files as NEXUS files using SeaView (should be four files total). Next, let’s estimate MP trees using PAUP*. Launch the PAUP* program and load your MUSCLE alignment generated using default parameters (“apoid_28S_muscle.nex”). If your alignments are saved as .nxs, make sure you select “All files (*.*)” as the file type.

Select **Data: Define Outgroup** and choose _Clypeadon taurulus_ and _Eucerceris tricolor_ (wasps closely related to the ingroup, bees). Run a boostrap analysis (**Analysis → Bootstrap/Jackknife Analysis**) with 100 replicates, random addition, and TBR branch swapping. 

:white_check_mark: __Compute a 50% Majority Rule consensus tree, and save the tree as LastName_apoid_28S_muscle_tree.pdf.__

Repeat this analysis for the MUSCLE alignment with custom parameters and for the MAFFT alignment. Remember to define your outgroup for each analysis! 

:white_check_mark: __Save these consensus trees as LastName_apoid_28S_muscle_custom_tree.pdf and LastName_apoid_28S_mafft_tree.pdf, respectively.__

:white_check_mark: __Are there any notable differences among the three trees? If so, what are they?__

:white_check_mark: __Are there any well-supported differences (bootstrap support >75)?__

:white_check_mark: __Based on your findings, which alignment algorithms produce the most similar trees?__

### Some other MSA software
[PASTA](https://github.com/smirarab/pasta)  
[PRANK](https://github.com/ariloytynoja/prank-msa/tree/master)  
[UPP](https://github.com/smirarab/sepp/tree/master)  
[T-Coffee](https://tcoffee.org/Projects/tcoffee/index.html)

### Other software for viewing and analyzing alignments
[AliView](https://ormbunkar.se/aliview/)  
[MEGA](https://www.megasoftware.net/)

### References

Gouy, M., S. Guindon, and O. Gascuel. 2010. SeaView version 4: a multiplatform graphical user interface for sequence alignment and phylogenetic tree building. Molecular Biology and Evolution. 27: 221-224.
