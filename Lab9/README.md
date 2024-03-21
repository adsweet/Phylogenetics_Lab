# Laboratory 8
## Gene Trees and Networks

### Introduction 

Throughout this semester, we have been estimating phylogenetic relationships using a single gene or by concatenating several genes together. However, there can be variation among different parts of the genome that result in variation in phylogenetic trees; in other words, there is not necessarily “one tree” that represents every gene or genomic region. Phenomena such as horizontal gene transfer (HGT), introgression through hybridization, and Incomplete Lineage Sorting (ILS) can generate this variation. One of the ways to account for this variation is by analyzing gene trees and/or through phylogenetic networks.

A gene tree approach can account for ILS by using the Multispecies Coalescent model. This is a population genetics model that accounts for variation in by using the probabilities that different alleles coalesce a certain amount of time in the past. Thus, it is a backwards-looking model. The multispecies coalescent seeks to use this model to summarize different gene trees into a single species tree. There are several approaches and algorithms for doing this, but some of the most commonly-used rely on quartets (groups of four tips in different genes trees) to summarize many gene trees. Because they are summary methods (and are not inferring the gene trees), these methods tend to be very fast, especially relative to other methods that simultaneously estimate the gene tree and species trees. 

Another way to account for variation in trees is to use phylogenetic method approaches. It could be that your group of interest does not evolve in a tree-like manner (because of HGT, introgression, etc.), so depicting evolution as a bifurcating tree is not appropriate. In these cases, it is more useful to show multiple alternative relationships at once, or to indicate if there is evidence for introgression between branches.

__Objectives__

In today’s lab, you will:

1. Gain experience working with approaches to estimate species trees.
2. Learn how to conduct phylogenetic network analyses.

### Exercise 9.1 

### Gene trees and species trees 

Before jumping into this week’s lab on gene trees, let’s check on the concatenated analysis you started at the end of last week’s lab. Log on to your Jetstream VM and launch your Web Shell. Reattach to your screen:

```
screen -r iqtree
```

You should see some information on the screen, including an indication that your run completed by listing the date and time. 

:white_check_mark: __How long did it take for your analysis to run?__

You should also see several files beginning with “hemipteroid_concat.” These are your output files from IQ-Tree. Download the file ending in “treefile” and open it in FigTree. 

:white_check_mark: __Save the tree (don’t worry about rooting) as a PDF (LastName_hemipteroid_concat.pdf).__

__ASTRAL__

Now let’s move on to our species tree analysis. We’ll use a genomic-level data set of “flowers of the sea” tube-dwelling worms from [Tilic et al. (2020)](https://www.sciencedirect.com/science/article/abs/pii/S1055790320301640). 

Navigate to your volume or home directory and create a directory called `Lab9_gene_trees`. Next, navigate to the `Phylogenetics_Lab` directory and update the directory with `git pull`. You should see a new directory called Lab9_gene_trees_networks. Copy the content of this directory to your `Lab9_gene_trees` directory and extract it using `tar`.

First, let’s run a concatenated analysis. We’ll just run an unpartitioned analysis. I also suggest running this analysis in a screen so you can work on other things. If you’re not still attached to the iqtree screen, reattach with `screen -r iqtree`. Start an IQ-Tree analysis using the following command:

```
iqtree -s Tilic_concat.fasta -m WAG -bb 1000
```

Because these data are amino acids, we need to use an amino acid substitution model (the WAG model). After launching your IQ-Tree run, detach from your screen (`CTRL+a d`). The analysis will likely take too long to complete during lab, so I’ve provided a finished tree file in your data folder (“Tilic_concat.tre”). 

:white_check_mark: __Open this file in FigTree, root on _Prionospio_, show the bootstrap values, and save as LastName_tilic_concat.pdf.__

While the concatenated analysis is running, let’s move on to the species tree analysis. For this type of analysis, you would generate separate trees for each of your gene alignments, likely using a for loop, like this:

```
for gene in *.fasta;do
	iqtree -s $gene -m MFP -bb 1000;
done
```

This would test for the best substitution model for each gene alignment and then estimate a tree. In the end, you will have X number of tree files for X number of gene alignments. You would then combine all of those gene tree files into a single file like this:

```
cat *treefile >gene_trees.trees
```

You are now ready for running a species tree analysis. For this lab, we’ll skip those previous steps and just start with a bunch of gene trees, but you should be familiar with the code above if you’re doing a species tree analysis for your project. 

We will run a species tree analysis in the program [ASTRAL](https://github.com/smirarab/ASTRAL), which uses gene trees as input and summarizes them to produce a single species tree. This approach accounts for Incomplete Lineage Sorting (ILS) through the Multispecies Coalescent model. To run ASTRAL, use the following command:

```
astral -i Tilic_gene_trees.trees -o Tilic_species_tree.tre
```

You should see a file called “Tilic_species_tree.tre.” Open this file in FigTree, root on Prionospio, and show the branch labels (“label”). These are your local posterior probabilities, which is a branch support measure from ASTRAL based on the frequency of quartets. Like Bayesian posterior probabilities, these values are between 0 and 1. 

:white_check_mark: __Save the tree as a PDF (LastName_tilic_species_tree.pdf).__

:white_check_mark: __How do the concatenated and coalescent trees compare to one another?__

Pay particular attention to the branch lengths. In the ASTRAL species tree, internal branches are meaningful (in coalescent units), but the lengths of terminal branches are meaningless.

__SVDquartets__

Another coalescent approach is to use each site of an alignment to estimate a species tree. This type of analysis is particularly useful for SNP data. We will run this analysis using [SVDquartets](https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/), which is implemented in PAUP*.

Download the file fernandez.nex and open it in a text editor. These are SNP data from a clade of flowering plants in the genus _Linaria_, and come from [Fernandez-Mazuecos et al. (2018)](https://academic.oup.com/sysbio/article/67/2/250/3953673).

:white_check_mark: __Briefly describe the structure and content of this file.__

:white_check_mark: __How many SNPs are in this file?__

Next, launch PAUP* and load/execute the NEXUS file. 

Go to Analysis → SVDquartets

Set the number of quartets to exhaustive, select the bootstrap option, and select the bootstrapping option (default of 100 bootstrap replicated). 

:white_check_mark: __Save the resulting majority-rule consensus tree as LastName_svd.pdf by using the Print Tree option in View trees, or by writing the tree to a .tre file and opening/saving in FigTree (remember to show the bootstrap values!).__ 

There is no specific outgroup in this data set, so either root on the ALGI group or use midpoint rooting (FigTree: Tree → Midpoint Root).

Let’s also run an ML analysis on all the data. Run IQ-Tree with the following command:

```
iqtree -s fernandez_trim_variable.nex -m GTR+ASC -bb 1000
```

The “ASC” addition to the GTR model is accounting for Ascertainment Bias, which you need to run with SNP data to avoid overestimating branch lengths and biasing the tree.

:white_check_mark: __How do the SVDquartets and IQ-Tree trees compare?__

### Exercise 9.2

__Networks__

It could be that your group of interest does not evolve in a tree-like manner. For example, population-level divergence, hybridization, and horizontal gene transfer (HGT) can lead to an evolutionary history that is not truly bifurcating. In these cases, it is not ideal to represent the evolutionary history as a tree, but as a phylogenetic network. They key for networks is looking for reticulations rather than nodes of two daughter lineages.

__SplitsTree__

Let’s make a network using the program SplitsTree . Find the program icon on the desktop and click to open. You should see a screen like this:




Download the file `symbiont.nex` from Jetstream and open it in SplitsTree (File → Open). This is an alignment of several genes from various genera of bacteria, including free-living (such as _E. coli_) and endosymbionts of insects (including _Sodalis_ spp.). You should immediately see an image pop up on your screen. This is produced from a NeighborNet algorithm, which is based on Neighbor Joining, but shows conflicts/alternative branching patterns rather than a strictly bifurcating tree. These alternatives are shown as reticulations. 

:white_check_mark: __Save the network as LastName_symbiont_network.pdf (File → Export image…).__

:white_check_mark: __What are some basic observations about the network?__

:white_check_mark: __What do you think it means if certain relationships are highly reticulated (have many lines connecting them)?__

Next, create an NJ tree in SplitsTree by going to Trees → BioNJ. A window will pop up; hit “Apply.” 

:white_check_mark: __How does the NJ tree compare to the NeighborNet network?__

__PhyloNet__

Another approach to producing phylogenetic networks is to use gene trees for estimating introgression while simultaneously accounting for ILS. To do this, we will use the pseudo-maximum likelihood approach in the program PhyloNet .

Download the file yeast_trees.nex and open it in a text editor. You should notice this is a NEXUS file that contains several gene trees and commands in a “PhyloNet block” at the end of the file. The gene trees are from several species of yeast .

:white_check_mark: __What commands are in the PhyloNet block?__

Run PhyloNet from the Web Shell using the command:

```
phylonet yeast_trees.nex
```

When the program finishes running, you should see five inferred networks. Copy the Newick format after the first “Visualize in Dendroscope.” Open the program Dendroscope  from your desktop and paste in this Newick tree (CTRL+v should work). NOTE: IF YOU CANNOT FIND DENDROSCOPE, I HAVE PROVIDED A PNG FILE OF THE TREE IN YOUR DATA FOLDER. DOWNLOAD THIS FILE AND OPEN IT TO VIEW.

:white_check_mark: __Export/save the file as LastName_phylonet.pdf.__

:white_check_mark: __What do you notice about this phylogeny? What does the red line indicate?__
