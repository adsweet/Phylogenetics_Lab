# Laboratory 6
## R and Distance Methods

### Introduction

The philosophy of phenetics (creating classifications based on overall similarity) has fallen out of favor in systematics, but some of the phylogenetic reconstruction methods used involving distance or similarity have persisted. Distance methods have an advantage over parsimony methods in some molecular contexts in that they allow the incorporation of statistical models of sequence evolution. Distance methods are very fast and can be used to analyze extremely large data sets in a reasonable period of time. The disadvantage of distance methods is that, depending on the choice of method, they can be significantly less accurate than parsimony, likelihood, and Bayesian methods in estimating the correct phylogeny. Also, because these methods convert character data into distance data, they are not useful to study character evolution.

Today, distance methods are used in several contexts:

1.	Neighbor joining (NJ), BIONJ[^1] , and unweighted pair group method with arithmetic mean (UPGMA)[^2] may be applied as a shortcut to obtain a roughly accurate tree. For instance, NJ can be used to obtain a starting tree for parsimony or likelihood analyses, to speed up the search. This basic strategy is also employed in bioinformatics search tools or sequence alignment tools to obtain “guide trees,” because it is fast for large data sets.

[^1]: BIONJ is a more recent version of NJ that can be more accurate (Gascuel 1997)
[^2]: UPGMA always produces ultrametric trees (i.e., tip to root length is equal for all tips). It is the simplest distance method (Sokal and Michener 1958)

3.	Employing sequence evolution models and converting molecular distances between taxa into estimates of time since divergence (“molecular clock analyses”). Likelihood and Bayesian methods are more appropriately used for this purpose, although they are more computationally intensive. Distance methods that apply “best fit” to the data, such as Fitch-Margoliash, are used in this context. We will talk more about the estimation of divergence times in a later lab.

4.	For handling continuous types of data for which character-based methods are inapplicable. This can include some morphometric data, immunological cross-reactions, and frequencies (for instance allele frequencies from allozyme analysis).  While generally less accurate in estimating phylogenies than character-based methods, using distance methods may avoid the loss of information that results when continuous data are converted to discrete characters. 

5.	To construct quick phylogenies for exceptionally large datasets, such as analyses of genomic data for multiple taxa or handling of datasets with several hundred to thousands of taxa. Such datasets may be too computationally demanding to be analyzed with other methods, though many advances have been made on this front. 

6.	Neighbor joining phylogenies and their corresponding distance matrices are often used for analyses at the species-population interface, where changes are rare and most often occur without homoplasy (making it less necessary to compensate for hidden changes). Here, overall distances can be used as a first line of evidence regarding species boundaries.

**Objectives**

In this lab you will:
- Be introduced to the basics of R
- Be introduced to using R for phylogenetic analysis
- Become familiar with running neighbor joining and UPGMA distance methods

### Phylogenetics in R

First, let’s become more familiar with the basics R. If you have used R before, this section will probably be review for you, but it’s always helpful to practice the basics. R is an environment that is designed for statistical computing (https://www.r-project.org/). R is free for anyone to download and use, allows users to develop their own publicly available “packages” for more specific tasks, and there are many (many!) helpful online resources and tutorials. Because of this, R has quickly become a popular (if not the most popular) statistics program in many scientific fields, including phylogenetics. There are also many packages that allow you to produce high-quality figures. 

The only “downside” to R is that it relies on its own syntax (i.e., programming language) to run analyses. There is a little bit of a learning curve to master the syntax, but once you have the basics down, you’ll be well on your way to using R on a regular basis. Let’s begin by going over some of the basic syntax. We’ll be using RStudio (https://www.rstudio.com/), which is a nice GUI program for running R.

First, log on to Jetstream and access your course Instance. Next, launch the Web Desktop. You should see something like this:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/jetstream_desktop.png" width=50% height=50%>

Notice the “RStudio” link (blue circle with an “R”) at the bottom of your screen. Double-click that icon to launch RStudio. You should see a screen like this:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/rstudio.png" width=50% height=50%>

Your screen should have four distinct sections. If you don’t see the upper-left section, click the green plus symbol in the upper left or go to File -> New File -> R Script. Here are brief descriptions of these sections, moving from the upper left: 1) A text file for writing/saving your code. You can run a line of code by clicking the “Run” button above this screen or by hitting CTRL (Command on Mac) + Enter. 2) Your environment variables. This will populate as you run code. 3) Interactive window for viewing plots, help files, etc. 4) Your console/terminal window. This is where you can run your code, either by typing in directly or running from the text window above.

As an introduction to the basic syntax of R, we’ll be writing our code in the text window (upper left window). First, type in the following line and hit “Run”:
```
X <- 50
```
This command assigns the number “50” to the variable “X.” Notice the variable now shows up in your upper right window. Next, type and run this command:
```
Y <- “practice”
```
Notice that we put the word “practice” in quotes. When you’re working in R, if a particular word (or set of characters, i.e. “strings”) is not a built-in variable or a variable you have assigned, you need to put it between quotation marks. Typically, numbers do not need to be in quotes (but this is not universally true). Finally, type and run this command:
```
m <- mean(c(10,30,40))
```
This is calling a command (`mean`) to calculate the average of the numbers 10, 30, and 40. The `c` is needed when you want to run a command on a combined set. Now, go to the console (bottom left), type m, and hit Enter. You should see a number printed to your console screen.

Ok, that was a very (very) brief overview of the basic syntax of R, but let’s try to apply some of this syntax to phylogenetic analysis, specifically for distance methods. 

### _Exercise 6.1_

Another nice feature of RStudio is that you can also interact with your computer using a terminal window, exactly like you have in previous labs through the Web Shell. In the bottom left window, select the “Terminal” tab. You should see a layout that you are very familiar with from previous labs. You can use bash commands in this window. Change to your mounted volume or home directory and create a new directory called “Lab6.” Next, change to the `Phylogenetics_Lab` directory and update the data with:
```
git pull
```
You should see a new `Lab6` directory. Change into this directory and copy the file (“Lab6_data.tar.gz”) into your own `Lab6` directory. Move to your `Lab6` directory and extract the files with the command:
```
tar -xvzf Lab6_data.tar.gz
```
You should see two new files in your `Lab6` directory. Download the file “beetle_cox1.fasta” and open it in SeaView. These data are mitochondrial sequences (from the _cox1_ gene) from several species of beetles, including several sequences that were obtained from specimens that are thousands of years old (the oldest over 44,000 years old)! The specimens were preserved in packrat middens in California and Mexico (Smith et al. 2021). Do these sequences need to be aligned? If so, align it using MAFFT or MUSCLE. Name your aligned file “beetle_cox1_aligned.fasta”.

:white_check_mark: __Did you align the sequences? If so, what program and parameters did you use?__

Download your aligned file and view it in SeaView to ensure the sequence are adequately aligned.

Before running distance methods in R, let’s generate trees with some other methods as a comparison.

First, estimate a ML tree with IQ-Tree. Be sure to run a model test and at least 1,000 ultrafast bootstrap replicates. 

:white_check_mark: __Download the resulting tree, open it in FigTree, root on the branch with “Maset” and “Masap”, add the bootstrap values as branch labels, ladderize the tree, and save as a PDF (LastName_beetle_ML.pdf).__

Second, estimate an MP tree in PAUP*. Use SeaView to convert your alignment to a NEXUS file and execute it in PAUP*. Run an MP heuristic search with 100 bootstrap replicates. Be sure to set your outgroups. 

:white_check_mark: __Save your majority-rule consensus tree as LastName_beetles_MP.pdf.__

:white_check_mark: __What search parameters did you use for your MP search?__

Ok, at long last, we are ready to run our distance analyses in R.

To begin, we need to load a couple packages for our analyses. Packages are pre-loaded sets of commands created by users for specific purposes. There are many packages for phylogenetics, but we will only be using a few in this course. Today we will be working with the packages APE and Phangorn. APE is one of the most widely used phylogenetics packages and contains many useful tools. To load the packages, add these lines to your text window in RStudio and run both of them:
```
library(ape)
library(phangorn)
```
Next, we need to change to the directory that contains our data. In R, the command to change directory is `setwd()`. Run this command in your console (be sure you are switched back from the terminal!) and put your path within quotations between the parentheses. For example:
```
setwd(“~/asweet/Lab6”)
```
Next, read in your FASTA alignment file using one of the functions from APE:
```
beetle <- read.FASTA(“beetle_cox1_aligned.fasta”)
```
You can see that this `beetle` object is a “List” that is 13 elements long (the 13 taxa). It is in the data class “DNAbin,” which is used by APE and other packages to run functions on DNA data. 

We will first use APE to run a neighbor-joining (NJ) distance analysis. More specifically, we will be running a BIONJ analysis. To run neighbor-joining, we first need to compute distance matrices. These can then be summarized as neighbor-joining trees. To compute a distance matrix from a DNA object, we can use the dist.dna function. Let’s start with computing raw distances.
```
beetle_dist_raw <- dist.dna(beetle, model=”raw”)
```
You will see `beetle_dist_raw` object appear in the Environment tab. In addition, it is a 'dist' object class (a distance matrix). To display the results, you can simply type out the object name in the console and hit `Return`. Because the matrix is pretty wide, it will wrap around in your console when displayed. What you will see are the pairwise distances between every specimen.
```
beetle_dist_raw
```
Now, let’s run the neighbor joining analysis! This will complete very quickly.
```
beetle_nj_raw <- bionj(beetle.dist.raw)
```
Let’s visualize the results:
```
plot(beetle.dist.raw)
```
You might see something is wrong: we did not root the tree. We could export this tree as a Newick (using the `write.tree()` command to write a Newick format file), open it FigTree, and then reroot, but we don’t have to do that since we can root it using APE. We use the outgroup argument to specify our outgroup (we’ll just use “Maset”), and also set `resolve.root = TRUE`. The default is to leave the root in an unresolved trichotomy with the ingroup. We’ll go ahead and ladderize it as well to order the nodes in displaying the tree. We will also add the scale bar so we have a scale bar for the distances on the branch lengths.
```
beetle_nj_raw <- root(beetle_nj_raw, outgroup = "Maset", resolve.root = TRUE)
beetle_nj_raw <- ladderize(beetle_nj_raw, right=F)
plot(beetle_nj_raw)
add.scale.bar(ask=T)
```
You will need to click on the figure to place your scale bar.
Now, let’s save this tree. Find the “Export” button within the Plots tab. Save this tree as a PDF. 

:white_check_mark: __Change the file name to LastName_NJ.pdf and save the file.__

Preview the resulting PDF file to make sure the export was successful. 

:white_check_mark: __How does this tree compare to the ML and parsimony trees in phylogenetic relationships and in general features?__

:white_check_mark: __Are there any polytomies in the Neighbor Joining tree?__

Models can also be applied to attempt to improve your distance analysis under a more realistic model of substitution. Distances calculated using models of sequence evolution are called corrected distances. Various models are implemented in the `dist.dna()` function. Let's run the default model, the Kimura-2P (or K80) model, which distinguishes transitions from transversions, and thus changes the branch lengths to represent this new substitution rate matrix. Models correct for potential hidden changes (e.g., from a C to a T and back to a C again), resulting in greater inferred distances for areas deeper in the tree. Correcting for these molecular rates is important when using distances for inferring divergence times and can impact relative branch lengths and topology. 

Run the neighbor joining analysis again, this time specifying "K80" for the model. Because it is the default model, it would be run if we did not specify a model, but let’s specify it anyway to be explicit.
```
beetle_dist_k2p <- dist.dna(beetle, model = "K80")
beetle_nj_k2p <- bionj(beetle_dist_k2p)
beetle_nj_k2p <- root(beetle_nj_k2p, outgroup = "Maset", resolve.root = TRUE)
beetle_nj_k2p <- ladderize(beetle_nj_k2p, right=F)
plot(beetle_nj_k2p)
add.scale.bar(ask=T)
```
:white_check_mark: __Again, save the tree (save as LastName_NJK2P.pdf).__

:white_check_mark: __How does this tree compare to your original neighbor joining tree? (Hint: Look at the scale bar.)__

The beetle alignment we just ran had relatively few taxa, but one advantage of distance methods is their speed. You won’t notice the difference on smaller files, but it can be quite noticeable with larger files. To demonstrate this, let’s work with the file “bee_csd.fasta,” which is an alignment of honey bee sequences from a portion of the CSD gene. Download the file and open it in SeaView.

:white_check_mark: __How many sequences are in this file?__

Convert the file to a NEXUS, execute it in PAUP*, and run a heuristic search with 100 boostrap replicates. 

:white_check_mark: __How long did the MP analysis take (or did it even finish)?__

Next, run a Maximum Likelihood analysis on the bee alignment in IQ-TREE with 100 ultrafast bootstrap replicates.

:white_check_mark: __How long did the ML analysis take?__

Now, let’s run the same data through a NJ analysis in R:
```
bee <- read.FASTA(“bee_csd.fasta”)
bee_dist <- dist.dna(bee)
bee_nj <- bionj(bee_dist)
```
:white_check_mark: __How long did the NJ analysis take with the bee data?__

Finally, let’s estimate a UPGMA tree (which is often used for making guide trees in alignment algorithms) on the beetle data.

First, load the R package phangorn:
```
library(phangorn)
```
:white_check_mark: __Root your tree on “Maset,” ladderize, plot, and save as a PDF (LastName_beetles_UPGMA.pdf).__

:white_check_mark: __How does this tree compare to the MP and ML trees?__

:white_check_mark: __Describe the main differences between the NJ and UPGMA trees.__
 
References

Gascuel, O. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular Biology and Evolution. 14: 685-695.  

Paradis, E., Claude, J., and Strimmer, K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics. 20: 289-290.
Schliep, K. 2011. Phangorn: phylogenetic analysis in R. Bioinformatics. 27: 592-593.  

Sokal, R. and Michener, C.B. 1958. A statistical method for evaluating systematic relationships. University of Kansas Science Bulletin. 38: 1409-1438.
