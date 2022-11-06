# Laboratory 5
## Substitution Models and Maximum Likelihood Methods

### Introduction

 
R.A. Fisher

Maximum likelihood was first applied to phylogeny estimation by Edwards and Cavalli-Sforza (1964), both students of the great statistician Ronald A. Fisher, who invented likelihood as a statistical method (Fisher 1912, 1921, 1922).  It was not until the work of Joe Felsenstein (Felsenstein 1981), which capitalized on earlier attempts by Jerzy Neyman (1971) and others, that a practical computerized method was found for calculating tree likelihoods from sequences. In the context of molecular phylogenetics, maximum likelihood attempts to find the tree that maximizes the likelihood that the sequence data we observe would have been produced, given a particular evolutionary model that specifies how the sequences are expected to evolve (“substitution models”). The use of an evolutionary model allows the explanation of a dataset in association with a tree to be more complete and precise (thus the trees that maximize the fit should in principle be more accurate).  Hence, model selection is important in maximum likelihood, and model misspecification can lead to an incorrect tree. Furthermore, the model may differ among different subsets of columns, or partitions. Testing whether or not partitions are evolving under similar models is another aspect of selecting the best model. With regard to selecting a realistic model given those that have been developed to date, objective methods using statistical criteria have been introduced. Distance and Bayesian phylogenetic methods also rely on models of sequence evolution, so future labs will continue to focus on model-based inference.

One program in particular, IQ-TREE  (Nguyen et al. 2015), was developed to overcome the hurdles associated with large datasets. While model selection, partition finding, and maximum likelihood estimation were previously implemented in separate software, IQ-TREE advances on previous methods by combining all into a single piece of software. 

**Objectives**

In today’s lab, we will explore determining the best fit model both with and without partitions, maximum likelihood estimation, and ultrafast bootstrapping.

•	To understand the basic principles of maximum likelihood (ML) inference and how ML analyses are conducted in IQ-TREE.

•	To become familiar with the ModelFinder program within IQ-TREE and to know how it selects from among a large number of possible models to find the models and partitioning scheme that best fits your data.  

•	To gain some insights into how results may vary across different methods of phylogenetic analysis.


#### _Laboratory Exercise 1: Testing for the best model_

Because we are using molecular data (DNA), we first have to consider which model will be appropriate to use. IQ-TREE determines the best-fit model using the ModelFinder program, which is implemented within IQ-TREE (Kalyaanamoorthy et al., 2017). Models can be specified to IQ-TREE using the -m flag, but you may also specify -m to test for the best model(s) using ModelFinder. For DNA substitution model, IQ-TREE implements all alternative time-reversible substitution models ranging from the simplest Jukes-Cantor model (JC), where all substitutions have the same rate, to the most complex General Time Reversible (GTR) model, where all base substitutions have different rates. Many other models are available, which you can read more about in the documentation: http://www.iqtree.org/doc/Substitution-Models.

Reversible nucleotide substitution models have three main components: 

1. The substitution rate matrix that determines categories of rate substitutions (e.g. single rate, multiple rates such as transitions vs. transversions, or all different rates)
2. The base frequencies (e.g. empirical base frequencies, equal base frequencies, or estimated base frequencies)
3. Rate heterogeneity among sites, indicated by either a gamma parameter that specifies the distribution of rates (+G (gamma) or +R (FreeRate) plus the number of rate categories), or invariant sites (+I) that specifies the proportion of sites that are not changing at all. ModelFinder also reports the type of frequencies (+F is the default for models with unequal frequencies, +FQ for equal frequencies).

To do determine the best model, ModelFinder evaluates the models that best match your data while adjusting for over-fitting to your data. To understand this, consider a bunch of data points in a plot of two variables (right). A line that perfectly explains the data would connect each dot to each other. Yet connecting the dots has no predictive value as to the general trend from which these points are naturally varying, usually represented by a line or curve. ModelFinder will effectively penalize overfitting the model to your data to allow the model to be a better predictor of the actual evolutionary trend. The program does this by testing using two different criteria, the Akaike Information Criteria (AIC) and the Bayesian Information Criteria (BIC) on many different models of sequence evolution. The Akaike information criterion compares the likelihoods of different models and applies a penalty to models that make more assumptions.  The Bayesian information criterion tends to penalize more for overparameterization. 

Now let’s practice searching for optimal models using IQ-Tree.

Logon to your Jetstream account and activate the Web Shell. Create a “Lab5” directory in your mounted volume. Next, navigate to the `A-State_Phylogenetics_Systematics_Data` directory and update the contents with git pull. Copy the content to your “Lab5” directory, change into the “Lab5” directory, and un-tar the file. You should see several files beginning with “Bombus.”

To only run ModelFinder in IQ-TREE, we will call IQ-TREE and specify the “Bombus.nex” file as the sequence file using the `-s` flag. This file contains mitochondrial 16S sequences for 15 species of bumble bees (genus _Bombus_) and one outgroup taxon, a stingless bee (Meliponini). Download the sequence file and open it in SeaView.

**Do the sequences look like they are aligned?

To tell IQ-TREE to not infer the phylogeny and only find the best model, we will specify the model -m as “MF”. By default, IQ-TREE will determine the type of data (e.g. DNA vs. amino acid), but you can specify this as well. We are going to specify to use only one thread (-nt):
```
iqtree -nt 1 -s Bombus.nex -m MF
```
As the program is running, you will see output in the main window like the following:

```
ModelFinder will test 286 DNA models (sample size: 556) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            2914.221     29  5886.443     5889.751     6011.745
  2  JC+I          2748.219     30  5556.439     5559.982     5686.062
  3  JC+G4         2729.101     30  5518.202     5521.744     5647.825
  4  JC+I+G4       2725.585     31  5513.169     5516.956     5647.113
  5  JC+R2         2724.508     31  5511.016     5514.802     5644.960
 ```

The models are listed in order of increasing complexity. Model 1 of 286 is the Jukes-Cantor (1969) model, which assumes that the nucleotides appear in equal frequencies and that transitions and transversions are equally likely. The models get more complex, making more assumptions, as the tests run towards the final model. As the models become more complex, they also include more parameters, including estimated proportion of invariable sites, transition/transversion ratios, gamma shape parameters, and estimated base frequencies. 

The Model column shows the names of the different models. The -LnL shows the negative log-likelihood of the models. Generally, a lower negative log-likelihood (which is actually a “higher” likelihood, since these are already NEGATIVE log-likelihood values) is preferred. df represents the degrees of freedom of the model. You will see the program also estimates and reports AIC, AICc, and BIC scores. AICc corrects for small sample sizes. As you can see, different models are not necessarily the best fit based on different information criteria. By default, IQ-TREE uses BIC to select the best model, but this can be changed by adding the -AIC or -AICc flag.

By default, IQ-TREE will also create several files, which it notes at the bottom of the standard out, an IQ-TREE report (.iqtree), a treefile in Newick format which is used in fitting the model (IQ-TREE uses a quick parsimony-based method for this), and the log file (which is a file that contains the output that was printed to screen).

Download the “Bombus.nex.iqtree” files and open it in a text editor. 

How many parsimony informative sites are in the alignment?

**What is the top model based on the BIC scores?

Now download the file “Bombus.nex.model.gz,” open it in a text editor, and scroll to the bottom.

**Do AIC and AICc have the same best model as BIC? If not, which model(s) do those metrics favor?

#### _Laboratory Exercise 2: ML Tree Inference

Now that you’ve done some model testing, let’s move to the next step: estimating phylogenetic trees with the ML approach. We will also be doing this step in IQ-TREE.

Before running an ML analysis, we should also generate an MP tree for comparions. Launch PAUP*, load in the Bombus.nex file, and run a bootstrapped parsimony analysis using the most appropriate parameters. Obtain a consensus of the resulting trees and save the consensus tree as LastName_parsimony.pdf.

**What conditions did you use for your parsimony analysis?

Next, let’s run an ML analysis on the same Bombus.nex file. We’re going to run both standard and ultrafast bootstrapping (Hoang et al. 2018). Ultrafast bootstrapping is a different method for assessing node support that is similar to bootstrapping, but is more computationally efficient for large, phylogemic datasets. Practically, the most important thing to note is that the range of ultrafast bootstraps differ, such that values above 95 are considered good support, compared to values above 70 for traditional bootstrapping. We are also setting the prefix this time for all of the output files. 

Traditional bootstrap:
```
iqtree -s Bombus.nex -m MFP -b 100 -pre Bombus_tree_bs
```
Ultrafast bootstrapping:
```
iqtree -s Bombus.nex -m MFP -bb 1000 -pre Bombus_tree_uf
```
It should only take a few minutes for the analyses to run. After the best-fit model is selected, IQ-TREE will generate numerous starting trees, pick the best 20, and perform branch-swapping heuristics (Nearest Neighbor Interchange - NNI) to find the maximum likelihood tree . It will compute the likelihood of trees as it searches through both tree space and parameter space as it tries to optimize the tree. Afterwards, it ends the tree search and optimizes the best-fit model on the best tree, reporting the likelihood of the best tree. Finally, it outputs all of the stats on the run and the analysis files.

IQ-TREE will also write the final trees to Bombus_tree_bs.treefile and Bombus_tree_uf.treefile. Download these files. We will open the files in the program FigTree, a GUI program for viewing and manipulating phylogenetic trees. Search for FigTree on your computer and launch the program. You should see a window like this: 

















Open your tree files from IQ-TREE with File --> Open. First, open Bombus_tree_bs.treefile. When prompted, type in "bootstrap" to label the values. This tree is technically unrooted (we could have specified an outgroup using the -o flag). Now, on the left hand side, check the box to the left of “Branch Labels” and click on the triangle. Select "bootstrap" to be displayed on the branch labels. You should now see the bootstrap values on the tree. Next, we want to root the tree. Recall that the outgroup is Meliponini. Click on the branch leading to Meliponini and click “Reroot” on the top menu. You can also organize the tree by using either Tree > Increasing Node Order or Decreasing Node Order. This rotates the nodes by node order, which can enhance legibility of the tree. I usually prefer showing branches in increasing order. Another nice feature of FigTree is the flexibility to customize the appearance of your tree. Play around with the options in “Appearance” and font size of Tip Labels and Branch Labels to maximize the clarity of your tree.

Follow the same steps to open your Bombus_tree_uf.treefile in FigTree. 

**Save both files (LastName_Lab5_iqtree_bs.pdf and LastName_Lab5_iqtree_uf.pdf). 

**How do these trees compare to one another?

**How do they compare to the parsimony tree?
	
#### _Laboratory Exercise 3: Partitioned Analysis

Next, let’s try running a partitioned analysis in IQ-TREE (If you use a partitioned analysis, cite Chernomor et al. 2016). Partitioning allows you to use different models for different subsets of your data. In practice, phylogeneticists will often concatenate their sequence alignments (combine all alignments into one big alignment, i.e. supermatrix) and partition this supermatrix by gene, codon position, or some other criteria. You can also partition individual gene alignments, but we will focus on working with concatenated data for this lab.

For this exercise, you will use alignments from two different genes in Bombus: EF-1α (“Bombus_ef1a.phy”) and PEPCK (“Boumbus_pepck.phy”). First, let’s concatenate these alignments. We will use FASconCAT-G for this purpose. This program will look for alignment files in your working directory (based on the filename extension) and concatenate them.  

Make a new directory (called “concat” or something similar) and copy the two Bombus gene files (“Bombus_ef1a.phy” and “Boumbus_pepck.phy”) into the new directory. Now change into your new directory. To launch FASconCAT-G, run the following command:

```
fasconcat-g -s -l
```

When the program finishes running, list your files. You should see several new files beginning with “FcC.” Rename these files:

```
mv FcC_supermatrix.fas Bombus_concat.fasta
mv FcC_supermatrix_partition.txt Bombus_concat_partitions.txt
```

Download your newly renamed files. Open the FASTA file in Seaview; it should be the combined length of your two gene alignments. Open the text file in a text editor; you should see two lines, both beginning with “DNA,” showing each partition name (in this case the two genes), and the range of sites for each partition.

Next, we will run IQ-TREE to estimate the best models and a ML tree. In this run, we will use the option “-p,” which will be our partition file. This tells IQ-TREE to 1) determine what the best partitioning scheme should be (Separate models? A single model for both partitions?), tests for the best model(s), and estimates the ML tree under these models. Run the analysis with following command:
```
iqtree -s Bombus_concat.fasta -spp Bombus_concat_partitions.txt -m MFP+MERGE -bb 1000 -pre Bombus_tree_concat
```

Download and open the *.iqtree file. 

**How many partitions are optimal? What was/were the best model(s)? Were there difference among AIC, AICc, and BIC (look in your *.model file)?

Next, download the tree file and open it in Figtree. Root the tree on Meliponini, ladderize, and show the bootstrap values as branch labels. Save the tree as LastName_Bombus_concat.pdf.

**How does this tree compare to your single-gene (16S) tree?

Keep in mind that it is a good idea to explore the impact of model choice on the results of your own analysis. For example, for your own project data you may wish to see if results differ between the AIC model and BIC model if these criteria differ in model choice, partitioning by codon position for protein-coding genes, etc.

As a final comparison, let’s run your concatenated alignment in the program RAxML, which is another commonly-used ML software. However, RAxML has fewer models compared to IQ-TREE and does not run a model search. We’ll just use a single partition. Run the following command:
```
raxml -s Bombus_concat.fasta -n Bombus_concat_raxml -m GTRGAMMAI -f a -x 1234 -p 4321 -# 200
```
The “-m” option specifies the model, in this case GTR + I + Γ. The “-f a” option specifies a rapid bootstrap test (“-#” is the number of bootstrap replicates), and the “-x” and “-p” are starting seeds (to run the exact same analysis again these numbers would need to be the same in each run, but usually you can put “random” numbers of your choosing. 

After the run finishes, your should see several files beginning with “RAxML.” Download file “RAxML_bipartitions.Bombus_concat_raxml.” This is your tree file, so go ahead and open it Figtree. Make sure to root, ladderize, and show bootstrap values on your tree. Save the tree as LastName_Bombus_concat_raxml.pdf. 

**How does this tree compare to the tree from IQ-TREE? 

**Based on your comparison, do you think partitioning and model testing is a crucial component of phylogenetic inference? If not, are there circumstances where it would be important? If yes, are there instances when it might not be as important?
