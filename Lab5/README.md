# Laboratory 5
## Substitution Models and Maximum Likelihood Methods

### Introduction

Maximum likelihood was first applied to phylogeny estimation by Edwards and Cavalli-Sforza (1964), both students of the great statistician Ronald A. Fisher, who invented likelihood as a statistical method (Fisher 1912, 1921, 1922).  It was not until the work of Joe Felsenstein (Felsenstein 1981), which capitalized on earlier attempts by Jerzy Neyman (1971) and others, that a practical computerized method was found for calculating tree likelihoods from sequences. In the context of molecular phylogenetics, maximum likelihood attempts to find the tree that maximizes the likelihood that the sequence data we observe would have been produced, given a particular evolutionary model that specifies how the sequences are expected to evolve (“substitution models”). The use of an evolutionary model allows the explanation of a dataset in association with a tree to be more complete and precise (thus the trees that maximize the fit should in principle be more accurate).  Hence, model selection is important in maximum likelihood, and model misspecification can lead to an incorrect tree. Furthermore, the model may differ among different subsets of columns, or partitions. Testing whether or not partitions are evolving under similar models is another aspect of selecting the best model. With regard to selecting a realistic model given those that have been developed to date, objective methods using statistical criteria have been introduced. Distance and Bayesian phylogenetic methods also rely on models of sequence evolution, so future labs will continue to focus on model-based inference.

One program in particular, [IQ-TREE](http://www.iqtree.org/), was developed to overcome the hurdles associated with large datasets. While model selection, partition finding, and maximum likelihood estimation were previously implemented in separate software, IQ-TREE advances on previous methods by combining all into a single piece of software. 

In today’s lab, we will explore determining the best fit model both with and without partitions, maximum likelihood estimation, and ultrafast bootstrapping.

**Objectives**
- To understand the basic principles of maximum likelihood (ML) inference and how ML analyses are conducted in IQ-TREE.
- To become familiar with the ModelFinder program within IQ-TREE and to know how it selects from among a large number of possible models to find the models and partitioning scheme that best fits your data.  
- To gain some insights into how results may vary across different methods of phylogenetic analysis.


### _Exercise 5.1: Testing for the best model_

Because we are using molecular data (DNA), we first have to consider which model will be appropriate to use. IQ-TREE determines the best-fit model using the ModelFinder program, which is implemented within IQ-TREE (Kalyaanamoorthy et al., 2017). Models can be specified to IQ-TREE using the -m flag, but you may also specify -m to test for the best model(s) using ModelFinder. For DNA substitution model, IQ-TREE implements all alternative time-reversible substitution models ranging from the simplest Jukes-Cantor model (JC), where all substitutions have the same rate, to the most complex General Time Reversible (GTR) model, where all base substitutions have different rates. Many other models are available, which you can read more about in the documentation: http://www.iqtree.org/doc/Substitution-Models.

Reversible nucleotide substitution models have three main components: 

1. The substitution rate matrix that determines categories of rate substitutions (e.g. single rate, multiple rates such as transitions vs. transversions, or all different rates)
2. The base frequencies (e.g. empirical base frequencies, equal base frequencies, or estimated base frequencies)
3. Rate heterogeneity among sites, indicated by either a gamma parameter that specifies the distribution of rates (+G (gamma) or +R (FreeRate) plus the number of rate categories), or invariant sites (+I) that specifies the proportion of sites that are not changing at all. ModelFinder also reports the type of frequencies (+F is the default for models with unequal frequencies, +FQ for equal frequencies).

To determine the best model, ModelFinder evaluates the models that best match your data while adjusting for over-fitting to your data. To understand this, consider a bunch of data points in a plot of two variables. A line that perfectly explains the data would connect each dot to each other. Yet connecting the dots has no predictive value as to the general trend from which these points are naturally varying, usually represented by a line or curve. ModelFinder will effectively penalize overfitting the model to your data to allow the model to be a better predictor of the actual evolutionary trend. The program does this by testing using two different criteria, the Akaike Information Criteria (AIC) and the Bayesian Information Criteria (BIC) on many different models of sequence evolution. The Akaike information criterion compares the likelihoods of different models and applies a penalty to models that make more assumptions.  The Bayesian information criterion tends to penalize more for overparameterization. 

Now let’s practice searching for optimal models using IQ-Tree.

Logon to your Jetstream account and activate the Web Shell. Create a “Lab5” directory in your mounted volume or in your home directory. Next, navigate to the `Phylogenetics_Lab` directory and update the contents with `git pull`. Copy the content to your `Lab5` directory, change into the `Lab5` directory, and  extract the tar file. You should see several files beginning with “Bombus.”

To only run ModelFinder in IQ-TREE, we will call IQ-TREE and specify the “Bombus.nex” file as the sequence file using the `-s` flag. This file contains mitochondrial _16S_ sequences for 15 species of bumble bees (genus _Bombus_) and one outgroup taxon, a stingless bee (from the tribe Meliponini). Download the sequence file and open it in SeaView.

:white_check_mark: __Do the sequences look like they are aligned?__

To tell IQ-TREE not to infer the phylogeny and only find the best model, we will specify the model ```-m``` as ```MF```. By default, IQ-TREE will determine the type of data (e.g. DNA vs. amino acid), but you can specify this as well.
```
iqtree -s Bombus.nex -m MF
```
As the program is running, you will see output in the main window like the following:

```
ModelFinder will test up to 286 DNA models (sample size: 556) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  GTR+F         2643.558     37  5361.116     5366.545     5520.984
  2  GTR+F+I       2459.816     38  4995.631     5001.364     5159.820
  3  GTR+F+G4      2422.722     38  4921.444     4927.177     5085.633
  4  GTR+F+I+G4    2417.280     39  4912.560     4918.606     5081.070
  5  GTR+F+R2      2416.197     39  4910.394     4916.441     5078.904
 ```

The models are grouped in order of decreasing complexity. The last models listed are Jukes-Cantor (JC) models, which assumes that the nucleotides appear in equal frequencies and that transitions and transversions are equally likely. The models get more complex, making more assumptions, as the tests run towards the final model. As the models become more complex, they also include more parameters, including estimated proportion of invariable sites, transition/transversion ratios, gamma shape parameters, and estimated base frequencies. 

The Model column shows the names of the different models. The -LnL shows the negative log-likelihood of the models. Generally, a lower negative log-likelihood (which is actually a “higher” likelihood, since these are already NEGATIVE log-likelihood values) is preferred. df represents the degrees of freedom of the model. You will see the program also estimates and reports AIC, AICc[^1], and BIC scores. AICc corrects for small sample sizes. As you can see, different models are not necessarily the best fit based on different information criteria. By default, IQ-TREE uses BIC to select the best model, but this can be changed by adding the `-AIC` or `-AICc` flag.

[^1]: AICc is AIC with a correction for small sample sizes. See [Burnham and Anderson (2004)](https://www.uvm.edu/~bmitchel/NR385/Burnham_Multimodel_inference.pdf) for recommendations on when to use AIC vs. AICc.

By default, IQ-TREE will also create several files, which it notes at the bottom of the standard out, an IQ-TREE report (.iqtree), a treefile in Newick format which is used in fitting the model (IQ-TREE uses a quick parsimony-based method for this), and the log file (which is a file that contains the output that was printed to screen).

Download the “Bombus.nex.iqtree” file and open it in a text editor. 

:white_check_mark: __How many parsimony informative sites are in the alignment?__

:white_check_mark: __What is the top model based on the BIC scores?__

Now download the file “Bombus.nex.model.gz,” open it in a text editor, and scroll to the bottom.

:white_check_mark: __Do AIC and AICc have the same best model as BIC? If not, which model(s) do those metrics favor?__

### _Exercise 5.2: ML Tree Inference_

Now that you’ve done some model testing, let’s move to the next step: estimating phylogenetic trees with the ML approach. We will also be doing this step in IQ-TREE.

Before running an ML analysis, we should also generate an MP tree for comparions. Launch PAUP*, load in the Bombus.nex file, and run a bootstrapped parsimony analysis using the most appropriate parameters. 

:white_check_mark: __Obtain a consensus of the resulting trees and save the consensus tree as LastName_parsimony.pdf.__

:white_check_mark: __What conditions did you use for your parsimony analysis?__

Next, let’s run an ML analysis on the same Bombus.nex file. We’re going to run both standard and ultrafast bootstrapping (Hoang et al. 2018). Ultrafast bootstrapping is a different method for assessing node support that is similar to bootstrapping, but is more computationally efficient for large, phylogemic datasets. Practically, the most important thing to note is that the range of ultrafast bootstraps differ, such that values above 95 are considered good support, compared to values above 70 for traditional bootstrapping. We will also specify the ```MFP``` (ModelFinder Plus) for the ```-m``` parameter, which tells IQ-TREE to first search for the best model before running a tree search. We are also setting the prefix this time for all of the output files.[^2]

[^2]: Note that if you try to rerun without specifying the prefix, IQ-TREE will first attempt to use the default (the alignment file name, or the partition file name if a partition file is given). However, since we’ve already used the default filename, IQ-TREE will alert you that the run has already completed and won’t start a new analysis. If you wanted to simply rerun an analysis, you can get overwrite previous files by specifying the `-redo` argument.

Traditional bootstrap:
```
iqtree -s Bombus.nex -m MFP -b 100 -pre Bombus_tree_bs
```
Ultrafast bootstrapping:
```
iqtree -s Bombus.nex -m MFP -B 1000 -pre Bombus_tree_uf
```
It should only take a few minutes for the analyses to run. After the best-fit model is selected, IQ-TREE will generate numerous starting trees, pick the best 20, and perform branch-swapping heuristics (Nearest Neighbor Interchange - NNI) to find the maximum likelihood tree.[^3] It will compute the likelihood of trees as it searches through both tree space and parameter space as it tries to optimize the tree. Afterwards, it ends the tree search and optimizes the best-fit model on the best tree, reporting the likelihood of the best tree. Finally, it outputs all of the stats on the run and the analysis files.

[^3]: You could also specify the best model from your previous ModelFinder analysis. To do this, you would specify `-m <model>` instead of `-m MFP`.

IQ-TREE will also write the final trees to Bombus_tree_bs.treefile and Bombus_tree_uf.treefile. Download these files. We will open the files in the program [Figtree](https://github.com/rambaut/figtree/releases), a GUI program for viewing and manipulating phylogenetic trees. Search for Figtree on your computer and launch the program. You should see a window like this: 

![Figtree image](../Images/figtree_screen.png)

Open your tree files from IQ-TREE with File --> Open. First, open the Bombus_tree_bs.treefile. When prompted, type in "bootstrap" to label the values. This tree is technically unrooted (we could have specified an outgroup using the `-o` flag in IQ-TREE). Now, on the left hand side, check the box to the left of “Branch Labels” and click on the triangle. Select "bootstrap" to be displayed on the branch labels. You should now see the bootstrap values on the tree. Next, we want to root the tree. Recall that the outgroup is Meliponini. Click on the branch leading to Meliponini and click “Reroot” on the top menu. You can also organize the tree by using either Tree > Increasing Node Order or Decreasing Node Order. This rotates the nodes by node order, which can enhance legibility of the tree. I usually prefer showing branches in increasing order. Another nice feature of FigTree is the flexibility to customize the appearance of your tree. Play around with the options in “Appearance” and font size of Tip Labels and Branch Labels to maximize the clarity of your tree.

Follow the same steps to open your Bombus_tree_uf.treefile in FigTree. 

:white_check_mark: __Save both files (LastName_Lab5_iqtree_bs.pdf and LastName_Lab5_iqtree_uf.pdf)__

:white_check_mark: __How do these trees compare to one another?__

:white_check_mark: __How do they compare to the parsimony tree?__
	
### _Exercise 5.3: Partitioned Analysis_

Next, let’s try running a partitioned analysis in IQ-TREE. Partitioning allows you to use different models for different subsets of your data. In practice, phylogeneticists will often concatenate their sequence alignments (combine all alignments into one big alignment, i.e. supermatrix) and partition this supermatrix by gene, codon position, or some other criteria. You can also partition individual gene alignments, but we will focus on working with concatenated data for this lab.

For this exercise, you will use alignments from two different genes in Bombus: EF-1α (“Bombus_ef1a.phy”) and PEPCK (“Boumbus_pepck.phy”). First, let’s concatenate these alignments. We will use [AMAS](https://github.com/marekborowiec/AMAS) for this purpose. The AMAS script will look for alignment files in your working directory based on a filename extension and concatenate them.  

Make a new directory (called “concat” or something similar) and copy the two Bombus gene files (“Bombus_ef1a.phy” and “Boumbus_pepck.phy”) into the new directory. Now change into your new directory. To launch AMAS, run the following command:

```
amas concat -f phylip -i *.phy -d dna -u fasta --part-format raxml
```
This command takes a PHYLIP file of DNA sequences as input, looks for any files ending in .phy, and outputs a concatenated file in FASTA format with a partition file formated for the program RAxML (more on this program in a minute).

When the program finishes running, list your files. You should see two new files: "concatenated.out" and "partitions.txt".

Rename these files by running:
```
mv concatenated.out Bombus_concat.fasta
mv partitions.txt Bombus_concat_partitions.txt
```

Download these files. Open the "concatenated.out" file in Seaview; it should be the combined length of your two gene alignments. Open the "partitions.txt" file in a text editor; you should see two lines, both beginning with “DNA,” showing each partition name (in this case the two genes), and the range of sites for each partition.

Next, we will run IQ-TREE to estimate the best models and a ML tree. In this run, we will use the option `-p`, which will be our partition file. This tells IQ-TREE to determine what the best partitioning scheme should be (Separate models? A single model for both partitions?) while allowing each partition to have its own rate of evolution, tests for the best model(s), and estimates the ML tree under these models. Run the analysis with following command:
```
iqtree -s Bombus_concat.fasta -p Bombus_concat_partitions.txt -m MFP+MERGE -B 1000 -pre Bombus_tree_concat
```

Download and open the *.iqtree and *.model files. 

:white_check_mark: __How many partitions are optimal? What were the best model(s)? Were there differences among AIC, AICc, and BIC (look in your *.model file)?__

Next, download the tree file and open it in Figtree. Root the tree on Meliponini, ladderize, and show the bootstrap values as branch labels. 

:white_check_mark: __Save the tree as LastName_Bombus_concat.pdf__

:white_check_mark: __How does this tree compare to your single-gene (16S) tree?__

Keep in mind that it is a good idea to explore the impact of model choice on the results of your own analysis. For example, for your own project data you may wish to see if results differ between the AIC model and BIC model if these criteria differ in model choice, partitioning by codon position for protein-coding genes, etc.

As a final comparison, let’s run your concatenated alignment in the program [RAxML](https://github.com/stamatak/standard-RAxML), which is another commonly-used ML software. However, RAxML has fewer models compared to IQ-TREE and does not run a model search. We’ll just use a single partition. Run the following command:
```
raxml -s Bombus_concat.fasta -n Bombus_concat_raxml -m GTRGAMMAI -f a -x 1234 -p 4321 -# 200
```
The “-m” option specifies the model, in this case GTR + I + Γ. The “-f a” option specifies a rapid bootstrap test (“-#” is the number of bootstrap replicates), and the “-x” and “-p” are starting seeds (to run the exact same analysis again these numbers would need to be the same in each run, but usually you can put “random” numbers of your choosing). 

After the run finishes, your should see several files beginning with “RAxML.” Download file “RAxML_bipartitions.Bombus_concat_raxml.” This is your tree file, so go ahead and open it Figtree. Make sure to root, ladderize, and show bootstrap values on your tree. 

:white_check_mark: __Save the tree as LastName_Bombus_concat_raxml.pdf.__ 

:white_check_mark: __How does this tree compare to the tree from IQ-TREE?__

:white_check_mark: __Based on your comparison, do you think partitioning and model testing is a crucial component of phylogenetic inference? If not, are there circumstances where it would be important? If yes, are there instances when it might not be as important?__

### References

Akaike, H.  1974.  A new look at the statistical model identification. IEEE  Trans. Autom. Contr. 19:716-723. 

Chernomor, O., von Haeseler A., Minh B.Q. 2016. Terrace aware data structure for phylogenomic inference from supermatrices. Systematic Biology 65:997-1008.

Edwards, A. W. F. and L. L. Cavalli-Sforza. 1964.  Reconstruction of evolutionary trees.  Pp. 67-76 in: Phenetic and Phylogenetic Classification, V. H. Heywood and J. McNeill, eds., Systematics Association Publ. 6, London,

Felsenstein, J. 1981.  Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution 17: 368-376.

Fisher, R. A.  1912.  On an absolute criterion for fitting frequency curves. Messenger of Mathematics 41: 155-160. 

Fisher, R. A. 1921. On the “probable error” of a coefficient of correlation deduced from a small sample. Metron 1: 3-32. 

Fisher, R. A. 1922. On the mathematical foundations of theoretical statistics. Philosophical Transactions of the Royal Society of London  A 222: 309-368. 

Hoang, D.T., Chernomor O., von Haeseler, Minh B.Q., Vinh L.S. 2018. UFBoot2: Improving the ultrafast bootstrap approximation. Molecular Biology and Evolution 35:518-522.

Kalyaanamoorthy, S., Minh B.Q., Wong T.K.F., von Haeseler A., Jermiin L.S. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14(6): 587-589.

Neyman, J. 1971.  Molecular studies of evolution: a source of novel statistical problems. Pp 1-27 in: Statistical Decision Theory and Related Topics, S. S. Gupta and J. Yackel, eds. Academic Press, New York. 

Nguyen, L.-T., Schmidt H.A., von Haeseler A., Minh B.Q. 2015. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution 32(1): 268-274.

