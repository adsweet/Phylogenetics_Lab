# Laboratory 7
## Bayesian Inference

### Introduction

The idea that we could infer phylogenies by implementing Bayes' theorem was introduced in the 1960s (Felsenstein's 1968 dissertation even touched on it), but only in the early 2000s, with the introduction of programs like MrBayes, has it been widely implemented.  Bayesian inference is similar to maximum likelihood because one of the central features of Bayes' theorem is the likelihood function. As in maximum likelihood, we use probabilistic models of sequence evolution. Unlike other approaches, Bayesian analysis provides us with estimates of branch support in the form of posterior probabilities. These branch supports, unlike bootstrap values (which technically represent a measure of repeatability of the data), represent the probability that a clade is correct given the model and the data. Bootstrap values do not represent probabilities, in this sense, and have been shown to be more conservative than posterior probabilities (however, ultrafast bootstrap values, which we covered in the last lab, are also less conservative than bootstrap values). Posterior probabilities tend to be inflated, such that values over 0.95 posterior probability are typically considered strong support.

Where did Bayes' theorem come from? Rev. Thomas Bayes was a minister in southern England during the 1700s, and after he retired from the church he focused his efforts on probability theory. In 1763 he developed his famous theorem, above. This has been adapted to phylogenetic inference with the following equation:

 ```math
P(\theta|D) = \frac{P(D|\theta)P(\theta)}{\int P(\theta)P(D|\theta)d\theta}
```

Where P(&theta;|D) is the posterior probability of the hypothesis, represented by &theta;. In phylogenetics, &theta; includes the topology and its associated branch lengths, and the substitution model parameters. Given that we are now working with multiple continuous parameters, the denominator is represented by a complex integral instead of a summation.

P(&theta;) is the prior and P(D|&theta;) is the likelihood function, which you have seen applied in the ML method. An appropriate model (&theta;) needs to be selected to evaluate the likelihood function, and prior distributions need to be chosen for all the parameters in the model. Most often, uniform (flat) priors are used, which place equal weight on all parameter values. Branch length and substitution model parameters often use continuous distributions, such as the uniform or gamma probability distributions.

Because of the complexity of Bayesian estimation using multiple models, it is not actually possible to estimate the denominator. To get around this, we commonly use the Markov Chain Monte Carlo random walk simulation method to obtain the estimate. This is described in a series of steps:
1. The Markov chain starts at some state &theta;<sub>i</sub> (this could be the initial prior state)
2. A new state, &theta;<sub>j</sub>, is proposed by taking a draw from a specified probability distribution. This is known as the proposal mechanism
3. The following formula, known as the acceptance ratio, is calculated:
   
 ```math
R = \frac{P(&theta;_j)P(D|&theta;_j)}{P(&theta;_i)P(D|&theta;_i)}
```
		 
4. If R $`\geq 1`$ then &theta;<sub>j</sub> is accepted
5. If R $`< 1`$ then a number is drawn from a uniform distribution (0,1)
6. If this number is $`<`$ R then &theta;<sub>j</sub> is accepted
7. If this number is $`>`$ R then the Markov chain stays in the same state

The MCMC results in many trees of similar likelihood. As implement, Bayesian analyses summarize the relationships found on these trees. The posterior probability for each node is, in practice, representative of how often each relationship occurred in the collection of sampled trees. In principle, the initial trees are discarded (called the "burn-in"), with only the trees that have reached stationarity being kept. Stationarity means that the likelihood of the trees stays within a very narrow range for a very many generations of the MCMC. Multiple runs of Bayesian analyses are usually performed to verify that the best point of stationarity has been reached, because the MCMC algorithm is stochastic (random), and may get stuck in local likelihood peaks.

Today we will estimate relationships between gall wasps (Cynipidae) and their close relatives (superfamily Cynipoidea) using Bayes’ theorem and data from four genes (cytochrome oxidase I (_cox1_), elongation factor 1-alpha (_EF-1&alpha;_), longwave rhodopsin (LW _Rh_), and 28S ribosomal RNA (_28S_)) and morphology. This lab was adapted from Jim Whitfield and Sydney Cameron's systematics lab taught at the University of Illinois, portions of which were adapted directly from Mike Alfaro, which he designed for the Workshop in Applied Phylogenetics (2003, in Bodega Bay, CA).

**Objectives**

In this lab you will:
- Be introduced to running Bayesian phylogenetic analyses
- Be introduced to assessing MCMC runs using R and other tools

### Exercise 7.1

Log on to Jetstream and launch the Web Shell on your Instance.

Navigate to your volume or home directory and create a new directory called `Lab7_bayes`

Next, change to the `Phylogenetics_Lab` directory and update the data with:

```
git pull
```

You should see a new directory in `Phylogenetics_Lab` called `Lab7` which contains another directory called `Data`. Copy the contents of the `Data` directory to your `Lab7_bayes` directory. Change back to your `Lab7_bayes` directory and extract the contents of this file using:

```
tar -xvzf Lab7_Bayes_data.tar.gz
```

Download the file "gall_wasp.nex" and open it in a text editor. It is a standard NEXUS file except for the MrBayes command block at the end. The block currently in this file divides the data matrix into character sets ("charset") and partitions so that they can be treated differently in future analyses. These are actually commands that you can run when you are running MrBayes in interactive mode. For the sake of learning MrBayes, we will start off running it interactive mode. To launch MrBayes, type the following command into your terminal:
```
 mb
```
After hitting Enter, you should see some information about MrBayes followed by the prompt:
```
MrBayes >
```
This tells indicates that you are successfully in the MrBayes interactive mode. You will enter in all your MrBayes commands after this prompt (similar to the `$` in the terminal window). Next, type the following command after the `MrBayes >` prompt and hit "Enter":
```
execute gall_wasp.nex
```
MrBayes will tell you what it reads in your file and whether it read successfully (look at this output carefully). For example, notice if it correctly reads in the taxa and characters, and identifies the character partitions for data types. Note that there are 25 morphological characters incompatible with the specified coding bias, since MrBayes assumes that all morphological characters are variable. We will let it behave as default for now (automatically removes the characters). Also check if the MrBayes block was read in correctly, including setting the outgroup, defining charsets, defining partitions from the charsets, setting the number of runs to 1. You should be able to see how all of these behaviors correspond to specific lines in the NEXUS file. Before we move on, type this command at the prompt:
```
help
```

You will get a list of commands that can be used in your analysis.  We’ll focus on these commands in the lab today: set, lset, prset, and mcmcp. You can find out more about individual commands by typing “help” before it.  Try this one:
```
help lset
```
MrBayes tells you that `This command sets the parameters of the likelihood model...`, and that `The correct usage is: lset <parameter>=<option> ... <parameter>=<option>`, etc.  You can always use help if you get confused about what you are doing and how the command should be written. So, now we can start analyzing our data. First, we will set up a partition that divides the characters by type (i.e., molecules and morphology; this is determined in our data block where it’s stated that "datatype=mixed"):
```
set partition=default
```
To see how the priors are currently defined for these two partitions, use this command:
```
help prset
```
You’ll see a long list of priors and what they refer to, and then at the end will be two tables (one for each partition) that list the current prior settings. Understanding these parameters is beyond the scope of this lab, but I encourage you to spend some time examining these priors and what they mean. Right now, they are set to default values (called flat or uninformative priors), and we will leave them as is. We are ready to begin analyzing our data, but we have to set the number of generations, etc. using the mcmcp command. Examine its format first:
```
help mcmcp
```
Let’s do a small run with 10,000 generations (ngen) and save the branch lengths (savebrlens). We’ll save trees every 10th generation (samplefreq) and print our progress to the screen (printfreq) every 10th generation. Keep in mind that a “real” MrBayes analysis should have 1,000,000 or many more generations (ngen=1000000) (more is better) with a sample frequency of at least 1000 (samplefreq=1000). This is because it can typically actually take this long to properly sample the posterior distribution. Note it can take even longer depending on the complexity of the model and the strength of signal in the data from which to estimate parameters. Your current command line should look like this:
```
mcmcp ngen=10000 printfreq=10 samplefreq=10 savebrlens=yes
```
To start the Markov chain (i.e., start your analysis) use the mcmc command:
```
mcmc 
```
Every 10th generation will be printed to the screen. In the MrBayes window you will see a summary of the chain conditions. The screen will be updated with six columns that track the status of the Markov chain. The first column shows the current generation number. The next four columns show the log likelihoods of each chain’s tree. This is technically a Metropolis-Coupled MCMC, so there is one cold chain and three hot chains. The hot chains move more rapidly through parameter space, allowing the MCMC to better escape local optima. Notice that one of the columns is surrounded by square brackets. This is the cold chain, the chain that is actually being recorded to the out file. Note that the square brackets move around between chains. This indicates that a search of tree space by a heated chain yielded a more likely tree topology and this heated chain subsequently became the new cold chain. The last column is the estimated time left in the analysis. When the analysis is finished, MrBayes asks you if you want to continue – type `n` and hit Enter.  

Now exit the interactive MrBayes program.
```
quit 
```
Now you should have several new files in the same directory as your original NEXUS file, but two in particular are important to note: "gall_wasp.nex.p" (list of recorded parameters during the analysis) and "gall_wasp.nex.t" (all the recorded trees visited). Download these two files and open them in a text editor and examine their structure. 

In a typical Bayesian analysis, we want to discard the early samples because the trees were not yet at stationarity (the chains were still searching and had not yet reached the set of most likely trees). These samples to be discarded are referred to burn-in. In order to assess the burn-in, we will use the program Tracer . 

Open the program [Tracer](https://github.com/beast-dev/tracer/releases/tag/v1.7.2). Click on the "+" sign below the Trace Files window to load your .p file. Make sure you tell the program to look for all types of files, as it defaults to looking for only .log files!

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/tracer_plus.png" width=25% height=25%>

Load the "gall_wasp.nex.p" file. By default, Tracer sets the burn-in at 10% (1000) of the number of generations. Set this value to 0 to see all of the generations by double-clicking on the value below the Burn-In column. To choose an appropriate burn-in value, we will first click on the Trace tab at the top of the right-hand window.

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/tracer_load.png" width=25% height=25%>

:white_check_mark: __What did you choose as a burn-in?__

:white_check_mark: __Did your analysis reach stationarity?__

You should see that the traces are still steadily increasing, rather than reaching a stable level. This means that we need to run the analysis for longer so that you can better sample the posterior distribution. Because we will be running it for longer, let us now use the non-interactive mode.

To run MrBayes commands in non-interactive mode, we simply need to add the commands to the MrBayes block in a NEXUS file. To modify the Bayes block, open the gall_wasp.nex file in a text editor. Delete the `mcmcp nruns=1` line. To keep things simple we set our first MrBayes analysis to only have one run, but it is a better practice to run multiple runs. The defaul number of runs for MrBayes is 2. After the partition command and before the end of the Bayes block, enter the following lines, or copy and paste them in.

```
set partition=default;
mcmcp ngen=1000000 printfreq=100 samplefreq=500 savebrlens=yes;
mcmc;
```

You’re now ready to run your analysis. 

Note that this run will take about 10 minutes. To ensure you can still work in the terminal (without opening another tab), type the following command in your terminal window to create a screen:
```
$ screen -S mrbayes
```
The screen command lets you open another terminal window, which allows you to start an analysis in the screen, detach from the screen, and continue doing other things while your first analysis runs in the background.

Once you’re in a new screen, run MrBayes by typing the following command:
```
mb gall_wasp.nex
```
This is generally how you’ll want to run MrBayes because saving your commands in a Bayes block in a NEXUS file helps keep a record of your commands.

After starting your run, detach from your screen by typing `CTRL + a`, then `d` after lifting up from `CTRL + a`. To check on your run, you can reattach to the screen with the command:
```
$ screen -r mrbayes
```
When you are out of your screen, you can list the files in your `Lab7_bayes` directory. Notice there are now two .p and .t files (run1 and run2), because we set the number of runs back to the default in the mcmcp command (nruns=2). Note that MrBayes uses the NEXUS file name by default for all output files.

Once your analysis finishes, download the .p files and take a look at them in Tracer. If you change the Burn-In to 0, you will see the entire posterior distribution. Here, we can see the likelihood values (plotted on the y-axis) vs. the generation number (plotted on the x-axis). The first step will be choosing a rough spot in the distribution as a preliminary burn-in. This will be the point where the graph seems to "level off." I have marked a red line at a rough point in the following example:

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/burnin.png" width=50% height=50%>

Estimate the generation number where the graph levels off and type this number into the Burn-In box. The Trace plot will change to reflect this burn-in, and you can then more closely investigate the data to find the appropriate burn-in value. An appropriate burn-in value will produce a roughly randomly distributed pattern across a horizontal line, the so-called “fuzzy caterpillar”. When the appropriate burn-in is selected, the plot will instead look like the following (the grayed out area is part of the burn-in, above right).

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/bayes_tracer.png" width=50% height=50%>

When you believe you have the correct burn-in, click to view the "Estimates" tab. This tab will show you a distribution of the likelihood values. If the burn-in was properly chosen and stationarity has been reached, the plot will converge on a normal distribution (compare this to when you have burn-in set to 0).

<img src="https://github.com/adsweet/Phylogenetics_Lab/blob/main/Images/bayes_estimates.png" width=50% height=50%>

In the column on the left-hand side, you can choose to look at other estimated parameters.  For example, specific rates of nucleotide change, base frequencies, and mutation rates may be estimated. The estimates that are present depend on your model. You will also notice a column labeled “ESS.” This is the estimated sample size of the parameter, and is another measure typically used for testing whether stationarity was reached. In general, values above 200 are excellent, and values over 100 are acceptable. If the ESS is below that for any parameter, you may have to run MrBayes for a longer number of generations. Tracer also allows you to compare plots from several independent analyses to determine whether you have converged on the same tree. 

:white_check_mark: __Did your run reach stationarity?__

:white_check_mark: __What burn-in value did you select for this run?__

:white_check_mark: __What was the mean likelihood of the post-burn-in data?__

Another way to assess your MCMC runs is to check if your trees have reached convergence between different runs. To do this, we will use the R package RWTY (R We There Yet).

Log out of your shell and load the Web Desktop. Open RStudio and load the package:
```
library(rwty)
```
Change to your Lab7_bayes directory with `setwd()`

Next, load in your MrBayes data with the command:
```
gall <- load.multi()
```
You should now see a large list object called “gall” that has trees, .p files, etc. This is your MrBayes output.

To check if your runs have converged, we will use the Average Standard Deviation of Split Frequencies (ASDSF), which basically looks at the variation in your trees among multiple runs. As a rule of thumb, multiple runs have reached stationarity when the ASDSF values go below 0.01. To check your ASDSF, run the command:
```
makeplot.asdsf(gall[c(1,2)])
```
:white_check_mark: __Save your ASDSF plot as LastName_ASDSF.pdf.__

:white_check_mark: __Based on the ASDSF plot, would you update your burn-in?__

:white_check_mark: __RWTY?__

After you have determined the burn-in, we will look at the parameters that MrBayes estimated for your run using the `sump` command ("summarize the parameters"). Note that it is necessary to divide the burn-in from Tracer by the sample frequency chosen in MrBayes (e.g., if Tracer burn-in is 100,000 generations and your samplefreq was 500, you would use `burnin=200`). 

Make sure you navigate back to the folder that contains your .p files, either through your RStudio terminal tab or by opening a Web Shell. Run `mb` to open up the interactive MrBayes prompt, and run the below:
```
sump burnin=200 filename=gall_wasp.nex;
```
sump will look for files in the current working directory that start with gall_wasp.nex and end with either .p or .run#.p. Here you can see all of the parameters and view a small, rough diagram plotting the likelihoods similar to the one you saw in Tracer. Once again, the parameter estimates will depend on your model. The sump command will give you slightly different information from the Tracer estimates – for example, the harmonic mean of likelihoods, which can be useful for comparisons between runs, is only given using sump.

Now that we know our burn-in value and have examined the parameters, we will obtain a consensus tree of your visited/recorded trees. We will use the `sumt` command ("summarize the trees"). Note that sumt requires the matrix to be loaded into memory first. The easiest way to do this is to execute the gall_wasp.nex file, but __before doing that, delete the mcmc command__ (otherwise your run will start again and overwrite your output files). Quit MrBayes and then edit the NEXUS file. After deleting the mcmc command, you can launch MrBayes again and run the following commands:
```
execute gall_wasp.nex 
sumt burnin=200 
```
This outputs a 50% majority rule consensus tree (saved into a tree file called "gall_wasp.nex.con") with posterior probability values and branch lengths (actually two trees are saved - one with and one without posterior probabilities), as well as two other files called "gall_wasp.nex.parts" (contains the posterior probability of all bipartitions and displays the majority rules consensus tree with posterior probabilities – we’ll ignore this in lab) and "gall_wasp.nex.trprobs" (contains the trees that were found during the MCMC search, sorted by posterior probability – we’ll also ignore this). The numbers on the branches are posterior probability values. As a rule of thumb, posterior probability values ≥0.95 are considered good support, roughly equivalent to a standard bootstrap of 70-75.

Open "gall_wasp.nex.con.tre" in FigTree, root on "Ibalia", turn on Branch Labels, and select Display: prob. 

:white_check_mark: __Export the file as LastName_MRC.pdf.__

### Exercise 7.2

We just completed an analysis where all the molecular and morphological characters were treated with the equivalent of the Jukes-Cantor model (the default model). Now let’s make our models and analyses a little more sophisticated. Unlike IQ-TREE, MrBayes does not automatically implement model and partition testing. As a reminder, you can run ModelFinder using the `-m MF` option in IQ-TREE. However, MrBayes only implements a few models, so when using ModelFinder to test models for MrBayes, you should set it to only test models available in MrBayes by providing the argument `-mset mrbayes`. 

Let’s say we ran ModelFinder and found that the general time reversible (GTR) model with a proportion of the sites invariable (I), and the rate for the remaining sites are drawn from a gamma distribution [G], referred to in MrBayes as invgamma, fits our _cox1_, _EF-1&alpha;_, and _28S_ datasets (GTR+I+G model).  The HKY + G model fits the LW _Rh_ data better. It can be hard to remember how to set models, but this guide is helpful: https://gist.github.com/brantfaircloth/895282

To set up the model, we can modify the Bayes block again. Make a copy of the gall_wasp.nex file called gall_wasp2.nex:
```
cp gall_wasp.nex gall_wasp2.nex
```

Make a new directory called “gall_wasp2” and move gall_wasp2.nex file into this directory.

Download the gall_wasp2.nex file and open it in a text editor. Delete the line:

```
set partition=default
```

and in place of that type in (or copy/paste) the following new command block. It should appear before the mcmcp command:
```
set partition=data;
ctype ordered: 1 3 5;
lset applyto=(2,3,5) nst=6 rates=invgamma;
lset applyto=(4) nst=2 rates=gamma;
unlink shape=(all) pinvar=(all) statefreq=(all) revmat=(all);
```

First, we set up the partition called “data” (it’s defined in the first MrBayes block) using the set command. We forced characters 1, 3, and 5 to be ordered using the `ctype` command. Then we used the `lset` command to assign likelihood models to each molecular character set. This is done using `applyto` and specifying the character set in the order in which they are listed in the partition “data” (e.g. _cox1_ is partition 2, _EF-1&alpha;_ is partition 3, etc.). The `unlink` command allows each character set to have its own parameters, free from other character sets, by specifying `all` partitions. Use the `help unlink` command to learn more. 

Remember, your file should still have the `mcmcp` command with the MCMC parameters. For now, delete the `mcmc` command. After making these edits, save the file and upload to your `gall_wasp2` folder.

Let’s check that all parameters were set correctly. Start up MrBayes in interactive mode:
```
mb
```
Then, execute the new file:
```
execute gall_wasp2.nex
```

It is a good idea to run the `execute` command before launcing an MCMC run to make sure your file was read correctly.

If all our parameters check out, exit MrBayes, go back to your gall_wasp2.nex, and add the ```mcmc;``` command (don’t forget the semicolon!) below the ```mcmcp``` line. Reattach to your screen:
```
screen -r mrbayes
```
and launch MrBayes
```
mb gall_wasp2.nex
```

Detach from your screen and wait for the analysis to finish.

While you’re waiting for your MrBayes run to finish, estimate an unpartitioned ML tree using IQ-Tree on the DNA data from the gall wasps (“gall_wasp_dna.fasta” from your data download). Be sure to run a model test and at least 1,000 ultrafast bootstrap replicates.

:white_check_mark: __Open the resulting tree in Figtree and save as LastName_gall_ML.pdf.__

Also compute an NJ distance tree in R (remember to set the correct working directory!):
```
library(ape)
gall_dna <- read.FASTA("gall_wasp_dna.fasta")
gall_dist <- dist.dna(gall_dna, model = "K80")
gall_dist_nj <- bionj(gall_dist)
gall_dist_nj <- root(gall_dist_nj, "Ibalia", resolve.root = T)
plot(gall_dist_nj)
```
:white_check_mark: __Save the resulting NJ tree as LastName_gall_NJ.pdf.__

When your MrBayes analysis has finished, download the .p files and examine both in Tracer and plot the ASDSF in `RWTY`. Once again, notice the estimates for different model parameters from the MCMC runs (including substitution rates, κ, etc.).

Remember to set the appropriate working directory (setwd()) before running RWTY.

:white_check_mark: __What is an appropriate burn-in value for this run?__

Now perform a sumt that incorporates your newly determined burn-in value. __Remember to delete the ```mcmc``` line from gall_wasp2.nex before executing.__
```
mb
```
```
execute gall_wasp2.nex
sumt burnin=<input your burn-in here>
```

When this is complete you can exit MrBayes by typing `quit`. Open your .con tree file in Figtree and show posterior probabilities (triangle next to Branch Labels, Display: posterior probabilities).

:white_check_mark: __Save this tree as LastName_MRC_Mixed.pdf.__

:white_check_mark: __How do the ML, NJ, and Bayesian trees compare?__
	

### References

Huelsenbeck, J. P., F. Ronquist, R. Nielsen, and J. P. Bollback.  2001.  Bayesian inference of phylogeny and its impact on evolutionary biology. Science 294: 2310-2314.

Huelsenbeck, J. P and F. Ronquist. 2001. MRBAYES: Bayesian inference of phylogeny. Bioinformatics 17: 754-755.

Huelsenbeck, J. P., B. Larget, R. E. Miller, F. Ronquist.  2002.  Potential applications and pitfalls of Bayesian inference of phylogeny.  Systematic Biology 51 (5): 673-688.

Ronquist, F., J. P. Huelsenbeck.  2003.  MrBayes3: Bayesian phylogenetic inference under mixed models.  Bioinformatics 19 (2): 1572-1574.

Yang, Z. 2006. Computational Molecular Evolution. Oxford, Oxford University Press.

Yang, Z. and B. Rannala.  1997.  Bayesian phylogenetic inference using DNA sequences: a Markov chain Monte Carlo method.  Molecular Biology and Evolution 14 (7): 717-724.
