# Placoderm Data Science Project

## Introduction
*This was a project performed by Charles Sansone with assistance from Dr. Lauren Sallan over Summer 2017*. Placoderms are an ancient extinct class or armored fish that existed from 430mya to 360mya. Notably, placoderms were the first fish to develop true teeth, and were among the first to develop jaws. Therefore, Placoderms are a sister group to living jawed vertebrates. This project hoped to reveal the ancestral habitats of placoderms and therefore the setting where the origination of jaws occurred. Whichever Benthic Assemblage (BA) zone it occurred in would inform us as to whether Placodermi were an apex predator in open water, in competition with invertebrates in reefs, etc. The origination environment would also tell us whether the diversity of placoderms, and thus jawed fishes, originated in a restricted environment, or as a result of transitions.

To get a good overview of the project, go to **Plots -> CURF2017_CHARLES_SANSONE**. Download and open the file. This is a scientific research poster I made and presented at the **CURF2017 Research Expo**. It details the steps I took, the plots I generated, and the results of my project, all in an accessible manner.

If you would like to view important code that I used, the bulk of R code is contained in **Tree Timescaling and ancThresh -> AncThreshNewScale.R**. If you would like to run my experiment yourself, continue on below.

### Prerequisites
Make sure that you have the following programs and installed on your computer: R-Studio, Mesquite, and MrBayes. Next begin by downloading the files in my repository, and place them into folders with the same names as they appear on here. In addition, you will need to install the packages "Phytools" and "Ape" onto R.

## Running the Experiment
Open MrBayes, then type **execute filename**. The filename should be the path to the .nex file titled MRP_matrix. Following this, you should add parameters to the MCMC(Markov Chain Monte Carlo) algorithm via the command **lset applyto = (1) coding = variable rates = gamma; [morphology] unlink statefreq = (all) revmat = (all) shape = (all); prset applyto = (all) ratepr = variable;** Next, type **mcmc** and MrBayes will run the MCMC algorithm (namely a Markov K model). When MrBayes prompts asking whether to continue the analysis, study if the average standard deviation of the split frequencies is below 0.01, and if not continue the analysis until it does. After this, note the total number of generations and the average standard deviation of split frequencies. Now, decide on the burn-in. This will be (total number of generations) / 100 * 0.25. Then run the following command, with XX being the burnin you found:
**sumt contype=halfcompat showtreeprobs=yes burnin=XX**
This will make MrBayes summarize the remaining trees and output a .tre file which you can then open in Mesquite. This .tre file will contain many different trees that have different branchings.

Next, you will want to create a majority consensus tree of all the trees in this .tre file. Mesquite is able to do this quite easily for you: **Taxa&Trees -> Make New Trees Block From -> Consensus Trees from Tree Blocks -> Stored Tree Blocks -> Majority-Rule Consensus -> OK** Next, you will supplement this tree with any occurences you can find of the same genus as the ones in this tree. Add them as polytomies to the nodes of that genus. When I did this I eventually ended up with the supplemented_supertree_edited.nex file, so you can use this one as well if that suits you.

Following this, go to the Tree Timescaling and ancThresh folder. Select AncThreshNewScale.R and run it in R-Studio. This script requires two inputs: 1. a PCM of the placoderm occurences that provides age ranges within 5 million years. 2. the phylogeny we made in the previous folder. One important caveat is that the PCM and Phylogeny occurence names need to be **identical**, both in order and in naming (see my two files for an example). In addition, eliminate all spaces from the Fauna that each occurence was found in, and abbreviate them as much as possible (see my file Placoderm_PCM_new.txt for an example).

Now you are ready to execute the R script. Go through the file, line by line. The comments should provide ample guidance as to how to proceed in Step 1. In Step 2, we use the ancThresh package and it gets more complicated. Basically, we will run 3 different models of ancThresh: "BM", "OU", and "Lambda". Then we compare the DIC(Deviance Information Criterion) scores and whichever model has the lowest DIC score is the one you continue with in Step 3 (most likely OU).

Next in Step 3 we will obtain the mean BA scores of the data, which will demonstrate how easily the fish moved from locale to locale. We will also plot the probability distributions for the parameters of our ancThresh model. Interpretations of the results are included in the Poster provided in the **Plots** folder. 

## Authors
* **Charles Sansone** - University of Pennsylvania ('20) - pursuing Bachelors of Arts and Sciences in Computer Science and Economics (Dual Major)
* **Dr. Lauren Sallan** - Martin Meyerson Assistant Professor in Interdisciplinary Studies

## Acknowledgments
* I would again like to acknowledge Dr. Lauren Sallan, my faculty mentor under Penn Earth & Environmental Science
* PURM for funding my research
* University of Pennsylvania College of Arts and Sciences
