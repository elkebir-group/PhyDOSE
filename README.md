# PhyDOSE - Phylogenetic Design Of Single-cell sequencing Experiments

The input to PhyDOSE is a set of candidate trees, a frequency matrix obtained from bulk data, and a confidence level. PhyDOSE provides a minimum number of single cells needed in a follow-up single-cell sequencing (SCS) experiment to determine the true phylogeny among the set of given phylogenies with at least the desired probability.

![Overview of PhyDOSE](overview.png)


### Step 1: Generate the minimal distinguishing feature family for each tree 

From the terminal, 
```
/../../designBatch.sh

```

### Step 2: Use PhyDOSE to calculate k^* 

```
Rscript PhyDOSE.r [path to set of candidate trees] [confidence level] [false negative rate] [path to distinguishing features]
```
By default, PhyDOSE uses a confidence level of 0.95, a false negatve rate of 0 and assumes that the distinguishing feature files are located in a subdirectory of thge candidate trees file named **distFeats** as this is where PhyDOSE would create them during processing.
Therefore, only the input file needs to be specified but the other arguments can be optionally supplied by the user if they differ from the default values. 


### Step 3[OPTIONAL]: Reconciling the SCS Experiment
If an SCS experiment has been conducted and a user wishes to reconcile the experiment against PhyDOSE's distinguishing features, the following command line script can be run.

```
Rscript PhyDOSE.r [path to SCS csv file][path to distinguishing features directory][path to mutation name mapping file]
```

