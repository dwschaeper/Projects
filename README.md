# Projects
The projects repository is a sampling of the work that I have completed.
Each folder contains an individual work. A description of each is as follows.

Bacteriophage-mutations: Collect many types of mutational data from bactiophages from phagesDB.

Predator/Prey: A graphical simulation of predator prey relationship given parameters than can be adjusted in the main function.

miRNA-classifier: A directory contains miRNA expression data from non-metastatic and metastatic tumors. This data is then processed and the top features are selected through information gain. Those top features are then used to create the classifier.

# Individual Works
DEG_analysis.R  
This R script creates a volcano plot and heatmap of DEGs utilizing the DESeq2 package and input files of gene counts, metadata, and a annotation file.  

predictive_models.R  
This R script generates a model to predict the mean RFFT value of individuals within a dataset. Also, different models are made to investigate the effect of female ornamentation in the ruby throat bird on their fitness.  
  
statistics_and_distributions.R  
This R script uses statistal test and correlation to investigate the relationship between variables in different datasets. It also utilizes Binomial and Poisson distributions to calculate disease probablities and expected prevalence.  
  
RNA_half_life.py  
This python 3+ script takes an input of 3 time course datasets from one file. It calculates the half life of each RNA based
upon the average slope across the 3 time courses. After that, it splits the data into the top and bottom 10% half life
duration and performs a simple functional enrichment analysis using GO.  
  
operon_prediction.py  
This python 3+ script takes an input of 3 ptt files and predicts the operons for the organism. It also takes a .gff file
and predicts the operons from the contigs.  
  
motif_find_and_search.py  
This python 3+ script takes a nucleotide count matrix and creates a PWM. It also takes a file of a list of
sequences, and uses the PWM to predict if the motif from the counts matrix is in the sequences. It outputs the top
30 hits to a csv file.  
  
simple_network_analysis.py  
This python 3+ script takes a PPI text file to build a graph from. It determines if the graph is scale free as well as
the calculates the average clustering coefficient, degree, and clustering coefficient for each node and outputs to a file. Given 2
protein lists it will calculate the shortest path length of all pairs and compare the distribution of the path lengths
using a Wilcox test.
