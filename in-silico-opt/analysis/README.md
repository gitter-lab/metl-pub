## Clustering and Analysis



***Note***: If you are hoping to replicate the results from the clustering portion of simulated annealing from the metl manuscript, feel free to contact Bryce Johnson (`bcjohnson7@wisc.edu`) before attempting to run. 

The code is not in an easy to run framework. The clustering (`cluster.py`) and analysis (`analysis.py`) are not meant for widespread adoption. 
This makes the code difficult to run, so please feel free to reach out if interested. 


### Repository structure

- `cluster.py` : Code for clustering the sequences that are found through simulated annealing. 
- `analysis.py` : contains code for running one off experiments of the data that was generated from the 10,000 experiments.
- `preprocesing_gfp_run.py` : Code to preprocess and understand the 64 training variants that were chosen for the GFP design scenario. 
- `comparison_stats.py` : Compares two clusters of sequences against each other for position, score, or amino acid distribution. 
- `general_stats.py`: Look at the statistics of a single cluster in terms of amino acid distribution, position distribution and score. 
- `figure_plots.py` : Make the figures in official plotting style in the METL manuscript 
- `training_set.py`: Further analysis of the training set. 
- `unit_tests.py` : Unit tests for clustering and final sequences that are chosen for experimental designs. 
- `distances.py`: Various distances that were tested for computing clusters, with varying compute times