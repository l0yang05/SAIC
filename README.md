# SAIC
Single cell RNA-seq data analysis algorithms
We have developed a novel algorithm, called SAIC (Single cell Analysis via Iterative Clustering), that identifies the optimal set of signature genes to separate single cells into distinct groups. Our method utilizes an iterative clustering approach to perform an exhaustive search for the best parameters within the search space, which is defined by a number of initial centers and P values. The end point is identification of a signature gene set that gives the best separation of the cell clusters. Using a simulated data set, we showed that SAIC can successfully identify the pre-defined signature gene sets that can correctly separated the cells into predefined clusters. We applied SAIC to two published single cell RNA-seq datasets. For both datasets, SAIC was able to identify a subset of signature genes that can cluster the single cells into groups that are consistent with the published results. The signature genes identified by SAIC resulted in better clusters of cells based on DB index score, and many genes also showed tissue specific expression. 
