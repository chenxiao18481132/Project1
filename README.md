# Project1
Single-cell RNA-seq data are difficult to reduce dimension by using classical methods because of the dropout characteristics. Dropout is an event in which a transcript is not detected in the seqyencing data owning to a failyre to capture it. Pierson and Yau (2015) present the Zero-Inflated Factor Analysis(ZIFA) to analyze zero-inflated single-cell gene expression. This pROJECT will derive an expectation-maximization algorithm to estimate model parameters and apply the algorithm to simulation study and real data to show the effectiveness of this algorithm.

ZIFA_EM.R includes the EM algorithm. zifa simulation.R is the code of simulation study. zifa real data.R is the code of the zifa model for CORTEX data.
ZIFA_VAE is the code for ZIFA VAE.
