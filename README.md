# Robust Principal Component Analysis and Related Methods
Performs robust principal component analysis (RPCA), computationally efficient truncated RPCA, and L2 regularization for noise reduction. RPCA introduced by Candes et al. (2011) recovers the low-rank and sparse components from the input data, using Principal Component Pursuit. We have proposed two extensions of RPCA that are computationally efficient and reduce noise inherent in high-dimensional genomic data. Truncated RPCA enables application of RPCA for large datasets by only using the adaptively selected top k singular vectors. We have also implemented L2 regularization to reduce noise with an application to single cell RNA-seq. This tRPCA with L2 regularization effectively decompose the input data into low-rank, sparse, and noise components. 

# Installation

To use a stable version from CRAN:
```R
install.packages("rpca")
```

To use a development version from GitHub:
```R
install.packages("devtools")
library("devtools")
install_github("macieksk/rpca")
```

# Reference
Gogolewski K.*, Sykulski M.*, Chung N.C., Gambin A. (2018) Truncated Robust Principal Component Analysis and Noise Reduction for Single Cell RNA-seq Data. ISBRA 2018. Lecture Notes in Computer Science. https://link.springer.com/chapter/10.1007/978-3-319-94968-0_32

Candes, E. J., Li, X., Ma, Y., Wright, J. (2011) Robust principal component analysis? JACM
