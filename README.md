# Robust Principal Component Analysis and Related Methods
Performs robust principal component analysis and related matrix decomposition techniques. Robust PCA (introduced by Candes, E. J., Li, X., Ma, Y., & Wright, J. (2011). Robust principal component analysis?. JACM) recovers the low-rank and sparse components from the input data, using Principal Component Pursuit. Truncated version of robust PCA further improves this matrix decomposition for high-dimensional data in a computation- and memory-efficient manner. Particularly, the top k singular vectors (using irlba proposed by J. Baglama and L. Reichel (2005) Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, SIAM J. Sci. Comput.) are used adaptively. In a case of highly noisy data, the noise reduction is provided via L2 regularization.

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
