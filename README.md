# InferN0
Inference of regulatory Network under L0 regularization

The analysis of the coexpression of genes in single cell can be used to retrieve potential regulatory elements; however, the preriquisite for this is to distinguish between intrinsic variability from either experimental varibility and/or sampling variability. This motivates the use of the L0 regulatization, as opposed to L1 or L2 regulazation, since it is capable of finding the most probable causal associations in highly sparse regulatory networks with the lowest false discovery rate:

![](man/PrecisionRecallBenchmark.png)
![](man/SingleCellNormalization.png)
![](man/ExampleOutput.png)

============
#Installation

Prerequisites:
  - gcc
  - Rpackages: Rcpp, RcppArmadillo, ggplot2

Alternatively, use provided script to install prerequisites:

  git clone https://github.com/lfhandfield/InferN0.git
  cd InferN0
  sudo sh install.sh
