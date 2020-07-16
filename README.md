# ALONE
Implements a genetic algorithm for learning the structure of a Bayesian network with the specific application
of inferring the temporal order of mutations in cancer. The algorithm is described in 
[this paper](https://www.biorxiv.org/content/10.1101/2020.04.13.040022v2.article-metrics).

## Usage
The algorithm is currently implemented as an R script, so you will have to source `main.R` and `helpers.R` before using the algorithm. The core function is `GA`, 
which takes as input a binary mutation matrix in which mutations are columns and samples are rows.

`demo.R` contains an example with colorectal cancer CNA data, which is included in the package. 
