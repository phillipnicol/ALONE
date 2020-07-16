###demo.R
###November 21, 2019
###Demo run of the genetic algorithm for colorectal cancer

#source the files
source('~/Documents/GitHub/ALONE/Implementation/main.R')
source('~/Documents/GitHub/ALONE/Implementation/helpers.R')

###Replace with your path if needed:
ColorectalCancer <- readRDS("~/Documents/GitHub/ALONE/Data/ColorectalCancer.rds")

###Preprocessing, add a column of all 1s (root node)
ColorectalCancer <- cbind(1, ColorectalCancer)

###Run genetic algorithm, print graph
Results = GA(ColorectalCancer, 100, 300, 0.5, 0.5, 0.01, 0.05, 0.01, rep(0.05, ncol(ColorectalCancer)),
             suppress = FALSE, method = "SA", penalty = "hard")
plot(Results$igraph_obj)
