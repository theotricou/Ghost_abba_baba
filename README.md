# Ghost lineages deceive introgression tests and call for a new null hypothesis

This repository contains the data and code associated with the manuscript entitled "Ghost lineages deceive introgression tests and call for a new null hypothesis" by T Tricou, E Tannier and DM de Vienne.



* The folder `approach1` contains all scripts required to run ms simulation and generate datasets to compute the proportion of erroneous interpretation of D-statistics. The python3 script converts a species tree, usually generated using Zombi (https://github.com/AADavin/Zombi) in a ms format. Precisely in a coala (https://github.com/statgenlmu/coala) format to use R for simulations. It also contains the R script to run ms simulation and compute all possible D-statistics.
* The folder `approach2` contains the R functions to generate the dataset that was used to estimate the proportion of erroneous interpretation of D-statistics with the approximated approach based on tree topology only (not using ms). The folder also contains the functions used for producing most of the plots of the manuscript.
