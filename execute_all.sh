#!/usr/bin/env bash

# This bash script can be used to run all python files
# sequentially for the BetweenNET algorithm. The required
# python libraries are numpy, scipy and networkx. The
# input files are provided within the data folder.
cd src


#########################################################
#
#	Compute Betweenness Values
#
#########################################################


printf "Compute Betweenness Values...\n\n"

./betweeness data/Betweenness

#########################################################
#
#	Generate Data for Random Walk
#
#########################################################


printf "Generate Data for Random Walk...\n\n"

python generate_outliers.py data Betweenness human_genes
python construct_bipartite_graph.py data patient_ids InfluenceMatrix MutationMatrix OutliersMatrix

#########################################################
#
#	Perform Random Walk
#
#########################################################


printf "\nRunning Random Walk...\n\n"

python random_walk.py


#########################################################
#
#	Rank Genes
#
#########################################################


printf "\nFinding Connected Components...\n\n"

python rank_genes.py
printf "\nDONE!!!\n"
