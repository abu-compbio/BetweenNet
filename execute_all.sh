#!/usr/bin/env bash

# This bash script can be used to run all python files
# sequentially for the BetweenNET algorithm. The required
# python libraries are numpy, scipy and networkx. The
# input files are provided within the data folder.
cd src
cancer_type="BRCA"
verbose=f #if "t", then all steps will be displayed in terminal
mkdir "../out/"$cancer_type

#########################################################
#
#	Compute Betweenness Values
#
#########################################################

printf "################################################\n"
printf "    1 - Betweenness Calculation...\n"
printf "################################################\n\n\n"
#g++-5 -std=c++0x -I$LEDAROOT/incl -L$LEDAROOT betweenness.cpp -lGeoW -lD3 -lW -lP -lG -lL -lX11 -lm -O1 -no-pie -o betweenness
#./betweeness cancer_type

#########################################################
#
#	Generate Data for Random Walk
#
#########################################################

printf "################################################\n"
printf "    2 - Selection of Outliers...\n"
printf "################################################\n\n\n"
python generate_outliers.py $cancer_type -v-$verbose

printf "################################################\n"
printf "    3 - Constructing Bipartite Graph ...\n"
printf "################################################\n\n\n"
python construct_bipartite_graph.py $cancer_type -v-$verbose

#########################################################
#
#	Perform Random Walk
#
#########################################################

printf "################################################\n"
printf "    4 - Random Walk Process...\n"
printf "################################################\n\n\n"

python random_walk.py $cancer_type -v-$verbose


#########################################################
#
#	Rank Genes
#
#########################################################


printf "################################################\n"
printf "    5 - Ranking Mutated Genes...\n"
printf "################################################\n\n\n"
python BetweenNet.py $cancer_type -v-$verbose
printf "\nDONE!!!\n"
