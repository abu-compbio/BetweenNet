# BetweenNET: Ranking Cancer Drivers via Betweenness-based Outlier Detection and Random Walks

### This is the original repository for the BetweenNet paper


**Installing the Dependencies**

```
pip install -r requirements.txt
```
## **Betweenness**

### **Input**

1. The PPI network file:

The file is located at data/IntAct_network

```
gene1 gene2 confidence value
MDM2  TP53  0.99
APP APP 0.99
MYC MAX 0.98
...
```
2. Patients Dataset (Normal and Tumor):

Files are located at data/[cancer]/betweenness_input/
```
gene_id        TCGA-A7-A0CE-01A    gene_id        TCGA-A7-A0CE-11A
ACIN1|22985    2916.8574           ACIN1|22985    3377.1429
ACLY|47        4599.7118           ACLY|47        3405.7143
removed        removed             ACMSD|130013   26.3492
ACN9|57001     11.1874             ACN9|57001     176.1905
...
```
Note: "removed"  refers to a non-expressed/mutated gene


### **Betweenness Calculation**

```
cd src
[1st] Compile betweenness.cpp using the following command line:
g++-5 -std=c++0x -I$LEDAROOT/incl -L$LEDAROOT betweenness.cpp -lGeoW -lD3 -lW -lP -lG -lL -lX11 -lm -O1 -no-pie -o betweenness

[2nd] Run the following command to compute betweenness values for samples in the selected [cacner]
./betweenness
```

####**Output**
Files will be located at out/[cancer]/Betweenness <br/>
Ex: TCGA-22-5478-01.txt
```
Gene         Betweenness value
A1BG    :    676.579
A2M     :    80833.4
A2ML1   :    78.5454
AAAS    :    1606.19
AACS    :    0
AADAC   :    0
AAGAB   :    904.682
...
```


## **Construction of the Bipartite Graph**

### **Input**
1. Mutation Data
The file is located at data/[ cancer ]/mutation_data.txt

```
p1      g1 g2 g6 .... gn
p2      g19 g2 g16 .... gn
p3      g3 g1 g16 .... gn
....
```


2. Outliers Data
The file contains a matrix of outlier genes.
The file is located at data/ [ cancer ] /outliers_data.csv

```
Genes   g1      g2      g3     ...  gn
p1      True    True    False  ...  True
p2      True    False   False  ...  False
p3      False   True    False  ...  False
p4      True    False   True   ...  True
....
```

### **Selection of Outlier Genes**
To generate the outliers matrix, run  "generate_outliers.py" script, as follows:
```
python generate_outliers.py [cancer] -v-[t/f]
```
### **Output**
Files will be located in out/ <br/>
1-out/  [ cancer ] /graph_nodes.txt <br/>
2-out/ [ cancer ] /graph_edges.txt <br/>
3-out/ [ cancer ] /graph_mut_freq.txt <br/>
4-out/ [ cancer ] /bipartite_graph.gml<br/>


## ** **
To construct the bipartite graph, run "construct_bipartite_graph.py" script, as follows:
```
python construct_bipartite_graph.py [cancer] -v-[t/f]
```

## **Random Walk**
**Note:** For the random walk process, we used the random walk python code for both breast and lung datasets, and the iterative version (described in our paper) for pan-cancer dataset for three iterations.
### **Input**

1 - The nodes to index file mapping:
The file is located at out/graph_nodes.txt
```
index GeneName
1 A1BG
2 A1CF
```
2 - Bipartite edges file:
The file is located at out/graph_edges.txt
```
node_i_Index node_j_Index weight
0 1 1
0 2 1
```
3 - The gene and its corresponding mutation frequency:
This file contains the mutation frequencies to be assigned as heats during the random walk.

```
A1BG 0.00353697749196
A2M 0.0128617363344
A4GALT 0.00064308681672
```

To apply random walk process on the constructed bipartite graph, run random_walk.py script, as follows:
```
python random_walk.py [cancer] -v-[t/f]
```


### **Output**
To rank mutated genes in the bipartite graphs, run:
```
cd src
python BetweenNet.py [cancer] -v-[t/f]
```
The ranking of mutated genes will be located in /out/ [cancer] /BetweenNet.txt
