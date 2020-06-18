# BetweenNET: Ranking Cancer Drivers via Betweenness-based Outlier Detection and Random Walks

### This is the original repository for the BetweenNet paper


**Installing the dependencies**

```
pip install -r requirements.txt
```
## **Betweenness**

### **Input**

1. The PPI network file:

The file is located at data/InfluenceMatrix

```
gene1 gene2 confidence value
MDM2  TP53  0.99
APP APP 0.99
MYC MAX 0.98
...
```
2. Patients data (Normal and Tumor):

Files are located at data/Betweenness
```
gene_id        TCGA-A7-A0CE-01A    gene_id        TCGA-A7-A0CE-11A
ACIN1|22985    2916.8574           ACIN1|22985    3377.1429
ACLY|47        4599.7118           ACLY|47        3405.7143
removed        removed             ACMSD|130013   26.3492
ACN9|57001     11.1874             ACN9|57001     176.1905
...
```
Note: "removed"  refers to a gene is not expressed or is a mutated gene

### **Run Betweenness calculation**
```
cd src
./betweenness ./data/Betweenness
```

### **Output**
Files will be located at out/Betweenness <br/>
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

## **Generate input for Random walk**

### **Input**
1. Mutation matrix
The file contains a matrix of mutated genes.
The file is located at data/MutationMatrix.csv

```
Genes   g1  g2  g3  g4 ... gn
p1      1   0   0   1  ... 1
p2      0   1   1   0  ... 1
p3      1   1   1   0  ... 0
p4      1   0   0   1  ... 1
....
```
2. Outliers matrix
The file contains a matrix of outlier genes.
The file is located at data/OutliersMatrix.csv

```
Genes   g1      g2      g3     ...  gn
p1      True    True    False  ...  True
p2      True    False   False  ...  False
p3      False   True    False  ...  False
p4      True    False   True   ...  True
....
```

### **Outliers matrix generation**
To generate the outliers matrix, run  generate_outliers.py script, as follow:

```
python generate_outliers.py [Input directory] [Betweenness Results Path] [Genes list]
```
### **Output**
Files will be located at out/ <br/>
1-out/graph_nodes.txt <br/>
2-out/graph_edges.txt <br/>
3-out/graph_mut_freq.txt <br/>


## **Random Walk**
### **Input**

1 - The nodes to index file mapping:
The file is located at out/graph_nodes.txt
```
index GeneName
1 A1BG
2 A1CF
```
2 - Edge file:
The file is located at out/graph_edges.txt
```
node_i_Index node_j_Index weight
0 1 1
0 2 1
```
3 - The gene and its corresponding mutation frequency:
This file contains the mutation frequncies, which are assigned as heats during the random walk.

```
A1BG 0.00353697749196
A2M 0.0128617363344
A4GALT 0.00064308681672
```

### **Output**
To rank genes, run:
```
cd src
python rank_genes.py
```
The file will be located at /out/ranking.txt

## **Run**
To run BetweenNET on the given input files:
```
cd src
sh execute_all.sh
```
