#***********************************
# Run the command 'python random_walk.py' in the Terminal to start the
# the random walk process on the bipartite graph
# Or directly run this file in your chosen IDE.
# This may take several hours.
# When is done, the files of prepared data are saved in the subfolders of out/[cacner]/.
NOTE: we used this code for both breast and lung cancer datasets
#***********************************


import scipy.sparse as sp
import numpy as np
import argparse
from scipy.sparse.linalg import inv
from tqdm import tqdm
import sys


#load gene ids
def load_gene_to_id():
    filename = "../out/"+cancer_type+"/graph_nodes.txt"
    #dictionary maps ids to gene names
    id_to_gene = {}
    with open(filename) as ifile:
        lines = ifile.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[line[1]] = int(line[0]) - GENE_ID_OFFSET
    return id_to_gene

def load_edge_list():
    filename = "../out/"+cancer_type+"/graph_edges.txt"

    edge_list = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            edge_list.append([int(line[0]) - GENE_ID_OFFSET, int(line[1]) - GENE_ID_OFFSET])
    return edge_list

############## Matrix W
def normalizeW(A):
    if verbose:
        print("1 - Computing W")
    n = len(A)
    W = np.zeros((n, n))
    for j in tqdm(range(n)):
        d_j = float(A[:,j].sum())
        if(d_j == 0):
            continue
        for i in range(n):
            W[i,j] = A[i,j]/d_j
    return W

def computeW():
    gene_to_id = load_gene_to_id()
    edge_list = load_edge_list()

    N = len(gene_to_id)
    w = np.zeros((N, N))

    for edge in edge_list:
        w[edge[0]][edge[1]] = 1
        w[edge[1]][edge[0]] = 1
    return normalizeW(w)

############## Matrix H
def computeH():
    if verbose:
        print("2 - Computing H")
    gene_to_id = load_gene_to_id()
    lines = []
    with open("../out/"+cancer_type+"/graph_mut_freq.txt") as f:
        lines = f.readlines()
    N = len(gene_to_id)
    h = np.zeros((N, N))
    for line in lines:
        line = line.split()
        if line[0] not in gene_to_id.keys():
            continue
        gene_id = gene_to_id[line[0]]
        h[gene_id][gene_id] = line[1]
    return h

############## Matrix F
def computeF(W, beta):
    n = W.shape[0]
    if verbose:
        print("3 - Computing F")
    return beta*np.linalg.inv(np.eye(n)-(1-beta)*W)

global verbose
network_beta=0.4

path="../out/"
if len(sys.argv) > 1:
    cancer_type=sys.argv[1]
    v=sys.argv[2].split("-v-")[1]
    if v=='t':
        verbose=True
    else:
        verbose=False
GENE_ID_OFFSET=1

#compute W
W = computeW()
sp_w = sp.csc_matrix(W)
del W # deleting the variable to save memory
sp.save_npz(path + cancer_type+"/random_walk_w_matrix.npz", sp_w)

#compute H
H = computeH()
sp_h = sp.csc_matrix(H)
sp.save_npz(path + cancer_type+"/random_walk_h_matrix.npz", sp_h)
if verbose:
    print("H")
    print(H)

#compute F
F = computeF(sp_w , network_beta)
sp_f = sp.csc_matrix(F)
sp.save_npz(path + cancer_type+"/random_walk_f_matrix.npz", sp_f)
if verbose:
    print("F")
    print(F)

#compute E
if verbose:
    print("4 - Compute E")
E = F.dot(H)
sp_e = sp.csc_matrix(E)
sp.save_npz(path + cancer_type+"/random_walk_e_matrix.npz", sp_e)
if verbose:
    print("E")
    print(E)
