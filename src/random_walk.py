import scipy.sparse as sp
import numpy as np
import argparse
from scipy.sparse.linalg import inv
from tqdm import tqdm
GENE_ID_OFFSET=1

def load_gene_to_id():
    filename = "../out/graph_nodes.txt"
    id_to_gene = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[line[1]] = int(line[0]) - GENE_ID_OFFSET
    return id_to_gene


def load_edge_list():
    filename = "../out/graph_edges.txt"

    edge_list = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            edge_list.append([int(line[0]) - GENE_ID_OFFSET, int(line[1]) - GENE_ID_OFFSET])
    return edge_list

def normalizeW(A):
    print("1 - Normalize W")
    n = len(A)
    W = np.zeros((n, n))
    for j in tqdm(range(n)):
        d_j = float(A[:,j].sum())
        if(d_j == 0):
            continue
        W[j,:]/= d_j
    return W

def computeF(W, beta):
    n = W.shape[0]
    print("2 - Computing F")
    return beta*np.linalg.inv(np.eye(n)-(1-beta)*W)


def computeH():
    print("3 - Computing H")
    gene_to_id = load_gene_to_id()
    lines = []
    with open("../out/graph_mut_freq.txt") as f:
        lines = f.readlines()
    N = len(gene_to_id)
    h = np.zeros((N, N))
    yes = 0
    no = 0

    for line in lines:
        line = line.split()
        if line[0] not in gene_to_id.keys():
            no += 1
            continue
        yes += 1
        gene_id = gene_to_id[line[0]]
        h[gene_id][gene_id] = line[1]
    print(yes)
    print(no)
    return h

def computeW():
    gene_to_id = load_gene_to_id()
    edge_list = load_edge_list()

    N = len(gene_to_id)
    w = np.zeros((N, N))

    for edge in edge_list:
        w[edge[0]][edge[1]] = 1
        w[edge[1]][edge[0]] = 1
    return normalizeW(w)


print()
network_beta=0.4
path="../out/"
W = computeW()
sp_w = sp.csc_matrix(W)
del W
sp.save_npz(path + "random_walk_w_matrix.npz", sp_w)

H = computeH()
sp_h = sp.csc_matrix(H)
sp.save_npz(path + "random_walk_h_matrix.npz", sp_h)
print("H")
print(H)

F = computeF(sp_w , network_beta)
sp_f = sp.csc_matrix(F)
sp.save_npz(path + "random_walk_f_matrix.npz", sp_f)
print("F")
print(F)




print("4 - Compute E")
E = F.dot(H)
sp_e = sp.csc_matrix(E)
sp.save_npz(path + "random_walk_e_matrix.npz", sp_e)
print("E")
print(E)
