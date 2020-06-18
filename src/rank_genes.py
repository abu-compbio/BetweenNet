@import scipy.sparse as sp
import numpy as np

def load_matrix():
    e=sp.load_npz("../out/random_walk_e_matrix.npz")
    return e


def load_gene_names():
    with open("../out/BRCA_index.txt","r") as findex:
        gene_index={}
        for line in findex.readlines():
            line=line.strip().split("\t")
            gene_index[line[1]]=int(line[0])-1
    return gene_index

def evaluate_matrix():
    final_score=np.sum(random_walk_matrix, axis=1)
    print(final_score)
    score2id={}
    gene_index=load_gene_names()
    for gene in gene_index:
        score2id[gene]=final_score[int(gene_index[gene])]

    rank=[]
    for g in sorted(score2id.items(), key=lambda x: x[1],reverse=True):
        if len(g[0].split("_"))==1:
            rank.append(g[0])

    return rank


def logger(rank):


    with open("../out/ranking.txt","w") as ofile:
        for gene in rank:

            ofile.write(gene)
            ofile.write(str("\n"))

def main():
    global random_walk_matrix
    random_walk_matrix=load_matrix()
    print(random_walk_matrix.shape)
    logger(evaluate_matrix())



if __name__ == "__main__":
    main()
