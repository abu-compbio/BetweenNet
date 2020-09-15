@import scipy.sparse as sp
import numpy as np
import networkx as nx

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
    global rw_scores
    rw_scores={}
    for g in sorted(score2id.items(), key=lambda x: x[1],reverse=True):
        if len(g[0].split("_"))==1:
            rw_scores[g[0]] = float(g[1])

    return rw_scores
def get_non_null_degrees(G):
    degrees = [(val,node) for (node, val) in G.degree()]

    null_degree=[]
    non_null_degree=[]
    for d in degrees:
        if d[0]==0:

            null_degree.append(d)
        else:
            non_null_degree.append(d)

    mutations=[out[1] for out in non_null_degree if len(out[1].split("_"))==1]

    return mutations



def logger(rank):

    G=nx.read_gml("betweennet_graph.gml")
    max_rw_score=rw_scores[max(mut_rw_score, key=mut_rw_score.get)]
    with open("../out/ranking.txt","w") as ofile:


        for alpha in [0.5]:


            drivers=[]
            bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}

            while (len(bipartite1_nodes)!=0):
                bipartite0_mutations= {n for n, d in G.nodes(data=True) if d['bipartite']==0}
                i+=
                all_mutation_degrees={}
                max_degree=-1
                for mutated_gene in bipartite0_mutations:
                    mutated_gene_degree=G.degree(mutated_gene)
                    if max_degree< mutated_gene_degree:
                        max_degree=mutated_gene_degree

                for mutated_gene in bipartite0_mutations:

                    mutated_gene_degree=G.degree(mutated_gene)
                    rw_score=rw_scores[mutated_gene]/max_rw_score

                    all_mutation_degrees[mutated_gene]=alpha*(mutated_gene_degree/max_degree)+(1-alpha)*rw_scores


                sorted_mutation_list= [k for k, v in sorted(all_mutation_degrees.items(), key=lambda kv: (-kv[1], kv[0]))]



                ##Get the first gene -mutated- to be removed
                node_to_remove=sorted_mutation_list[0]

                drivers.append(node_to_remove)
                outliers_to_delete= list(G.neighbors(node_to_remove))
                for outlier_to_rmv in outliers_to_delete:
                    G.remove_node(outlier_to_rmv)
                G.remove_node(node_to_remove)

                bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}

                ofile.write(node_to_remove)
                ofile.write(str("\n"))




def main():
    global random_walk_matrix
    random_walk_matrix=load_matrix()
    print(random_walk_matrix.shape)
    logger(evaluate_matrix())



if __name__ == "__main__":
    main()
