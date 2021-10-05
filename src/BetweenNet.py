#***********************************
# Run the command 'python BetweenNet.py [cacner] -v-[t/f]' in the Terminal to start the
# the process of ranking mutated genes in the bipartite graph
#
# This may take several minutes.
# When is done, the files of prepared data are saved in the subfolders of data/[cacner]/BetweenNet.txt.
#***********************************


import scipy.sparse as sp
import numpy as np
import networkx as nx

def load_matrix():
    e=sp.load_npz("../out/"+cancer_type+"/random_walk_e_matrix.npz")
    return e

def load_gene_names():
    with open("../out/"+cancer_type+"/graph_nodes.txt","r") as findex:
        gene_index={}
        for line in findex.readlines():
            line=line.strip().split("\t")
            gene_index[line[1]]=int(line[0])-1
    return gene_index

def evaluate_matrix(random_walk_matrix):

    #sum weights of the incoming edges (sum of entire row of the corresponding gene)
    final_score=np.sum(random_walk_matrix, axis=1)

    #map scores to gene names
    score2id={}
    gene_index=load_gene_names()
    for gene in gene_index:
        score2id[gene]=final_score[int(gene_index[gene])]
    global rw_scores

    #select only scores of the mutated genes
    rw_scores={}
    for g in sorted(score2id.items(), key=lambda x: x[1],reverse=True):
        if len(g[0].split("_"))==1:
            rw_scores[g[0]] = float(g[1])

    return rw_scores

def rank_mutated_genes(rw_scores):

    G=nx.read_gml("../out/"+cancer_type+"/bipartite_graph.gml")
    max_rw_score=rw_scores[max(rw_scores, key=rw_scores.get)]
    with open("../out/"+cancer_type+"/BetweenNet.txt","w") as ofile:


        for alpha in [0.5]:

            drivers=[]
            bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}

            while (len(bipartite1_nodes)!=0):
                bipartite0_mutations= {n for n, d in G.nodes(data=True) if d['bipartite']==0}

                m_scores={}
                # compute max degree
                max_degree=max([G.degree(mutated_gene) for mutated_gene in  bipartite0_mutations])

                for mutated_gene in bipartite0_mutations:
                    #compute the degree of the selected mutated gene
                    mutated_gene_degree=G.degree(mutated_gene)
                    # get the edge weight of the mutated gene m after random walk
                    rw_score=rw_scores[mutated_gene]/max_rw_score
                    #calculate the score of the mutated gene
                    m_scores[mutated_gene]=alpha*(mutated_gene_degree/max_degree)+(1.0-alpha)*rw_score

                #sort mutated genes based on their scores
                sorted_mutation_list= [k for k, v in sorted(m_scores.items(), key=lambda kv: (-kv[1], kv[0]))]



                ##Get the first gene -mutated- to be removed
                node_to_remove=sorted_mutation_list[0]
                #remove the mutated gene with highest score and its neighbors
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
    if len(sys.argv) > 1:
        cancer_type=sys.argv[1]
    random_walk_matrix=load_matrix()
    rw_scores = evaluate_matrix(random_walk_matrix)
    
    _mutated_genes(rw_scores)


if __name__ == "__main__":
    main()
