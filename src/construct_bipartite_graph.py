#***********************************
# Run the command 'python construct_bipartite_graph.py' in the Terminal to start the
# the process of constructing the bipartite network and generate input data for random walk process
# This may take several minutes.
# When is done, the files of prepared data are saved in the subfolders of out/[cacner]/.
#***********************************


import networkx as nx
import pandas as pd
import os,sys
import operator
from tqdm import tqdm
from networkx.algorithms import bipartite

class DataPreprocessing:
    @staticmethod
    def load_inf_graph(file):
        inf_graph={}
        with open("../"+input_directory+"/"+file+".txt") as ifile:
            for line in ifile.readlines():
                line=line.strip().split("\t")
                if line[0] not in inf_graph:
                    inf_graph[line[0]]=[line[1]]
                else:
                    inf_graph[line[0]].append(line[1])

                if line[1] not in inf_graph:
                    inf_graph[line[1]]=[line[0]]
                else:
                    inf_graph[line[1]].append(line[0])
        return inf_graph

    @staticmethod
    def load_mutations(file,cancer_type):

        #load mutation data
        with open("../"+input_directory+"/"+cancer_type+"/"+file+".txt","r") as mutation_data:

            #dictionary to store patient and its set of mutated genes
            patients_vs_mutations={}

            #dictionary to store mutated genes and set of samples, it is used to calculate mutation frequancies
            mutations_vs_patients={}

            # set of all mutations to construc bipartite graph
            all_mutations=set()

            for line in mutation_data.readlines():
                patient_id=line.strip().split("\t")[0]
                patients.append(patient_id)

                #list of mutated genes2
                genes=line.strip().split("\t")[1].split(",")

                for gene in genes:
                    all_mutations.add(gene)
                    if gene not in mutations_vs_patients:
                        mutations_vs_patients[gene]=set()
                        mutations_vs_patients[gene].add(patient_id)
                    else:
                        mutations_vs_patients[gene].add(patient_id)
                patients_vs_mutations[patient_id]=genes
        return all_mutations,patients_vs_mutations,mutations_vs_patients,patients

    @staticmethod
    def load_outliers(file,cancer_type):

        #load outliers matrix
        outliers_matrix=pd.read_csv("../"+input_directory+"/"+cancer_type+"/"+file+".csv",index_col=0)
        outliers_matrix=outliers_matrix.transpose()

        # patients in outliers matrix
        p_out_matrix=list(outliers_matrix.columns)

        #set of outliers to be used in constructing the bipartite grpah
        all_outliers=set()

        patient_vs_outliers={}
        outliers_vs_patient={}

        for p in patients:
            outliers_vs_patient[p]=[]
            if p in p_out_matrix:

                #select only a set of dysregulated genes in p
                outliers_for_spec_patient=outliers_matrix.index[outliers_matrix[p] == True].tolist()
                patient_vs_outliers[p]=outliers_for_spec_patient
                for g in outliers_for_spec_patient:
                    if g not in outliers_vs_patient:
                        outliers_vs_patient[g]=[p]
                    else:
                        outliers_vs_patient[g].append([p])

        for p in patient_vs_outliers:
            for outlier in patient_vs_outliers[p]:
                new_name=outlier+str("_")+p
                all_outliers.add(new_name)
        return all_outliers,patient_vs_outliers,outliers_vs_patient

class Graph:
    @staticmethod
    def construct_bipartite_graph(patient_mutations,patient_outliers):

        G = nx.Graph()

        #set of mutated geens (partition 1)
        G.add_nodes_from(mutated_genes, bipartite=0)

        #set of mutated geens (partition 2)
        G.add_nodes_from(outlier_genes, bipartite=1)
        if verbose:
            pbar = tqdm(range(len(patients)))
        k=0

        #for each p in P
        for p in patients:
            if verbose:
                pbar.update(1)
            #check that p exists in both mutation and dysrgulated genes data
            if p in patient_mutations and p in patient_outliers:
                # for o dysrgulated in p
                for outlier in patient_outliers[p]:
                    #check if the node exist in the PPI
                    if outlier in inf_graph:
                        # if the outlier gene is also mutated then skip.
                        if outlier not in patient_mutations[p]:
                            #for m mutated in p
                            for mutation in patient_mutations[p]:
                                if mutation in inf_graph:
                                    if mutation in inf_graph[outlier]:
                                        outlier_=outlier+str("_")+p
                                        G.add_edges_from([(mutation,outlier_ )])

            else:
                print("Errorr Missing Patient: ",p)
        bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}
        bipartite0_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}

        #deleting nodes with 0 degree
        for node in bipartite1_nodes:
            if int(G.degree(node))==0:
                G.remove_node(node)
        for node in bipartite0_nodes:
            if int(G.degree(node))==0:
                G.remove_node(node)
        if verbose:
            print("Graph successfully generated with a size of: ",len(G.nodes())," nodes, and ",len(G.edges())," edges")


        return G

def main():
    global patients,mutated_genes,outlier_genes,inf_graph,input_directory,verbose,patients
    patients=[]

    # make sure that the files declared below are correct
    input_directory="data"
    if len(sys.argv) > 1:
        cancer_type=sys.argv[1]
        v=sys.argv[2].split("-v-")[1]
        if v=='t':
            verbose=True
        else:
            verbose=False


    influence_matrix="IntAct_network"
    mutations_matrix="mutation_data"
    outliers_matrix="outliers_data"


    #load data
    inf_graph=DataPreprocessing.load_inf_graph(influence_matrix)
    mutated_genes,patient_mutations,mutation_patients,patients=DataPreprocessing.load_mutations(mutations_matrix,cancer_type)
    outlier_genes,patient_outliers,outlier_patients=DataPreprocessing.load_outliers(outliers_matrix,cancer_type)

    #construct the bipartite graph
    G=Graph.construct_bipartite_graph(patient_mutations,patient_outliers)
    nx.write_gml(G, "../out/"+cancer_type+"/bipartite_graph.gml")


    #generate input for random walk step
    bipartite0_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}
    bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}


    #genereate index of gene for random walk algorithm
    gene_id_map={}
    id=1
    output_index_file=open("../out/"+cancer_type+"/graph_nodes.txt","w")
    for gene in bipartite0_nodes:
        gene_id_map[gene]=id
        row=str(gene_id_map[gene])+str("\t")+gene+str("\t")+str(0)+str("\n")
        output_index_file.write(row)
        id+=1
    for gene in bipartite1_nodes:
        gene_id_map[gene]=id
        row=str(gene_id_map[gene])+str("\t")+gene+str("\t")+str(0)+str("\n")
        output_index_file.write(row)
        id+=1
    output_index_file.close()


    #genereate index of edges for random walk algorithm
    output_edge_file=open("../out/"+cancer_type+"/graph_edges.txt","w")
    edges_set=set()
    edges=G.edges()
    for edge in edges:
        nodeA=edge[0]
        nodeB=edge[1]
        if (nodeA,nodeB) not in edges_set:
            edges_set.add((nodeA,nodeB))
            row=str(gene_id_map[nodeA])+str("\t")+str(gene_id_map[nodeB])+str("\t")+str(1)+str("\n")
            output_edge_file.write(row)
    output_edge_file.close()


    #compute mutation frequencies
    output_mut_freqe=open("../out/"+cancer_type+"/graph_mut_freq.txt","w")
    visited=[]

    for gene in bipartite0_nodes:
        if gene in gene_id_map:
            mut_occ=len(mutation_patients[gene])/len(patients)
            row = gene+str("\t")+str(mut_occ)+str("\n")
            output_mut_freqe.write(row)
        else:
            print("Error, Missing mutated gene in the index data")
            continue

    for gene in bipartite1_nodes:
        gene_=gene.split("_")[0]
        if gene in gene_id_map:
            mut_occ=float(len(outlier_patients[gene_])/len(patients))
            row = gene+str("\t")+str(mut_occ)+str("\n")
            output_mut_freqe.write(row)
        else:
            print("Error, Missing dysregulated gene in the index data")
            continue
    output_mut_freqe.close()




if __name__ == "__main__":
    main()
