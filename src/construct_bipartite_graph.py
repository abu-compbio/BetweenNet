import networkx as nx
import pandas as pd
import os,sys
import operator
from tqdm import tqdm
from networkx.algorithms import bipartite

class Data_Preprocessing:
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

    def load_mutations_matrix(file):
        mutation_data=pd.read_csv("../"+input_directory+"/"+file+".csv",sep="\t")
        patients_mut=[pat for pat in mutation_data["Genes"]]
        patients_vs_mutations={}
        gene_vs_patients={}
        all_mutations=set()
        i=0
        for pat in patients:
            if pat in patients_mut:
                patient_row = mutation_data.loc[mutation_data['Genes'] == pat]
                idx=list(patient_row.index)[0]
                mutations_for_spec_patient=[cols for cols in patient_row if str(patient_row.get_value(idx, cols, takeable=False)) == "1"]
                for g in mutations_for_spec_patient:
                    all_mutations.add(g)
                    if g not in gene_vs_patients:
                        gene_vs_patients[g]=set()
                        gene_vs_patients[g].add(pat)
                    else:
                        gene_vs_patients[g].add(pat)
                patients_vs_mutations[pat]=mutations_for_spec_patient
        return all_mutations,patients_vs_mutations,gene_vs_patients


    def load_outliers_matrix(file):
        outliers_data=pd.read_csv("../"+input_directory+"/"+file+".csv")
        outliers_data=outliers_data.set_index(outliers_data["Genes"])
        patients_in_outliers_data=[col for col in outliers_data.index]
        outliers_data=outliers_data.transpose()
        all_outliers=set()
        all_outliers_=set()
        outliers_vs_patients={}
        gene_vs_patients_outliers={}

        for pat in patients:
            gene_vs_patients_outliers[pat]=[]
            if pat in patients_in_outliers_data:
                outliers_for_spec_patient=outliers_data.index[outliers_data[pat] == True].tolist()
                outliers_vs_patients[pat]=outliers_for_spec_patient
                for g in outliers_for_spec_patient:
                    all_outliers_.add(g)
                    if g not in gene_vs_patients_outliers:
                        gene_vs_patients_outliers[g]=[pat]
                    else:
                        gene_vs_patients_outliers[g].append([pat])


        outliers_vs_patients_={}
        i=0
        for pat in outliers_vs_patients:
            outliers_vs_patients_[pat]=[]
            for outlier in outliers_vs_patients[pat]:
                new_name=outlier+str("_")+pat
                all_outliers.add(new_name)
                outliers_vs_patients_[pat].append(new_name)
        return all_outliers,outliers_vs_patients,gene_vs_patients_outliers


class Graph:
    def construct_bipartite_graph(patient_mutations,patient_outliers):
    
        G = nx.Graph()
        G.add_nodes_from(mutated_genes, bipartite=0)
        G.add_nodes_from(outlier_genes, bipartite=1)
        pbar = tqdm(range(len(patients)))
        k=0
        for pat in patients:
            pbar.update(1)
            if pat in patient_mutations and pat in patient_outliers:

                for outlier in patient_outliers[pat]:
                    if outlier in inf_graph:
                        for mutation in patient_mutations[pat]:
                            if mutation in inf_graph:
                                
                                if mutation in inf_graph[outlier]:
                                    outlier_=outlier+str("_")+pat
                                    G.add_edges_from([(mutation,outlier_ )])


        bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}
        bipartite0_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}

        #deleting nodes with 0 degree
        for node in bipartite1_nodes:
            if int(G.degree(node))==0:
                G.remove_node(node)
        for node in bipartite0_nodes:
            if int(G.degree(node))==0:
                G.remove_node(node)
        print("Graph successfully generated with a size of: ",len(G.nodes())," nodes, and ",len(G.edges())," edges")


        return G

    def map_nodes2ids(self,G):
        edges=G.edges()
        nodes_set=set()
        for edge in edges:
            nodeA=edge[0]
            nodeB=edge[1]
            nodes_set.add(nodeA)
            nodes_set.add(nodeB)
        gene_id_map={}
        id=1
        for gene in bipartite0_nodes:
            gene_id_map[gene]=id
            id+=1
        for gene in bipartite1_nodes:
            gene_id_map[gene]=id
            id+=1







def main():
    print("3 - Construct Bipartite Graph")
    global patients,mutated_genes,outlier_genes,inf_graph
    global input_directory
    # count the arguments
    arguments = len(sys.argv) - 1
    if arguments < 4:
        print("________________________________________________________________________________________")
        print('Please run the code using the following command line and Arguments:  ')
        print("python construct_bipartite_graph.py [Input Directory] [Patient Ids] [Influence Matrix] [Mutations Matrix] [Outliers Matrix]")
        print("________________________________________________________________________________________\n\n")
        sys.exit("")


    input_directory=sys.argv[1]
    patients_ids=sys.argv[2]
    influence_matrix=sys.argv[3]
    mutations_matrix=sys.argv[4]
    outliers_matrix=sys.argv[5]


    with open("../"+input_directory+"/"+patients_ids+".txt") as ifile:
        patients=[pat.strip() for pat in ifile.readlines()]


    #load data
    inf_graph=Data_Preprocessing.load_inf_graph(influence_matrix)
    mutated_genes,patient_mutations,mutation_patients=Data_Preprocessing.load_mutations_matrix(mutations_matrix)

    outlier_genes,patient_outliers,outlier_patients=Data_Preprocessing.load_outliers_matrix(outliers_matrix)
    G=Graph.construct_bipartite_graph(patient_mutations,patient_outliers)
    nx.write_gml(G, "../out/bipartite_graph.gml")


    bipartite1_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}
    bipartite0_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}


    #prepare data for random-walk
    #nodes=map_nodes2ids(G)
    gene_id_map={}
    id=1
    for gene in bipartite0_nodes:
        gene_id_map[gene]=id
        id+=1
    for gene in bipartite1_nodes:
        gene_id_map[gene]=id
        id+=1

    #genereate index of gene for random walk algorithm
    output_index_file=open("../out/graph_nodes.txt","w")
    for gene in gene_id_map:
        row=str(gene_id_map[gene])+str("\t")+gene+str("\t")+str(0)+str("\n")
        output_index_file.write(row)
    output_index_file.close()
    #genereate index of edges for random walk algorithm
    output_edge_file=open("../out/graph_edges.txt","w")
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

    output_mut_freqe=open("../out/graph_mut_freq.txt","w")
    visited=[]
    for gene in bipartite1_nodes:
        gene_=gene.split("_")[0]
        if gene in gene_id_map:
            mut_occ=len(outlier_patients[gene_])
            row = gene+str("\t")+str(mut_occ/len(patients))+str("\n")
            output_mut_freqe.write(row)

    for gene in bipartite0_nodes:

        if gene in gene_id_map:
            mut_occ=len(mutation_patients[gene])
            row = gene+str("\t")+str(mut_occ/len(patients))+str("\n")
            output_mut_freqe.write(row)
        else:
            continue

    output_mut_freqe.close()





    print("*************************************************")
    print("*                Successfuly Done               *")
    print("*************************************************")


if __name__ == "__main__":
    main()
