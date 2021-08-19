#***********************************
# Run the command 'python generate_outliers.py [cancer] -v-[t/f]' in the Terminal to start the
# the process of evaluating the betweenness values and generate the set of
# dysrgulated genes for each patient by fitting a truncated normal distribution
# This may take several minutes.
# When is done, the files of prepared data are saved in the subfolders of data/[cacner]/.
#***********************************

#generate outliers
import sys
import pandas as pd
import os
import copy
from tqdm import tqdm
import numpy as np
from scipy.stats import truncnorm



def calculate_BW_difference(genes,directory):
    patients = os.listdir(directory)
    if verbose:
        print("     Computing BW Differences: STARTED")

    # list of already treated samples
    patient_treated = []

    #dictionary to store the betweenness diff. values
    bw_diff_values = {}
    #first key is the patient ids
    bw_diff_values['patients']=[]

    #add genes in TCGA as keys to the dictionary
    for gene in genes:
        bw_diff_values[gene]=[]

    if verbose:
        pbar=tqdm(total=len(patients))
    for patient in patients:
        if verbose:
            pbar.update(1)
        if ".txt" not in patient:
            continue

        #retreive patient id from the file name
        patient_id=patient.split("_bw_")[1].split(".txt")[0]

        #check if patient already treated
        if patient_id in patient_treated:
            continue

        #append the patient name to the treated list
        patient_treated.append(patient_id)

        #create a dictionary to store patient p's betweenness values
        normal_bw_values = {g:-1 for g in genes}
        tumor_bw_values = {g:-1 for g in genes}



        #reading data from the first file (Normal graph)
        with open(directory + "/normal_bw_" + patient_id+".txt") as file1:
            for line in file1.readlines():
                if len(line.strip().split(':')[0]) == 0:
                    continue
                normal_bw_values[line.strip().split('\t:\t')[0]] = float(line.strip().split('\t:\t')[1])


        #reading data from the first file (Tumor graph)
        try:
            file2 = open(directory + "/tumor_bw_" + patient_id+".txt")
        except:

            print("ERROR missing tumor file for patient: ",patient_id)
            continue
        for line in file2.readlines():
            if len(line.strip().split(':')[0]) == 0:
                continue
            tumor_bw_values[line.strip().split('\t:\t')[0]] = float(line.strip().split('\t:\t')[1])


        #store the data in the created dictionary
        bw_diff_values['patients'].append(patient_id)
        # calculate betweenness difference
        for gene in genes:
            if normal_bw_values[gene]==-1 or tumor_bw_values[gene] ==-1:
                normal_bw_values[gene]=0
                tumor_bw_values[gene]=0
            bw_diff_values[gene].append(abs(normal_bw_values[gene] - tumor_bw_values[gene]))
    if verbose:
        print("     Computing BW Differences: DONE")
    return bw_diff_values

def load_gene_list(gene_list_input):
    human_genes_list=[]
    full_path="../"+input_directory+"/"+gene_list_input
    with open(full_path) as gene_list_file:
        for gene in gene_list_file.readlines():
            gene=gene.strip()
            human_genes_list.append(gene)
    gene_list_file.close()
    return human_genes_list

def calculate_mean_std_values(bw_diff_matrix):
    if verbose:
        print("\n     Computing Mean and Standard deviations: STARTED")
    bw_diff_matrix=bw_diff_matrix.set_index("patients")
    genes_standard_deviation_mean={}

    if verbose:
        pbar = tqdm(range(len([gene for gene in bw_diff_matrix])))

    for gene in bw_diff_matrix:
        if verbose:
            pbar.update(1)
        data_specefic_gene=bw_diff_matrix[gene].tolist()
        data_specefic_gene_sorted=sorted([x for x in data_specefic_gene])

        truncnorm_fit=truncnorm.fit(data_specefic_gene_sorted, fix_a=0, fix_b=max(data_specefic_gene), loc=0, scale=1)
        mean_value=truncnorm_fit[2]
        std_value=truncnorm_fit[3]

        genes_standard_deviation_mean[gene]=(std_value,mean_value)
    if verbose:
        print("\n     Computing  Mean and Standard deviations: DONE")
    return genes_standard_deviation_mean

def generate_outliers_matrix(bw_diff_matrix,genes_standard_deviation_mean):
    if verbose:
        print("\n     Generating Outliers matrix: STARTED")
    patients=[pat for pat in bw_diff_matrix["patients"]]
    bw_diff_matrix=bw_diff_matrix.set_index("patients")

    patient_outliers_dic={}
    for gene in bw_diff_matrix:
        if gene == "patients":
            continue
        patient_outliers_dic[gene]={}
        for patient in patients:
            if bw_diff_matrix[gene][patient]<=(0.5*genes_standard_deviation_mean[gene][0])+genes_standard_deviation_mean[gene][1]:
                patient_outliers_dic[gene][patient]=False
            else:
                patient_outliers_dic[gene][patient]=True
    if verbose:
        print("\n     Generating Outliers matrix: DONE")
    return patient_outliers_dic


def main():
    global input_directory,verbose
    input_directory="data"
    if len(sys.argv) > 1:
        cancer_type=sys.argv[1]
        v=sys.argv[2].split("-v-")[1]
        if v=='t':
            verbose=True
        else:
            verbose=False

    btw_results_path="Betweenness"
    gene_list_input="human_genes"

    #list files
    file_list=[]
    for file in os.listdir("../out/"+cancer_type+"/"+btw_results_path):
        file_list.append("../out/"+btw_results_path+"/"+file)

    genes_list=load_gene_list(gene_list_input)
    bw_diff=calculate_BW_difference(genes_list,"../out/"+cancer_type+"/"+btw_results_path)

    #convert dictionary to pandas matrix
    bw_diff_matrix = pd.DataFrame(bw_diff)

    #calculate mean and standard deviation
    genes_standard_deviation_mean=calculate_mean_std_values(bw_diff_matrix)

    #generate outlier list
    list_outliers=generate_outliers_matrix(bw_diff_matrix,genes_standard_deviation_mean)

    outliers_matrix=pd.DataFrame.from_dict(list_outliers)
    outliers_matrix.to_csv("../data/"+cancer_type+"/outliers_data.csv")




if __name__ == "__main__":
    main()
