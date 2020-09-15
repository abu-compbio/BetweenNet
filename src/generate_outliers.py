#generate outliers
import sys
import pandas as pd
import os
import copy
from tqdm import tqdm
import numpy as np
from scipy.stats import truncnorm
from scipy.optimize import fmin_slsqp
np.seterr(divide = 'ignore')


def tumorOrNormal(s):
    if s == "tumor":
        return "normal_bw_"
    else:
        return "tumor_bw_"


def calculate_BW_difference(genes,directory):
    patients = os.listdir(directory)
    print("calculate BW differences started")


    patient_treated = []

    bw_diff_values = {}
    bw_diff_values['patients']=[]
    for gene in genes:
        bw_diff_values[gene]=[]

    count_patient=0
    for patient in patients:

        if "Icon" in patient or "DS_Store" in patient:
            continue
        #check if patient already treated
        print(patient)
        if patient.split("_bw_")[1].split(".tgenet")[0] in patient_treated:
            continue

        normal_bw_values = copy.deepcopy(bw_diff_values)
        tumor_bw_values = copy.deepcopy(bw_diff_values)
        #append the patient name to the treated list
        patient_treated.append(patient.split("_bw_")[1].split(".tgenet")[0])

        #reading data from the first file
        with open(directory + "/" + patient) as file1:
            list1 = []
            for line in file1.readlines():

                if len(line.strip().split(':')[0]) == 0:
                    continue
                normal_bw_values[line.strip().split('\t:\t')[0]] = float(line.strip().split('\t:\t')[1])


        #reading data from the second file by checking if the first one was normal or tumor
        try:
            file2 = open(directory + "/" + tumorOrNormal(patient.split('_bw_')[0])+patient.split('_bw_')[1])
        except:
            continue
        list2 = []
        for line in file2.readlines():
            if len(line.strip().split(':')[0]) == 0:
                continue
            tumor_bw_values[line.strip().split('\t:\t')[0]] = float(line.strip().split('\t:\t')[1])



        bw_diff_values['patients'].append(patient.split('_bw_')[1].split('.')[0])

        for gene in genes:


            if (type(normal_bw_values[gene]) is list) or type(tumor_bw_values[gene]) is list:
                normal_bw_values[gene]=0
                tumor_bw_values[gene]=0



            bw_diff_values[gene].append(abs(normal_bw_values[gene] - tumor_bw_values[gene]))

        count_patient+=1






    print("calculate BwC differences finished")
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

    print("\nCalculate Meand and Standard deviations: started")
    bw_diff_matrix=bw_diff_matrix.set_index("patients")
    genes_standard_deviation_mean={}

    pbar = tqdm(range(len([gene for gene in bw_diff_matrix])))

    for gene in bw_diff_matrix:
        pbar.update(1)
        data_specefic_gene=bw_diff_matrix[gene].tolist()
        data_specefic_gene_sorted=sorted([x for x in data_specefic_gene])

        truncnorm_fit=truncnorm.fit(data_specefic_gene_sorted, fix_a=0, fix_b=max(data_specefic_gene), loc=0, scale=1)
        std_value=truncnorm_fit[3]
        mean_value=truncnorm_fit[2]

        genes_standard_deviation_mean[gene]=(std_value,mean_value)

    print("\nCalculate Meand and Standard deviations: done")
    print("*************************************************\n\n\n")
    return genes_standard_deviation_mean

def generate_outliers_matrix(bw_diff_matrix,genes_standard_deviation_mean):
    print("\nGenerating outliers matrix: started")
    patients=[pat for pat in bw_diff_matrix["patients"]]
    bw_diff_matrix=bw_diff_matrix.set_index("patients")



    print(bw_diff_matrix.head())
    patient_outliers_dic={}
    for gene in bw_diff_matrix:
        if gene == "patients":
            continue
        patient_outliers_dic[gene]={}
        for patient in patients:
            #if gene in genes_standard_deviation_mean:
            if bw_diff_matrix[gene][patient]<=(0.5*genes_standard_deviation_mean[gene][0])+genes_standard_deviation_mean[gene][1]:
                patient_outliers_dic[gene][patient]=False
            else:
                patient_outliers_dic[gene][patient]=True
    print("\nGenerating outliers matrix: done")
    print("*************************************************\n\n\n")
    return patient_outliers_dic


def main():
    global input_directory
    # count the arguments
    arguments = len(sys.argv) - 1
    if arguments < 3:
        print("________________________________________________________________________________________")
        print('Please run the code using the following command line and Arguments:  ')
        print("python generate_outliers.py [Input directory] [Betweenness Results Path] [Genes list]")
        print("________________________________________________________________________________________\n\n")
        sys.exit("")


    input_directory=sys.argv[1]
    btw_results_path=sys.argv[2]
    gene_list_input=sys.argv[3]

    #list files
    file_list=[]
    for file in os.listdir("../out/"+btw_results_path):
        file_list.append("../out/"+btw_results_path+"/"+file)

    genes_list=load_gene_list(gene_list_input)
    bw_diff=calculate_BW_difference(genes_list,"../out/"+btw_results_path)

    #convert dictionary to pandas matrix
    bw_diff_matrix = pd.DataFrame(bw_diff)

    #calculate mean and standard deviation
    genes_standard_deviation_mean=calculate_mean_std_values(bw_diff_matrix)

    #generate outlier list
    list_outliers=generate_outliers_matrix(bw_diff_matrix,genes_standard_deviation_mean)

    outliers_matrix=pd.DataFrame.from_dict(list_outliers)
    outliers_matrix.to_csv("../data/OutliersMatrix.csv")

    print("*************************************************")
    print("*                Successfuly Done               *")
    print("*************************************************")


if __name__ == "__main__":
    main()