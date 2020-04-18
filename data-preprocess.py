#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:30:41 2019

@author: Murat Cihan Sorkun
"""

import numpy as np
import pandas as pd
import preprocess


def main():
    
    data_path_a="data/raw-dataset-A.csv"
    data_path_h="data/raw-dataset-H.csv"
    
#raw dataset various solubility metrics (g/L, mg/L..) with Name and CAS Number 
    preprocess_dataset_a(data_path_a)
    
#raw dataset has solubility values(LogS) with SLN representations    
    preprocess_dataset_h(data_path_h)





def preprocess_dataset_a(data_path):

    print("\n\n**\nPreprocess of raw-dataset-A is started\n**")
    data = pd.read_csv(data_path, header=0)     
    smiles_list=preprocess.collect_smiles_from_web(data['Name'].values,data['CAS'].values,verbose=0)                
    
    inchi_list,inchikey_list=preprocess.smiles_to_inchi_inchikey(smiles_list,verbose=0)
    id_list=preprocess.generate_id_list("A",len(smiles_list))
    name_list=preprocess.collect_names_from_web(inchikey_list,smiles_list,verbose=0)
    logs_list=preprocess.convert_to_logs(smiles_list,data['Metric'].values,data['OriginalSol'].values,verbose=0)
    prediction_list=preprocess.collect_predictions_from_web(smiles_list,verbose=0)
    
    dataset_a_df = pd.DataFrame(np.column_stack([ id_list, name_list, inchi_list, inchikey_list, smiles_list, logs_list, prediction_list]), 
                                   columns=[ 'ID', 'Name', 'InChI', 'InChIKey', 'SMILES', 'Solubility', 'Prediction'])
    
    
    #filter dataset by removing missing information. (for strings: "XXX" and for numeric: "999")
    dataset_a_df_clean = preprocess.clean(dataset_a_df) 
    
    
    #update ID after filtering
    id_list=preprocess.generate_id_list("A",len(dataset_a_df_clean.index))
    id_clean_df = pd.DataFrame({'ID': id_list})
    dataset_a_df_clean.update(id_clean_df)
    
    #write dataset into CSV file
    dataset_a_df_clean.to_csv('../results/dataset-A.csv', index=False)
    print("**\nPreprocessed dataset-A is written into dataset-A.csv\n**")
    
    return



def preprocess_dataset_h(data_path):

    
    print("\n\n**\nPreprocess of raw-dataset-H is started\n**")
    
    
    data = pd.read_csv(data_path, header=0)     
    smiles_list=preprocess.sln_to_smiles(data['SLN'].values,verbose=0)              
    
    inchi_list,inchikey_list=preprocess.smiles_to_inchi_inchikey(smiles_list,verbose=0)
    id_list=preprocess.generate_id_list("H",len(smiles_list))
    name_list=preprocess.collect_names_from_web(inchikey_list,smiles_list,verbose=0)
    logs_list=data['Solubility'].values
    prediction_list=preprocess.collect_predictions_from_web(smiles_list,verbose=0)
    
    dataset_h_df = pd.DataFrame(np.column_stack([ id_list, name_list, inchi_list, inchikey_list, smiles_list, logs_list, prediction_list]), 
                                   columns=[ 'ID', 'Name', 'InChI', 'InChIKey', 'SMILES', 'Solubility', 'Prediction'])
    
    
    #filter dataset by removing missing information. (for strings: "XXX" and for numeric: "999")
    dataset_h_df_clean = preprocess.clean(dataset_h_df) 
    
    
    #update ID after filtering
    id_list=preprocess.generate_id_list("H",len(dataset_h_df_clean.index))
    id_clean_df = pd.DataFrame({'ID': id_list})
    dataset_h_df_clean.update(id_clean_df)
    
    #write dataset into CSV file
    dataset_h_df_clean.to_csv('../results/dataset-H.csv', index=False)
    print("**\nPreprocessed dataset-H is written into dataset-H.csv\n**")
    
    return


if __name__== "__main__":
    main()