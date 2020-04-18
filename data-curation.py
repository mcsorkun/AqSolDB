#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:31:39 2019

@author: Murat Cihan Sorkun
"""

import pandas as pd
import merge
import descriptors


data_a = pd.read_csv('data/dataset-A.csv', header=0) 
data_b = pd.read_csv('data/dataset-B.csv', header=0) 
data_c = pd.read_csv('data/dataset-C.csv', header=0) 
data_d = pd.read_csv('data/dataset-D.csv', header=0) 
data_e = pd.read_csv('data/dataset-E.csv', header=0) 
data_f = pd.read_csv('data/dataset-F.csv', header=0) 
data_g = pd.read_csv('data/dataset-G.csv', header=0) 
data_h = pd.read_csv('data/dataset-H.csv', header=0) 
data_i = pd.read_csv('data/dataset-I.csv', header=0) 


# first merge data directly    
data_merged= pd.concat([data_a,data_b,data_c,data_d,data_e,data_f,data_g,data_h,data_i])
data_merged=data_merged.reset_index(drop=True)

# curate the data    
data_curated = merge.curate(data_merged)


# generate the descriptors    
descriptors = descriptors.generate(data_curated['SMILES'].values)

# combine descriptors with 
data_curated=pd.concat([data_curated, descriptors], axis=1)

# write into CSV file
data_curated.to_csv("results/data_curated.csv", sep=',', encoding='utf-8', index=False)


