#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 13:52:29 2019

@author: Murat Cihan Sorkun
"""

import numpy as np
import pandas as pd
import math


# Selects most reliable values among the multiple occurences
def curate(data_merged):

    smiles_list= data_merged['SMILES'].tolist()
    solubility_list= data_merged['Solubility'].tolist()
    id_list= data_merged['ID'].tolist()
    inchi_list= data_merged['InChI'].tolist()
    name_list= data_merged['Name'].tolist()
    prediction_list= data_merged['Prediction'].tolist()    


# define variables and assign default values    
    dif_val_list = []
    dif_sol_list = []
    
    same_value_counter=0  #same molecules with same values
    different_value_counter_2times=0 #same molecules with different values (2 occurences)
    different_value_counter_mutiple=0 #same molecules with different values (more than 2 occurences)
    
    ocurrence_count=[-999]*len(id_list)
    SD=[-999]*len(id_list)
    reliability_group=["-"]*len(id_list)
    selected_list=[0]*len(id_list)


# First step: Remove same molecules with same solubility values (change their SMILES into "XXX")
    for i in range(0,len(id_list)):        
        same_value_List=[]       
        if(smiles_list[i] != "XXX" ):
            same_value_List.append(i)
# collect same molecules with in range of 0.01 solubility value          
            for j in range(i+1,len(id_list)):  
                if(inchi_list[i]==inchi_list[j]):
                    if(math.fabs(solubility_list[i]-solubility_list[j])<=0.01):
                        same_value_List.append(j)

# select the best source according to: 1:name existance 2: size of the dataset (already in ordered according to size)                            
            if(len(same_value_List)>1):
                bestId=same_value_List[0]
                for sameId in same_value_List:
                    if((pd.isnull(name_list[bestId]) or name_list[bestId]=="-") and ( not pd.isnull(name_list[sameId]) and name_list[sameId]!="-")):
                        bestId=sameId
                same_value_List.remove(bestId)
                for sameId in same_value_List:
                    smiles_list[sameId]="XXX"                    
                    same_value_counter=same_value_counter+1                    
    print ("Total removed same molecule with same value: "+str(same_value_counter))
           

# Second step: Select the most reliable solubility value among the same molecules (change unselected SMILES into XXX)    
    for i in range(0,len(id_list)):       
        same_molecule_List=[]       
        
# collect same molecules with different solubility value   
        if(smiles_list[i] != "XXX" and selected_list[i]==0):    
            same_molecule_List.append(i)
            for j in range(i+1,len(id_list)):  
                if(smiles_list[j] != "XXX" and inchi_list[i]==inchi_list[j]):
                    same_molecule_List.append(j)

# if occurrence count=1 (Group:G1)        
            if(len(same_molecule_List)==1):
                selected_list[i]=1
                reliability_group[i]="G1"
                SD[i]=0
                ocurrence_count[i]=1

# if occurrence count = 2  (closest to reference (prediction) method )                  
            elif(len(same_molecule_List)==2):

# calculate difference betweeen values and prediction (tie breaker)            
                diff1=math.fabs(solubility_list[same_molecule_List[0]]-prediction_list[same_molecule_List[0]])   
                diff2=math.fabs(solubility_list[same_molecule_List[1]]-prediction_list[same_molecule_List[1]])

                bestId=same_molecule_List[0]
                if(diff1<=diff2):
                    smiles_list[same_molecule_List[1]]="XXX"
                    different_value_counter_2times=different_value_counter_2times+1
                    bestId=same_molecule_List[0]
                    selected_list[bestId]=1
                else:
                    smiles_list[same_molecule_List[0]]="XXX"
                    different_value_counter_2times=different_value_counter_2times+1 
                    bestId=same_molecule_List[1]
                    selected_list[bestId]=1
                    
# decide reliability group (if SD>0.5 Group:G2, else Group:G3)   
                diff=math.fabs(solubility_list[same_molecule_List[0]]-solubility_list[same_molecule_List[1]]) 
                if(diff>1):                                                
                    reliability_group[bestId]="G2"
                else:
                    reliability_group[bestId]="G3"                    

# store differences and SD and occurrence count            

                
                SD[bestId]=diff/2
                ocurrence_count[bestId]=2               
                                                      
# if occurrence count > 2 (closest to mean method )                        
            elif(len(same_molecule_List)>2):
                total=0
                different_solubility_values_list=[]
                
                for sameId in same_molecule_List:
                    total=total+solubility_list[sameId]
                    different_solubility_values_list.append(solubility_list[sameId])
                    
                mean=total / len(same_molecule_List)    
                bestId=same_molecule_List[0]
                bestDiff=999
                
                for sameId in same_molecule_List:
                    diff=math.fabs(solubility_list[sameId]-mean) 
                    if(diff<bestDiff):
                        bestId=sameId
                        bestDiff=diff
                selected_list[bestId]=1                        
               
                std=np.std(different_solubility_values_list, axis=0)
                SD[bestId]=std
                ocurrence_count[bestId]=len(same_molecule_List)
                
# decide reliability group (if SD>0.5 Group:G4, else Group:G5)                                
                if(std>0.5):                                                
                    reliability_group[bestId]="G4"
                else:
                    reliability_group[bestId]="G5"                 
                
                same_molecule_List.remove(bestId)
                for sameId in same_molecule_List:
                    smiles_list[sameId]="XXX"
                    different_value_counter_mutiple=different_value_counter_mutiple+1                                        

    
# add reliability information to curated dataset and filter duplicates    
    data_merged['SD']=pd.Series(SD)    
    data_merged['Occurrences']=pd.Series(ocurrence_count)
    data_merged['Group']=pd.Series(reliability_group)
    data_merged=data_merged.drop(columns=['Prediction'])
    data_filtered=data_merged[data_merged['Group'] !="-"]
    data_filtered=data_filtered.reset_index(drop=True)    
    
                    
    return data_filtered


