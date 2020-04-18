#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:23:41 2019

@authors: Murat Cihan Sorkun, Elif Sorkun
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdSLNParse
from rdkit.Chem import rdinchi
from rdkit.Chem import Descriptors
import urllib.request as urllib2  
from time import sleep


#Convert SLN to SMILES 
def sln_to_smiles(sln_list,verbose=0):      
    
    error_counter=0
    warning_counter=0
    smiles_list= []
    inchis_from_sln=[]
    inchis_from_smiles=[]

    for elem in sln_list:
        mol=rdSLNParse.MolFromQuerySLN(elem) 

        if(mol!=None):
            smiles=Chem.MolToSmiles(mol)
            smiles_list.append(smiles)
            inchi, retcode, message, logs, aux =rdinchi.MolToInchi(mol)  
            inchis_from_sln.append(inchi)              
            mol=Chem.MolFromSmiles(smiles) 
            
            if(mol!=None):                
                inchi, retcode, message, logs, aux =rdinchi.MolToInchi(mol)
                inchis_from_smiles.append(inchi)  
            else: 
                inchis_from_smiles.append("XXX")
                if(verbose!=0):
                    print(elem+ ": SMILES cannot converted to mol object, added XXX instead!")

        else:        
            smiles_list.append("XXX")
            print(elem+ ": SLN can not converted to mol object, added XXX instead!")
            error_counter=error_counter+1
        
    for i1,i2 in zip(inchis_from_sln,inchis_from_smiles):
        if(i1!=i2):
            if(verbose!=0):
                print("Warning:"+i1+" - "+i2)
            warning_counter=warning_counter+1
     
    print("\nConversion from SLN to SMILES is completed.")    
    print("Total errors:"+str(error_counter))
    print("Total warnings:"+str(warning_counter)+"\n")    
              
    return smiles_list

#Returns InChI InChIKey from SMILES
def smiles_to_inchi_inchikey(smiles,verbose=0):
          
    error_counter=0
    warning_counter=0
    inchis=[]
    inchis2=[]
    inchiKeys=[]

    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        if(mol!=None):
            inchi, retcode, message, logs, aux =rdinchi.MolToInchi(mol)
            inchiKey=rdinchi.InchiToInchiKey(inchi)   
            inchis.append(inchi)  
            inchiKeys.append(inchiKey)   
            try:
                mol, retcode, message, logs=rdinchi.InchiToMol(inchi)
                if(mol!=None):                
                    inchi2, retcode, message, logs, aux =rdinchi.MolToInchi(mol)
                    inchis2.append(inchi2)  
                else: 
                    inchis2.append("XXX") 
                    if(verbose!=0):
                        print(elem+ ": InChI cannot converted to mol object, added XXX instead!")
            except:
                if(verbose!=0):
                    print(retcode)
                    print(message)
                    print("Smiles:"+elem)
                inchis2.append("XXX")  
            
        else:  
            inchis.append("XXX")
            inchiKeys.append("XXX")
            inchis2.append("XXX")            
            if(verbose!=0):
                print(elem+ ": can not converted added XXX instead! ")            
            error_counter=error_counter+1

    for i1,i2 in zip(inchis,inchis2):
        if(i1!=i2):
            if(verbose!=0):
                print("Warning:"+i1+" - "+i2)
            warning_counter=warning_counter+1

    print("\nGeneration of InChI and InChIKey from SMILES is completed.")    
    print("Total errors:"+str(error_counter))
    print("Total warnings:"+str(warning_counter)+"\n")    
                                          
    return inchis,inchiKeys

#Returns ID list by dataset name and size
def generate_id_list(name,size):
          
    id_list= []

    for i in range(1, size+1):
        id_list.append(name+"-"+str(i))
        
    return id_list

#Converts different solubility units into the LogS
def convert_to_logs(smiles_list,unit_list,solubility_list,verbose=0):

    error_counter=0
    logs_list=[]
    
    for i in range(len(smiles_list)):  
        try:            
            mw=Descriptors.MolWt(Chem.MolFromSmiles(smiles_list[i]))                        
            
            try:
                if(unit_list[i]=="g/L"):
                    logs_list.append(np.log10(float(solubility_list[i])/mw))
                
                elif(unit_list[i]=="mg/L"):
                    logs_list.append(np.log10(float(solubility_list[i])/(mw*1000)))
                
                elif(unit_list[i]=="µg/L"):
                    logs_list.append(np.log10(float(solubility_list[i])/(mw*1000000)))
    
                elif(unit_list[i]=="M"):
                    logs_list.append(np.log10(float(solubility_list[i])))
                    
                elif(unit_list[i]=="µM"):
                    logs_list.append(np.log10(float(solubility_list[i])/1000000))
                
                elif(unit_list[i]=="LogS"):
                    logs_list.append(solubility_list[i])
        
                else:
                    if(verbose!=0):
                        print("Unknown metric:"+unit_list[i])  
                    error_counter=error_counter+1
                    logs_list.append(999)
            except:
                if(verbose!=0):
                    print("Error during the calculation! Solubility:" +str(solubility_list[i]) +" unit:"+unit_list[i])      
                error_counter=error_counter+1
                logs_list.append(999)
        except:
            if(verbose!=0):
                print("Error during the creation of mol object! Smiles:" +smiles_list[i])      
            error_counter=error_counter+1
            logs_list.append(999)
            
    print("\nConversion of units into LogS is completed.")   
    print("Total instances:"+str(len(smiles_list)))
    print("Total errors:"+str(error_counter)+"\n")
        
    return logs_list


#Collects SMILES from the Chemical Identifier Resolver web service of the National Cancer Institute using Name and CAS number
def collect_smiles_from_web(name_list,cas_list,verbose=0):

    error_counter=0
    smiles_list = []
       
    for i in range(len(name_list)):
        if(verbose!=0):
            print(name_list[i])
        try:
            sleep(0.04)
            target_url="https://cactus.nci.nih.gov/chemical/structure/"+cas_list[i].lstrip("0").rstrip(" ")+"/smiles"    
            context = urllib2.urlopen(target_url) 
            smiles=context.read().decode('utf-8')               
            smiles_list.append(smiles)      
            
        except:
            if(verbose!=0):
                print("Cannot collected SMILES using CAS:"+cas_list[i])    
                print("Trying to collect using Name:"+name_list[i])  
            try:
                sleep(0.04)
                target_url="https://cactus.nci.nih.gov/chemical/structure/"+name_list[i].rstrip(" ").replace(" ", "%20")+"/smiles"
                context = urllib2.urlopen(target_url) 
                smiles=context.read().decode('utf-8')
                smiles_list.append(smiles)            
            except:     
                if(verbose!=0):
                    print("Cannot collected SMILES using Name:"+name_list[i])
                smiles_list.append("XXX")
                error_counter=error_counter+1
                
    print("\nCollecting SMILES from webservice is completed.")    
    print("Total instances:"+str(len(name_list)))
    print("Total errors:"+str(error_counter)+"\n")
    
    return smiles_list           


#Collects names from the Chemical Identifier Resolver web service of the National Cancer Institute using InChIKey and SMILES
def collect_names_from_web(incikey_list,smiles_list,verbose=0):

    error_counter=0
    name_list = []
       
    for i in range(len(incikey_list)):
        try:
            sleep(0.04)
            target_url="https://cactus.nci.nih.gov/chemical/structure/"+incikey_list[i]+"/iupac_name"    
            context = urllib2.urlopen(target_url) 
            name=context.read().decode('utf-8')               
            name_list.append(name)      
            
        except:
            if(verbose!=0):
                print("Cannot collected name using InChIKey:"+incikey_list[i])    
                print("Trying to collect using SMILES:"+smiles_list[i])  
            try:
                sleep(0.04)
                target_url="https://cactus.nci.nih.gov/chemical/structure/"+smiles_list[i].replace("#","%23")+"/iupac_name"
                context = urllib2.urlopen(target_url) 
                name=context.read().decode('utf-8')
                name_list.append(name)             
            except:     
                if(verbose!=0):
                    print("Cannot collected name using SMILES:"+smiles_list[i])
                name_list.append("XXX")
                error_counter=error_counter+1
                
    print("\nCollecting names from webservice is completed.")    
    print("Total instances:"+str(len(incikey_list)))
    print("Total errors:"+str(error_counter)+"\n")
    
    return name_list     


#Collects solubility predictions from the ALOGPS 2.1 using InChIKey and SMILES
def collect_predictions_from_web(smiles_list,verbose=0):

    error_counter=0
    prediction_list = []
       
    for i in range(len(smiles_list)):
        try:
            sleep(0.04)
            target_url="http://146.107.217.178/web/alogps/calc?SMILES="+smiles_list[i]    
            context = urllib2.urlopen(target_url) 
            text=context.read().decode('utf-8')
            if("could not" not in text):
                prediction=text.split(" ")[5]
                prediction_list.append(prediction)   
            else:
                raise Exception("Cannot collected prediction for:"+smiles_list[i])
            
        except:
            if(verbose!=0):
                print("Cannot collected prediction for:"+smiles_list[i])    
            prediction_list.append(999)
            error_counter=error_counter+1
                
    print("\nCollecting predictions from webservice is completed.")    
    print("Total instances:"+str(len(smiles_list)))
    print("Total errors:"+str(error_counter)+"\n")
    
    return prediction_list   

#filter dataset by removing missing information. (for strings: "XXX" and for numeric: "999")
def clean(dataset,verbose=1):

    dataset_clean = dataset[(dataset.InChIKey != "XXX") & (dataset.Solubility != "999") ].reset_index(drop=True)
    
    print("Total removed during the cleaning:"+str(len(dataset)-len(dataset_clean))+"\n")
    
    return dataset_clean   



