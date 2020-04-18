#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:38:55 2019

@author: Murat Cihan Sorkun
"""

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd


def generate(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolMR = Descriptors.MolMR(mol)
        desc_HeavyAtomCount = Descriptors.HeavyAtomCount(mol)
        desc_NumHAcceptors = Descriptors.NumHAcceptors(mol)
        desc_NumHDonors = Descriptors.NumHDonors(mol)
        desc_NumHeteroatoms = Descriptors.NumHeteroatoms(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        desc_NumValenceElectrons = Descriptors.NumValenceElectrons(mol)              
        desc_NumAromaticRings = Descriptors.NumAromaticRings(mol)        
        desc_NumSaturatedRings = Descriptors.NumSaturatedRings(mol)        
        desc_NumAliphaticRings = Descriptors.NumAliphaticRings(mol)
        desc_RingCount = Descriptors.RingCount(mol)        
        desc_TPSA = Descriptors.TPSA(mol)
        desc_LabuteASA = Descriptors.LabuteASA(mol)        
        desc_BalabanJ = Descriptors.BalabanJ(mol)
        desc_BertzCT = Descriptors.BertzCT(mol) 
                
        row = np.array([desc_MolWt,desc_MolLogP,desc_MolMR,desc_HeavyAtomCount,desc_NumHAcceptors,desc_NumHDonors,desc_NumHeteroatoms,
                        desc_NumRotatableBonds,desc_NumValenceElectrons,desc_NumAromaticRings,desc_NumSaturatedRings,
                        desc_NumAliphaticRings,desc_RingCount,desc_TPSA,desc_LabuteASA,desc_BalabanJ,desc_BertzCT])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MolWt","MolLogP","MolMR","HeavyAtomCount","NumHAcceptors","NumHDonors","NumHeteroatoms",
                 "NumRotatableBonds","NumValenceElectrons","NumAromaticRings","NumSaturatedRings",
                 "NumAliphaticRings","RingCount","TPSA","LabuteASA","BalabanJ","BertzCT"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    print("Total descriptors generated: 17 x "+str(len(smiles)))
    return descriptors