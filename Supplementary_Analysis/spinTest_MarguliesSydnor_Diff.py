#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  21 17:32:52 2024

@author: f_meck01
"""

import numpy as np
import pandas as pd
import os
from netneurotools import datasets as nntdata
from neuromaps import nulls, stats, datasets, parcellate
import h5py
from scipy import io
import time

studyDir = '.../'
permutations = 100000
data4spinTestPath = studyDir + '/data4Spintest_allDiff.csv'

#read text file into NumPy array
data4spinTest = pd.read_csv(data4spinTestPath)
varNames = data4spinTest.columns[1:]
data4spinTest = data4spinTest.to_numpy()

orderFile = studyDir + 'order.txt'
order = np.loadtxt(orderFile,  dtype='int')
order = order-1 # make it python indeces

parcels = order.shape[0]

possibleCols = range(2,52)
numCols = len(possibleCols)

data4spinTestOrder = data4spinTest[order,] # reorder data

dataRS = data4spinTestOrder[:,possibleCols]

lausanne250 = nntdata.fetch_cammoun2012(version='fslr32k')['scale125']
parcL = parcellate.Parcellater(lausanne250[0], 'fsLR', hemi = 'L') # fit the parcellator
parcR = parcellate.Parcellater(lausanne250[1], 'fsLR', hemi = 'R') # fit the parcellator

if not os.path.exists(studyDir + '_NullsDiff'):
    os.mkdir(studyDir + '_NullsDiff')
    
## MARGULIES
for annotation in datasets.available_annotations():
    print(annotation)
margulies_1 = datasets.fetch_annotation(source = 'margulies2016', desc = 'fcgradient01', space = 'fsLR', den = '32k')
margulies_2 = datasets.fetch_annotation(source = 'margulies2016', desc = 'fcgradient02', space = 'fsLR', den = '32k')
sydnor = datasets.fetch_annotation(source = 'sydnor2021', desc = 'SAaxis', space = 'fsLR', den = '32k')

margulies_1_parc_L = parcL.fit_transform(margulies_1[0], 'fsLR', hemi = 'L') # parcellate the data (archetypal gradient)
margulies_1_parc_R = parcR.fit_transform(margulies_1[1], 'fsLR', hemi = 'R') # parcellate the data (archetypal gradient)
margulies_1_parc_L = np.delete(margulies_1_parc_L, (3), axis=0)
margulies_1_parc_R = np.delete(margulies_1_parc_R, (3), axis=0)
margulies_1_parc = np.concatenate((margulies_1_parc_L, margulies_1_parc_R), axis=0)

margulies_2_parc_L = parcL.fit_transform(margulies_2[0], 'fsLR', hemi = 'L') # parcellate the data (archetypal gradient)
margulies_2_parc_R = parcR.fit_transform(margulies_2[1], 'fsLR', hemi = 'R') # parcellate the data (archetypal gradient)
margulies_2_parc_L = np.delete(margulies_2_parc_L, (3), axis=0)
margulies_2_parc_R = np.delete(margulies_2_parc_R, (3), axis=0)
margulies_2_parc = np.concatenate((margulies_2_parc_L, margulies_2_parc_R), axis=0)

sydnor_parc_L = parcL.fit_transform(sydnor[0], 'fsLR', hemi = 'L') # parcellate the data (archetypal gradient)
sydnor_parc_R = parcR.fit_transform(sydnor[1], 'fsLR', hemi = 'R') # parcellate the data (archetypal gradient)
sydnor_parc_L = np.delete(sydnor_parc_L, (3), axis=0)
sydnor_parc_R = np.delete(sydnor_parc_R, (3), axis=0)
sydnor_parc = np.concatenate((sydnor_parc_L, sydnor_parc_R), axis=0)

# combine all
data4spinTestOrderMaps = np.c_[data4spinTestOrder, margulies_1_parc,margulies_2_parc,sydnor_parc]
varNames = varNames.tolist()
varNames.extend(["Margulies_1", "Margulies_2", "Sydnor"])

for iVar in range(48,51): #possibleCols: 
    variable = varNames[iVar-1]
    rotated = nulls.alexander_bloch(pd.to_numeric(data4spinTestOrderMaps[:,iVar]), atlas='fsLR', density='32k', parcellation=lausanne250, n_perm=permutations, seed=420)
    outputFileName = studyDir + '_NullsDiff/' + 'Spin_' + variable +'_rotated.h5'
    h5f = h5py.File(outputFileName, 'w')
    h5f.create_dataset('dataset_1', data=rotated)
    h5f.close()
    print('DONE: ' + variable)
    time.sleep(1)
    
permutations = 100000
   
possibleCols = range(48,55)
numCols = len(possibleCols)

correlation_results = np.ones((numCols,numCols,2))
for i in range(numCols):
    variable = varNames[possibleCols[i]-1]
    print(variable)
    dataFramePath = studyDir + '_NullsDiff/' + 'Spin_' + variable +'_rotated.h5'
    h5f = h5py.File(dataFramePath,'r')
    rotated = h5f['dataset_1'][:]
    h5f.close()
    pos_i = possibleCols[i]
    for j in range(numCols):
        pos_j = possibleCols[j]
        if not i == j:
            correlation_results[i,j,0], correlation_results[i,j,1] = stats.compare_images(pd.to_numeric(data4spinTestOrderMaps[:,pos_i]),
                                                                                          pd.to_numeric(data4spinTestOrderMaps[:,pos_j]), nulls=rotated)
        time.sleep(1)

outputMatFile = studyDir + 'spinTest_Sub_results_allDiff.mat'
io.savemat(outputMatFile, {'results': correlation_results})

