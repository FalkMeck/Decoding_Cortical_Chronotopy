#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:13:35 2024

@author: fehringj
"""

import numpy as np
import pandas as pd
# import os
from netneurotools import datasets as nntdata
from neuromaps import stats as nmStats
import h5py
from scipy import io
import time
from scipy import stats

studyDir = '.../'
permutations = 100000
data4spinTestPath = studyDir + '/data4spinTest.csv'

#read text file into NumPy array
data4spinTest = pd.read_csv(data4spinTestPath)
varNames = data4spinTest.columns[1:]
data4spinTestNP = data4spinTest.to_numpy()

orderFile = studyDir + 'orderSpin.txt' # provided in folder, gives transforms order of parcels to correct order for spin test
order = np.loadtxt(orderFile,  dtype='int')
order = order-1 # make it python indeces

parcels = order.shape[0]

dataAllCont = data4spinTestNP[order,] # reorder data

if not os.path.exists(studyDir + '_Nulls'):
    os.mkdir(studyDir + '_Nulls')
    
possibleCols = [2,4,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21]

dataAllCont[:,4] = np.where(dataAllCont[:,4] == "RC_Club", 3, 
                   np.where(dataAllCont[:,4] == "RC_Feeder", 2,1))
dataAllCont[:,6] = np.where(dataAllCont[:,6] == "DC_Club", 3, 
                   np.where(dataAllCont[:,6] == "DC_Feeder", 2,1))

lausanne250 = nntdata.fetch_cammoun2012(version='fslr32k')['scale125']

for iVar in possibleCols: #possibleCols: 
    variable = varNamesAll[iVar-1]
    rotated = nulls.alexander_bloch(pd.to_numeric(dataAllCont[:,iVar]), atlas='fsLR', density='32k', parcellation=lausanne250, n_perm=permutations, seed=420)
    outputFileName = studyDir + '_Nulls/' + 'Mean_' + variable +'_rotated.h5'
    h5f = h5py.File(outputFileName, 'w')
    h5f.create_dataset('dataset_1', data=rotated)
    h5f.close()
    print('DONE: ' + variable)
    time.sleep(1)
    
numCols = len(possibleCols)

correlation_results = np.ones((numCols,numCols,2))
mi = 0
for i in possibleCols:
    mj = 0
    variable = varNamesAll[i-1]
    print(variable)
    dataFramePath = studyDir + '_Nulls/' + 'Mean_' + variable +'_rotated.h5'
    h5f = h5py.File(dataFramePath,'r')
    rotated = h5f['dataset_1'][:]
    h5f.close()
    for j in possibleCols:
        if not i == j:
            
            if (i == 4 or i == 6) and not (j == 4 or j == 6): # if first variable is either RC or DC
                Group1 = dataAllCont[dataAllCont[:,i] == 1,j]
                Group2 = dataAllCont[dataAllCont[:,i] == 2,j]
                Group3 = dataAllCont[dataAllCont[:,i] == 3,j]
                correlation_results[mi,mj,0] = stats.kruskal(Group3, Group2, Group1)[0]
                KW_results = np.zeros((permutations,1))
                for k in range(permutations):
                    Group1_k = dataAllCont[rotated[:,k] == 1,j]
                    Group2_k =  dataAllCont[rotated[:,k] == 2,j]
                    Group3_k =  dataAllCont[rotated[:,k] == 3,j]
                    KW_results[k,0] = stats.kruskal(Group3_k, Group2_k, Group1_k)[0]
                
                KWboot = KW_results[:,0] >= correlation_results[mi,mj,0]
                correlation_results[mi,mj,1] = KWboot.sum()/permutations
           
            if (j == 4 or j == 6) and not (i == 4 or i == 6):  # if second variable is either RC or DC
                Group1 = dataAllCont[dataAllCont[:,j] == 1,i]
                Group2 = dataAllCont[dataAllCont[:,j] == 2,i]
                Group3 = dataAllCont[dataAllCont[:,j] == 3,i]
                correlation_results[mi,mj,0] = stats.kruskal(Group3, Group2, Group1)[0]
                KW_results = np.zeros((permutations,1))
                for k in range(permutations):
                    Group1_k = rotated[dataAllCont[:,j] == 1,k]
                    Group2_k =  rotated[dataAllCont[:,j] == 2,k]
                    Group3_k =  rotated[dataAllCont[:,j] == 3,k]
                    KW_results[k,0] = stats.kruskal(Group3_k, Group2_k, Group1_k)[0]
                KWboot = KW_results[:,0] >= correlation_results[mi,mj,0]
                correlation_results[mi,mj,1] = KWboot.sum()/permutations
            
            if not (i == 4 or i == 6) and not (j == 4 or j == 6):  # if the variable is neither RC nor DC
                correlation_results[mi,mj,0], correlation_results[mi,mj,1] = nmStats.compare_images(pd.to_numeric(dataAllCont[:,i]),
                                                                                          pd.to_numeric(dataAllCont[:,j]), nulls=rotated)
            
            if (i == 4 or i == 6) and (j == 4 or j == 6):  # if both variables are either RC or DC
                ctTabTrue = pd.crosstab(dataAllCont[:,i], dataAllCont[:,j])
                correlation_results[mi,mj,0] =  stats.chi2_contingency(ctTabTrue)[0]
                Chi_results = np.zeros((permutations,1))
                for k in range(permutations):
                    ctTab_k = pd.crosstab(rotated[:,k], dataAllCont[:,j])
                    Chi_results[k,0] =  stats.chi2_contingency(ctTab_k)[0]
                Chiboot = Chi_results[:,0] >= correlation_results[mi,mj,0]
                correlation_results[mi,mj,1] = Chiboot.sum()/permutations
            time.sleep(1)
        mj = mj +1   
    mi = mi +1

outputMatFile = studyDir + 'spinTest_MeanClass_results.mat'
io.savemat(outputMatFile, {'results': correlation_results})
