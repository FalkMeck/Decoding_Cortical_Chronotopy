#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 09:00:03 2023

Run the nusicanceRegrssionPipeline_FM for all subjects, without the multiprocessing part

@author: f_meck01
"""
import os
os.chdir('.../Analysis')

import _Code.nuisanceRegressionPipeline_FM as nrp

nrp.step1_createNuisanceRegressors()
nrp.step2_nuisanceRegression(model='24pXaCompCorXVolterra',spikeReg=True,zscore=False)
