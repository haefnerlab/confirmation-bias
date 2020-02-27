#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:27:08 2020

@author: liushizhao
"""

from mat4py import loadmat
import os
import numpy as np
import copy

datafolder = '../dscData'
dataName = 'syntheticData_priorC9.mat'

syntheticData = loadmat(os.path.join(datafolder,dataName))
choice = syntheticData['sim_results']
signals = syntheticData['signals']

params = syntheticData['params']

paramsNew = copy.deepcopy(params)
fieldstofit = ['gamma']
xval = [0.05]
for i,f in enumerate(fieldstofit):
    paramsNew[f] = xval[i]

