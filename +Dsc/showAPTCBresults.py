#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:33:13 2020

@author: liushizhao
"""

import pickle
import os
import numpy as np
from scipy.io import loadmat
from APTCB_utils import makeFigures
if __name__ == "__main__":
    C = 2
    g = 2
    l = 2
    s = 1
    results_folder = '../dscData/resultsSyntheticDataCB_version3'
    results_name = 'lshc_resTrial1priorC%dgamma%dlapse%dsamples%d.pkl' %(C,g,l,s)
    with open(os.path.join(results_folder,results_name),'rb') as res:
        results = pickle.load(res)
    
    log = results[0]
    posterior = results[1]
    
    # load parameters information
    paraInfo = loadmat(os.path.join('../dscData/SyntheticDataCB','genDataInfo.mat'))
    prior_C_list = paraInfo['genDataInfo']['prior_C_list'].tolist()[0][0][0]
    gamma_list = paraInfo['genDataInfo']['gamma_list'].tolist()[0][0][0]
    lapse_list = paraInfo['genDataInfo']['lapse_list'].tolist()[0][0][0]
    samples_list = paraInfo['genDataInfo']['samples_list'].tolist()[0][0][0]
    true_params = [prior_C_list[C-1],gamma_list[g-1],lapse_list[l-1]]
                   
    labels_params = ['prior_C','gamma','lapse']
                       
    prior_min = np.array([0,0,0])
    prior_max = np.array([1,0.5,0.5])               
    makeFigures(log,posterior,true_params,labels_params,prior_min,prior_max)
