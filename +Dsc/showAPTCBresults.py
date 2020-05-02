#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:33:13 2020

@author: liushizhao
"""

import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from delfi.utils.viz import samples_nd
from scipy.io import loadmat
def makeFigures(posterior,true_params):
    
    fig = plt.figure(figsize=(15,5))
    
    plt.plot(log[1]['loss'],lw=2)
    plt.xlabel('iteration')
    plt.ylabel('loss');
    
    
    
    labels_params = ['prior_C','gamma','lapse','samples']
    prior_min = np.array([0,0,0,1])
    prior_max = np.array([1,0.5,0.5,50])

    prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)
    
    posterior_samples = posterior[1].gen(10000)
    
    ###################
    # colors
    hex2rgb = lambda h: tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    
    # RGB colors in [0, 255]
    col = {}
    col['GT']      = hex2rgb('30C05D')
    col['SNPE']    = hex2rgb('2E7FE8')
    col['SAMPLE1'] = hex2rgb('8D62BC')
    col['SAMPLE2'] = hex2rgb('AF99EF')
    
    # convert to RGB colors in [0, 1]
    for k, v in col.items():
        col[k] = tuple([i/255 for i in v])
    
    ###################
    # posterior
    fig, axes = samples_nd(posterior_samples,
                           limits=prior_lims,
                           ticks=prior_lims,
                           labels=labels_params,
                           fig_size=(5,5),
                           diag='kde',
                           upper='kde',
                           hist_diag={'bins': 50},
                           hist_offdiag={'bins': 50},
                           kde_diag={'bins': 50, 'color': col['SNPE']},
                           kde_offdiag={'bins': 50},
                           points=[true_params],
                           points_offdiag={'markersize': 5},
                           points_colors=[col['GT']],
                           title='');
if __name__ == "__main__":
    C = 2
    g = 2
    l = 2
    s = 1
    results_folder = '../dscData/resultsSyntheticDataCB'
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
    true_params = [prior_C_list[C-1],gamma_list[g-1],lapse_list[l-1],samples_list[s-1]]
                   
                   
    makeFigures(posterior,true_params)
