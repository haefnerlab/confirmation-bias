#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 09:41:32 2020

@author: liushizhao
"""

import numpy as np
import pickle as pk
import os
import matplotlib.pyplot as plt
from scipy import stats

# load analytical results
name = 'LineFitting/analytical.pk'
with open(name, 'rb') as handle:
    file  = pk.load(handle)
analytical_posterior = file[0]
alpha_space = file[1]
beta_space = file[2]
# normalize the posterior, from log space to probabilty space summing up to 1
analytical_posterior = np.exp(analytical_posterior)
analytical_posterior = analytical_posterior /np.sum(analytical_posterior[:])
# draw samples from analytical posterior
# create a population of alph&beta
alphav, betav = np.meshgrid(alpha_space , beta_space)

k = 10000
analytical_posterior_flat = analytical_posterior.flatten()
index = np.random.choice(np.arange(len(analytical_posterior_flat)),size = k,replace = True,p = analytical_posterior_flat)
analytical_samples = np.zeros((k,2))
analytical_samples[:,0] = alphav.flatten()[index]
analytical_samples[:,1] = betav.flatten()[index]
#%%
samples_list = [2000,3000,5000,10000]

resultsFolder = 'LineFitting/paramDraws/'
ks_value_apt = np.zeros((len(samples_list),1))
ks_value_mh = np.zeros((len(samples_list),1))
for i,N in enumerate(samples_list):
    # load apt results
    name = 'posterior_apt_1000_%d_3.pk' %N
    with open(os.path.join(resultsFolder,name), 'rb') as handle:
        file = pk.load(handle)
    
    n_rounds = len(file)
    posterior = file[n_rounds-1]
    apt_samples = posterior.gen(k)
    
    # calculate the ks-statistics between analytical samples and apt samples
    ks_alpha,_ = stats.ks_2samp(analytical_samples[:,0], apt_samples[:,0])
    ks_beta,_ = stats.ks_2samp(analytical_samples[:,1], apt_samples[:,1])
    
    ks_value_apt[i] = np.mean([ks_alpha,ks_beta])

    # load sampling results
    name = 'mh_samples_1000_%d_3.pk' %N
    with open(os.path.join(resultsFolder,name), 'rb') as handle:
        file = pk.load(handle)
    mh_samples = file
    # calculate the ks-statistics between analytical samples and apt samples
    ks_alpha,_ = stats.ks_2samp(analytical_samples[:,0], mh_samples[:,0])
    ks_beta,_ = stats.ks_2samp(analytical_samples[:,1], mh_samples[:,1])
    
    ks_value_mh[i] = np.mean([ks_alpha,ks_beta])
    
plt.figure()
plt.plot(ks_value_apt)
plt.plot(ks_value_mh)