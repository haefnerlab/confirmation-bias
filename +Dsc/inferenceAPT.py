#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 15:26:46 2020

@author: liushizhao
"""

#%% 
import matlab.engine
from delfi.simulator.BaseSimulator import BaseSimulator
import os

import delfi.distribution as dd
import delfi.generator as dg
import numpy as np
import delfi.inference as infer
from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
import copy
from mat4py import loadmat
#%%
def ConfirmationBiasSimulator(xval,fieldstofit,params, signals, engine):
    # run matlab code run.vectorized and read results
    print("Hello I am in Simulator")
    # print(type(signals))
    # print(type(params[0]))
    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        paramsNew[f] = float(xval[i])

    #signals = np.ndarray.tolist(signals)
    sim_results = engine.Model.runVectorized(paramsNew, signals)
    return sim_results
# define a simulator class linking to the real simulator

class ConfirmationBias(BaseSimulator):
    
    def __init__(self, signals, fieldstofit,params,engine, seed = None):
        
        dim_param = 5
        super().__init__(dim_param = dim_param,seed = seed)
        self.signals = signals
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        sim_results = self.simulate(xval,self.fieldstofit,self.params,self.signals, self.engine)
        return sim_results
#%%

class ConfirmationBiasStats(BaseSummaryStats):
    def __init__(self,choice):
        self.choice = choice
    def calc(self):
        stats = self.choice
        return stats
#%%
# load the MATLAB engine
engine = matlab.engine.start_matlab()

datafolder = '../dscData'
filename = 'syntheticData_priorC5.mat'
syntheticData = loadmat(os.path.join(datafolder, filename))
# load data
choice = syntheticData['sim_results']
signals = syntheticData['signals']
params = syntheticData['params']

# define prior over model parameters
fieldstofit = ['prior_C','gamma','samples','lapse','sensor_noise']
seed_p = 2
prior_min = np.array([0,0,1,0,0])
prior_max = np.array([1,1,100,1,5])
prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)

# inference parameters
seed_inf = 1

pilot_samples = 2000

# training schedule
n_train = 2000
n_rounds = 1

# fitting setup
minibatch = 100
epochs = 100
val_frac = 0.05

# network setup
n_hiddens = [50,50]

# convenience
prior_norm = True

# MAF parameters
density = 'maf'
n_mades = 5         # number of MADES


# define generator class
s = ConfirmationBiasStats(choice)
# use multiple processes of simulators
n_processes = 1
seeds_m = np.arange(1,n_processes+1,1)
m = []
for i in range(n_processes):
    m.append(ConfirmationBias(signals, fieldstofit,params,engine, seeds_m[i]))
g = dg.MPGenerator(models=m, prior=prior, summary=s)
#%%
# do inference
# inference object
res = infer.SNPEC(g,
                obs=choice,
                n_hiddens=n_hiddens,
                seed=seed_inf,
                pilot_samples=pilot_samples,
                n_mades=n_mades,
                prior_norm=prior_norm,
                density=density)

# train
log, _, posterior = res.run(
                    n_train=n_train,
                    n_rounds=n_rounds,
                    minibatch=minibatch,
                epochs=epochs,
                silent_fail=False,
                proposal='prior',
                val_frac=val_frac,
                verbose=True,)
