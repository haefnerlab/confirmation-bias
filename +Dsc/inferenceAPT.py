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
def ConfirmationBiasSimulator(xval,fieldstofit,params, dataPath, engine):
    # run matlab code run.vectorized and read results
   # print("Hello I am in Simulator")

    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        if type(params[f]) is float:
            paramsNew[f] = float(xval[i])
        elif type(params[f]) is int:
            paramsNew[f] = int(xval[i])


    sim_results = engine.Model.runVectorizedAPT(paramsNew, dataPath)
    #print("Hello I am out of simulator")
    
    sim_choice = np.squeeze(np.asarray(sim_results['choices']))



    return sim_choice
# define a simulator class linking to the real simulator

class ConfirmationBias(BaseSimulator):
    
    def __init__(self, dataPath, fieldstofit,params,engine, seed = None):
        
        dim_param = 5
        super().__init__(dim_param = dim_param,seed = seed)
        self.dataPath = dataPath
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        sim_results = self.simulate(xval,self.fieldstofit,self.params,self.dataPath, self.engine)
        return sim_results
#%%

class ConfirmationBiasStats(BaseSummaryStats):
    def __init__(self):
        self.n_summary = 1
    def calc(self,choice):
        stats = choice

     
        return stats
#%%
# load the MATLAB engine
if __name__ == "__main__":
    engine = matlab.engine.start_matlab()
    
    datafolder = '../dscData'
    filename = 'syntheticData_priorC5.mat'
    dataPath = os.path.join(datafolder,filename)
    syntheticData = loadmat(dataPath)
    # load data
    choice = np.squeeze(np.asarray(syntheticData['sim_results']['choices']))

    params = syntheticData['params']
    
    # define prior over model parameters
    fieldstofit = ['prior_C','gamma','samples','lapse']
    seed_p = 2
    prior_min = np.array([0,0,1,0])
    prior_max = np.array([1,1,100,1])
    prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    
    
    #sim_resultsTest = ConfirmationBiasSimulator([0.5,0.01,5,0,0],fieldstofit,params, dataPath, engine)
    
    
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
    s = ConfirmationBiasStats()
    # use multiple processes of simulators
    n_processes = 1
    seeds_m = np.arange(1,n_processes+1,1)
    m = ConfirmationBias(dataPath, fieldstofit,params,engine, seeds_m[0])
    
    g = dg.Default(model=m, prior=prior, summary=s)
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
