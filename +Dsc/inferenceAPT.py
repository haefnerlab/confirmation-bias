#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 15:26:46 2020

@author: liushizhao
"""

#%% 
# define simulator
import matlab.engine
from delfi.simulator.BaseSimulator import BaseSimulator
import scipy.io
import os


def ConfirmationBiasSimulator(params,signals):
    # run matlab code run.vectorized and read results
    
    sim_results = eng.Model.runVectorized(params, signals)
    return sim_results
# define a simulator class linking to the real simulator

def ConfirmationBias(BaseSimulator):
    def __init__(self,signal,seed = None):
        dim_param = 5
        super().__init__(dim_param = dim_param,seed = seed)
        self.signal = signal
        self.simulate = ConfirmationBiasSimulator
    def gen_single(self,params):
        sim_results = self.simulate(params,self.signal)
        return sim_results

#%%
# summary statics of observed data
# but since our data is only 1/-1 choice we don't need to define summary statistics
# How about also use lpo?

from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
class ConfirmationBiasStats(BaseSummaryStats):
    def __init__(self,choice):
        self.choice = choice
    def calc(self):
        stats = self.choice
        return stats
#%%

import delfi.distribution as dd
import delfi.generator as dg
import numpy as np
import delfi.inference as infer

if __name__ == "__main__":

    # load the MATLAB engine
    eng = matlab.engine.start_matlab()

    datafolder = '../dscData'
    filename = 'syntheticData_priorC5.mat'
    mat = scipy.io.loadmat(os.path.join(datafolder, filename))
    # load data


    
    # define prior over model parameters
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
    n_processes = 10
    seeds_m = np.arange(1,n_processes+1,1)
    m = []
    for i in range(n_processes):
        m.append(ConfirmationBias(signal),seed=seeds_m[i])
    g = dg.MPGenerator(models=m, prior=prior, summary=s)

 
    # do inference
    # inference object
    res = infer.SNPEC(g,
                    obs=obs_stats,
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