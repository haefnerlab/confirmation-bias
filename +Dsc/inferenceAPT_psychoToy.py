#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:21:47 2020

@author: liushizhao
"""

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
import logging
import datetime
from scipy.io import savemat
import matplotlib.pyplot as plt
from delfi.utils.viz import samples_nd


def psychoToySimulator(x,x_center,params, engine,seed=None):
    choice = engine.Dsc.psychoBernoulliSimulator(x.tolist(),x_center,params.tolist())
    sim_choice = np.squeeze(np.asarray(choice))
    return sim_choice
class psychoToyModel(BaseSimulator):
    def __init__(self, x, x_center,engine,dim_param,seed = None):


        super().__init__(dim_param=dim_param, seed=seed)
        self.x = x
        self.x_center = x_center
        self.simulate = psychoToySimulator
        self.engine = engine
    def gen_single(self,params):

        params = np.asarray(params)

        assert params.ndim == 1, 'params.ndim must be 1'

        hh_seed = self.gen_newseed()

        choice = self.simulate(self.x, self.x_center,params, self.engine, seed=hh_seed)

        return {'choice': choice}
    
class psychoToyStats(BaseSummaryStats):
    
    def __init__(self,x):
        self.x = x
    def calc(self,repetition_list):
        seg = 20
        stats = []
    
        for r in range(len(repetition_list)):

            data = repetition_list[r]
            x = self.x
            choice = np.asarray(data['choice'])
            N = x.shape[0]
            x_sorted = np.asarray([x for x,choice in sorted(zip(x,choice))])
            choice_sorted = np.asarray([choice for x,choice in sorted(zip(x,choice))])
            sum_stats_vec = []
            for i in np.arange(seg):
                ind = np.arange(i*int(N/seg),(i+1)*int(N/seg)).astype(int)
                p_empirical = np.mean(choice_sorted[ind])

                sum_stats_vec.extend([p_empirical])
            stats.append(sum_stats_vec)
        return stats
#%%
if __name__ == "__main__":
    
    sensitivity = 1.0
    bias = 3.0

    true_params = [sensitivity ,bias]
    labels_params = ['sensitivity' ,'bias']
    # generate observation
    N = 1000
    x0 = np.random.uniform(0, 20, N)
    x_center = 10.0
    engine = matlab.engine.start_matlab()
    
    # define linear noise model
    m = psychoToyModel( x0, x_center,engine,dim_param = 2)
    
    # generate and show observation
    obs0 = m.gen_single(true_params)
    plt.figure()
    plt.plot(x0,obs0['choice'],'bo')
    # run multiple times
    # check the match between empirical p and theorical p
    nRun = 20
    obs0Multi = np.zeros((N,nRun))
   
    for i in range(nRun):
        obs0Multi[:,i] = m.gen_single(true_params)['choice']
    plt.figure()
    plt.plot(x0,np.mean(obs0Multi,axis = 1),'bo')    
    # define prior
    prior_min = np.array([0.1,-5])
    prior_max = np.array([2,5])
    seed_p = 2
    prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    
    s = psychoToyStats(x0)
    g = dg.Default(model=m, prior=prior, summary=s)
    
    # define statsitics summary of observation 
    obs_stats = s.calc([obs0])
    
    #%%
    
    # set hyparameters 
    # training schedule
    n_train = 5000
    n_rounds = 2
    seed_inf = 1
    pilot_samples = 2000


    val_frac = 0.05
    # network setup
    n_hiddens = [50,50]
    minibatch = 500
    epochs = 100
    
    prior_norm = True
    
    # MAF parameters
    density = 'mog'
    n_mades = 5         # number of MADES
    



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
                    proposal='gaussian',
                    val_frac=val_frac,
                    verbose=True,)
    
    
    #%%
    fig = plt.figure(figsize=(15,5))

    plt.plot(log[0]['loss'],lw=2)
    plt.xlabel('iteration')
    plt.ylabel('loss');
    
    
    
    


    prior_min = g.prior.lower
    prior_max = g.prior.upper
    prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)
    
    posterior_samples = posterior[0].gen(10000)
    
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
        