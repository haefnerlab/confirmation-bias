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
    
    def __init__(self,x_list,x0,option = 'choice_prob'):
        self.x_list = x_list
        self.x0 = x0
        self.option = option
    def calc(self,repetition_list):
       
        stats = []
    
        for r in range(len(repetition_list)):
            data = repetition_list[r]
            signals = self.x0
            x_list = self.x_list
            choices = np.asarray(data['choice'])
            
            if self.option == 'choice_prob':
                sum_stats_vec = []
                for x in x_list:
                    ind = np.where(x0 == x)[0]
                    sum_stats_vec.extend([np.mean(choices[ind])])
            elif self.option == 'choice_triggered_signal':
                    ind1 = np.where(choices == 0)[0]
                    ind2 = np.where(choices == 1)[0]
                    sum_stats_vec = np.hstack((np.mean(signals[ind1],axis = 0),np.mean(signals[ind2],axis = 0),\
                               np.std(signals[ind1],axis = 0)/np.sqrt(len(ind1)),\
                               np.std(signals[ind2],axis = 0)/np.sqrt(len(ind2)))) 
            

            stats.append(sum_stats_vec)
        return stats
#%%
if __name__ == "__main__":
    engine = matlab.engine.start_matlab()
    
    sensitivity_list = [0.5,0.7,1.0,1.2,1.5]
    bias_list = [-3.0,-2.0,-1.0,1.0,2.0,3.0]
    signal_0_mean = np.zeros((len(sensitivity_list),len(bias_list)))
    signal_1_mean = np.zeros((len(sensitivity_list),len(bias_list)))
    signal_0_sem = np.zeros((len(sensitivity_list),len(bias_list)))
    signal_1_sem = np.zeros((len(sensitivity_list),len(bias_list)))
    for i,sensitivity in enumerate(sensitivity_list):
        for j, bias in enumerate(bias_list):
    

            true_params = [sensitivity ,bias]
            labels_params = ['sensitivity' ,'bias']
            # generate observation
            N = 1000
            min_x = -10
            max_x = 10
            x_list = np.linspace(min_x, max_x, 21)
            x_center = (min_x + max_x) / 2 
            
            x0 = np.repeat(x_list, 50)
        
           
            
            # define linear noise model
            m = psychoToyModel( x0, x_center,engine,dim_param = 2)
            
            # generate and show observation
            obs0 = m.gen_single(true_params)
            plt.figure()
            plt.plot(x0,obs0['choice'],'bo')
           
            # define prior
            prior_min = np.array([0.1,-5])
            prior_max = np.array([2,5])
            seed_p = 2
            prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
            
            s1 = psychoToyStats(x_list,x0,option = 'choice_prob')
            s2 = psychoToyStats(x_list,x0,option = 'choice_triggered_signal')
            #g = dg.Default(model=m, prior=prior, summary=s)
            
            # define statsitics summary of observation 
            obs_stats1 = s1.calc([obs0])
            obs_stats2 = s2.calc([obs0])
            signal_0_mean[i,j] = obs_stats2[0][0]
            signal_1_mean[i,j] = obs_stats2[0][1]
            signal_0_sem[i,j] = obs_stats2[0][2]
            signal_1_sem[i,j] = obs_stats2[0][3]
            
    #%%
    plt.figure()
    plt.subplot(221)
    plt.imshow(signal_0_mean)
    plt.xlabel('bias')
    plt.xticks(np.arange(6),bias_list)
    plt.ylabel('sensitivity')
    plt.yticks(np.arange(5),sensitivity_list)
    plt.colorbar()
    plt.title('signal_0_mean')
    
    plt.subplot(222)
    plt.imshow(signal_1_mean)
    plt.xlabel('bias')
    plt.xticks(np.arange(6),bias_list)
    plt.ylabel('sensitivity')
    plt.yticks(np.arange(5),sensitivity_list)
    plt.colorbar()
    plt.title('signal_1_mean')
    
    plt.subplot(223)
    plt.imshow(signal_0_sem)
    plt.xlabel('bias')
    plt.xticks(np.arange(6),bias_list)
    plt.ylabel('sensitivity')
    plt.yticks(np.arange(5),sensitivity_list)
    plt.colorbar()
    plt.title('signal_0_sem')
    
    plt.subplot(224)
    plt.imshow(signal_1_sem)
    plt.xlabel('bias')
    plt.xticks(np.arange(6),bias_list)
    plt.ylabel('sensitivity')
    plt.yticks(np.arange(5),sensitivity_list)
    plt.colorbar()
    plt.title('signal_1_sem')
    
    plt.tight_layout()
    # plt.figure()
    # plt.subplot(121)
    # plt.plot(obs_stats[0])
    # plt.subplot(122)
    
    # #%%
    
    # # set hyparameters 
    # # training schedule
    # n_train = 5000
    # n_rounds = 2
    # seed_inf = 1
    # pilot_samples = 2000


    # val_frac = 0.05
    # # network setup
    # n_hiddens = [50,50]
    # minibatch = 500
    # epochs = 100
    
    # prior_norm = True
    
    # # MAF parameters
    # density = 'mog'
    # n_mades = 5         # number of MADES
    



    # # inference object
    # res = infer.SNPEC(g,
    #                 obs=obs_stats,
    #                 n_hiddens=n_hiddens,
    #                 seed=seed_inf,
    #                 pilot_samples=pilot_samples,
    #                 n_mades=n_mades,
    #                 prior_norm=prior_norm,
    #                 density=density)
    
    # # train
    # log, _, posterior = res.run(
    #                     n_train=n_train,
    #                     n_rounds=n_rounds,
    #                     minibatch=minibatch,
    #                 epochs=epochs,
    #                 silent_fail=False,
    #                 proposal='gaussian',
    #                 val_frac=val_frac,
    #                 verbose=True,)
    
    
    # #%%
    # fig = plt.figure(figsize=(15,5))

    # plt.plot(log[0]['loss'],lw=2)
    # plt.xlabel('iteration')
    # plt.ylabel('loss');
    
    
    
    


    # prior_min = g.prior.lower
    # prior_max = g.prior.upper
    # prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)
    
    # posterior_samples = posterior[0].gen(10000)
    
    # ###################
    # # colors
    # hex2rgb = lambda h: tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    
    # # RGB colors in [0, 255]
    # col = {}
    # col['GT']      = hex2rgb('30C05D')
    # col['SNPE']    = hex2rgb('2E7FE8')
    # col['SAMPLE1'] = hex2rgb('8D62BC')
    # col['SAMPLE2'] = hex2rgb('AF99EF')
    
    # # convert to RGB colors in [0, 1]
    # for k, v in col.items():
    #     col[k] = tuple([i/255 for i in v])
    
    # ###################
    # # posterior
    # fig, axes = samples_nd(posterior_samples,
    #                        limits=prior_lims,
    #                        ticks=prior_lims,
    #                        labels=labels_params,
    #                        fig_size=(5,5),
    #                        diag='kde',
    #                        upper='kde',
    #                        hist_diag={'bins': 50},
    #                        hist_offdiag={'bins': 50},
    #                        kde_diag={'bins': 50, 'color': col['SNPE']},
    #                        kde_offdiag={'bins': 50},
    #                        points=[true_params],
    #                        points_offdiag={'markersize': 5},
    #                        points_colors=[col['GT']],
    #                        title='');
        