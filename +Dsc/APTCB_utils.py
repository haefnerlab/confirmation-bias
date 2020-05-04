#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 09:44:29 2020

@author: liushizhao
"""


from delfi.simulator.BaseSimulator import BaseSimulator


from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
import copy
import matplotlib.pyplot as plt
import numpy as np
from delfi.utils.viz import samples_nd
#%%
def ConfirmationBiasSimulator(xval,fieldstofit,params, signals, engine):
    # run matlab code run.vectorized and read results
   

    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        if f == 'samples':
            paramsNew[f] = int(xval[i])
        elif type(params[f]) is int:
            paramsNew[f] = float(xval[i])


    sim_results = engine.Model.runVectorizedAPT(paramsNew, signals.tolist())
    sim_choice = np.squeeze(np.asarray(sim_results['choices']))



    return sim_choice

# define a simulator class linking to the real simulator
class ConfirmationBias(BaseSimulator):
    
    def __init__(self, signals,fieldstofit,params,engine, seed = None):
        
        dim_param = 4
        super().__init__(dim_param = dim_param,seed = seed)
        self.signals = signals
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        choices = self.simulate(xval,self.fieldstofit,self.params,self.signals, self.engine)
        return {'choices':choices}
    
class ConfirmationBiasStats(BaseSummaryStats):
    
    def __init__(self,signals,infoType,condition = 'choices_infoType',useStd = True):
        self.signals = signals
        self.infoType = infoType
        self.condition = condition
        self.useStd = useStd
    def calc(self,repetition_list):

        stats = []
    
        for r in range(len(repetition_list)):

            data = repetition_list[r]
            choices = np.asarray(data['choices'])
            signals = np.asarray(self.signals)
            infoType = np.asarray(self.infoType)
            
            if self.condition == 'choices_infoType':
                # group trials according to both choices and catergory information (4 groups)
                ind1 = np.where( (choices == 1) & (infoType == 1))[0]
                ind2 = np.where( (choices == 1) & (infoType == -1))[0]
                ind3 = np.where( (choices == -1) & (infoType == 1))[0]
                ind4 = np.where( (choices == -1) & (infoType == -1))[0]
                if self.useStd: 
                    # include standard error mean in summary statistics
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.mean(signals[ind3,:],axis = 0),np.mean(signals[ind4,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)),\
                                               np.std(signals[ind3,:],axis = 0)/np.sqrt(len(ind3)),\
                                               np.std(signals[ind4,:],axis = 0)/np.sqrt(len(ind4)))) 
                else:
                    # only use average of signals as summary statistics 
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                   np.mean(signals[ind3,:],axis = 0),np.mean(signals[ind4,:],axis = 0)))    
            elif self.condition == 'choices':
                # group trials according to just choices (2 groups)
                ind1 = np.where(choices == 1)[0]
                ind2 = np.where(choices == -1)[0]
                if self.useStd:
                    # include standard error mean in summary statistics
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    # only use average of signals as summary statistics 
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                              np.mean(np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1))),\
                                                  np.mean(np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))))
            elif self.condition == 'infoType':
                # group trials according to just category information (2 groups)
                ind1 = np.where(infoType == 1)[0]
                ind2 = np.where(infoType == -1)[0]
                if self.useStd:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                              np.mean(np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1))),\
                                                  np.mean(np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))))
                    
                    
            stats.append(sum_stats_vec)
        return stats
def makeFigures(log,posterior,true_params,labels_params,prior_min,prior_max):
    
    fig = plt.figure(figsize=(15,5))
    
    plt.plot(log[-1]['loss'],lw=2)
    plt.xlabel('iteration')
    plt.ylabel('loss');
    
    

    prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)
    
    posterior_samples = posterior[-1].gen(10000)
    
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