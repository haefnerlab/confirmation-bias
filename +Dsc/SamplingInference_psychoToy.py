#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:41:35 2020

@author: liushizhao
"""

import numpy as np
import matplotlib.pyplot as plt
import pymc3 as pm


min_x = -10
max_x = 10
N = 20 # 10 trials at each orientation
x = np.linspace(min_x, max_x, 21) # the experiment is conducted at 21 x's

x_center = (min_x + max_x)/ 2
# t = np.linspace(-10, 10, 101)
bias_true = .8
sensitivity_true = 1.2

true_params = [sensitivity_true, bias_true]
labels_params = ['sensitivity', 'bias']
    
p_true = (1/(1 + np.exp(-(bias_true+sensitivity_true*(x-x_center)) )))
np.random.seed(1234)
choices = (np.random.uniform(size=[N, np.size(p_true)]) < p_true).T * 1
observed_choices = np.sum(choices, 1)
p_empi = np.mean(choices,axis = 1)

plt.figure()
plt.plot(p_true)
plt.plot(p_empi)
#%%
psychoToyModel = pm.Model()
with psychoToyModel:
    sensitivity = pm.Uniform('sensitivity', lower = .2, upper = 5)
    bias = pm.Uniform('bias', lower = -5, upper = 5)
    
    p = (1/(1 + np.exp(-(bias+sensitivity*(x-x_center)) )))
    
    Y_obs = pm.Binomial('choice',n = N,p = p,observed = observed_choices)
#%%
N = 5000
with psychoToyModel:
    # draw posterior samples
    trace = pm.sample(N)
pm.traceplot(trace)