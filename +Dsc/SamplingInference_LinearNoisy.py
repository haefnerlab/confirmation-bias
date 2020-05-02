#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:22:47 2020

@author: liushizhao
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

# slope (alpha), intercept (beta), variance (sigma)
true_params = np.array([0.5, -2, 0.2])

N = 1000
t = np.linspace(0, 10, 2)
x = np.random.uniform(0, 10, N) # we are limiting the data generation from 0 to 10
y = x * true_params[0] + true_params[1]
y_obs = y + true_params[-1] * np.random.randn(N)

plt.plot(x, y_obs, ".k", label="observations")
plt.plot(t, true_params[0]*t + true_params[1], label="truth")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(fontsize=14);
#%%
# define linear 
import pymc3 as pm

linearNoisyModel = pm.Model()

with linearNoisyModel:
    
    alpha = pm.Normal('alpha', mu=0, sigma=10)
    beta = pm.Normal('beta', mu=0, sigma=10)
    sigma = pm.HalfNormal('sigma', sigma=1)
    
    
    mu = alpha * x + beta
#     mu = 
    # Likelihood (sampling distribution) of observations
    Y_obs = pm.Normal('y_obs', mu = mu , sigma=sigma, observed=y_obs)
#%%
N = 1000
with linearNoisyModel:
    # draw posterior samples
    trace = pm.sample(N)
pm.traceplot(trace)