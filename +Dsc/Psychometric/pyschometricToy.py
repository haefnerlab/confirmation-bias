import numpy as np
import scipy.stats as st

class analyticalPsychometric(object):
    
    def __init__(self, x, N, x_center, init_params, prior):
        
        self.x = x
        self.N = N
        self.x_center = x_center
        self.init_params = init_params
        self.y = self.psychometricSimulator(self.init_params)
        self.prior = prior
        
    def psychometricSimulator2(self, params, seed=None):
        
        plogistic = 1/(1 + np.exp(-(params[1]+params[0]*(self.x - self.x_center)) ))
        choices = (np.random.uniform(size=[self.N, np.size(plogistic)]) < plogistic).T * 1
        observed_choices = np.sum(choices, 1)

        return observed_choices
    
    def psychometricSimulator(self, params, seed=None):
        """
        This returns the binary choices (N x size(number of stimulus values))
        """
        plogistic = 1/(1 + np.exp(-(params[1]+params[0]*(self.x - self.x_center)) ))
        choices = (np.random.uniform(size=[self.N, np.size(plogistic)]) < plogistic).T * 1 # N x size(x) array
        # observed_choices = np.sum(choices, 1)

        return choices
    
    def posteriorPsychometric(self, params, observed_choices):
        
        sensitivity = params[0]
        bias = params[1]
        
        plogistic = 1/(1 + np.exp(-( bias + sensitivity * (self.x - self.x_center) ) ))
        
        # choices = (np.random.uniform(size=[self.N, np.size(plogistic)]) < plogistic).T * 1
        # observed_choices = np.sum(choices, 1)
        lnlikeli = np.sum(st.binom.logpmf(observed_choices, self.N, plogistic))
        
        if sensitivity <= 0:
            return -np.inf # sensitivity cannot be negative
        else: 
            lnprior_s = -.5 * (sensitivity - self.prior['mu_sensitivity'])**2 / self.prior['var_sensitivity']**2
        
        lnprior_b = -.5 * (bias - self.prior['mu_bias'])**2  / self.prior['var_bias']**2
        
        lnposterior = lnlikeli + lnprior_s + lnprior_b
        
        if np.isnan(lnposterior):
            lnposterior = -np.inf
            
        return lnposterior
    
    def __call__(self, params):
        return self.posteriorPsychometric(params)