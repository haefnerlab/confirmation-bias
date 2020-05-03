import numpy as np

class analyticalLinearNoise(object):
    
    def __init__(self, x, err, init_params, prior):
        self.x = x
        self.init_params = init_params
        self.N = x.shape[0]
        self.err = err
        self.y = self.x * init_params[0] + init_params[1] + self.err * np.random.randn(self.N)
        self.prior = prior
        
    def linearNoiseSimulator(self, params, seed=None):

        alpha = params[0]
        beta = params[1]
        
        y = self.x * alpha + beta + self.err * np.random.randn(self.N)
        
        return y
    
    def posteriorLinearNoise(self, params):
        a = params[0]
        b = params[1]
        
        dy = self.y - a * self.x - b
        ivar = 1 / self.err**2
        
        lnlikeli =  -0.5 * (self.N* np.log(2*np.pi) + np.sum(2*np.log(self.err)) + np.sum(dy**2 * ivar))
        lnprior_a = -.5 * (a - self.prior['mu_alpha'])**2 / self.prior['var_alpha']**2
        lnprior_b = -.5 * (b - self.prior['mu_beta'])**2  / self.prior['var_beta']**2
        lnposterior = lnlikeli + lnprior_a + lnprior_b
        
        if np.isnan(lnposterior):
            lnposterior = -np.inf
            
        return lnposterior
    
    def __call__(self, params):
        return self.posteriorLinearNoise(params)