# -*- coding: utf-8 -*-
"""


@author: modified from jerometubiana
"""

import numpy as np
from scipy.special import erf,erfinv
from scipy.sparse import csr_matrix
import cy_utilities

from utilities import check_random_state,average,average_product



def logistic(x, out=None): ## from SKlearn
    if out==None:
        out = np.empty(np.atleast_1d(x).shape, dtype=np.float64)
    out[:] = x

    out *= .5
    np.tanh(out, out)
    out += 1
    out *= .5

    return out.reshape(np.shape(x))


#%% Layer class

class Layer():
    def __init__(self, N = 100, nature = 'Bernoulli', position = 'visible', n_c = 1, gauge = 'zerosum', random_state = None,zero_field = False):
        self.N = N
        self.nature = nature
        self.random_state = check_random_state(random_state)
        self.position = position
        self.tmp = 0
        self.mu_psi = np.zeros(N) # For batch_norm.

        if self.nature in  ['Bernoulli']:
            self.n_c = 1
            self.fields = np.zeros(N)
            self.fields0 = np.zeros(N) # useful for PT.

        if self.nature in ['Gaussian','ReLU','dReLU','ReLU+']:
            self.type=np.float
        else:
            self.type = np.int

        self.zero_field = zero_field

    def mean_from_inputs(self,psi, beta = 1):
        if self.nature == 'Bernoulli':
            if beta == 1:
                return logistic(psi + self.fields[np.newaxis,:])
            else:
                return logistic(beta* psi + self.fields0[np.newaxis,:] + beta*(self.fields[np.newaxis,:] - self.fields0[np.newaxis,:]) )

        

    def mean2_from_inputs(self,psi, beta = 1):
        if self.nature in ['Bernoulli','Potts']:
            return self.mean_from_inputs(psi,beta=beta)
        elif self.nature =='Spin':
            return np.ones(psi.shape)
       







    



    def var_from_inputs(self,psi,beta =1):
        if self.nature in ['Bernoulli','Potts']:
            mean = self.mean_from_inputs(psi,beta=beta)
            return mean * (1-mean)
        elif self.nature =='Spin':
            return 1- self.mean_from_inputs(psi,beta=beta)**2
       


    def sample_from_inputs(self,psi,beta=1,previous=(None,None)):
        if self.nature == 'Bernoulli':
            return (self.random_state.random_sample(size = psi.shape) < self.mean_from_inputs(psi,beta))
          
        elif self.nature == 'Bernoulli_coupled':
            (x,fields_eff) = previous

            if x is None:
                x = np.random.randint(0,high=2,size=[psi.shape[0],self.N])
            if fields_eff is None:
                fields_eff = self.fields[np.newaxis] + self.compute_output(x,self.couplings)

            if psi is not None:
                x,fields_eff=cy_utilities.Bernoulli_Gibbs_input_C(x.shape[0], self.N,beta,x,fields_eff,psi, self.fields0, self.couplings, self.random_state.randint(0,high=self.N,size=[x.shape[0],self.N]), self.random_state.rand(x.shape[0],self.N) )
            else:
                x,fields_eff=cy_utilities.Bernoulli_Gibbs_free_C(x.shape[0], self.N,beta,x,fields_eff, self.fields0,self.couplings, self.random_state.randint(0,high=self.N,size=[x.shape[0],self.N]), self.random_state.rand(x.shape[0],self.N) )

            return (x,fields_eff)


    def compute_output(self,config, couplings, direction='up'):

        if config.ndim == 1: config = config[np.newaxis, :] # ensure that the config data is a batch, at leas of just one vector
        n_data = config.shape[0]

        if direction == 'up':
            N_output_layer = couplings.shape[0]
            return np.tensordot(config, couplings, axes = (1,1))
        elif direction == 'down':
            N_output_layer = couplings.shape[1]
            return np.tensordot(config, couplings, axes = (1,0))


    def energy(self,config,remove_init = False):
        if config.ndim == 1: config = config[np.newaxis, :] # ensure that the config data is a batch, at leas of just one vector
        if self.nature in ['Bernoulli','Spin']:
            if remove_init:
                return -np.dot(config,self.fields - self.fields0)
            else:
                return -np.dot(config,self.fields)
        
        if self.nature in ['Bernoulli_coupled','Spin_coupled']:
            if remove_init:
                return - np.dot(config,self.fields-self.fields0) - 0.5* (np.dot(config , self.couplings) * config).sum(1)
            else:
                return - np.dot(config,self.fields) - 0.5* (np.dot(config , self.couplings) * config).sum(1)
        else:
            print 'nature not supported, energy'


    def logpartition(self,inputs,beta = 1):
        if inputs is None:
            inputs = np.zeros([1,self.N])
        if self.nature == 'Bernoulli':
            if beta == 1:
                return np.log(1+ np.exp(self.fields[np.newaxis,:] + inputs )).sum(1)
            else:
                return np.log(1+ np.exp( (beta* self.fields + (1-beta) * self.fields0)[np.newaxis,:] + beta * inputs )).sum(1)
        
        elif self.nature == 'Bernoulli_coupled':
            if beta ==0:
                return np.log(1+ np.exp(self.fields0)[np.newaxis,:] + 0 * inputs ).sum(1)
            else:
                print 'mean field equations not implemented for Bernoulli_coupled'
        else:
            print 'hidden type not supported'
            return


    def transform(self,psi):
        if self.nature == 'Bernoulli':
            return (psi+self.fields)>0
       

    def random_init_config(self,n_samples,N_PT=1):
        if self.nature in ['Bernoulli','Spin']:
            if N_PT>1:
                return self.sample_from_inputs( np.zeros([N_PT*n_samples, self.N]) ,beta = 0).reshape([N_PT,n_samples,self.N])
            else:
                return self.sample_from_inputs( np.zeros([n_samples, self.N]) ,beta = 0)

        elif self.nature in ['Bernoulli_coupled']:
            if N_PT>1:
                (x,fields_eff) = self.sample_from_inputs( np.zeros([N_PT*n_samples, self.N]) ,beta = 0)
                x = x.reshape([N_PT,n_samples,self.N])
                fields_eff = fields_eff.reshape([N_PT,n_samples,self.N])
            else:
                (x,fields_eff) = self.sample_from_inputs( np.zeros([N_PT*n_samples, self.N]) ,beta = 0)
            return (x,fields_eff)

    def internal_gradients(self,data_pos,data_neg, l1 = None, l2 = None,weights = None,weights_neg=None,value='data'):
        gradients = {}
        if self.nature in ['Bernoulli']:
            if value == 'data': # data_pos, data_neg are configurations
                mu_pos = average(data_pos,c=self.n_c,weights=weights)
                mu_neg = average(data_neg,c=self.n_c,weights=weights_neg)
            elif value == 'mean':
                mu_pos = average(data_pos,weights=weights)
                mu_neg = average(data_neg,weights=weights_neg)
            elif value == 'input':
                mu_pos = average( self.mean_from_inputs(data_pos),weights=weights)
                mu_neg = average( self.mean_from_inputs(data_neg),weights=weights_neg)

            gradients['fields'] =   mu_pos - mu_neg
            if weights is not None:
                 gradients['fields'] *= weights.sum()/data_pos.shape[0]
        

        elif self.nature in ['Bernoulli_coupled']:
            if value =='data':
                mu_pos = average(data_pos,weights=weights)
                mu_neg = average(data_neg,weights=weights_neg)
                comu_pos = average_product(data_pos, data_pos,weights=weights)
                comu_neg = average_product(data_neg, data_neg,weights=weights_neg)


            elif value == 'mean':
                mu_pos = data_pos[0]
                mu_neg = average(data_neg,weights=weights_neg)
                comu_pos = data_pos[1]
                comu_neg = average_product(data_neg, data_neg,weights=weights_neg)
            elif value == 'input':
                print 'not supported'
            if self.batch_norm:
                gradients['couplings'] = comu_pos - comu_neg - mu_pos[:,np.newaxis] * (mu_pos-mu_neg)[np.newaxis,:] - mu_pos[np.newaxis,:] * (mu_pos-mu_neg)[:,np.newaxis]
                gradients['fields'] = mu_pos-mu_neg - np.dot(gradient['couplings'],mu_pos)
            else:
                gradients['fields'] = mu_pos-mu_neg
                gradients['couplings'] = comu_pos - comu_neg
            if weights is not None:
                gradients['fields'] *= weights.sum()/data_pos.shape[0]
                gradients['couplings'] *= weights.sum()/data_pos.shape[0]

            if l2 is not None:
                gradients['couplings'] -= l2 * self.couplings
            if l1 is not None:
                gradients['couplings'] -= l1 * np.sign(self.couplings)

       




        if (self.position== 'hidden'):
            if self.nature in ['Bernoulli','Bernoulli_coupled']:
                mu_neg0 = self.mean_from_inputs( np.zeros([1,self.N]),beta=0)[0]
                gradients['fields0'] = mu_pos - mu_neg0
                if weights is not None:
                    gradients['fields0'] *= weights.sum()/data_pos.shape[0]
          


        if self.zero_field:
            gradients['fields'] *= 0
        return gradients



    def init_params_from_data(self,X,eps=1e-6,mean=False,weights=None):
        if self.nature in ['Bernoulli','Bernoulli_coupled']:
            if mean:
                mu = X
            else:
                mu = average(X,weights=weights)
            self.fields = np.log((mu+ eps)/(1-mu + eps))
            self.fields0 = self.fields.copy()
        
        if self.nature in ['Bernoulli_coupled']:
            self.couplings *=0
            self.couplings0 *=0

        if self.zero_field:
            self.fields *=0
            self.fields0 *=0


    def sample_and_energy_from_inputs(self,psi,beta=1,previous=None):
        if self.nature in ['Bernoulli_coupled']:
            (x,fields_eff) = self.sample_from_inputs(psi,beta=beta,previous=previous)
            if psi is None:
                    f = 0.5* (self.fields[np.newaxis] + fields_eff) - self.fields0[np.newaxis]
            else:
                    f = 0.5* (self.fields[np.newaxis] + fields_eff) + psi-self.fields0[np.newaxis]

        else:
            config = self.sample_from_inputs(psi,beta=beta)
            energy = self.energy(config,remove_init=True) - (config*psi).sum(1)
        return config,energy
    
