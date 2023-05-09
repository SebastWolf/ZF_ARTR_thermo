#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: modified from jerometubiana
"""

import numpy as np
import layer
import pgm
import utilities
from utilities import check_random_state
import cy_utilities
import time,copy
import itertools


class BM(pgm.PGM):
    def __init__(self,N = 100,nature = 'Bernoulli', n_c =1, random_state = None, gauge='zerosum',zero_field=False):
        self.N = N
        self.nature = nature
        self.random_state = check_random_state(random_state)
        if self.nature == 'Potts':
            self.n_c = n_c            
        else:
            self.n_c = 1
        self.zero_field = zero_field

            
        super(BM, self).__init__(n_layers = 1, layers_size = [self.N],layers_nature = [self.nature + '_coupled'], layers_n_c = [self.n_c] , layers_name = ['layer'])
        
                
        self.gauge = gauge

        self.layer = layer.Layer(N= self.N, nature = self.nature + '_coupled', position = 'visible', n_c = self.n_c, random_state = self.random_state, zero_field = self.zero_field,gauge=self.gauge)
        self.init_couplings(0.01)
        
    def init_couplings(self,amplitude):
        if (self.n_c >1):
            self.layer.couplings = amplitude * self.random_state.randn(self.N, self.N,self.n_c, self.n_c)
            self.layer.couplings = pgm.gauge_adjust_couplings(self.layer.couplings,self.n_c,self.n_c,gauge=self.gauge)
        else:
            self.layer.couplings = amplitude * self.random_state.randn(self.N,self.N)
            
            
    def markov_step_PT(self,(x,fields_eff)):
        for i,beta in zip(np.arange(self.N_PT), self.betas):
            x[i],fields_eff[i] = self.markov_step((x[i],fields_eff[i]),beta =beta)
        return (x,fields_eff)
    
    def markov_step_PT2(self,(x,fields_eff),E):
        for i,beta in zip(np.arange(self.N_PT), self.betas):
            (x[i],fields_eff[i]),E[i] = self.layer.sample_and_energy_from_inputs(None,beta=beta,previous=(x[i],fields_eff[i]))
        return (x,fields_eff),E



    
        
    def exchange_step_PT(self,(x,fields_eff),E,record_acceptance=True,compute_energy=True):
        if compute_energy:
            for i in np.arange(self.N_PT):
                E[i,:] = self.energy(x[i,:,:],remove_init=True)
        
        if self.record_swaps:
            particle_id = self.particle_id[-1].copy()
        for i in np.arange(self.count_swaps%2,self.N_PT-1,2):  
            proba = np.minimum( 1,  np.exp( (self.betas[i+1]-self.betas[i]) * (E[i+1,:]-E[i,:])   ) )
            swap = np.random.rand(proba.shape[0]) < proba
            if i>0:
                x[i:i+2,swap] = x[i+1:i-1:-1 ,swap]
                fields_eff[i:i+2,swap] = fields_eff[i+1:i-1:-1 ,swap]
                E[i:i+2,swap] = E[i+1:i-1:-1,swap]
                if self.record_swaps:
                    particle_id[i:i+2,swap] = particle_id[i+1:i-1:-1,swap]
                    
            else:
                x[i:i+2,swap] = x[i+1::-1,swap]
                fields_eff[i:i+2,swap] = fields_eff[i+1::-1 ,swap]
                E[i:i+2,swap] = E[i+1::-1,swap]
                if self.record_swaps:
                    particle_id[i:i+2,swap] = particle_id[i+1::-1,swap]

            
            if record_acceptance:
                self.acceptance_rates[i] = swap.mean()
                self.mav_acceptance_rates[i] = self.mavar_gamma * self.mav_acceptance_rates[i] +  self.acceptance_rates[i]*(1-self.mavar_gamma)
                
        if self.record_swaps:
            self.particle_id.append(particle_id)

        self.count_swaps +=1                
        return (x,fields_eff),E
    
    
    def update_betas(self,beta=1):
        super(BM,self).update_betas(beta=beta)
    
    
    def markov_step(self,(x,fields_eff),beta =1):
         (x,fields_eff) = self.layer.sample_from_inputs( None,previous = (x,fields_eff), beta = beta)
         return (x,fields_eff)

    def markov_step_and_energy(self, (x,fields_eff),E, beta=1):
        (x,fields_eff),E = self.layer.sample_and_energy_from_inputs(None,beta=beta,previous=(x,fields_eff))
        return (x,fields_eff),E         
    
    def compute_fields_eff(self,x,beta=1):
        if x.ndim ==1: x = x[np.newaxis,:]
        return self.layer.compute_output(x,self.layer.couplings) + self.layer.fields[np.newaxis]

    
    
    def energy(self,x,remove_init = False):
        if x.ndim ==1: x = x[np.newaxis,:]
        return self.layer.energy(x,remove_init = remove_init)
    
    def free_energy(self,x):      
        return self.energy(x)

    def compute_all_moments(self):
       if self.nature == 'Bernoulli':
           string = ','.join(['[0,1]' for _ in range(self.N)])
           exec('iter_configurations=itertools.product(%s)'%string)
       else:
           print 'no supported'
       configurations = np.array([config for config in iter_configurations])
       weights = -self.free_energy(configurations)
       maxi = weights.max()
       weights -= maxi
       weights = np.exp(weights)
       Z = weights.sum()
       if self.nature == 'Potts':
           mean = np.zeros([self.N,self.n_c] )
           for color in range(self.n_c):
               mean[:,color] = ((configurations == color) * weights[:,np.newaxis]).sum(0)/Z
           covariance = np.zeros([self.N,self.N,self.n_c,self.n_c] )
           for color1 in range(self.n_c):
               for color2 in range(self.n_c):
                   covariance[:,:,color1,color2] = ( (configurations[:,np.newaxis,:]==color2) * (configurations[:,:,np.newaxis]==color1) * weights[:,np.newaxis,np.newaxis]).sum(0)/Z
       else:
           mean = (configurations * weights[:,np.newaxis]).sum(0)/Z
           covariance = ( configurations[:,np.newaxis,:] * configurations[:,:,np.newaxis] * weights[:,np.newaxis,np.newaxis]).sum(0)/Z

       Z = Z * np.exp(maxi)
       return Z,mean,covariance

                
    def gen_data(self, Nchains = 10, Lchains = 100, Nthermalize = 0 ,Nstep = 1, N_PT =1, config_init = [], beta = 1,batches = None,reshape = True,record_replica = False, record_acceptance=None, update_betas = None,record_swaps=False):
        return super(BM,self).gen_data(Nchains = Nchains,Lchains=Lchains,Nthermalize=Nthermalize,Nstep=Nstep, N_PT = N_PT, config_init=config_init, beta =beta, batches = batches,reshape=reshape,record_replica=record_replica,record_acceptance=record_acceptance,update_betas=update_betas,record_swaps=record_swaps)


    
        