#%% Useful functions
import numpy as np
import numbers
from scipy.special import erf,erfinv
from scipy.sparse import csr_matrix
import cy_utilities
import itertools

def check_random_state(seed):
    if seed==None or seed==np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)





def cumulative_probabilities(X,maxi=1e9):
    out = np.exp(X)
    out[out>maxi] = maxi # For numerical stability.
    out = np.cumsum(out,axis=-1)
    out/= out[:,:,-1][:,:,np.newaxis]
    return out



def kronecker(X1,c1):
    l = []
    for color in range(c1):
        l.append( np.asarray(X1 == color,dtype=float))
    return np.rollaxis(np.array(l),0, X1.ndim+1)

def saturate(x,xmax):
    return np.sign(x) * np.minimum(np.abs(x),xmax)


#def average(X,c=1):
#    if (c==1):
#        return X.mean(0)
#    else:
#        out = np.zeros([X.shape[1],c])
#        for color in range(c):
#            out[:,color] = csr_matrix(X==color).mean(0)
#        return out

def average(X,c=1,weights=None):
    if (c==1):
        if weights is None:
            return X.mean(0)
        else:
            if X.ndim ==1:
                return (X * weights).sum(0)/weights.sum()
            elif X.ndim ==2:
                return (X * weights[:,np.newaxis]).sum(0)/weights.sum()
            elif X.ndim ==3:
                return (X * weights[:,np.newaxis,np.newaxis]).sum(0)/weights.sum()
    else:
        if weights is None:
            return cy_utilities.average_C(X,c)
        else:
            return cy_utilities.weighted_average_C(X,weights,c)



def average_product(X1, X2, c1 = 1,c2=1, mean1 = False, mean2= False,weights =None):
    if (c1 == 1) & (c2 == 1):
        if weights is None:
            return np.dot(X1.T,np.asarray(X2,dtype=float))/float(X1.shape[0])
        else:
            return (X1[:,:,np.newaxis] * X2[:,np.newaxis,:] * weights[:,np.newaxis,np.newaxis]).sum(0)/weights.sum()
    elif (c1 == 1) & (c2 <> 1):
        if mean2: # X2 in format [data,site,color]; each conditional mean for each color.
            if weights is None:
                return np.tensordot(X1,X2,axes=([0],[0]))/float(X1.shape[0])
            else:
                return np.tensordot(X1 * weights[:,np.newaxis],X2,axes=([0],[0]))/weights.sum()
        else:
            if weights is None:
                return cy_utilities.average_product_FxP_C(np.asarray(X1,dtype=float),X2,c2)
            else:
                return cy_utilities.average_product_FxP_C(np.asarray(X1*weights[:,np.newaxis],dtype=float),X2,c2) * weights.shape[0]/weights.sum()
    elif (c1 <>1) & (c2 ==1):
        if mean1: # X1 in format [data,site,color]; each conditional mean for each color.
            if weights is None:
                return np.swapaxes( np.tensordot(X1,X2,axes=([0],[0])),1,2)/float(X1.shape[0])
            else:
                return np.swapaxes( np.tensordot(X1,weights[:,np.newaxis,:]*X2,axes=([0],[0])),1,2)/weights.sum()
        else:
            if weights is None:
                return np.swapaxes(cy_utilities.average_product_FxP_C(np.asarray(X2,dtype=float),X1,c1),0,1)
            else:
                return np.swapaxes(cy_utilities.average_product_FxP_C(np.asarray(X2 * weights[:,np.newaxis],dtype=float),X1,c1),0,1) * weights.shape[0]/weights.sum()
    elif (c1 <>1) & (c2 <>1):
        if mean1 & mean2:
            if weights is None:
                return np.swapaxes(np.tensordot(X1,X2,axes = ([0],[0])), 1,2)/float(X1.shape[0])
            else:
                return np.swapaxes(np.tensordot(X1 * weights[:,np.newaxis,np.newaxis],X2,axes = ([0],[0])), 1,2)/weights.sum()
        elif mean1 & (~mean2):
            out = np.zeros([X1.shape[1], X2.shape[1],c1,c2])
            for color2 in range(c2):
                if weights is None:
                    out[:,:,:,color2] = np.swapaxes(np.tensordot(X1, (X2 == color2), axes = ([0],[0]) ), 1,2 )/float(X1.shape[0])
                else:
                    out[:,:,:,color2] = np.swapaxes(np.tensordot(X1 * weights[:,np.newaxis,np.newaxis], (X2 == color2), axes = ([0],[0]) ), 1,2 )/weights.sum()
            return out
        elif (~mean1) & mean2:
            out = np.zeros([X1.shape[1], X2.shape[1],c1,c2])
            for color1 in range(c1):
                if weights is None:
                    out[:,:,color1,:] = np.tensordot( (X1 == color1), X2, axes = ([0],[0]))/float(X1.shape[0])
                else:
                    out[:,:,color1,:] = np.tensordot( (X1 == color1), weights[:,np.newaxis,np.newaxis]* X2, axes = ([0],[0]))/weights.sum()
            return out

        else:
            if weights is None:
                return cy_utilities.average_product_PxP_C(X1,X2,c1,c2)
            else:
                return cy_utilities.weighted_average_product_PxP_C(X1,X2,weights,c1,c2)


def covariance(X1, X2, c1 = 1,c2=1, mean1 = False, mean2= False,weights =None):
    if mean1:
        mu1 = average(X1,weights=weights)
    else:
        mu1 = average(X1,c=c1,weights=weights)
    if mean2:
        mu2 = average(X2,weights=weights)
    else:
        mu2 = average(X2,c=c2,weights=weights)

    prod = average_product(X1,X2,c1=c1,c2=c2,mean1=mean1,mean2=mean2,weights=weights)

    if (c1>1) & (c2>1):
        covariance = prod - mu1[:,np.newaxis,:,np.newaxis] * mu2[np.newaxis,:,np.newaxis,:]
    elif (c1>1) & (c2==1):
        covariance = prod - mu1[:,np.newaxis,:] * mu2[np.newaxis,:,np.newaxis]
    elif (c1==1) & (c2>1):
        covariance = prod - mu1[:,np.newaxis,np.newaxis] * mu2[np.newaxis,:,:]
    else:
        covariance = prod - mu1[:,np.newaxis] * mu2[np.newaxis,:]
    return covariance




def bilinear_form(W, X1,X2,c1=1,c2=1):
    if (c1==1) & (c2==1):
        return np.sum( X1 * np.dot(W,X2.T).T,1)
    elif (c1 ==1) & (c2>1):
        return np.sum(X1 * cy_utilities.compute_output_C(X2.shape[1],c2,X2,W) ,1)
    elif (c1>1) & (c2 ==1):
        return cy_utilities.dot_Potts2_C(X1.shape[1],c1,X1, np.tensordot(X2,W,(1,1) ) )
    elif (c1>1) & (c2>1):
        return cy_utilities.dot_Potts2_C(X1.shape[1], c1, X1, cy_utilities.compute_output_Potts_C(X2.shape[1],c2,c1,X2,W) )


def copy_config(config,N_PT=1,record_replica=False):
    if type(config)==tuple:
        if N_PT>1:
            if record_replica:
                return config[0].copy()
            else:
                return config[0][0].copy()
        else:
            return config[0].copy()
    else:
        if N_PT>1:
            if record_replica:
                return config.copy()
            else:
                return config[0].copy()
        else:
            return config.copy()


def make_all_discrete_configs(N,nature,c=1):
    if nature == 'Bernoulli':
        string = ','.join(['[0,1]' for _ in range(N)])
        exec('iter_configurations=itertools.product(%s)'%string)
    else:
        print 'no supported'
    configurations = np.array([config for config in iter_configurations])
    return configurations
