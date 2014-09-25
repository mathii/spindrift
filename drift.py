# Functions for doing drift-related calculations. 
from __future__ import division, print_function
import numpy as np
import pdb

###########################################################################

def estimate_covariance(freq_data, T):
    """
    Using the notation of the Berg and Coop paper
    """
    G=np.array(freq_data, dtype=float)
    G=G.T

    M,K = G.shape

    eps=np.mean(G, axis=0)

    var2=np.expand_dims(1/eps*(1-eps),axis=0)
    TGs=T.dot(G*var2)
    F=TGs.dot(TGs.T)/(K-1)
    return F
            
###########################################################################

def Fst_matrix(freq, pops):
    """
    Method by Nick Patterson. Should be equivalent to the inbreed=True
    Option in Eigenstrat/smartpca
    """
    npops=len(pops)
    Fst=np.zeros( (npops, npops), dtype="float")
    
    for i in range(npops-1):
        for k in range(i+1, npops):
            s=freq["total"][:,i]
            t=freq["total"][:,k]
            u=freq["count"][:,i]
            v=freq["count"][:,k]

            EX=np.mean(np.square(np.true_divide(u,s)-np.true_divide(v,t)))
            Eh1=np.mean(np.true_divide(u*(s-u), s*(s-1)))
            Eh2=np.mean(np.true_divide(v*(t-v), t*(t-1)))
            Eh1s=np.mean(np.true_divide(u*(s-u), s*s*(s-1)))
            Eh2t=np.mean(np.true_divide(v*(t-v), t*t*(t-1)))

            Nhat=EX-Eh1s-Eh2t
            Dhat=Nhat+Eh1+Eh2

            # if Nhat/Dhat<0:
            #     pdb.set_trace()
            
            Fst[i,k]=Fst[k,i]=Nhat/Dhat

    return {"Fst":Fst, "pops":pops}

###########################################################################
