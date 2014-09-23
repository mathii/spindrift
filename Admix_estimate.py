# Classes for making inference about admixture dates and proportions

from __future__ import division, print_function
import numpy as np
import scipy as sp
import scipy.optimize as opt
import sys
import pdb

###########################################################################

def admix_objective(par, x,y,z):
    """
    par = (alpha, tau1, phi1)
    """
    a, p, t = par

    t1=1+(1-a)*(1-np.exp(-t))-2*x
    t2=1+a*(1-np.exp(-p))-2*y
    t3=5-a*np.exp(-2*p)-(1-a)*np.exp(-2*t)-10*z
    return np.array([t1,t2,t3])

###########################################################################
    
class patterson_estimator:
    """
    This implements the estimator described by Nick Patterson in his 
    technical note. Data should have three populations and we are 
    modelling the third as an admixture of the first two. We are 
    assuming that the reference alelle is the ancestral allele. 
    """

    def __init__(self, data, order):
        """
        Not much
        """
        data.check()
        if not len(data.pops)==3:
            raise Exception("Patterson estimator needs exactly three populations")
        self.data=data
        self.order=order

    def estimate_admix(self):
        """
        Estimate the admixture proportions and times. 
        this isn't exactly what Nick proposed, but
        I think it's in the spirit
        """

        p=1-self.data.freq                # Derived allele frequency
        pops=self.data.pops

        idx=[]
        for ip in pops: idx.append(np.where(self.order==ip)[0][0])
        idx=np.array(idx)

        pops=pops[idx]
        p=p[:,idx]
                
        w1=2*p[:,0]*(1-p[:,0])*p[:,1]     # 01||1
        w2=2*p[:,1]*(1-p[:,1])*p[:,0]     # 1||01
        w3=4*p[:,0]*(1-p[:,0])*p[:,1]*(1-p[:,1])     # 01||01

        x=np.sum(p[:,2]*w1)/np.sum(w1)
        y=np.sum(p[:,2]*w2)/np.sum(w2)
        z=np.sum(p[:,2]*w3)/np.sum(w3)
        
        soln=opt.root(admix_objective, (0.2,0.1,0.1), args=(x,y,z))
        if not soln["success"]:
            print( "Root finder failed to converge; solving least squares", file=sys.stderr)
            soln=opt.root(admix_objective, (0.2,0.1,0.1), args=(x,y,z), method="lm")

        if soln["success"]:
            pars=soln["x"]
        else:
            raise Exception("Root finder and least square fit both failed to converge")

        print( "\nFor %s = a %s + (1-a) %s\n"%(pops[2], pops[0],
                                               pops[1]), file=sys.stderr)
        print("a = %1.4f\np1= %1.4f\nt1= %1.4f\n"%tuple(pars), file=sys.stderr)

        pdb.set_trace()
                
        return pars
        
###########################################################################
# END CLASS
    
