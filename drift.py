# Functions for doing drift-related calculations. 
from __future__ import division, print_function
import numpy as np
import sys
import pdb

###########################################################################

def population_count(data, inbred=[]):
    """
    Compute population counts. 
    If the population is inbred, we pick a random allele.  
    """
    npop=len(data["POPS"])
    pops=data["POPS"]
    
    count=np.zeros(shape=(data["GT"].shape[0], npop), dtype=np.dtype('i4'))
    total=np.zeros(shape=(data["GT"].shape[0], npop), dtype=np.dtype('i4'))

    for j, pop in enumerate(data["POPS"]):
        sub_data=data["GT"][:,data["IPOP"]==pop]
        count_0=np.sum(sub_data==0, axis=1)
        count_1=np.sum(sub_data==1, axis=1)
        count_2=np.sum(sub_data==2, axis=1)
        if pop in inbred:
            count_ps_1=np.array([np.random.binomial(c1, 0.5) if c1>0 else 0 for c1 in count_1])
            count[:,j]=count_ps_1+count_2
            total[:,j]=count_0+count_1+count_2
        else:
            count[:,j]=count_1+2*count_2
            total[:,j]=2*(count_0+count_1+count_2)
            
    bad=(total<2).any(axis=1)
    count=count[~bad]
    total=total[~bad]
    print("Removed " +str(sum(bad)) + " SNPs with <2 alleles in one population", file=sys.stderr)
    data["COUNT"]=count
    data["TOTAL"]=total

###########################################################################

def remove_monomorphic(data):
    """
    Remove monomorphic snps - have to have compuated counts
    """
    monomorphic = np.logical_or(np.all(data["COUNT"]==0, axis=1),   
                   np.all(data["COUNT"]==data["TOTAL"], axis=1))
    
    for what in ["GT", "COUNT", "TOTAL", "FREQ", "TFREQ"]:
        if what in data.keys():
            data[what]=data[what][~monomorphic,:]

    for what in ["REF", "ALT", "POS", "CHR", "SNPID"]:
        if what in data.keys():
            data[what]=data[what][~monomorphic]
            
    print("Removed "+str(sum(monomorphic))+" monomorphic SNPs", file=sys.stderr)
    
###########################################################################

def population_freq(data, transform=True):
    """
    Compute population frequencies - can add adjsted frequencies normalised
    to mean. Need to call population_count first. 
    """
    data["FREQ"]=data["COUNT"]/data["TOTAL"]
    if transform:
        mean=np.sum(data["COUNT"], axis=1)/np.sum(data["TOTAL"], axis=1)
        data["TFREQ"]=data["FREQ"].copy()
        for i in range(len(data["POPS"])):
            data["TFREQ"][:,i]=(data["FREQ"][:,i]-mean)/(mean*(1-mean))

###########################################################################

def estimate_covariance(freq_data):
    """
    Using the notation of the Berg and Coop paper
    """
    G=np.array(freq_data, dtype=float)
    G=G.T

    M,K = G.shape

    T=np.zeros((M-1,M), dtype=float)-1/M
    np.fill_diagonal(T, (M-1)/M)


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

