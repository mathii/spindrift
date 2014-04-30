# Generic i/o functions. 
from __future__ import division, print_function
import numpy as np
import pdb

###########################################################################

def read_eigenstrat(root, pops=None):
    """
    Read a genotype file in  ancestrymap format into a numpy array.
    Assume that the geno, ind and snp file are named root.geno root.ind
    and root.snp respectively. 
    """

    sample_pop=[]
    include=[]

    INDS=[]
    SNPS=[]
    CHR=[]
    POS=[]
    REF=[]
    ALT=[]
    
    # Load individuals, only those in pops if it is specified
    ind=open(root+".ind", "r")
    n_ind=0
    for i, line in enumerate(ind):
        ind, sex, this_pop=line.split()

        if not pops or this_pop in pops:
            n_ind+=1
            include.append(i)
            sample_pop.append(this_pop)
            INDS.append(ind)

    # Read all snps
    snp=open(root+".snp", "r")
    n_snp=0

    alleles=0                             # +1 - include, -1: not included
    for line in snp:
        if not alleles:
            if len(line.split())==6:
                alleles=1
            elif len(line.split())==4:
                alleles=-1
        if  alleles==1:
            name, chrom, gmap, pos, ref, alt = line.split()
            REF.append(ref)
            ALT.append(alt)
        elif alleles==-1:
            name, chrom, gmap, pos = line.split()
        else:
            raise Exception("Neither 4 nor 6 columns in .snp file")
            
        SNPS.append(name)
        CHR.append(chrom)
        POS.append(int(pos))
        n_snp+=1

    # Specified genotypes. 
    geno=np.genfromtxt(root+".geno", dtype='i1', delimiter=1, usecols=include)

    return {"GT":geno, "POPS":pops, "INDS":INDS, "SNPID":np.array(SNPS), 
            "CHR":np.array(CHR), "POS":np.array(POS), "IPOP":np.array(sample_pop),
            "REF":np.array(REF), "ALT":np.array(ALT)}

###########################################################################


###########################################################################
