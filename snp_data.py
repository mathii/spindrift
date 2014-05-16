# Class for loading and manipulating SNP data. 

from __future__ import division, print_function
import numpy as np
import sys
import pdb

###########################################################################
# datatype definitions
dt_snp1=np.dtype([("ID", np.str_, 16), ("CHR", np.str_, 2), ("POS", np.int32)])
dt_snp2=np.dtype([("ID", np.str_, 16), ("CHR", np.str_, 2), ("POS", np.int32), 
                  ("REF", np.str_, 1), ("ALT", np.str_, 1)])
dt_ind=np.dtype([("IND", np.str_, 32), ("POP", np.str_, 32)])

###########################################################################

class snp_data:
    """
    Class for loading and manipulating snp data
    virtual base class - derived classes need to implement 
    the load function. 
    """
    def __init__(self, file_root, pops=None, no_gt=True, inbred=[]):
        """
        File root is the root of an eigenstrat file. 
        """
        ind, snp, geno, pops=self.load(file_root, pops)

        self.ind=ind
        self.snp=snp
        self.pops=pops
        self.geno=geno

        self.add_population_counts(inbred=inbred)
        if no_gt:                     # Only keep the genotypes if we want them
            del self.geno
        self.remove_monomorphic()
        self.add_population_freqs()

    def load(self, file_root, pops):
        """
        Virtual
        """
        raise Exception("Instantiate using a derived class")

    def add_population_counts(self, inbred=[]):
        """
        Compute population counts. 
        If the population is inbred, we pick a random allele.  
        """
        npop=len(self.pops)
        
        count=np.zeros(shape=(self.geno.shape[0], npop), dtype=np.dtype('i4'))
        total=np.zeros(shape=(self.geno.shape[0], npop), dtype=np.dtype('i4'))
        
        for j, pop in enumerate(self.pops):
            sub_data=self.geno[:,self.ind["POP"]==pop]
            count_0=np.sum(sub_data==0, axis=1)
            count_1=np.sum(sub_data==1, axis=1)
            count_2=np.sum(sub_data==2, axis=1)
            if pop in inbred:
                count_ps_1=np.array([np.random.binomial(c1, 0.5) if 
                                     c1>0 else 0 for c1 in count_1])
                count[:,j]=count_ps_1+count_2
                total[:,j]=count_0+count_1+count_2
            else:
                count[:,j]=count_1+2*count_2
                total[:,j]=2*(count_0+count_1+count_2)
            
        bad=(total<2).any(axis=1)
        count=count[~bad]
        total=total[~bad]
        print("Removed " +str(sum(bad)) + 
              " SNPs with <2 alleles in one population", 
              file=sys.stderr)

        self.count=count
        self.total=total
        
    def remove_monomorphic(self):
        """
        Remove monomorphic snps - have to have compuated counts
        """
        monomorphic = np.logical_or(np.all(self.count==0, axis=1),   
                                    np.all(self.count==self.total, axis=1))

        self.snp=self.snp[~monomorphic]
        self.count=self.count[~monomorphic,:]
        self.total=self.total[~monomorphic,:]
        if hasattr(self, "geno"):
            self.geno=self.geno[~monomorphic,:]
    
        print("Removed "+str(sum(monomorphic))+" monomorphic SNPs", file=sys.stderr)


    def add_population_freqs(self):
        """
        Compute population frequencies - can add adjsted frequencies normalised
        to mean. Need to call population_count first. 
        """
        self.freq=self.count/self.total

        # transform
        mean=np.sum(self.count, axis=1)/np.sum(self.total, axis=1)
        self.tfreq=self.freq.copy()
        for i in range(len(self.pops)):
            self.tfreq[:,i]=(self.freq[:,i]-mean)/(mean*(1-mean))
            
    
###########################################################################
# END CLASS

class eigenstrat_data(snp_data):
    """
    Unpacked Eigenstrat data
    """
        
    def load(self, file_root, pops=None):
        """
        Load from an unpacked eigenstrat file into our internal format.
        If pops is specified then load only the specified populations. 
        """
        # Read .ind file
        ind=np.genfromtxt(file_root+".ind", dtype=dt_ind, usecols=(0,2))   # ignore sex

        include=np.ones(len(ind), dtype=bool)
        if pops:
            include=np.in1d(ind["POP"], pops)
            ind=ind[include]
        else:
            pops=np.unique(ind["POP"])
            
        # Read .snp file
        snp_file=open(file_root+".snp", "r")
        line=snp_file.readline()
        bits=line.split()
        snpdt=dt_snp1                     # does the snp file have the alleles in?
        snpcol=(0,1,3)
        if len(bits) not in [4,6]:
            raise Exception("Could not read snp file")
        elif len(bits)==6:
            snpdt=dt_snp2
            snpcol=(0,1,3,4,5)

        snp_file.seek(0)
        snp=np.genfromtxt(snp_file, dtype=snpdt, usecols=snpcol)
        snp_file.close()

        # Finally, read .geno file
        geno=np.genfromtxt(file_root+".geno", dtype='i1', delimiter=1, 
                           usecols=np.where(include)[0])

        return ind,snp,geno,pops
        
###########################################################################
# END CLASS
