# Class for loading and manipulating SNP data. 

from __future__ import division, print_function
import numpy as np
import sys
import pyEigenstrat as pE
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
    def __init__(self, file_root, pops=None, no_gt=True, inbred=[], sparse=2, 
                 snps=None, freqs=True, inds=None, exclude_inds=None, remove_monomorphic=True):
        """
        File root is the root of an eigenstrat file. 
        """
        self.load(file_root, pops=pops, inds=inds, exclude_inds=exclude_inds, snps=snps)
        self.add_population_counts(inbred=inbred)
        if hasattr(self, "geno") and no_gt:  # Only keep the genotypes if we want them
            del self.geno
        if remove_monomorphic:
            self.remove_monomorphic()
        if sparse:
            self.remove_sparse(n=sparse)
        self.add_population_freqs()

    def load(self, file_root, pops, snps, inds):
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
            
        self.count=count
        self.total=total

    def check(self):
        """
        Check that everything is the right shape and so on
        """
        NSNP=self.snp.shape[0]
        NPOP=len(self.pops)
        NIND=len(self.ind)

        if hasattr(self, "count"):
            if self.count.shape[0]!=NSNP:
                raise Exception("Count has wrong number of snps")
            if self.count.shape[1]!=NPOP:
                raise Exception("Count has wrong number of populations")
            if self.count.shape!=self.total.shape:
                raise Exception("Count and total out of shape")

        for what in ["geno", "alt_count", "ref_count"]:
            if hasattr(self, what):
                if self.geno.shape[0]!=NSNP:
                    raise Exception("%s has wrong number of snps"%(what,))
                if self.geno.shape[1]!=NIND:
                    raise Exception("%s has wrong number of individuals"%(what,))
            
        if hasattr(self, "freq"):
            if self.freq.shape[0]!=NSNP:
                raise Exception("Frequency has wrong number of snps")
            if self.freq.shape[1]!=NPOP:
                raise Exception("Frequency has wrong number of populations")

    def filter_snps(self, include):
        """
        Include should be a numpy vector of length equal to 
        the number of snps - boolean, defines whether the 
        SNPs are retained. 
        """
        self.snp=self.snp[include]
        self.count=self.count[include,:]
        self.total=self.total[include,:]
        for what in ["geno", "freq", "ref_count", "alt_count"]:
            if hasattr(self, what):
                setattr(self, what, getattr(self, what)[include,:])

        self.check()
        
    def remove_sparse(self, n=2):
        """
        Remove any snps with less that n alleles in any population
        i.e too much missing data. 
        """
        bad=(self.total<n).any(axis=1)
        self.filter_snps(~bad)
        
        print("Removed " +str(sum(bad)) + 
              " SNPs with <" + str(n) + " alleles in one population", 
              file=sys.stderr)

    def remove_monomorphic(self):
        """
        Remove monomorphic snps - have to have compuated counts
        """
        monomorphic = np.logical_or(np.all(self.count==0, axis=1),   
                                    np.all(self.count==self.total, axis=1))

        self.filter_snps(~monomorphic)
        print("Removed "+str(sum(monomorphic))+" monomorphic SNPs", file=sys.stderr)


    def add_population_freqs(self):
        """
        Compute population frequencies - can add adjsted frequencies normalised
        to mean. Need to call population_count first. 
        """
        old_settings=np.seterr(divide='ignore', invalid='ignore')
        self.freq=self.count/self.total

        # transform
        mean=np.sum(self.count, axis=1)/np.sum(self.total, axis=1)
        self.tfreq=self.freq.copy()
        for i in range(len(self.pops)):
            self.tfreq[:,i]=(self.freq[:,i]-mean)/(mean*(1-mean))
            
        np.seterr(**old_settings)
        
###########################################################################
# END CLASS

class eigenstrat_data(snp_data):
    """
    Either packed or unpacked Eigenstrat data using the pyEigenstrat module. 
    """
        
    def load(self, file_root, pops=None, inds=None, exclude_inds=None, snps=None):
        """
        Load from an unpacked eigenstrat file into our internal format.
        If pops is specified then load only the specified populations. 
        """
        data=pE.load(file_root, pops=pops, inds=inds, exclude_inds=exclude_inds, snps=snps)
        self.ind=data.ind
        self.snp=data.snp
        self.pops=np.unique(data.ind["POP"])
        self.geno=data.geno()

###########################################################################
# END CLASS

class frequency_data(snp_data):
    """
    This class looks like a snp_data class but actually just loads
    a file of frequencies. This has m allele frequencies in n populations.
    n+3 columns total, m+1 rows with a header row.
    The first three columns are snpID, chr, pos in that order
    The remaining N columns are populations
    Population names taken from the column headings
    Entries are converted to floats. 

    Um, in hindisght, maybe this is a case for duck typing. 
    """

    def load(self, file_root, pops=None, inds=None, exclude_inds=None, snps=None):
        """
        load the frequency file
        """

        if inds or exclude_inds:
            raise Exception("Neither inds nor exclude_inds do anything here")
        
        self.ind=[]
        freq_file=open(file_root, "r")
        header=freq_file.next()
        file_pops=np.array(header.split()[3:])
        which_pops=np.in1d(file_pops, np.array(pops))
        self.pops=list(file_pops[which_pops])
        self.freq=np.loadtxt(freq_file, dtype=float, usecols=5+np.where(which_pops)[0])

        freq_file.seek(0)
        freq_file.next()
        snpdt=pE.dt_snp2 
        snpcol=(0,1,2,3,4)
        self.snp=np.genfromtxt(freq_file, dtype=snpdt, usecols=snpcol)
        freq_file.close()

        if snps:
            include_snps=np.in1d(self.snp["ID"], np.array(snps))
            self.snp=self.snp[include_snps,]
            self.freq=self.freq[include_snps,]
        
    def add_population_counts(self, inbred):
        """
        Do nothing
        """
        pass

    def remove_monomorphic(self):
        """
        TODO: Implement this
        """
        pass

    def remove_sparse(self, n=None):
        """
        TODO: Implement this?
        """
        pass

    def add_population_freqs(self):
        """
        We already did this! 
        """
        pass
        
        
        
###########################################################################
# END CLASS

class read_data(snp_data):
    """
    This class has read-level data, and generates genotypes from that. 
    """

    def load(self, file_root, pops=None, inds=None, exclude_inds=None, snps=None):
        """
        This loads in data in Nick Patternson's snp .ans format.
        Load, but look for a .ans file instead of a .geno file.
        9 Columns: SNPID, chr, pos, ref, alt, "::"? SampleID ref alt reads
        """
        # Read .ind file
        ind,pops,_include=pE.load_ind_file(file_root, pops, inds, exclude_inds)    
        # Read .snp file
        snp,_snp_include=pE.load_snp_file(file_root, pops, snps)
        # Read .ans file (nickformat)
        ans=load_ans_file(file_root, snp, ind)

        alt_count=np.zeros((len(snp), len(ind)), dtype=int)
        ref_count=alt_count.copy()

        print("Bypassing allele checks", file=sys.stderr)
        print("Warning: Flipping all alleles to .ans", file=sys.stderr)

        for i_s,s in enumerate(snp):
            for i_i,i in enumerate(ind["IND"]):
                # if ans[s["ID"]][i]["ref"]!=s["REF"]:
                #     pdb.set_trace()
                #     raise Exception("Reference allele mismatch for %s"%(s["ID"]))
                # if ans[s["ID"]][i]["alt"]!=s["ALT"]:
                #     raise Exception("Alternative allele mismatch for %s"%(s["ID"]))
                s["REF"]=ans[s["ID"]][i]["ref"]
                s["ALT"]=ans[s["ID"]][i]["alt"]
                alt_count[i_s,i_i]=ans[s["ID"]][i]["alt_ct"]
                ref_count[i_s,i_i]=ans[s["ID"]][i]["ref_ct"]

        self.ind=ind
        self.snp=snp
        self.pops=pops
        self.alt_count=alt_count
        self.ref_count=ref_count 

    # TODO: add count/freq handlers to this class. 
    def add_population_freqs(self):
        self.freq=np.true_divide(1.0*self.count, self.total)

    def add_population_counts(self, inbred=[]):
        """
        Add population counts based on the reads. methods should be 
        max: maximise diversity. het if there are both alleles, hom otherwise.
        TODO: EM: Infer using EM algorithm. 
        """
        print("Adding counts with max diversity method", file=sys.stderr)
        npop=len(self.pops)
        count=np.zeros(shape=(self.alt_count.shape[0], npop), dtype=np.dtype('i4'))
        total=np.zeros(shape=(self.alt_count.shape[0], npop), dtype=np.dtype('i4'))

        for j, pop in enumerate(self.pops):
            sub_ac=self.alt_count[:,self.ind["POP"]==pop]
            sub_rc=self.ref_count[:,self.ind["POP"]==pop]
            if not sub_ac.shape[1]:
                pass
            else:
                this_count=0*sub_ac
                this_total=0*sub_ac
                for ki in range(this_count.shape[0]):
                    for kj in range(this_count.shape[1]):
                        if sub_ac[ki,kj]+sub_rc[ki,kj]==0:
                            pass
                        elif sub_ac[ki,kj]==0:
                            this_count[ki,kj]=1
                            this_total[ki,kj]=1
                        elif sub_rc[ki,kj]==0:
                            this_count[ki,kj]=0
                            this_total[ki,kj]=1
                        else:
                            this_count[ki,kj]=1
                            this_total[ki,kj]=2

                count[:,j]=np.sum(this_count, axis=1)
                total[:,j]=np.sum(this_total, axis=1)

        self.count=count
        self.total=total            
        
    def remove_monomorphic(self):
        pass
    
###########################################################################
# END CLASS

def load_ans_file(file_root, snp, ind):
    """
    Load the ans file into a dict data structure for fast lookup. Load
    only those snps and individuals included
    """
    snp_inc=set(snp["ID"])
    ind_inc=set(ind["IND"])
    ans={}
    ans_file=open(file_root+".ans", "r")
    for line in ans_file:
        snp,chrom,pos,ref,alt,dd,ind,ref_ct,alt_ct=line.split()  # @UnusedVariable
        if snp not in snp_inc:
            continue
        if snp not in ans:
            ans[snp]={}
        if ind in ind_inc:   
            ans[snp][ind]={"ref":ref, "alt":alt, "ref_ct":int(ref_ct), "alt_ct":int(alt_ct)}

    ans_file.close()
    return ans
        
###########################################################################
