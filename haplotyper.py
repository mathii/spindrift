# Infer hapltoypes from reads or genotypes. 
from __future__ import division, print_function
import numpy as np
from scipy import stats
from copy import deepcopy
import pdb

# Mutation/error probability
MT_PROB=0.001

class read_haplotyper():
    """
    Infer haplotypes from read level data. 
    """

    def __init__(self, data, hap_snps, hap_alleles):
        self.data=data
        self.hap_snps=hap_snps
        self.hap_alleles=hap_alleles 

    def haplotype_probs(self):
        """
        Estimate haplotype probabilities. For each individual, am hxh matrix
        where h is the number of haplotypes. 
        """

        n_ind=len(self.data.ind)
        n_hap=len(self.hap_alleles)

        haps=self.hap_alleles.keys()
        haps.sort()

        snps=np.array(self.hap_snps)
        
        results={}
        p_alleles=self.polarise_alleles_to_ref()
        reordering = np.array([np.where(s==snps)[0][0] for s in self.data.snp["ID"]])

        for i,indi in enumerate(self.data.ind["IND"]):
            this_prob=np.zeros((n_hap,n_hap),dtype=np.float64)
            for x in range(n_hap):
                for y in range(n_hap):
                    hap0=haps[x]
                    hap1=haps[y]

                    gt=[a+b for a,b in zip(p_alleles[hap0], p_alleles[hap1])]
                    p=[MT_PROB+(0.5-MT_PROB)*g for g in gt]
                    
                    ref_reads=self.data.ref_count[reordering,i]
                    alt_reads=self.data.alt_count[reordering,i]
                    
                    probs=[np.log(stats.binom.pmf(a,b,c)) for a,b,c in 
                           zip(alt_reads, alt_reads+ref_reads, p)]
                    this_prob[x,y]=np.exp(sum(probs))

            this_prob=this_prob/np.sum(this_prob)
            # if indi=="S0271":
            #     pdb.set_trace()
            results[indi]=this_prob

        return results,haps
        
    def polarise_alleles_to_ref(self):
        """
        Polarise the hap alleles to the reference data snps.
        Return 0 if ref 1 if alt. 
        """

        new_alleles=deepcopy(self.hap_alleles)
        for i, snp in enumerate(self.hap_snps):
            which_snp=np.where(self.data.snp["ID"]==snp)[0][0]
            ref=self.data.snp["REF"][which_snp]
            alt=self.data.snp["ALT"][which_snp]
            for k,v in new_alleles.iteritems():
                if v[i]==ref:
                    v[i]=0
                elif v[i]==alt:
                    v[i]=1
                else:
                    raise Exception("Allele mismatch at hap %s, snp %s (%s)"%(k,i,snp))

        return new_alleles
        
