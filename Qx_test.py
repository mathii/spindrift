from __future__ import division, print_function
import numpy as np
import numpy.lib.recfunctions as nprec
import sys, drift, copy
import pdb

class Qx_test:
    """
    Class that implements the Qx test
    Generally effects can change, but data should 
    be constant. 
    """

    N_COVARIANCE_SNPS=10000          # Number of snps to sample to estimate the covariance 
    N_FREQ_BINS=20                   # Number of bins to sample from.  
    
    def __init__(self, effects, data):
        self.effects=effects.effects.copy()   # copied
        self.data=data                        # Not copied - TODO make this a class

        self.filter_effects_against_data()
        self.augment_effects() 

        self.data_bins=self.data_frequency_bins()
        
    def copy(self):
        new=copy.copy(self)
        new.effects=self.effects.copy()
        return new
        
###########################################################################
        
    def filter_effects_against_data(self):
        """
        Take a dataset (Todo: change the dataset object to a class) and
        filter out all the snps that are not in the dataset. Also flip the
        alleles so that the EFFECT alelles is the REF allele 
        """

        # First, filter alleles: 
        data=self.data
        dt=np.dtype([("CHR", np.str_, 2), ("POS", np.int32), ("REF", np.str_, 1),
                      ("ALT", np.str_, 1)])
        data_positions=np.array(zip(data["CHR"], data["POS"], data["REF"],
                                    data["ALT"]), dtype=dt)

        new_effects=nprec.join_by(["CHR", "POS"], self.effects, data_positions,
                                       usemask=False, jointype="inner")

        print( "Removed "+str(len(self.effects)-len(new_effects))+
               " effects not in data",file=sys.stderr)

        flipped=0
        removed=0
        for rec in new_effects:
            if rec["EFFECT"]==rec["REF"] and rec["OTHER"]==rec["ALT"]:
                pass
            elif rec["OTHER"]==rec["REF"] and rec["EFFECT"]==rec["ALT"]:
                flipped+=1
                rec["OTHER"]=rec["ALT"]
                rec["EFFECT"]=rec["REF"]
                rec["BETA"]=-rec["BETA"]
            else:
                removed+=1
                rec["EFFECT"]=rec["OTHER"]="N"

        new_effects=new_effects[new_effects["EFFECT"]!="N"]
        print( "Removed "+str(removed)+" non-matching alleles",file=sys.stderr)
        print( "Flipped "+str(flipped)+" alleles",file=sys.stderr)

        self.effects=new_effects[["CHR", "POS", "EFFECT", "OTHER", "BETA"]]

###########################################################################

    def augment_effects(self):
        """
        Add the population frequency information to the effects. 
        """
        data=self.data
        npop=len(data["POPS"])
        dt=np.dtype([("CHR", np.str_, 2), ("POS", np.int32), ("FREQ", np.float, (npop,))])
        data_positions=np.array(zip(data["CHR"], data["POS"], data["FREQ"]), dtype=dt)
        
        new_effects=nprec.join_by(["CHR", "POS"], self.effects, data_positions,
                                  usemask=False, jointype="inner")

        self.effects=new_effects

###########################################################################

    def genetic_values(self):
        """
        Mean genetic values for each population
        """
        Z=2*np.sum(np.expand_dims(self.effects["BETA"],axis=1)*self.effects["FREQ"], axis=0)
        return Z

###########################################################################

    def scaling_factor(self):
        """
        Total amount of phenotypic variation
        """
        mean_freq=np.mean(self.effects["FREQ"], axis=1)
        VA=4*np.sum(self.effects["BETA"]*self.effects["BETA"]*mean_freq*(1-mean_freq))
        return VA

###########################################################################
    
    def T_matrix(self):
        """
        The matrix which mean centers and scales the allele frequencies
        """
        M=len(self.data["POPS"])
        T=np.zeros((M-1,M), dtype=float)-1/M
        np.fill_diagonal(T, (M-1)/M)
        return T

    
###########################################################################

    def Qx(self):
        """
        Calculate the test statistic
        """
        T=self.T_matrix()
        Z=self.genetic_values()
        Zp=T.dot(np.expand_dims(Z, axis=1))
        VA=self.scaling_factor()
        # F=drift.estimate_covariance(self.data["TFREQ"])
        F=drift.estimate_covariance(self.sample_snp_freq(10000))
        C=np.linalg.cholesky(F)

        X=np.linalg.inv(C).dot(Zp)/np.sqrt(VA)
        Qx=sum(X*X)[0]
        return Qx
        # print(Qx)
        
###########################################################################

    def test(self):
        """
        Actually run the test, compute p-values. 
        TODO: add option method=exact/bootstrap
        """
        pass

###########################################################################

    def effect_frequency_bin_counts(self):
        """
        Compute the proportion of effect snps in frequency bins
        """
        bins=self.N_FREQ_BINS
        breaks=np.arange(0,bins+1)/bins
        hist=np.histogram(self.effects["FREQ"][:,0], breaks)[0]
        return hist/sum(hist)

###########################################################################

    def data_frequency_bins(self):
        """
        Compute the IDs of data points in effect bins. 
        """
        bins=self.N_FREQ_BINS
        breaks=np.arange(0,bins+1)/bins
        data_freq=np.array(self.data["FREQ"][:,0])
        IDS=[None]*bins
        for i in range(bins):
            IDS[i]=np.where(np.logical_and(data_freq>breaks[i], data_freq<breaks[i+1]))[0]
        return IDS
        
###########################################################################

    def sample_snp_freq(self, K, match=True):
        """
        Sample a random set of K snps. If match=True then 
        match the frequencies to the distribution effect frequencies 
        defined by the data_frequency_bins.
        """
        data=self.data["FREQ"]
        G=np.array(data, dtype=float)

        if match:
            eff_bins=self.effect_frequency_bin_counts()
            if len(eff_bins)!=self.N_FREQ_BINS or len(self.data_bins)!=self.N_FREQ_BINS:
                raise Exception("Length of effect bins does not match length of data bins")
            
            new_ids=np.zeros(K, dtype=int)
            new_bins=np.random.choice(self.N_FREQ_BINS, K, p=eff_bins, replace=True)
            for k in range(K):
                new_ids[k]=np.random.choice(self.data_bins[new_bins[k]], 1)
            sample=G[new_ids,:]
            return(sample)
        else:
            sample=G[np.random.randint(G.shape[0],size=K),:]
            return(sample)

###########################################################################

    def bootstrap_statistic(self, nboot):
        data=self.data["FREQ"]
        G=np.array(data, dtype=float)

        for i in range(nboot):
            boot=self.copy()
            sample=(np.random.randint(G.shape[0],size=boot.effects.shape[0]))
            boot.effects["FREQ"]=G[sample,:]
            stat=boot.Qx()
            print(stat)
            
###########################################################################
# END CLASS


