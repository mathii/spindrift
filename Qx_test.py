from __future__ import division, print_function
import numpy as np
import scipy as sp
import numpy.lib.recfunctions as nprec
import sys, drift, copy
import pdb
MATPLOTLIB_AVAILABLE=True
try:
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB_AVAILABLE=False
    
class Qx_test:
    """
    Class that implements the Qx test
    Generally effects can change, but data should 
    be constant. 
    """

    N_COVARIANCE_SNPS=10000      # Number of snps to sample to estimate the covariance but... 
    MAX_FRAC_RESAMPLE=10           # ...sample at most 1/this proportion of snps
    N_FREQ_BINS=20                 # Number of bins to sample from.  
    
    def __init__(self, effects, data):
        self.effects=effects.effects.copy()   # copied to make bootstrapping easier
        self.data=data                        # Not copied
        self._n_cov_snps=max(self.N_COVARIANCE_SNPS, 
                             len(data.snp)//self.MAX_FRAC_RESAMPLE)
        
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
        Take a dataset and
        filter out all the snps that are not in the dataset. Also flip the
        alleles so that the EFFECT alelles is the REF allele 
        """

        # First, filter alleles: 
        new_effects=nprec.join_by(["CHR", "POS"], self.effects, 
                                  self.data.snp[["CHR", "POS", "REF", "ALT"]],
                                  usemask=False, jointype="inner")

        print( "Removed "+str(len(self.effects)-len(new_effects))+
               " effect SNPS not in data",file=sys.stderr)

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

        tmp_snp=nprec.append_fields(self.data.snp, "FREQ", None, 
                                    dtypes=[(float,self.data.freq.shape[1])],
                                    usemask=False)
        tmp_snp["FREQ"]=self.data.freq
        
        new_effects=nprec.join_by(["CHR", "POS"], self.effects, tmp_snp,
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
        M=len(self.data.pops)
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
        F=drift.estimate_covariance(self.sample_snp_freq(self._n_cov_snps, match=False))
        C=np.linalg.cholesky(F)

        X=np.linalg.inv(C).dot(Zp)/np.sqrt(VA)
        Qx=sum(X*X)[0]
        return Qx
        
###########################################################################

    def test(self, output_root=None, nboot=None):
        """
        Actually run the test, compute p-values. 
        TODO: add option method=exact/bootstrap
        """
        print("Running Qx test", file=sys.stderr)
        Z=self.genetic_values()
        print("\nGenetic values:\n"+ 
              "\n".join([x+" : "+str(y) for x,y in zip(self.data.pops, Z)]), file=sys.stderr)

        Qx=self.Qx()
        print("\nQx = %2.3f"%(Qx), file=sys.stdout)
        
        X2df = len(Z)-1
        X2p=1 - sp.stats.chi2.cdf(Qx, X2df)
        print("X^2_%d p = %1.4f"%(X2df, X2p), file=sys.stdout)

        bootp=None
        if nboot:
            print("Bootstrapping", file=sys.stderr)
            boots=self.bootstrap_statistic(nboot)
            bootp=np.mean(boots>Qx)
            print("Bootstrap p = %1.4f"%(bootp), file=sys.stdout)

        if output_root:
            with open(output_root + ".Qx.stat.txt", "w") as f:
                f.write(str(Qx))
            with open(output_root + ".Qx.boot.txt", "w") as f:
                np.savetxt(f, boots)
            if MATPLOTLIB_AVAILABLE:
                self.plot_bootstrap_hist(boots, X2df, Qx, bootp, output_root+".Qx.boot.pdf",
                                         ",".join(self.data.pops))
                
###########################################################################

    def plot_bootstrap_hist(self, stats, df, Qx, bootp, out_file, title=""):
        """
        Plot a histogram of the bootsrap values
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for item in [fig, ax]:
            item.patch.set_visible(False)
            
        N = max(10, len(stats)/20)
        n, bins, patches = ax.hist(stats, N, normed=1, facecolor='#377EBA', 
                                   alpha=0.75, label="Bootstrap")
        xs = np.arange(bins[1], bins[-1]+1,0.1)
        ys = sp.stats.chi2.pdf( xs, df)
        ax.plot(xs, ys, 'b', linewidth=2, label="Theoretical")
        ax.set_xlabel('Qx')
        ax.set_ylabel('Density')
        ax.set_title(title)
        plt.axvline(Qx, color="r", linewidth=2, label="Observed")
        plt.text(Qx+1,  sp.stats.chi2.pdf( Qx, df)+0.01
                 , "p = %1.4f"%(np.mean(stats>Qx)), color="red")
        plt.legend(frameon=False)
        plt.savefig(out_file)
        
###########################################################################

    def effect_frequency_bin_counts(self):
        """
        Compute the proportion of effect snps in frequency bins
        """
        bins=self.N_FREQ_BINS
        breaks=np.arange(0,bins+1)/bins
        hist=np.histogram(np.mean(self.effects["FREQ"], axis=1), breaks)[0]
        return hist/sum(hist)

###########################################################################

    def data_frequency_bins(self):
        """
        Compute the IDs of data points in effect bins, averaged over all populations 
        """
        bins=self.N_FREQ_BINS
        breaks=np.arange(0,bins+1.)/bins
        data_freq=np.mean(self.data.freq, axis=1)
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
        G=self.data.freq

        if match:
            eff_bins=self.effect_frequency_bin_counts()
            if len(eff_bins)!=self.N_FREQ_BINS or len(self.data_bins)!=self.N_FREQ_BINS:
                raise Exception("Length of effect bins does not match length of data bins")
            new_ids=np.zeros(K, dtype=int)
            new_bins=np.random.choice(self.N_FREQ_BINS, K, p=eff_bins, replace=True)
            for k in range(K):
                new_ids[k]=np.random.choice(self.data_bins[new_bins[k]], 1)
                sample=G[new_ids,:]
            return sample
        else:
            sample=G[np.random.randint(G.shape[0],size=K),:]
            return sample

###########################################################################

    def bootstrap_statistic(self, nboot):
        stats=np.zeros(nboot, dtype=float)
        N_effects=self.effects.shape[0]

        for i in range(nboot):
            boot=self.copy()
            boot.effects["FREQ"]=self.sample_snp_freq(N_effects)
            stats[i]=boot.Qx()

        return stats 

###########################################################################

    def output_effects(self, file_root):
        """
        Write the used effects to a file. 
        """
        np.savetxt(file_root+".Qx.used.txt",
                   self.effects[["CHR", "POS", "EFFECT", "OTHER", "BETA"]],
                   fmt="%s\t%d\t%s\t%s\t%f")
        
###########################################################################
# END CLASS


