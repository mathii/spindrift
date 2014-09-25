# object for predicting phenotypes based on either quantitative or discrete trait sets

from __future__ import division, print_function
import numpy as np
import numpy.lib.recfunctions as nprec
import sys, copy
import pdb


###########################################################################

class predictor():
    """
    Predicts individual level phenotypes
    Modifies data and effects.
    """

    def __init__(self, data, effects):
        self.effects=effects.effects
        self.data=data
        self.filter_effects()
        
    def filter_effects(self):
        """
        Merge effects and data, and flip effect alleles 
        """
        effect_positions=self.effects[["CHR", "POS"]]
        data_positions=self.data.snp[["CHR", "POS"]]

        effect_include=np.in1d(effect_positions, data_positions)
        data_include=np.in1d(data_positions, effect_positions)

        self.data.filter_snps(data_include)
        self.effects=self.effects[effect_include]
        # Just give up and convert to float. I have no idea why int doesn't work here
        # but it's something to do with the fact that you can't have None as a numpy int
        # wheras float gets converted to nan. 
        tmp_data=nprec.append_fields(self.data.snp, "GENO", None, dtypes=[(float,self.data.geno.shape[1])],usemask=False)
        tmp_data["GENO"]=self.data.geno
        self.effects=nprec.join_by(["CHR", "POS"], self.effects, tmp_data, usemask=False, jointype="inner")
        flipped=0
        removed=0
        for rec in self.effects:
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

        self.effects=self.effects[self.effects["EFFECT"]!="N"]
        print( "Removed "+str(removed)+" non-matching alleles",file=sys.stderr)
        print( "Flipped "+str(flipped)+" alleles",file=sys.stderr)

    def predict_effects(self):
        """
        Predict the trait values for each individual. 
        Replace missing data (gt=9) with the expected value from 
        other individuals. 
        """

        # Replace missing data with mean
        g_mat=copy.copy(self.effects["GENO"])
        for i in range(g_mat.shape[0]):
            g_mat[i,g_mat[i,:]==9.]

        # Compute genetic values
        values=np.sum(g_mat*np.expand_dims(self.effects["BETA"],axis=1), axis=0)
        return values
        
    def print_values(self, outfile=sys.stdout):
        """
        Just output the values
        """
        values=self.predict_effects()

        print("\nIndividual\tPopulation\tValue", file=outfile)
        for i in range(len(self.data.ind)):
            val=(self.data.ind["IND"][i], self.data.ind["POP"][i], values[i])
            print("%s\t%s\t%1.4f" % val, file=outfile)
        
            
###########################################################################

