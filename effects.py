# Class for loading and manipulating SNP effects

from __future__ import division, print_function
import numpy as np
import numpy.lib.recfunctions as nprec
import sys
import pdb

###########################################################################

class effects:
    """
    Class for loading and manipulating SNP effects
    """
    
    def __init__(self, effects_file):
        """
        file should have 5 columns, corresponding to CHR, POS, EFFECT, 
        OTHER and BETA
        """
        num_lines = sum(1 for line in open(effects_file))
        
        dt=np.dtype([("CHR", np.str_, 2), ("POS", np.int32), ("EFFECT", np.str_, 1),
                     ("OTHER", np.str_, 1), ("BETA", np.float)])
        self.effects=np.zeros(num_lines, dtype=dt)
                     
        data=open(effects_file)
        for i, line in enumerate(data):
            chrom, pos, effect, other, beta = line.split()
            self.effects[i]=((chrom, int(pos), effect, other, float(beta)))
        data.close()

        if i==num_lines-1:
            print( "Read "+str(num_lines)+ " effect lines", file=sys.stderr)
        else:
            raise Exception("Unexpected number of lines")

###########################################################################
# END CLASS
