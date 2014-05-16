# Implementing the QST-like test from Berg & Coop

from __future__ import division, print_function
import numpy as np
import sys, getopt
import drift, effects, Qx_test, snp_data
from scipy import stats
import pdb

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    gwas: Gwas data. 3-col: CHR, POS, EFFECT
    pops: Comma separated list of populations to include
    inbred: Comma separated list of pops that might be inbred
    nboot: Number of bootstrap replicates. 
    """
    options ={ "data":"", "gwas":"", "pops":[], "inbred":[], "nboot":1000, "out":None}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:g:p:n:o:i:",
                                   ["data", "gwas", "pops", "nboot", "out", "inbred"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-g","--gwas"]:       options["gwas"] = a
        elif o in ["-p","--pops"]:       options["pops"] = a.split(",")
        elif o in ["-i","--inbred"]:     options["inbred"] = a.split(",")
        elif o in ["-n","--nboot"]:      options["nboot"] = int(a)

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options
   
###########################################################################

def main(options):
    # Load population data - TODO: Move all this to an eigenstrat class
    data=snp_data.eigenstrat_data(options["data"], options["pops"], True, options["inbred"])

    # Load gwas data
    gwas=effects.effects(options["gwas"])

    test=Qx_test.Qx_test(gwas, data)
    Qx=test.test( output_root=options["out"], nboot=options["nboot"])
   
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
