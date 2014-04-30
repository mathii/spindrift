# Implementing the QST-like test from Berg & Coop

from __future__ import division, print_function
import numpy as np
import pdb
import sys, getopt
import drift, effects, Qx_test
import spindrift_io as io
from scipy import stats

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    gwas: Gwas data. 3-col: CHR, POS, EFFECT
    pops: Comma separated list of populations to include
    nboot: Number of bootstrap replicates. 
    """
    options ={ "data":"", "gwas":"", "pops":[], "nboot":1000}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:g:p:n:", ["data", "gwas", "pops", "nboot"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-g","--gwas"]:       options["gwas"] = a
        elif o in ["-p","--pops"]:       options["pops"] = a.split(",")
        elif o in ["-n","--nboot"]:      options["nboot"] = int(a)

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options

###########################################################################

   
###########################################################################

def main(options):
    # Load population data - TODO: Move all this to an eigenstrat class
    data=io.read_eigenstrat(options["data"], options["pops"])
    drift.population_count(data)
    drift.remove_monomorphic(data)
    drift.population_freq(data)
    del(data["GT"])                       # Save memory

    # Load gwas data
    gwas=effects.effects(options["gwas"])

    test=Qx_test.Qx_test(gwas, data)

    Qx=test.Qx()
    print(str(Qx), file=sys.stderr)
    test.bootstrap_statistic(10)
    # pdb.set_trace()
    
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
