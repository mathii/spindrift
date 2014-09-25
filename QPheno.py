# Infer quantitative trait values for each individual using genotypes or read counts. 

from __future__ import division, print_function
import numpy as np
import sys, getopt
import snp_data, effects, predictor
from scipy import stats
from parse import parse_pops
import pdb

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.ans .snp .ind}
    pops: Comma separated list of populations to include
    inbred: Comma separated list of pops to treat as inbred (i.e. pick random allele)
    out: root of output file. 
    read: Use read level data. 
    """
    options ={ "data":"", "pops":[], "inbred":[], "out":"", "read":False, "effects":""}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:p:o:i:e:r",
                                   ["data", "pops",  "out", "inbred", "effects", "read"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-e","--effects"]:    options["effects"] = a
        elif o in ["-p","--pops"]:       options["pops"] = parse_pops(a)
        elif o in ["-i","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-r","--read"]:       options["read"] = True

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options

###########################################################################

###########################################################################

def main(options):
    """
    Run Phenotype predictor. 
    """

    eff=effects.effects(options["effects"])
    
    if options["read"]:    
        raise Exception("Not implemented yet - todo")
    else: 
        data=snp_data.eigenstrat_data(options["data"], options["pops"], 
                                False, options["inbred"], sparse=0)
    
    pred=predictor.predictor(data, eff)
    pred.print_values()
    
###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)


