# Tests for admixture times and proportions. 

from __future__ import division, print_function
import numpy as np
import sys, getopt
import Admix_estimate, snp_data
import pdb

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    1,2, and 3: Corresponding to populations A,B and C
    inbred: Comma separated list of pops that might be inbred
    """
    options ={ "data":"", "1":"", "2":"", "3":"", "inbred":[],
               "out":"out"}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:1:2:3:i:o:",
                                   ["data", "1", "2", "3", "inbred", "out"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-i","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-1","--1"]:          options["1"] = a
        elif o in ["-2","--2"]:          options["2"] = a
        elif o in ["-3","--3"]:          options["3"] = a

    options["pops"]=[options[x] for x in ["1", "2", "3"]]
            
    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options
   
###########################################################################

def main(options):
    # Load population data - TODO: Move all this to an eigenstrat class
    data=snp_data.eigenstrat_data(options["data"], options["pops"], 
                                  True, options["inbred"], 0)

    test=Admix_estimate.patterson_estimator(data, np.array([options[x] for x in ["1", "2", "3"]]))
    test.estimate_admix()
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    

