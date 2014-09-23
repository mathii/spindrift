# Implementing the QST-like test from Berg & Coop

from __future__ import division, print_function
import sys, getopt
import effects, Qx_test, snp_data
from parse import parse_list_of_pops, parse_pops

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    gwas: Gwas data. 3-col: CHR, POS, EFFECT OTHER BETA
    pops: Comma separated list of populations to include
    inbred: Comma separated list of pops that might be inbred
    nboot: Number of bootstrap replicates. 
    """
    options ={ "data":"", "gwas":"", "pops":[[]], "inbred":[], "full":False,
               "nboot":1000, "out":None, "used":False, "match":True}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:g:p:n:o:i:umf",  # @UnusedVariable
                                   ["data", "gwas", "pops", "nboot", "out", 
                                    "inbred", "used", "nomatch", "full"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-g","--gwas"]:       options["gwas"] = a
        elif o in ["-p","--pops"]:       options["pops"] = parse_list_of_pops(a)
        elif o in ["-i","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-n","--nboot"]:      options["nboot"] = int(a)
        elif o in ["-u","--used"]:       options["used"] = True
        elif o in ["-m","--nomatch"]:    options["match"] = False
        elif o in ["-f","--full"]:       options["full"] = True

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options
   
###########################################################################

def main(options):
    # Load population data 
    
    output_file=open(options["out"]+"results.txt", "w")
    output_file.write("Pops\tQx\tP.X2\tP.boot\n")
    
    output_root=None
    full_results=len(options["pops"])==1 or options["full"]

    gwas=effects.effects(options["gwas"])
    
    for pops in options["pops"]:
        print("\nLoading: "+",".join(pops), file=sys.stderr)

        data=snp_data.eigenstrat_data(options["data"], pops, True, options["inbred"])
        test=Qx_test.Qx_test(gwas, data, match=options["match"])

        if full_results:
            output_root=options["out"]+"_".join(pops)
        
        test.test( output_file=output_file, output_root=output_root, nboot=options["nboot"])
        if options["used"]:
            test.output_effects(options["out"])
    
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
