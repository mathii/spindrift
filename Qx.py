# Implementing the QST-like test from Berg & Coop

from __future__ import division, print_function
import sys, getopt
import effects, Qx_test, snp_data
from parse import parse_list_of_pops, parse_pops
import pdb

###########################################################################

def parse_options():
    """
    -d data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    -q frequency: data file in frequency format 
    -g gwas: Gwas data. 3-col: CHR, POS, EFFECT OTHER BETA
    -p pops: Comma separated list of populations to include
    [-i] inbred: Comma separated list of pops that might be inbred
    [-c] center: Populations to center allele frequency around. 
    [-n] nboot: Number of bootstrap replicates. 
    [-x] exclude: Exclude these individuals 
    """
    options ={ "data":"", "gwas":"", "pops":[[]], "inbred":[], "full":False,
               "nboot":1000, "out":None, "used":False, "match":True, "center":[],
               "exclude":[], "values":False, "frequency":""}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:q:g:p:n:o:i:c:x:umfv",  # @UnusedVariable
                                   ["data", "frequency", "gwas", "pops", "nboot", "out",
                                    "inbred", "center", "exclude", "used", "nomatch",
                                    "full", "values"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-q","--frequency"]:  options["frequency"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-g","--gwas"]:       options["gwas"] = a
        elif o in ["-c","--center"]:     options["center"] = parse_pops(a)
        elif o in ["-p","--pops"]:       options["pops"] = parse_list_of_pops(a)
        elif o in ["-i","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-x","--exclude"]:    options["exclude"] = parse_pops(a) #actually individuals not populations
        elif o in ["-n","--nboot"]:      options["nboot"] = int(a)
        elif o in ["-u","--used"]:       options["used"] = True
        elif o in ["-m","--nomatch"]:    options["match"] = False
        elif o in ["-f","--full"]:       options["full"] = True
        elif o in ["-v","--values"]:     options["values"] = True

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    if options["data"] and options["frequency"]:
        raise Exception("Specify only one of data (-d) and frequency (-q)")
    
    return options
   
###########################################################################

def main(options):
    # Load population data 
    
    output_file=open(options["out"]+".results.txt", "w")
    headings=["Pops", "Qx", "P.X2"]
    if options["nboot"]:
        headings.append("P.boot")
    if options["values"]:
        headings.append("Genetic.values")

    output_file.write("\t".join(headings)+"\n")
    
    output_root=None
    full_results=options["full"]

    gwas=effects.effects(options["gwas"])

    for pops in options["pops"]:
        print("\nLoading: "+",".join(pops), file=sys.stderr)

        if options["data"]:
            data=snp_data.eigenstrat_data(options["data"], pops, True, options["inbred"])
        elif options["frequency"]:
            data=snp_data.frequency_data(options["frequency"], pops)
        else:
            raise Exception("No data or frequency file specified")

        test=Qx_test.Qx_test(gwas, data, center=options["center"], match=options["match"])

        if full_results:
            output_root=options["out"]+"_".join(pops)
        
        test.test( output_file=output_file, output_root=output_root, nboot=options["nboot"], values=options["values"])
        if options["used"]:
            test.output_effects(options["out"])
    
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
