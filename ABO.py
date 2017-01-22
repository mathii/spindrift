# Try and infer ABO blood group from genotypes or read counts. 

from __future__ import division, print_function
import numpy as np
import sys, getopt
import snp_data, haplotyper
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
    options ={ "data":"", "pops":[], "inbred":[], "out":"", "read":True}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:p:o:i:r",
                                   ["data", "pops",  "out", "inbred", "read"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-p","--pops"]:       options["pops"] = parse_pops(a)
        elif o in ["-i","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-r","--read"]:       options["read"] = True

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options

###########################################################################

# These alleles *have* to be on the +strand
# The last one is suspect. 
# ABO_SNPS=["rs8176720", "rs1053878", "rs8176741", "rs514659"]
# ABO_ALLELES={"A1":["T","G","G", "C"], "A1wA2":["T", "A", "G", "C"],
#              "B":["C","G","A", "C"], "O1":["T","G","G", "A"],
#              "O1w":["C","G","G", "A"], "O2":["C","G","G", "C"]}
# ABO_TYPES=["A", "B", "AB", "O"]

ABO_SNPS=["rs8176719", "rs8176746", "rs8176747" ]
ABO_ALLELES={"A":["I", "G", "C"], "B":["I", "T", "G"], "O":["D", "G", "C"]}
ABO_TYPES=["A", "B", "AB", "O"]

###########################################################################

def type_from_haps(hap1, hap2):
    """
    Return blood type given two haplotypes
    """
    haps=[hap1, hap2]
    haps.sort()
    
    if haps==["A", "A"]:
        return "A"
    if haps==["B", "B"]:
        return "B"
    if haps==["A", "O"]:
        return "A"
    if haps==["B", "O"]:
        return "B"
    if haps==["A", "B"]:
        return "AB"
    if haps==["O", "O"]:
        return "O"
    
    raise Exception("Can't type haps %s and %s"%(hap1, hap2))

    
###########################################################################

def collapse_types(results, haps):
    """
    Collapse the sub-alleles and get a type (A,B,AB or O)
    Haps should be sorted ABO_ALLELES, but just in case
    """    
    types=np.array([h[0] for h in haps])
    probs=np.zeros(len(ABO_TYPES), dtype=np.float64)
    for x in range(len(haps)):
        for y in range(len(haps)):
            abo_type=type_from_haps(types[x], types[y])
            wh=np.where(np.array(ABO_TYPES)==abo_type)[0][0]
            probs[wh]=max(probs[wh],results[x,y])

    # Hack: Fix this. 
    probs=probs/sum(probs)
    
    return probs

###########################################################################

def output(data, type_info, options):
    """
    For each sample output the blood type info
    """
    abo_lower=[x.lower() for x in ABO_TYPES]
    abo_upper=[x.upper() for x in ABO_TYPES]
    
    outfile=open(options["out"]+".abo.txt", "w")
    for indi in data.ind:
        sample=indi["IND"]
        pop=indi["POP"]
        # Pick the most likely haplotype and compute approximate posteriors. 
        best=np.argmax(type_info[sample])
        log_lik=np.log(type_info[sample]/np.max(type_info[sample]))
        best_p=type_info[sample][best]
        best_type=abo_lower[best]
        if np.max(np.delete(log_lik, best))<=-2:
            best_type=abo_upper[best]

        info=(sample, pop, best_type)+tuple(log_lik)
        outfile.write("%s\t%s\t%s\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n"%info)

    outfile.close()
        
###########################################################################

def main(options):
    """
    Run ABO haplotyper. 
    """
    if not options["read"]:
        raise Exception("Genotype-based typing not implemented")
    
    data=snp_data.read_data(options["data"], options["pops"], 
                                  True, options["inbred"], sparse=0, 
                                  snps=ABO_SNPS)

    typer=haplotyper.read_haplotyper(data, ABO_SNPS, ABO_ALLELES)
    results,haps=typer.haplotype_probs()

    type_info={}
    for k,v in results.iteritems():
        type_info[k]=collapse_types(v, haps)

    output(data, type_info, options)
        
###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
