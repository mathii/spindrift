#Extract and print allele frequencies for a set of populations. 

from __future__ import division, print_function
import numpy as np
import sys, getopt
import snp_data
from scipy import stats
from parse import parse_pops
import pdb

###########################################################################

def parse_options():
    """
    data: Root of genotype data in eigenstrat format, i.e. root{.geno .snp .ind}
    pops: Comma separated list of populations to include
    inbred: Comma separated list of pops to treat as inbred (i.e. pick random allele)
    snps: List of snps (in .snp format)
    out: root of output file. 
    """
    options ={ "data":"", "pops":None, "inds":None, "inbred":[], "snps":"", "out":"", "genotypes":False}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "d:p:i:o:n:s:g",
                                   ["data", "pops", "inds", "out", "inbred", "snps", "genotypes"])
    except Exception as err:
        print(str(err))
        sys.exit()

    for o, a in opts:
        if o in ["-d","--data"]:         options["data"] = a
        elif o in ["-o","--out"]:        options["out"] = a
        elif o in ["-s","--snps"]:       options["snps"] = a
        elif o in ["-p","--pops"]:       options["pops"] = parse_pops(a)
        elif o in ["-n","--inbred"]:     options["inbred"] = parse_pops(a)
        elif o in ["-i","--inds"]:       options["inds"] = parse_pops(a)
        elif o in ["-g","--genotypes"]:  options["genotypes"] = True

    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr)

    return options

###########################################################################

def output_frequencies(data, out, genotypes=False):
    """
    Just write out the frequency and count information for each snp.
    """
    freq_file=open(out+".freq", "w")
    count_file=open(out+".count", "w")
    total_file=open(out+".total", "w")
    gt_file=dgtstr=sgtstr=None
    if genotypes:
        gt_file=open(out+".gt", "w")
        dgtstr="\t".join(["%d"]*len(data.ind))
        sgtstr="\t".join(["%s"]*len(data.ind))

    fstr="\t".join(["%1.4f"]*len(data.pops))
    dstr="\t".join(["%d"]*len(data.pops))
    sstr="\t".join(["%s"]*len(data.pops))
    
    for this_file in [freq_file, count_file, total_file]:
        this_file.write("%s\t%s\t%s\t%s\t%s\t"%("ID", "CHR", "POS", "REF", "ALT"))
        this_file.write(sstr%tuple(data.pops))
        this_file.write("\n")

    if genotypes:
        gt_file.write("%s\t%s\t%s\t%s\t%s\t"%("ID", "CHR", "POS", "REF", "ALT"))
        gt_file.write(sgtstr%tuple(data.ind["IND"]))
        gt_file.write("\n")
        
    for i, snp in enumerate(data.snp):
        for this_file in [freq_file, count_file, total_file, gt_file]:
            if this_file:
                this_file.write("%s\t%s\t%d\t%s\t%s\t"%tuple(snp))
        freq_file.write(fstr%tuple(data.freq[i,]))
        count_file.write(dstr%tuple(data.count[i,]))
        total_file.write(dstr%tuple(data.total[i,]))
        if genotypes:
            gt_file.write(dgtstr%tuple(data.geno[i,:]))

        for this_file in [freq_file, count_file, total_file, gt_file]:
            if this_file:
                this_file.write("\n")
    
    for this_file in [freq_file, count_file, total_file, gt_file]:
        if this_file:
            this_file.close()
            
###########################################################################

def main(options):
    # Load population data - TODO: Move all this to an eigenstrat class
    snps,snp_include=snp_data.load_snp_file(options["snps"])
    data=snp_data.eigenstrat_data(options["data"], options["pops"], 
                                  not options["genotypes"], options["inbred"],
                                  sparse=0, snps=snps["ID"], inds=options["inds"])
    output_frequencies(data, options["out"], options["genotypes"])
    return

###########################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

        
        
    
