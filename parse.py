# Helper functions for parsing input and whatnot.  

from __future__ import division, print_function
import re

###########################################################################

def parse_pops(instr):
    """
    Pasrse a list of populations. Can either be a comma separated string, 
    or a file containing either whitespace or comma separated names. 
    """
    data=None
    try: 
        datafile=open(instr, "r")
        data=datafile.read().splitlines()
        data=[x for x in data if x!='']   # trailing newlines leave a blank
        if len(data)==1:
            data=re.findall(r"\w+", data[0])            
    except IOError:
        data=instr.split(",")

    return data
    
###########################################################################

def parse_list_of_pops(instr):
    '''
    parse a list of pops - return a list of lists, just one list if it's a single 
    set of pops or a mutltiple ones if we pass a file 
    '''
    data=None
    try:
        datafile=open(instr, "r")
        data=datafile.read().splitlines()
        data=[x for x in data if x!='']   # trailing newlines leave a blank
        data=[x.split(",") for x in data]
    except IOError:
        data=[instr.split(",")]
        
    return data

###########################################################################
