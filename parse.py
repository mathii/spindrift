# Helper functions for parsing input and whatnot.  

from __future__ import division, print_function
import re,pdb

###########################################################################

def parse_pops(input):
    """
    Pasrse a list of populations. Can either be a comma separated string, 
    or a file containing either whitespace or comma separated names. 
    """
    data=None
    try: 
        file=open(input, "r")
        data=file.read().splitlines()
        data=[x for x in data if x!='']   # trailing newlines leave a blank
        if len(data)==1:
            data=re.findall(r"\w+", data[0])            
    except IOError:
        data=input.split(",")

    return data
    
###########################################################################        
