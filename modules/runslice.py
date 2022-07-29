#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/


'''
A Convergent Amino Acid Substitution identification 
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).

MODULE NAME:    runslice.py
DESCRIPTION:    The slicer function.
DEPENDENCIES:   TBD
'''
from modules.alimport import *

### Function runslice (collects the )
def runslice(options_object):

    # Inputs (transferring parsed options_object to variables)
    the_alignment = options_object.single_alignment
    alignment_format = options_object.ali_format

    # Alignment slice: 1- Calculate column treshold

    with open(options_object.config_file) as cfg_handle:
        cfg_list = cfg_handle.read().splitlines()
    
    values = []

    for x in cfg_list:
        try:
            c = x.split("\t")
            values.append(c[1])
        except:
            pass
    
    fg_species = values.count("1")
    bg_species = values.count("0")


    # Alignment slicing: sum the null values (allowed_gaps + allowed_missing_species)

    sum_nulls_fg = 0

    for x in (options_object.max_fg_gaps_string,options_object.max_fg_miss_string):
        try:
            sum_nulls_fg = sum_nulls_fg + int(x)
        except:
            pass

    sum_nulls_bg = 0

    for x in (options_object.max_bg_gaps_string,options_object.max_bg_miss_string):
        try:
            sum_nulls_bg = sum_nulls_bg + int(x)
        except:
            pass



    fg_threshold = fg_species - sum_nulls_fg
    bg_threshold = bg_species - sum_nulls_bg

    c_threshold = min(fg_threshold, bg_threshold)


    # Alignment slice: 2- Filter positions (slice alignment)

    out = slice(the_alignment, alignment_format, c_threshold, float(options_object.max_gaps_pos_string))

    return out
