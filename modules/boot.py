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


MODULE NAME: boot.py
DESCRIPTION: bootstrap function
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.init_bootstrap import *
from modules.disco import process_position
from modules.caas_id import iscaas
from modules.alimport import *

from os.path import exists
import functools


# FUNCTION fetch_caas():
# fetches caas per each thing

# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg
    
    if max_all != "NO" and all_g > int(max_all):
        out = False

    elif max_fg != "NO" and gfg > int(max_fg):
        out = False

    elif max_bg != "NO" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for its gaps

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "NO" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "NO" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "NO" and mbg > int(max_m_bg):
        out = False

    return out



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, cycles, admitted_patterns):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    def filter_trait(trait, the_processed_position, max_bg_gaps, max_fg_gaps, max_all_gaps, max_bg_miss, max_fg_miss, max_all_miss):

        the_gfg = the_processed_position.trait2gaps_fg[trait]
        the_gbg = processed_position.trait2gaps_bg[trait]

        the_mfg = the_processed_position.trait2miss_fg[trait]
        the_mbg = processed_position.trait2miss_bg[trait]
    

        ### GAPS filtering.

        if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
            return False
        if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_fg):
            return False
        if maxgaps_all != "NO" and processed_position.trait2gaps_fg[trait] + processed_position.trait2gaps_bg[trait] > int(maxgaps_all):
           return False

        ### Missings filtering.
        
        if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
            return False
        if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_fg):
            return False
        if maxmiss_all != "NO" and processed_position.trait2miss_fg[trait] + processed_position.trait2miss_bg[trait] > int(maxmiss_all):
            return False

        else:
            aa_tag_fg_list = the_processed_position.trait2aas_fg[trait]
            aa_tag_fg_list.sort()
            aa_tag_fg = "".join(aa_tag_fg_list)

            aa_tag_bg_list = the_processed_position.trait2aas_bg[trait]
            aa_tag_bg_list.sort()
            aa_tag_bg = "".join(aa_tag_bg_list)

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag)

            if check.caas == True and check.pattern in admitted_patterns:
                return True
    
    output_traits = filter(functools.partial(
        filter_trait,
        the_processed_position = processed_position,
        max_bg_gaps = maxgaps_bg,
        max_fg_gaps = maxgaps_fg,
        max_all_gaps = maxgaps_all,

        max_bg_miss = maxmiss_bg,
        max_fg_miss = maxmiss_fg,
        max_all_miss = maxmiss_all,
    ),
                    valid_traits)



    output_traits = list(output_traits)

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))

    traitline = ",".join(output_traits)
    empval = str(int(count)/cycles)

    outline = "\t".join([position_name, count, str(cycles), empval, traitline])

    return outline
        

# FUNCTION disco_bootstrap()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file):


    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), sliced_object.d)
    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "resamplings")

    # Step 4: extract the raw caas

    output_lines = map(
        functools.partial(
            caasboot,
            list_of_traits = resampled_traits.alltraits,
            genename = the_genename,
            maxgaps_fg = max_fg_gaps,
            maxgaps_bg = max_bg_gaps,
            maxgaps_all = max_overall_gaps,

            maxmiss_fg = max_fg_miss,
            maxmiss_bg = max_bg_miss,
            maxmiss_all = max_overall_miss,

            admitted_patterns = the_admitted_patterns,
            cycles = resampled_traits.cycles) ,processed_positions
    )

    output_lines = list(output_lines)

    ooout = open(output_file, "w")

    for line in output_lines:
        print(line + "\t" + trait_config_file, file=ooout)
    
    ooout.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d