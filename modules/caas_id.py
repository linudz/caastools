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

MODULE NAME:    caas_id
DESCRIPTION:    Identification of caas mutation from MSA.
DEPENDENCIES:   pindex, alimport


TABLE OF CONTENTS
------------------------------------------

process_position()          processes a position from an imported alignment.
                            The output will be a dictionary that points each 
                            aminoacid (gaps included) to the
                            species sharing it.

check_scenario()            checks the scenario

fetch_caas()                fetches caas per each thing
'''                                                       

from modules.pindex import *
from modules.alimport import *
from os.path import exists

# Function process_position()
# processes a position from an imported alignment. The output will be
# a dictionary that points each aminoacid (gaps included) to the
# species sharing it.

def process_position(position, multiconfig, species_in_alignment):

    class processed_position():
        def __init__(self):
            self.position = ""
            self.aas2species = {}
            self.aas2traits = {}
            self.trait2aas_fg = {}
            self.trait2aas_bg = {}

            self.trait2ungapped_fg = {}
            self.trait2ungapped_bg = {}

            self.trait2gaps_fg = {}
            self.trait2gaps_bg = {}
            self.trait2gaps_all = {}

            self.trait2miss_fg = {}
            self.trait2miss_bg = {}
            self.trait2miss_all = {}

            self.trait2missings = {}

            self.gapped = []
            self.missing = []
        
    z = processed_position()

    # Load missing species

    confirmed_species = set(multiconfig.s2t.keys()).intersection(set(species_in_alignment))
    z.missing = list(set(multiconfig.s2t.keys()) -  confirmed_species)

    # Load aas2species
    for x in position.keys():
        z.position = position[x].split("@")[1]
        try:
            z.aas2species[position[x].split("@")[0]].append(x)
        except:
            z.aas2species[position[x].split("@")[0]] = [x]

    # Load aas2traits    
    for key in z.aas2species.keys():
        traits = []
        for species in z.aas2species[key]:
            try:
                for v in multiconfig.s2t[species]:
                    if v not in traits:
                        traits.append(v)
            except:
                pass
        
        z.aas2traits[key] = traits

        for t in traits:
            if t[-2:] == "_1" and key != "-":
                try:
                    z.trait2aas_fg[t[:-2]].append(key)
                except:
                    z.trait2aas_fg[t[:-2]] = [key]
                    pass
            if t[-2:] == "_0" and key != "-":
                try:
                    z.trait2aas_bg[t[:-2]].append(key)
                except:
                    z.trait2aas_bg[t[:-2]] = [key]
                    pass

    try:
        z.gapped = z.aas2species["-"]
    except:
        pass

    # Determine Ungapped Species

    for trait in z.trait2aas_bg.keys():
        
        # Present species (ungapped or missing)

        pfg = list(set(multiconfig.trait2fg[trait]) - set(z.gapped + z.missing))
        z.trait2ungapped_fg[trait] = pfg

        pbg = list(set(multiconfig.trait2bg[trait]) - set(z.gapped + z.missing))
        z.trait2ungapped_bg[trait] = pbg

        # Missing in alignment

        miss_fg = list(set(multiconfig.trait2fg[trait]).intersection(set(z.missing)))
        miss_bg = list(set(multiconfig.trait2bg[trait]).intersection(set(z.missing)))

        # Number of gaps
        gfg = len(set(multiconfig.trait2fg[trait]).intersection(set(z.gapped)))
        gbg = len(set(multiconfig.trait2bg[trait]).intersection(set(z.gapped)))

        # Number of missings
        mfg = len(miss_fg)
        mbg = len(miss_bg)

        gall = gfg + gbg
        mall = mfg + mbg

        z.trait2gaps_fg[trait] = gfg
        z.trait2gaps_bg[trait] = gbg
        z.trait2gaps_all[trait] = gall

        z.trait2miss_fg[trait] = mfg
        z.trait2miss_bg[trait] = mbg
        z.trait2miss_all[trait] = mall

        try:
            z.trait2missings[trait] = (miss_fg + miss_bg)
        except:
            z.trait2missings[trait] = "none"

    return z


# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg

    print("in filter_for_gaps:", max_bg, max_fg, max_all)
    print(all_g)

    
    
    if max_all != "nofilter" and all_g > int(max_all):
        out = False

    elif max_fg != "nofilter" and gfg > int(max_fg):
        out = False

    elif max_bg != "nofilter" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for its gaps

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "nofilter" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "nofilter" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "nofilter" and mbg > int(max_m_bg):
        out = False

    return out



# FUNCTION check_scenario()
# checks the scenario

def iscaas(input_string):
    
    class caaspositive():
        def __init__(self):
            self.caas = True
            self.scenario = "4"
    
    z = caaspositive()

    twosides = input_string.split("/")
    fg = list(twosides[0])
    bg = list(twosides[1])

    # Is this a CAAS?
    for x in fg:
        if x in bg:
            z.caas = False
            break

    # What is the scenario?
    if len(fg) == 1 and len(bg) == 1:
        z.scenario = "1"
    
    elif len(fg) == 1:
        z.scenario = "2"
    elif len(bg) == 1:
        z.scenario = "3"
    
    if len(fg) == 0 or len(bg) == 0:
        z.scenario = "null"
    
    return z

# FUNCTION fetch_caas():
# fetches caas per each thing

def fetch_caas(genename, processed_position, list_of_traits, output_file, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, admitted_scenarios = ["1","2","3"]):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    # Filter for the number of gaps

    for trait in valid_traits:
        print("in fetch_caas:", maxgaps_bg, maxgaps_fg, maxgaps_all)

        
        if filter_for_gaps(max_bg = maxgaps_bg, max_fg = maxgaps_fg, max_all = maxgaps_all, gfg = processed_position.trait2gaps_fg[trait], gbg = processed_position.trait2gaps_bg[trait]) == False:
            valid_traits.remove(trait)

        
        elif filter_for_missings(max_m_bg = maxmiss_bg, max_m_fg = maxmiss_fg, max_m_all = maxmiss_all, mfg = processed_position.trait2miss_fg[trait], mbg = processed_position.trait2miss_bg[trait]) == False:
            valid_traits.remove(trait)


    output_traits = []

    # Filter for scenario
    for x in valid_traits:
        aa_tag_fg_list = processed_position.trait2aas_fg[x]
        aa_tag_fg_list.sort()
        aa_tag_fg = "".join(aa_tag_fg_list)

        aa_tag_bg_list = processed_position.trait2aas_bg[x]
        aa_tag_bg_list.sort()
        aa_tag_bg = "".join(aa_tag_bg_list)

        tag = "/".join([aa_tag_fg, aa_tag_bg])

        check = iscaas(tag)

        if check.caas == True and check.scenario in admitted_scenarios:
            output_traits.append(x + "@" + tag + "@scenario" + check.scenario)
    

    # Print the output

    header = "\t".join([
        "Gene",
        "Trait",
        "Position",
        "Substitution",
        "Pvalue",
        "Scenario",
        "FFGN",
        "FBGN",
        "GFG",
        "GBG",
        "MFG",
        "MBG",
        "FFG",
        "FBG",
        "MS"

    ])

    if len(output_traits) > 0:

        if exists(output_file):
            out = open(output_file, "a")
        else:
            out = open(output_file, "w")
            print(header, file=out)

        for trait in output_traits:

            traitname = trait.split("@")[0]            
            change = trait.split("@")[1]
            thescenario = trait.split("@")[2]

            fg_species_number = str(len(processed_position.trait2ungapped_fg[traitname]))
            bg_species_number = str(len(processed_position.trait2ungapped_bg[traitname]))

            fg_ungapped = processed_position.trait2ungapped_fg[traitname]
            bg_ungapped = processed_position.trait2ungapped_bg[traitname]

            fg_ungapped.sort()
            bg_ungapped.sort()

            missings = "-"

            if len(processed_position.trait2missings[traitname]) > 0:
                missings = ",".join(processed_position.trait2missings[traitname])


            #pvalue_string = pvdict[genename + "@" + processed_position.position]
            print(  "\t".join(
                [genename,
                    traitname,
                    processed_position.position,
                    change,
                    "pvalue_missing",
                    thescenario,
                    fg_species_number,
                    bg_species_number,
                    str(processed_position.trait2gaps_fg[traitname]),
                    str(processed_position.trait2gaps_bg[traitname]),
                    str(processed_position.trait2miss_fg[traitname]),
                    str(processed_position.trait2miss_bg[traitname]),
                    ",".join(fg_ungapped),
                    ",".join(bg_ungapped),
                    missings]
            ), file = out)

        out.close()
