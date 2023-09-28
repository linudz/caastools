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

check_pattern()            checks the pattern

fetch_caas()                fetches caas per each thing
'''                                                       

from modules_py.pindex import *
from modules_py.alimport import *
from modules_py.hyper import *
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

            self.d = {}
        
    z = processed_position()
    z.d = position

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


# FUNCTION check_pattern()
# checks the pattern

def iscaas(input_string):
    
    class caaspositive():
        def __init__(self):
            self.caas = True
            self.pattern = "4"
    
    z = caaspositive()

    twosides = input_string.split("/")
    fg = list(twosides[0])
    bg = list(twosides[1])

    # Is this a CAAS?
    for x in fg:
        if x in bg:
            z.caas = False
            break

    # What is the pattern?
    if len(fg) == 1 and len(bg) == 1:
        z.pattern = "1"
    
    elif len(fg) == 1:
        z.pattern = "2"
    elif len(bg) == 1:
        z.pattern = "3"
    
    if len(fg) == 0 or len(bg) == 0:
        z.pattern = "null"
    
    return z

# FUNCTION fetch_caas():
# fetches caas per each thing

def fetch_caas(genename, processed_position, list_of_traits, output_file, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, admitted_patterns = ["1","2","3"]):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    # Filter for the number of gaps and missing species

    if len(valid_traits) > 0:

        for trait in valid_traits:

            ### GAPS filtering.

            if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxgaps_all != "NO" and processed_position.trait2gaps_fg[trait] + processed_position.trait2gaps_bg[trait] > int(maxgaps_all):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)


            ### Missing filtering.

            if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxmiss_all != "NO" and processed_position.trait2miss_fg[trait] + processed_position.trait2miss_bg[trait] > int(maxmiss_all):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)


            '''
            gaps_condition = missgaps(
                max_bg = maxgaps_bg,
                max_fg = maxgaps_fg,
                max_all = maxgaps_all,
                nfg = processed_position.trait2gaps_fg[trait],
                nbg = processed_position.trait2gaps_bg[trait])
            
            missings_condition = missgaps(
                max_bg = maxmiss_bg,
                max_fg = maxmiss_fg,
                max_all = maxmiss_all,
                nfg = processed_position.trait2miss_fg[trait],
                nbg = processed_position.trait2miss_bg[trait])
       
            
            if gaps_condition == False or missings_condition == False:
                valid_traits.remove(trait)
            '''

        output_traits = []

        # Filter for pattern
        for x in valid_traits:
            aa_tag_fg_list = processed_position.trait2aas_fg[x]
            aa_tag_fg_list.sort()
            aa_tag_fg = "".join(aa_tag_fg_list)

            aa_tag_bg_list = processed_position.trait2aas_bg[x]
            aa_tag_bg_list.sort()
            aa_tag_bg = "".join(aa_tag_bg_list)

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag)

            if check.caas == True and check.pattern in admitted_patterns:
                output_traits.append(x + "@" + tag + "@pattern" + check.pattern)
        

        # Print the output

        header = "\t".join([
            "Gene",
            "Trait",
            "Position",
            "Substitution",
            "Pvalue",
            "Pattern",
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
                thepattern = trait.split("@")[2]

                fg_species_number = str(len(processed_position.trait2ungapped_fg[traitname]))
                bg_species_number = str(len(processed_position.trait2ungapped_bg[traitname]))

                fg_ungapped = processed_position.trait2ungapped_fg[traitname]
                bg_ungapped = processed_position.trait2ungapped_bg[traitname]

                fg_ungapped.sort()
                bg_ungapped.sort()

                missings = "-"

                if len(processed_position.trait2missings[traitname]) > 0:
                    missings = ",".join(processed_position.trait2missings[traitname])
                
                # Starting the pvalue determination


                pv = calcpval_random(processed_position.d, genename, int(fg_species_number), int(bg_species_number))
                pvalue_string = str(pv)


                print("CAAS found in alignment", genename, "on position", processed_position.position, "with pvalue", pvalue_string)

                # End of the pvalue determination


                #pvalue_string = pvdict[genename + "@" + processed_position.position]
                print(  "\t".join(
                    [genename,
                        traitname,
                        processed_position.position,
                        change,
                        pvalue_string,
                        thepattern,
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