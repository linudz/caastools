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

MODULE NAME: disco.py
DESCRIPTION: runs the caas discovery on one single alignment. Returns non-validated caas candidate positions.
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: CT.

'''


from modules.caas_id import *
from modules.alimport import *
from modules.pindex import *
import os
from os.path import exists

### FUNCTION discovery()
### Scans one single alignment to identify the CAAS

def discovery(input_cfg, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, admitted_scenarios, output_file):

    # Step 1: import the trait into a trait object (load_cfg from pindex.py)
    trait_object = load_cfg(input_cfg)

    # Step 2: import the alignment int a processed position object (slice from alimport.py)
    p = sliced_object

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = trait_object, species_in_alignment = p.species), p.d)

    # Step 5: Overwrite the output file

    if exists(output_file):
        os.system("rm -r " + output_file)

    # Step 5: extract the raw caas
    for position in processed_positions:
        print("in ct script:", max_bg_gaps, max_fg_gaps, max_overall_gaps)
        fetch_caas( p.genename,
                    position,
                    trait_object.alltraits,

                    maxgaps_bg= max_bg_gaps,
                    maxgaps_fg= max_fg_gaps,
                    maxgaps_all= max_overall_gaps,

                    maxmiss_bg= max_bg_miss,
                    maxmiss_fg= max_fg_miss,
                    maxmiss_all= max_overall_miss,

                    admitted_scenarios=admitted_scenarios,
                    output_file = output_file
                    )