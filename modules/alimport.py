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


MODULE NAME: alimport.py
DESCRIPTION: MSA importation from various format through BioPython
INPUTS:      Input MSAs
CALLED BY:   caas_id.py, fastcaas_core.py, disco.py

TABLE OF CONTENTS
------------------------------------------
import_position()           Imports a position from a BioPython imported alignment

filter_position()           This function is designed to exclude those positions that are so conserved
                            that it is impossible (or unlikely) for them to return a CAAS.

slice()                     Filters the alignment and returns the

'''                                                       


from Bio import AlignIO
import functools


# FUNCTION import_position()
# Imports a position from a BioPython imported alignment

def import_position(position, imported_alignment):
    position_dictionary = {}

    for record in imported_alignment:
        position_dictionary[record.id] = record.seq[position] + "@" + str(position)
    
    return position_dictionary



# FUNCTION filter_position()                       #devnote TO BE EXPORTED IN FILE, GAPSRATIO INCLUDED
# This function is designed to exclude those positions that are so conserved
# that it is impossible (or unlikely) for them to return a CAAS.

def filter_position(imported_position, changes_threshold, max_gaps_ratio):

    aas = map(lambda x : x.split("@")[0], [x for x in imported_position.values()])

    seq = "".join(list(aas))

    outflag = True

    # Filter per gaps
    gaps_ratio = 0
    gaps_ratio = seq.count("-")/float(len(seq))

    if gaps_ratio > max_gaps_ratio:
        outflag = False
        return outflag
 
    # Filter per amino acid diversity (minimum changes)
    single_symbols = list(set(seq))

    try:
        single_symbols.remove("-")
    except:
        pass

    all_symbols = list(seq)

    if len(single_symbols) == 1:
        seconds = 0
    elif len(single_symbols) > 1:
        counts = []
        for x in single_symbols:
            counts.append(all_symbols.count(x))
        
        counts.remove(max(counts))        
        seconds = sum(counts)
    
    if seconds < changes_threshold:
        outflag = False

    return outflag


# FUNCTION slice()
# Generates a key file per each gene
 
def slice(alignment_file, alignment_format, column_threshold, max_gaps = 0.5):

    class slice_object():
        def __init__(self):
            self.d = []
            self.genename = ""
            self.species = []
    
    z = slice_object()
    imported_alignment = AlignIO.read(alignment_file, alignment_format)
    z.genename = alignment_file.split("/")[-1].split(".")[0]

    # SPECIES IN THE ALIGNMENT

    for x in imported_alignment:
        z.species.append(x.id)
    
    z.species = list(set(z.species))

    # IMPORTING POSITIONS

    imported_positions = list(map(functools.partial(import_position, imported_alignment = imported_alignment), [position for position in range(0,imported_alignment.get_alignment_length())]))

    # FILTERING POSITIONS
    z.d = list(filter(functools.partial(filter_position, changes_threshold = column_threshold, max_gaps_ratio = max_gaps), imported_positions))
    return z