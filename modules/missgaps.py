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

MODULE NAME:    missgaps
DESCRIPTION:    function to filter missing species and gaps.
DEPENDENCIES:   none
CALLED BY: caas_id.py

'''                                                   

# FUNCTION missgaps()
# filters a trait for its gaps

def missgaps(max_bg, max_fg, max_all, nfg, nbg):

    out = True

    all_g = nfg + nbg

    
    if max_all != "nofilter" and all_g > int(max_all):
        out = False

    if max_fg != "nofilter" and nfg > int(max_fg):       
            out = False

    if max_bg != "nofilter" and nbg > int(max_bg):       
            out = False


    return out

