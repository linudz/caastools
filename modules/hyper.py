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

MODULE NAME:    hyper.py
DESCRIPTION:    Pvalue assignment to CAAS prediction based on hypergeometric probability function.
DEPENDENCIES:   TBD
'''


from itertools import combinations
import functools
from scipy import stats as ss
import glob


# FUNCTION count_symbols() - Counts the symbols (AAs) 
def count_symbols(list1, list2):
    outdict = {}
    for x in list1:
        outdict[x] = list2.count(x)
    
    return outdict

# FUNCTION pstate() - Calculate the statistics of a symbol combination ()

def pstate(iterable, freq_dictionary, set_size, contrast_size):

    '''
    Hypergeometric variables

    M = IL MAZZO DI CARTE
    Mn = IL MAZZO DI CARTE SAPENDO CHE L'AVVERSARIO (BACKGROUND) HA COLTO LE SUE (M - CONTRAST_SIZE)
    N = IL NUMERO DI CARTE CHE VUOI NEL MAZZO = sum of the frequences of the iterable
    n = LE CARTE CHE HAI IN MANO (SET_SIZE)
    k LE CARTE CHE HAI IN MANO E CHE COINCIDONO CON (N) -----> n = k (SET_SIZE)


    '''
    
    M = functools.reduce(lambda a, b: a+b, freq_dictionary.values())
    Mn = M - contrast_size

    it_freqs = map(lambda x : freq_dictionary[x], iterable)
    N = functools.reduce(lambda a, b: a+b, it_freqs)
    n = set_size
    k = set_size

    if N > Mn:
        p = 0
    else:
        l = ss.hypergeom(Mn, n, N)
        p = l.pmf(k)

    return p

def sstate(combination, freq_dictionary, fg_size, bg_size):
    f = pstate(list(combination[0]), freq_dictionary, fg_size, bg_size)
    b = pstate(list(combination[1]), freq_dictionary, bg_size, fg_size)

    return f * b

# FUNCTION calcpval_random() - fetches a pvalue from a line, given foreground and background sizes
def calcpval_random(line_dictionary, genename, fg_size, bg_size):

    class lstats():
        def __init__(self):
            self.p = ""

            # Symbols
            self.s =[]                      # Symbols. All the amino acids in the position.
            self.sfreq = {}                 # Frequency per each symbol in th position
            self.c = []                     # One to many combinations
            self.uc = []                    # One to one combinations

            self.fgsize = fg_size
            self.bgsize = bg_size

            self.pvalue = 0

            self.out_tuple = ()

        def combine(self):
            for x in self.s:
                bg = list(set(self.s))
                bg.remove(x)
                self.c.append([[x], bg])
                self.c.append([bg, [x]])
            
            self.uc = list(combinations(self.s, 2))
            rec = map(lambda x : (x[1], x[0]), self.uc)

            self.uc = self.uc + list(rec)
        
        def dostat(self):

            combresults = list(map(functools.partial(sstate, freq_dictionary = self.sfreq, fg_size = self.fgsize, bg_size = self.bgsize), self.c))
            p_comb = functools.reduce(lambda a, b: a+b, combresults)
                        
            if len(self.s) > 2:                     # You will subtract one to one combinations only w 3+ AAs
                
                u_combresults = list(map(functools.partial(sstate, freq_dictionary = self.sfreq, fg_size = self.fgsize, bg_size = self.bgsize), self.uc))
                up_comb = functools.reduce(lambda a, b: a+b, u_combresults)
                self.pvalue = p_comb - up_comb

            else:
                self.pvalue = p_comb
            
            self.out_tuple = (genename + "@" + str(self.p), self.pvalue)

    z = lstats()

    symbols = list(map(lambda x : x.split("@")[0], line_dictionary.values()))
    z.p = list(map(lambda x : x.split("@")[1], line_dictionary.values()))[0]
    
    try:
        symbols.remove("-") # Clean the gaps...
    except:
        pass                # ...if any.

    # Load the z object

    z.s = set(symbols)
    z.sfreq = count_symbols(z.s, symbols)
    z.combine()
    z.dostat()

    return z.out_tuple

### FUNCTION genepval() - iterates the pvalue calculation over a whole gene

def genepval(imported_alignment, gpv_fg_size, gpv_bg_size, mode = "random"):
    pval_dictionary = {}

    if mode == "random":
        pvals = map(functools.partial(calcpval_random, genename = imported_alignment.genename, fg_size = gpv_fg_size, bg_size = gpv_bg_size), imported_alignment.d)
        pval_dictionary = dict(pvals)

    return pval_dictionary

