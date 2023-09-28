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
                # Esta modificación es mía jejej

MODULE NAME: fmtparser.py
DESCRIPTION: detects input alignment format and states the fmt variable accordingly. 
DEPENDENCIES: TBD.
CALLED BY: CT.

'''

from Bio import AlignIO
import sys

def detect_format(file_path):
    formats = ["clustal", "emboss", "fasta", "fasta-m10", "ig", "maf", "mauve", "msf", "nexus", "phylip", "phylip-sequential", "phylip-relaxed", "stockholm"]
    
    for fmt in formats:
        try:
            # Try to read the file with the current format
            AlignIO.read(file_path, fmt)
            return fmt
        except:
            continue
    return None

if __name__ == "__main__":
    file_path = sys.argv[1]
    fmt = detect_format(file_path)
    print(fmt)