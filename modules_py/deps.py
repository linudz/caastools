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

MODULE NAME:    dep.py
DESCRIPTION:    checks the dependencies for a specific tools (packages and inputs)
DEPENDENCIES:   none
'''

import pkg_resources

def check_dependencies(tool, list_of_deps):
    
    missing_dependencies = []

    for package in list_of_deps:
        try:
            dist = pkg_resources.get_distribution(package)
        except pkg_resources.DistributionNotFound:
            missing_dependencies.append(package)
    
    if len(missing_dependencies) > 0:
        print("\n\n****ERROR: ct", tool, "requires the following python packages:")
        print(", ".join(missing_dependencies))
        print("\nPlease, install the packages.")
        print("")
        print("")
        exit()

def check_inputs(tool, list_of_inputs):
    return 0