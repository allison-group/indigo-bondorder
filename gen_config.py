# Do things to generate a local config file
"""from pathlib import Path

# Directory structure setup
SOURCE_DIR = Path('/Users/ares/GitHub/indigo-bondorder/indigo')
WORK_DIR = Path('/Users/ares/tmp/')

# General BO things
INFINITY = 9e9
RUN_QBND = True
BASIS_LEVEL = 'def2svpd'  
COUNTERPOISE_CORRECTED = False
HYPERVALENT = True
ELECTRON_PAIRS = True
HYPERPENALTY = True
PREFILL_LOCATIONS = False
SUPPORTED_ELEMENTS = {'H','C','N','O','F','P','S','Cl','Br'}
TIMEOUT = 2.5
NUM_PROCESSES = 4

# GA things
POP_SIZE = 50
MUTATE_PROB = 0.1
MIN_GENERATIONS = 25
MAX_GENERATIONS = 100
BRUTEFORCE_CUTOFF = POP_SIZE * MIN_GENERATIONS 
CONVERGENCE = 20
ELITEISM_SIZE = 0.25
BREEDING_ELITEISM = 0.25
SEED_COUNT = 25

# FPT things
JAVA_PATH = Path('/usr/bin/java')
LIBTW_PATH = SOURCE_DIR / 'external'
TD_TIMEOUT = 60
ALLOW_FALLBACK = False
MAX_TREEWIDTH = 5

"""
