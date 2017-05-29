import os
import sys
import time
from setuptools import setup, find_packages
from setuptools.command.install import install as _install

if sys.version_info[:2] < (3,4):
    raise RuntimeError("Python version >= 3.4 required.")

version = "0.1"

class install(_install):
    def run(self):
        self.execute(pre_setup, ('src/indigox',),
                     msg="Generating config.py file")
        _install.do_egg_install(self)

def cpu_count():
    import multiprocessing
    return multiprocessing.cpu_count()

def find_java():
    # basically same as unix which command
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
    program = 'java'
    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    print('WARNING: unable to locate java executable.', file=sys.stderr)
    return ''

def pre_setup(in_dir):
    config_data = {'num_processes':cpu_count(),
                   'javapath':find_java(),
                   'installpath':in_dir,
                   'workpath':os.path.expanduser("~")+'/tmp',
                   'date':time.strftime('%c'),
                   'version':version,
               }

    config_template = """# config.py autogenerated by setup.py at {date}

from pathlib import Path
import os

# Directory structure setup
SOURCE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
WORK_DIR = Path('{workpath}')

# General BO things
INFINITY = 9e9
RUN_QBND = True
BASIS_LEVEL = 'def2svpd'  
COUNTERPOISE_CORRECTED = False
ALLOW_HYPERVALENT = True
ELECTRON_PAIRS = True
HYPERPENALTY = True
PREFILL_LOCATIONS = False
SUPPORTED_ELEMENTS = {{'H','C','N','O','F','P','S','Cl','Br'}}
TIMEOUT = 2.5
NUM_PROCESSES = {num_processes}

# A* things
HEURISTIC = 'tight'
INITIAL_LO_ENERGY = False

# BALL things
try:
    import BALLCore
except ImportError:
    BALL_AVAILABLE = False
else:
    BALL_AVAILABLE = True
MAX_SOLUTIONS = 50
BALL_DATA_FILE = SOURCE_DIR/'external/OriginalBO.xml'

# FPT things
JAVA_PATH = Path('{javapath}')
LIBTW_PATH = SOURCE_DIR / 'external'
TD_TIMEOUT = 60
ALLOW_FALLBACK = False
MAX_TREEWIDTH = 5

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

# LO things
INIT_WITH_GA = False

"""

    # write config file
    with open(in_dir+'/config.py','w') as f:
        f.write(config_template.format(**config_data))


# do the setup stuff
packages = find_packages('src')

package_data = {
    'indigox':['external/*.xml','external/*.jar'],
    }


if __name__ == "__main__":
    setup(  cmdclass = {'install':install},
            name="indigox",
            version=version,
            packages=packages,
            install_requires=['scipy>=0.18',
                              'networkx>=1.10',
                              'bitarray>=0.8',
                              'openbabel>=2.3',
                              'numpy>=1.10',],
            package_data = package_data,
            package_dir={'':'src'},
            zip_safe = False,)
