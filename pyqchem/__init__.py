# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (c) 2014, Andreas W. Hauser
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

#####################################################################
#                                                                   #
#               pyQchem - Input/Output-Tools for Q-Chem             #
#                                                                   #
#                           Version 0.8                             #
#                                                                   #
#####################################################################

# AWH, Jan 2014

############################## RULES ################################

# Object properties are only accessible via methods. Exceptions are 
# either info-objects (mostly in outputfile parsing) or cases where it
# really makes sense for the user and Python datatypes come in handy.
# In these rare exceptions, the property has to contain the word "list"
# if it is a list and "dict" if it is a dictionary.

# Indentation is 4 spaces. Hide subroutines and classes which are not
# needed by the user. This keeps the tab completion tidy and efficient.

# Note that "remove" functions start counting at 1, in contrast to
# the Python tradition. 

# Current addition: The jobfile class allows direct access to certain
# array objects (rem, basis, molecule ...). 

######################### STANDARD MODULES ##########################

from copy import deepcopy

############################ CONSTANTS ##############################

from . import constants

############################# MODULES ###############################

# First we need all visible inputfile classes. 

from .input_classes import *

# Then we add some hidden outputfile classes... 

from .output_classes import _outputfile
from .output_classes import _multioutput

# and one hidden inputfile class for users who need to create their own
# unsupported input arrays

from .input_classes import _unsupported_array

# Now we load some subroutines for filereading.

from .utilities import _readzmat
from .utilities import _readcartesian
from .utilities import _readtinker
from .utilities import _readpybel
from .utilities import _readinput

# Lastly, let's load some scripts for running files

from .running_scripts import *


########################### FILEHANDLING ############################

# This is the main filereading method. It reads all types of supported
# files. All other reading methods are hidden from the user. 

def read(filename, silent=False):
    """
    This method reads Q-Chem input files, output files, and coordinate files (xyz,zmat,txyz). 
    """

    extension = (filename.split("."))[-1]

    # catchall pybel list
    pybel_inlist = ['txyz', 'text', 'alc', 'castep', 'nwo', 'cdx', 'xml',
                    'pwscf', 'rsmi', 'xtc', 'g09', 'pcm', 'mopin', 'mopcrt',
                    'xyz', 'fchk', 'g03', 'cube', 'axsf', 'mpc', 'mpo', 'mop',
                    'pos', 'dat', 'moo', 'dx', 'mol', 'inchi', 'hin', 'cml',
                    'outmol', 'xsf', 'qcout', 'output', 'mdl', 'unixyz',
                    'pdbqt', 'gzmat', 'arc', 'out', 'c09out', 'feat', 'crk3d',
                    'got', 'mopout', 'tdd', 'mmod', 'bs', 'mmd', 'box', 'bgf',
                    'vmol', 'acr', 'pqs', 'crk2d', 'CONFIG', 'pdb', 'ck',
                    'c3d2', 't41', 'c3d1', 'CONTCAR', 'gamout', 'mmcif', 'txt',
                    'ct', 'therm', 'log', 'pc', 'dmol', 'molden', 'ml2',
                    'fract', 'msi', 'cdxml', 'g98', 'prep', 'gpr', 'cub', 'gam',
                    'gukin', 'cmlr', 'abinit', 'POSCAR', 'ins', 'tmol', 'png',
                    'cif', 'gamess', 'car', 'mcif', 'smi', 'can', 'caccrt',
                    'fhiaims', 'inp', 'gukout', 'sy2', 'fasta', 'mpqc', 'mold',
                    'molf', 'jout', 'yob', 'mcdl', 'ent', 'adfout', 'gro',
                    'smiles', 'fs', 'mol2', 'fa', 'pqr', 'g94', 'g92', 'fch',
                    'VASP', 'fck', 'HISTORY', 'fsa', 'gamin', 'rxn', 'mrv',
                    'sdf', 'gal', 'res', 'sd', 'ccc', 'acesout']

    # Do we have an inputfile?
    if extension in ("inp", "in", "IN", "INP", "qcin", "QCIN"):

        seperator = []
        with open(filename) as infile:
            for num, line in enumerate(infile):
                if "@@@" == line.strip() or "@@@@" == line.strip():
                    seperator.append(num)

        N_jobs = len(seperator)

        # Does the file contain multiple jobs? 
        if N_jobs > 0:
            if not silent:
                print("Batch Jobfile detected.")
            joblist = []
            infile = open(filename, "r")
            content = infile.readlines()
            seperator.insert(0, -1)
            seperator.append(len(content))

            # Create batch jobfile object
            re_file = multifile()

            # Populate it with jobs
            for k in range(N_jobs + 1):
                start = seperator[k] + 1
                end = seperator[k + 1]
                dummylines = content[start:end]
                dummy = _readinput(dummylines, silent)
                re_file.add(dummy)
            return re_file

        # No, it's a single job file    
        else:
            if not silent:
                print("Jobfile detected.")
            return _readinput(filename, silent)

    # Is it a z-matrix file?
    elif extension in ("zmat", "ZMAT", "Z", "z"):
        if not silent:
            print("Z-matrix file detected.")
        return _readzmat(filename)

    # Is it a cartesian coordinates file?
    elif extension in ("xyz", "XYZ"):
        if not silent:
            print("Cartesian coordinates file detected.")
        return _readcartesian(filename)

    # Is it a tinker coordinates file?
    elif extension in ("txyz", "TXYZ"):
        if not silent:
            print("Tinker coordinates file detected.")
        return _readtinker(filename)

    # Do we have a Q-Chem outputfile?
    if extension.lower() in ("out", "qcout", "qchem"):

        seperator = []
        with open(filename) as infile:
            for num, line in enumerate(infile):
                if "Welcome to Q-Chem" in line:
                    seperator.append(num)

        N_jobs = len(seperator)

        # Does the file contain multiple jobs? 
        if N_jobs > 1:
            if not silent:
                print("Batch-Outputfile detected.")
            joblist = []
            infile = open(filename, "r")
            content = infile.readlines()
            # seperator.insert(0,-1)
            seperator.append(len(content))

            # Create batch jobfile object
            re_file = _multioutput()

            # Populate it with jobs
            for k in range(N_jobs):
                start = seperator[k] + 1
                end = seperator[k + 1]
                dummylines = content[start:end]
                dummy = _outputfile(dummylines, silent)
                re_file.add(dummy)
            return re_file

        # No, it's a single job file    
        else:
            if not silent:
                print("Outputfile detected.")
            return _outputfile(filename, silent)

    # Is it anything pybel can read?
    elif extension in pybel_inlist:
        try:
            import pybel
        except ImportError:
            print("Error: File type not recognized.")
            return
        if not silent:
            print(pybel.informats[extension], " detected")
        return _readpybel(extension, filename)

        # What the heck? This is not a valid file.
    else:
        if not silent:
            print("Error: File type not recognized.")


if __name__ == "__main__":
    print("This file is supposed to be loaded as a module, not as main.")
