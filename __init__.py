# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (C) 2014  Andreas W. Hauser

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#####################################################################
#                                                                   #
#               pyQchem - Input/Output-Tools for Q-Chem             #
#                                                                   #
#                           Version 0.4                             #
#                                                                   #
#####################################################################

# AWH, Jan 2014

############################## RULES ################################

# Object properties are only accessible via methods. Exceptions are 
# valid in cases where it really makes sense for the user and Python
# datatypes come in handy. In these rare exceptions, the property has
# to contain the word "list" if it is a list and "dict" if it is
# a dictionary. Indentation is 4 spaces.

# Current addition: The jobfile class allows direct access to certain
# fragment objects (rem, basis, molecule ...). 

######################### STANDARD MODULES ##########################

from copy import deepcopy

############################ CONSTANTS ##############################

import constants

############################# MODULES ###############################

# First we need all visible inputfile classes. 

from input_classes import *

# Then we add some hidden outputfile classes... 

from output_classes import _outputfile
from output_classes import _multioutput

# and one hidden inputfile class for users who need to create their own
# unsupported input arrays

from input_classes import _unsupported_fragment

# Now we load some subroutines for filereading.

from utilities import _readzmat
from utilities import _readcartesian
from utilities import _readtinker
from utilities import _readinput

########################### FILEHANDLING ############################
        
# This is the main filereading method. It reads all types of supported
# files. All other reading methods are hidden from the user. 

def read(filename):
    extension = (filename.split("."))[-1]

    # Do we have an inputfile?
    if  extension in ("inp","in","IN","INP","qcin","QCIN"):

        seperator = []
        with open(filename) as infile:
            for num, line in enumerate(infile):
                if "@@@"==line.strip() or "@@@@"==line.strip():
                    seperator.append(num)

        N_jobs = len(seperator)

        # Does the file contain multiple jobs? 
        if N_jobs>0:
            print "Batch Jobfile detected."
            joblist = []
            infile = open(filename,"r")
            content = infile.readlines()
            seperator.insert(0,-1)
            seperator.append(len(content))
    
            # Create batch jobfile object
            re_file = multifile()
    
            # Populate it with jobs
            for k in range(N_jobs+1):
                start=seperator[k]+1
                end=seperator[k+1]
                dummylines = content[start:end]
                dummy = _readinput(dummylines)
                re_file.add(dummy)
            return re_file
    
        # No, it's a single job file    
        else:
            print "Jobfile detected."
            return _readinput(filename)
        
    # Is it a z-matrix file?
    elif extension in ("zmat","ZMAT","Z","z"):
        print "Z-matrix file detected."
        return _readzmat(filename)
          
    # Is it a cartesian coordinates file?
    elif extension in ("xyz","XYZ"):
        print "Cartesian coordinates file detected."
        return _readcartesian(filename)
        
    # Is it a tinker coordinates file?
    elif extension in ("txyz","TXYZ"):
        print "Tinker coordinates file detected."
        return _readtinker(filename)
    

    # Do we have a Q-Chem outputfile?
    if  extension in ("out","OUT","qcout","QCOUT"):

        seperator = []
        with open(filename) as infile:
            for num, line in enumerate(infile):
                if "Welcome to Q-Chem" in line:
                    seperator.append(num)

        N_jobs = len(seperator)

        # Does the file contain multiple jobs? 
        if N_jobs>1:
            print "Batch-Outputfile detected."
            joblist = []
            infile = open(filename,"r")
            content = infile.readlines()
            #seperator.insert(0,-1)
            seperator.append(len(content))
    
            # Create batch jobfile object
            re_file = _multioutput()
    
            # Populate it with jobs
            for k in range(N_jobs):
                start=seperator[k]+1
                end=seperator[k+1]
                dummylines = content[start:end]
                dummy = _outputfile(dummylines)
                re_file.add(dummy)
            return re_file
    
        # No, it's a single job file    
        else:
            print "Outputfile detected."
            return _outputfile(filename)

    # What the heck? This is not a valid file.
    else:
        print "Error: File type not recognized."

if __name__ == "__main__":
    print "This file is supposed to be loaded as a module, not as main."
   
   