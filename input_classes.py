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
#                      pyQchem - Input Classes                      #
#                                                                   #
#####################################################################

# This file contains all input file classes and their methods.

import numpy as _np
from . import constants
from . import running_scripts


############################# RUNDATA  ###############################

class _rundata(object):
    def __init__(self):
        self.name = ''
        self.loc53 = ''
        self.qchem = ''
        self.nt = 1
        self.np = 1
        self.timestamp = False

    def __str__(self):
        ret_str = "Submission status summary:\n" + 26 * "-" + "\n\n"
        if self.name != '':
            ret_str += "Filename is " + self.name + "\n"
        else:
            ret_str += "No filename provided, will use timestamp instead\n"
        if self.loc53 != '':
            ret_str += "53.0 is stored at \'" + self.loc53 + "\'\n"
        if self.qchem != '':
            ret_str += "Q-Chem version is " + self.qchem + "\n"
        if self.nt != 1:
            ret_str += "Using " + str(self.nt) + " threads\n"
        if self.np != 1:
            ret_str += "Using " + str(self.np) + " processors\n"
        return ret_str

    def info(self):
        print(self)


########################### MULTIFILE  ##############################

class multifile(object):

    def __init__(self, jobs=[]):
        self.list_of_jobs = []
        self.list_of_content = []
        for k in jobs:
            self.add(k)

    def add(self, new_job):
        ''' Adds an inputfile to your batch object.'''
        if type(new_job) == type(inputfile()):
            self.list_of_jobs.append(new_job)
            self.list_of_content.append(new_job._jtype)
        else:
            print("Can only add inputfiles.")

    def remove(self, position=0):  # if not specified delete last
        ''' Removes an inputfile from your batch object. If no other specified the last is removed.'''
        del self.list_of_content[position]
        del self.list_of_jobs[position]

    def __str__(self):
        if self.list_of_jobs == []:
            ret_str = "empty"
        else:
            ret_str = self.list_of_jobs[0].__str__()
        if len(self.list_of_jobs) > 1:
            for k in self.list_of_jobs[1:]:
                ret_str += "@@@\n\n" + k.__str__()
        return ret_str

    def write(self, filename):
        '''Writes the batch jobfile to disk.'''
        f = open(filename, 'w')
        str_ret = self.__str__()
        print(str_ret, file=f)
        f.close()

    def run(self, name='', loc53='', qchem='', nt=1, np=1, timestamp=False):
        '''Makes Q-Chem process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        loc53 ...... 53.0 file location
        nt ......... number of threads
        np ......... number of processors.
        timestamp... adds a timestamp to input and output if set True.

        If nothing specified, pyQChem will fall back on information in the corresponding runinfo object.'''

        running_scripts._run(self, name, loc53, qchem, nt, np, timestamp)


########################### INPUTFILE  ##############################

class inputfile(object):

    def __init__(self, arrays=[]):
        self.list_of_arrays = []
        self.list_of_content = []
        self.runinfo = _rundata()
        self.__jtype = "undef"
        for k in arrays:
            self.add(k)

    def add(self, new_array):
        ''' Adds an array to your inputfile object.'''
        if type(new_array) == type(rem_array()):
            self.rem = new_array
            if "rem" in self.list_of_content:
                index = self.list_of_content.index("rem")
                self.list_of_arrays[index] = new_array
            else:
                self.list_of_content.append("rem")
                self.list_of_arrays.append(new_array)
            self._jtype = new_array.jobtype()  # rem variable "jobtype" defines type

        elif type(new_array) == type(rem_frgm_array()):
            self.rem_frgm = new_array
            if "rem_frgm" in self.list_of_content:
                index = self.list_of_content.index("rem_frgm")
                self.list_of_arrays[index] = new_array
            else:
                self.list_of_content.append("rem_frgm")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(mol_array()):
            self.molecule = new_array
            if "molecule" in self.list_of_content:
                index = self.list_of_content.index("molecule")
                self.list_of_arrays[index] = new_array
            else:
                self.list_of_content.append("molecule")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(comment_array()):
            self.list_of_content.append("comment")
            self.list_of_arrays.append(new_array)

        elif type(new_array) == type(basis_array()):
            self.basis = new_array
            if "basis" in self.list_of_content:
                index = self.list_of_content.index("basis")
                self.list_of_arrays[index] = new_array
            else:
                self.list_of_content.append("basis")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(ecp_array()):
            self.ecp = new_array
            if "ecp" in self.list_of_content:
                index = self.list_of_content.index("ecp")
                self.list_of_arrays[index] = new_array
            else:
                self.list_of_content.append("ecp")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(_unsupported_array()):
            self.list_of_content.append(str(new_array.type))
            self.list_of_arrays.append(new_array)

        # poor man's typecasting because these feel equivalent to users
        elif (type(new_array) == cartesian or type(new_array) == zmat or type(
                new_array) == tinker):
            self.add(mol_array(new_array))

        else:
            print("Array type unknown.")

    def remove(self, position=0):  # if not specified delete last
        ''' Removes an array from your inputfile object. If no other specified the last is removed.'''
        del self.list_of_content[position]
        del self.list_of_arrays[position]

    def __str__(self):
        ret_str = ""
        for k in self.list_of_arrays:
            ret_str += k.__str__() + "\n"
        return ret_str

    def write(self, filename):
        f = open(filename, 'w')
        str_ret = self.__str__()
        print(str_ret, file=f)
        f.close()

    def info(self):
        '''A quick overview of your inputfile.'''  # Health check could be put here
        if "rem" and "molecule" in self.list_of_content:
            status = "valid"
        else:
            status = "invalid"

        print("Type: inputfile")
        print("Status: " + status)

    def run(self, name='', loc53='', qchem='', nt=1, np=1, timestamp=False):
        '''Makes Q-Chem process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        loc53 ...... 53.0 file location
        nt ......... number of threads
        np ......... number of processors.
        timestamp... adds a timestamp to input and output if set True.

        If nothing specified, pyQChem will fall back on information in the corresponding runinfo object.'''

        running_scripts._run(self, name, loc53, qchem, nt, np, timestamp)

    def __add__(self, other):
        # autoadd subarrays - works
        import copy
        new = copy.deepcopy(self)
        new.add(other)
        return new


#    def __radd__(self,other):
#	#not working right, but general form should be this
#	import copy
#	new=copy.deepcopy(self)
#	new.add(other)
#	return new

######################## INPUT FRAGMENTS ############################

class _array(object):

    def __init__(self, content="undef"):
        self.content = content

    def __str__(self):
        ret_str = "$undef\n"
        ret_str += str(self.content) + "\n"
        ret_str += "$end\n"
        return ret_str

    def write(self, filename):
        f = open(filename, 'w')
        str_ret = self.__str__()
        print(str_ret, file=f)
        f.close()

    def __add__(self, other):
        if isinstance(other, _array):
            a = inputfile()
            a.add(self)
            a.add(other)
            return a
        if isinstance(other, inputfile):
            return other + self

    def __radd__(self, other):
        if isinstance(other, _array):
            a = inputfile()
            a.add(self)
            a.add(other)
            return a


##################### UNSUPPORTED FRAGMENT ##########################

class _unsupported_array(_array):

    def __init__(self, arraytype="undef"):
        self.content = []
        self.type = arraytype

    def add_line(self, line):
        self.content.append(line.strip())

    def __str__(self):
        ret_str = "$" + str(self.type) + "\n"
        for k in self.content:
            ret_str += k + "\n"
        ret_str += "$end\n"
        return ret_str


######################### ZMAT FRAGMENT #############################

class zmat(_array):
    __tabstop = 10

    def __init__(self):
        self.__Natoms = 0
        self.__Nvariables = 0
        self.__lines = []
        self.__variables = {}

    def add_atom(self, line, position=0):  # if not specified add at the end
        '''Adds an atom to your Z-Matrix.'''
        if position == 0:
            self.__lines.append(line)
        else:
            self.__lines.insert(position, line)
        self.__Natoms += 1

    def remove_atom(self, position=0):  # if not specified delete last
        '''Removes an atom from your Z-Matrix. Takes the last if no other specified.'''
        del self.__lines[position]
        self.__Natoms -= 1

    def variable(self, key="variable", value="show"):
        '''Adds, changes or removes variable definitions.'''
        if value == "" and key in self.__variables:
            del self.__variables[key]
        elif value == "show":
            return self.__variables[key]
        else:
            self.__variables[key] = value

    def __str__(self):
        ret_str = ""
        for k in self.__lines:
            ret_str += k + "\n"
        ret_str += "\n"
        for key, value in self.__variables.items():
            ret_str += key + " " * (zmat.__tabstop - len(key)) + value + "\n"
        return ret_str


####################### CARTESIAN FRAGMENT ##########################

class cartesian(_array):
    def __init__(self, title="", atom_list=[]):
        import copy
        self.__title = title
        self.__Natoms = 0
        self.xyzs = []
        self.com = _np.array([0.0, 0.0, 0.0])
        self.centroid = _np.array([0.0, 0.0, 0.0])
        self.list_of_atoms = copy.deepcopy(atom_list)
        if atom_list != []:
            self.__Natoms = len(atom_list)
            for i in range(self.__Natoms):
                x = self.list_of_atoms[i][1]
                y = self.list_of_atoms[i][2]
                z = self.list_of_atoms[i][3]
                self.xyzs.append(_np.array([float(x), float(y), float(z)]))
            self.xyzs = _np.array(self.xyzs)
            self.__center_of_mass()

    def fix(self):
        """This fixes any odd errors resulting from modifying the number of atoms"""
        self.__Natoms = len(self.list_of_atoms)
        if self.__Natoms == 0:
            return
        self.xyzs = []
        for i in range(self.__Natoms):
            x = self.list_of_atoms[i][1]
            y = self.list_of_atoms[i][2]
            z = self.list_of_atoms[i][3]
            self.xyzs.append(_np.array([float(x), float(y), float(z)]))
        self.xyzs = _np.array(self.xyzs)
        self.__center_of_mass()

    def __center_of_mass(self):
        """This computes the centroid and center of mass using standard atomic masses"""
        # print self.xyzs, self.__Natoms
        self.com = _np.array([0.0, 0.0, 0.0])
        self.centroid = _np.array([0.0, 0.0, 0.0])
        if len(self.xyzs) == 0:
            return
        total_mass = 0.0
        self.centroid = sum(self.xyzs) / len(self.xyzs)
        wts = [constants.dict_of_atomic_masses[
                   self.list_of_atoms[i][0].replace("@", "")] for i in
               range(self.__Natoms)]
        for i, atom in enumerate(self.xyzs):
            wt = wts[i]
            total_mass = total_mass + wt
            self.com = self.com + atom * wt
        self.centroid = _np.array([i / self.__Natoms for i in self.centroid])
        self.com = _np.array([i / total_mass for i in self.com])

    def title(self, title="show"):
        if title == "show":
            return self.__title
        else:
            self.__title = title

    def add_atom(self, name="H", x="0", y="0", z="0"):
        self.list_of_atoms.append([name, x, y, z])
        self.fix()
        self.__center_of_mass()

    def remove_atom(self, position=0):
        del self.list_of_atoms[position]  # First atom is atom 1
        self.fix()
        self.__center_of_mass()

    def ghost(self):
        atoms = []
        for i in range(self.__Natoms):
            atoms.append(
                ['@' + self.list_of_atoms[i][0], self.list_of_atoms[i][1],
                 self.list_of_atoms[i][2], self.list_of_atoms[i][3]])
        return atoms

    def atoms(self):
        for i, k in enumerate(self.list_of_atoms):
            print(str(i + 1) + ":\t" + k[0] + "\t" + k[1] + "\t" + k[2] + "\t" +
                  k[3])

    def atomic_distance(self, a, b):
        """Gives the pair-wise distance between two atoms (counting from 0)"""
        from math import sqrt
        d = self.xyzs[b] - self.xyzs[a]
        return sqrt(d.dot(d))

    def print_centroid(self):
        print(str(self.centroid[0]) + '\t' + str(self.centroid[1]) + '\t' + str(
            self.centroid[2]))

    def print_center_of_mass(self):
        print(str(self.com[0]) + '\t' + str(self.com[1]) + '\t' + str(
            self.com[2]))

    def move(self, dir, amt=1.0):
        dir = _np.array(dir)
        for i in range(self.__Natoms):
            self.xyzs[i] = self.xyzs[i] + dir * amt
            self.list_of_atoms[i][1] = str(self.xyzs[i][0])
            self.list_of_atoms[i][2] = str(self.xyzs[i][1])
            self.list_of_atoms[i][3] = str(self.xyzs[i][2])
        self.__center_of_mass()

    def __str__(self):
        str_ret = str(self.__Natoms) + "\n" + self.__title + "\n"
        for k in self.list_of_atoms:
            str_ret += k[0] + "    " + k[1] + "    " + k[2] + "    " + k[
                3] + "\n"
        return str_ret

    def __sub__(self, other):
        if type(other) == type([]):  # let's move the atoms
            self.move(other, -1.0)
        if type(other) == type(
                self.com):  # let's move the atoms using a numpy array
            self.move(other, -1.0)

    def __add__(self, other):
        if type(other) == type([]):  # let's move the atoms
            self.move(other, 1.0)
        if type(other) == type(
                self.com):  # let's move the atoms using a numpy array
            self.move(other, 1.0)
        if type(other) == type(self):  # merge two cartesians
            atoms = self.list_of_atoms + other.list_of_atoms
            return cartesian(atom_list=atoms)
        if isinstance(other, _array):
            return other + mol_array(self)

    def __radd__(self, other):  # reverse of above
        if type(other) == type([]):
            self.move(other)
        if type(other) == type(self):
            atoms = self.list_of_atoms + other.list_of_atoms
            return cartesian(atom_list=atoms)
        if isinstance(other, _array):
            return other + mol_array(self)


####################### FRAGMENT FRAGMENT ##########################

class fragment(cartesian):
    def __init__(self, title="", fragment_list=[], atom_list=[]):
        self.__title = title
        if len(fragment_list) == 0:
            self.__Natoms = 0
            self.list_of_atoms = []
            self.number_of_fragments = 0
            self.fragment_list = []
        else:
            self.fragment_list = fragment_list
            self.__Natoms = 0
            for i in self.fragment_list:
                self.__Natoms += i.__Natoms
                self.__Natoms = 0

        if (atom_list != []):
            xyzs = []
            self.__Natoms = len(atom_list)
            for i in range(self.__Natoms):
                x = atom_list[i][1]
                y = atom_list[i][2]
                z = atom_list[i][3]
                xyzs.append(_np.array([float(x), float(y), float(z)]))
            self.xyzs = _np.array(xyzs)
            from .constants import covalent_radii, dict_of_atomic_numbers
            from math import sqrt
            nleft = len(atom_list)
            ifrag = cartesian(atom_list=[atom_list[0]])
            # print atom_list[0]
            # ifrag.add_atom(self,name=atom_list[0][0],x=atom_list[0][1],y=atom_list[0][2],z=atom_list[0][3])
            self.fragment_list.append(ifrag)
            for i in range(1, self.__Natoms):
                myxyz = self.xyzs[i]
                included = 0
                for myfrag in self.fragment_list:
                    for k in range(len(myfrag.xyzs)):
                        d = (myfrag.xyzs[k] - myxyz)
                        d = sqrt(d.dot(d))
                        try:
                            myA = dict_of_atomic_numbers[atom_list[i][0]]
                        except:
                            myA = atom_list[i][0]
                        try:
                            myB = dict_of_atomic_numbers[
                                myfrag.list_of_atoms[k][0]]
                        except:
                            myB = myfrag.list_of_atoms[k][0]
                        maxd = covalent_radii[myA] + covalent_radii[myB]
                        maxd = maxd * 1.3
                        if d < maxd:
                            myfrag.add_atom(name=atom_list[i][0],
                                            x=atom_list[i][1],
                                            y=atom_list[i][2],
                                            z=atom_list[i][3])
                            included = 1
                            break
                if included == 0:
                    ifrag = cartesian(atom_list=[atom_list[i]])
                    self.fragment_list.append(ifrag)

    def __str__(self):
        str_ret = str(self.__Natoms) + "\n" + self.__title + "\n"
        for l in self.fragment_list:
            for k in l.list_of_atoms:
                str_ret += k[0] + "    " + k[1] + "    " + k[2] + "    " + k[
                    3] + "\n"
        return str_ret


####################### TINKER FRAGMENT ##########################

class tinker(cartesian):

    def __init__(self, title=""):
        self.__title = title
        self.__Natoms = 0
        self.list_of_atoms = []
        self.dict_of_types = {}  # dictionary of atom types

    def title(self, title="show"):
        if title == "show":
            return self.__title
        else:
            self.__title = title

    def remove_atom(self, position):
        del self.list_of_atoms[position]  # First atom is atom 1
        self.__Natoms -= 1

    def add_atom(self, name="H", x="0", y="0", z="0", atomtype="0", con1="0",
                 con2="0", con3="0", con4="0"):
        if name not in self.dict_of_types:
            self.dict_of_types[name] = atomtype
        self.list_of_atoms.append(
            [name, x, y, z, atomtype, con1, con2, con3, con4])
        self.__Natoms += 1

    def change_type(self, atomtype="none", value="show"):
        '''Changes the definition number of atom "atomtype" to "value".'''
        if value == "" and atomtype in self.dict_of_types:
            del self.dict_of_types[atomtype]
        elif value == "show":
            return self.dict_of_types
        else:
            self.dict_of_types[atomtype] = value
            # atomtype definition has changed, list_of_atoms needs to be updated:
            Nchanges = 0
            for k in self.list_of_atoms:
                if k[0] == atomtype:
                    k[4] = value
                    Nchanges += 1
            print("Atomtype definition has changed, " + str(
                Nchanges) + " atoms updated in list_of_atoms.")

    def atoms(self):
        for i, k in enumerate(self.list_of_atoms):
            print(str(i + 1) + ":\t" + k[0] + "    " + k[1] + "    " + k[2] + \
                  "    " + k[3] + "    " + k[4] + "    " + k[5] + "    " + k[
                      6] + "    " + k[7] + "    " + k[8])

    def __str__(self):
        str_ret = str(self.__Natoms) + "\t" + self.__title + "\n"
        for i, k in enumerate(self.list_of_atoms):
            str_ret += str(i + 1) + "\t" + k[0] + "    " + k[1] + "    " + k[
                2] + \
                       "    " + k[3] + "    " + k[4] + "    " + k[5] + "    " + \
                       k[6] + "    " + \
                       k[7] + "    " + k[8] + "\n"
        return str_ret


######################### MOL FRAGMENT ##############################

class mol_array(_array):

    def __init__(self, geometry=""):
        if geometry == "":
            geometry = cartesian()
        self.content = {"CHARGE": "0", "MULTIPLICITY": "1",
                        "GEOMETRY": geometry}

    def charge(self, value="show"):
        '''Total charge of the molecule.'''
        if value == "show":
            return self.content["CHARGE"]
        else:
            self.content["CHARGE"] = value

    def multiplicity(self, value="show"):
        '''Spin multiplicity of the molecule.'''
        if value == "show":
            return self.content["MULTIPLICITY"]
        else:
            self.content["MULTIPLICITY"] = value

    def geometry(self, value="show"):
        '''Reads xyz, txyz or zmat coordinate array.'''
        if value == "show":
            return self.content["GEOMETRY"]
        elif value == "read":
            self.content["GEOMETRY"] = value
        else:
            if type(value) == type(cartesian()) or type(value) == type(
                    zmat()) or type(value) == type(tinker()):
                self.content["GEOMETRY"] = value
            else:
                print(
                    "Only cartesian, tinker or zmat arrays can be added here.")

    def clear(self):
        self.content = {"CHARGE": "0", "MULTIPLICITY": "1", "GEOMETRY": ""}

    def __str__(self):
        if self.content["GEOMETRY"] == "read":
            str_ret = "$molecule\nread\n$end\n"
        else:
            str_ret = "$molecule\n" + self.content["CHARGE"] + " " + \
                      self.content["MULTIPLICITY"] + "\n"
            if type(self.content["GEOMETRY"]) == type(cartesian()):
                for k in (self.content["GEOMETRY"]).list_of_atoms:
                    str_ret += k[0] + "    " + k[1] + "    " + k[2] + "    " + \
                               k[3] + "\n"
            elif type(self.content["GEOMETRY"]) == type(fragment()):
                for l in self.content["GEOMETRY"].fragment_list:
                    str_ret += "--\n0 1\n"
                    for k in l.list_of_atoms:
                        str_ret += k[0] + "    " + k[1] + "    " + k[
                            2] + "    " + k[3] + "\n"
            elif type(self.content["GEOMETRY"]) == type(zmat()):
                str_ret += (self.content["GEOMETRY"]).__str__()
            elif type(self.content["GEOMETRY"]) == type(tinker()):
                for k in (self.content["GEOMETRY"]).list_of_atoms:
                    str_ret += k[0] + "    " + k[1] + "   " + k[2] + \
                               "    " + k[3] + "    " + k[4] + "    " + k[5] + \
                               "    " + k[6] + "    " + k[7] + "    " + k[
                                   8] + "\n"
            str_ret += "$end\n"
        return str_ret

    def info(self):
        switch = 0
        if type(self.content["GEOMETRY"]) == type(cartesian()):
            coor_type = "cartesian coordinates"
            switch = 1
        elif type(self.content["GEOMETRY"]) == type(zmat()):
            coor_type = "Z-Matrix"
            switch = 1
        elif type(self.content["GEOMETRY"]) == type(tinker()):
            coor_type = "Tinker"
            switch = 1
        else:
            coor_type = "empty"
        print("Type: molecule array, " + coor_type)
        if switch == 1:
            print("Number of atoms: " + str(
                len((self.content["GEOMETRY"]).list_of_atoms)))


######################### BASIS FRAGMENT ############################

class basis_array(_array):

    def __init__(self):
        self.dict_of_atoms = {}

    def add(self, atom, line):
        if atom in self.dict_of_atoms:
            self.dict_of_atoms[atom].append(line)
        else:
            self.dict_of_atoms[atom] = list()
            self.dict_of_atoms[atom].append(line)

    def remove(self, atom):
        if atom in self.dict_of_atoms:
            del self.dict_of_atoms[atom]

    def __str__(self):
        ret_str = "$basis\n"
        for key in self.dict_of_atoms:
            ret_str += str(key) + "\n"
            for line in self.dict_of_atoms[key]:
                ret_str += line + "\n"
            ret_str += "****\n"
        ret_str += "$end\n"
        return ret_str


########################## ECP FRAGMENT #############################

class ecp_array(basis_array):

    def __str__(self):
        ret_str = "$ecp\n"
        for key in self.dict_of_atoms:
            ret_str += str(key) + "\n"
            for line in self.dict_of_atoms[key]:
                ret_str += line + "\n"
            ret_str += "****\n"
        ret_str += "$end\n"
        return ret_str


########################## OPT FRAGMENT #############################

######################### ALIST FRAGMENT ############################

####################### QM_ATOMS FRAGMENT ###########################

################## FORCE_FIELD_PARAMS FRAGMENT ######################

####################### COMMENT FRAGMENT ############################

class comment_array(_array):

    def __init__(self, content=""):
        self.content = content

    def __str__(self):
        ret_str = "$comment\n"
        lines = list(
            filter(bool, self.content.split("\n")))  # remove empty list entries
        for line in lines:
            ret_str += line.strip() + "\n"
        ret_str += "$end\n"
        return ret_str


######################### REM FRAGMENT ##############################

class rem_array(_array):
    __tabstop = 30

    def __init__(self, rem_init=""):
        self.dict_of_keywords = {"JOBTYPE": "sp"}
        rem_init = rem_init.splitlines()
        if len(rem_init) != 0:
            for i in rem_init:
                i = i.split(" ")
                if len(i) == 0:
                    i = i.split("=")
                if i[0].startswith("$"):
                    continue
                self.add(i[0], i[1])

    # -------------- Computer-generated List of REM keywords  -----------------

    def cc_dip(self, value="show"):
        '''
Name: CC_DIP
Type: INTEGER

Description: Initializes a EOM-DIP-CCSD calculation
    '''
        if value == "":
            if "CC_DIP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIP"] = value.lower()

    def cc_do_triples(self, value="show"):
        '''
Name: CC_DO_TRIPLES
Type: INTEGER

Description: This keyword initializes a EOM-CC(2,3) calculation. If {CC_IP_PROPER} is set then EOM-IP-CC(2,3) states are calculated.
    '''
        if value == "":
            if "CC_DO_TRIPLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_TRIPLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DO_TRIPLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_TRIPLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DO_TRIPLES"] = value.lower()

    def cc_iterate_on(self, value="show"):
        '''
Name: CC_ITERATE_ON
Type: INTEGER

Description: In active space calculations, use a mixed iteration procedure if the value is greater than 0.  Then if the RMS orbital gradient is larger than the value of CC_THETA_GRAD_THRESH, micro-iterations will be performed to converge the occupied-virtual mixing angles for the current active-space. The maximum number of space iterations is given by this option.
Recommendation: : Can be useful for non-convergent active space calculations    '''
        if value == "":
            if "CC_ITERATE_ON" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ITERATE_ON"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_ITERATE_ON" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ITERATE_ON"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_ITERATE_ON"] = value.lower()

    def cc_orbs_per_block(self, value="show"):
        '''
Name: CC_ORBS_PER_BLOCK
Type: INTEGER

Description: Specifies target (and maximum) size of blocks in orbital space.
    '''
        if value == "":
            if "CC_ORBS_PER_BLOCK" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ORBS_PER_BLOCK"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_ORBS_PER_BLOCK" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ORBS_PER_BLOCK"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_ORBS_PER_BLOCK"] = value.lower()

    def cc_preconv_sd(self, value="show"):
        '''
Name: CC_PRECONV_SD
Type: INTEGER

Description: Solves the EOM-CCSD equations, prints energies, then uses EOM-CCSD vectors  as initial vectors in EOM-CC(2,3). Very convenient for calculations using energy additivity schemes.
Recommendation: : Turning this option on is recommended    '''
        if value == "":
            if "CC_PRECONV_SD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_SD"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONV_SD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_SD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONV_SD"] = value.lower()

    def cc_reset_theta(self, value="show"):
        '''
Name: CC_RESET_THETA
Type: INTEGER

Description: The reference MO coefficient matrix is reset every n iterations to help   overcome problems associated with the theta metric as theta becomes large. 
    '''
        if value == "":
            if "CC_RESET_THETA" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESET_THETA"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_RESET_THETA" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESET_THETA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_RESET_THETA"] = value.lower()

    def cc_restr_ampl(self, value="show"):
        '''
Name: CC_RESTR_AMPL
Type: INTEGER

Description: Controls the restriction on amplitudes is there are restricted orbitals
    '''
        if value == "":
            if "CC_RESTR_AMPL" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTR_AMPL"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_RESTR_AMPL" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTR_AMPL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_RESTR_AMPL"] = value.lower()

    def cc_restr_triples(self, value="show"):
        '''
Name: CC_RESTR_TRIPLES
Type: INTEGER

Description: Controls which space the triples correction is computed in
    '''
        if value == "":
            if "CC_RESTR_TRIPLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTR_TRIPLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_RESTR_TRIPLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTR_TRIPLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_RESTR_TRIPLES"] = value.lower()

    def cc_rest_occ(self, value="show"):
        '''
Name: CC_REST_OCC
Type: INTEGER

Description: Sets the number of restricted occupied orbitals including frozen occupied  orbitals.
    '''
        if value == "":
            if "CC_REST_OCC" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REST_OCC"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REST_OCC" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REST_OCC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REST_OCC"] = value.lower()

    def cc_rest_vir(self, value="show"):
        '''
Name: CC_REST_VIR
Type: INTEGER

Description: Sets the number of restricted virtual orbitals including frozen virtual  orbitals.
    '''
        if value == "":
            if "CC_REST_VIR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REST_VIR"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REST_VIR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REST_VIR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REST_VIR"] = value.lower()

    def cc_theta_grad_thresh(self, value="show"):
        '''
Name: CC_THETA_GRAD_THRESH
Type: INTEGER

Description: RMS orbital gradient threshold [10-n] above which mixed iterations  are performed in active space calculations if CC_ITERATE_OV is  TRUE.
Recommendation: : Can be made smaller if convergence difficulties are encountered.    '''
        if value == "":
            if "CC_THETA_GRAD_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_GRAD_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_THETA_GRAD_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_GRAD_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_THETA_GRAD_THRESH"] = value.lower()

    def cc_tmpbuffsize(self, value="show"):
        '''
Name: CC_TMPBUFFSIZE
Type: INTEGER

Description: Maximum size, in Mb, of additional buffers for temporary arrays used to work with individual blocks or matrices.
Recommendation: : Should not be smaller than the size of the largest possible block.    '''
        if value == "":
            if "CC_TMPBUFFSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TMPBUFFSIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_TMPBUFFSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TMPBUFFSIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_TMPBUFFSIZE"] = value.lower()

    def gvb_orb_conv(self, value="show"):
        '''
Name: GVB_ORB_CONV
Type: INTEGER

Description: The GVB-CC wave function is considered converged when the root-mean-square  orbital gradient and orbital step sizes are less than  10-GVB_ORB_CONV. Adjust THRESH simultaneously.
Recommendation: : Use 6 for PP(2) jobs or geometry optimizations. Tighter convergence (i.e. 7 or higher) cannot always be reliably achieved.    '''
        if value == "":
            if "GVB_ORB_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_ORB_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_ORB_CONV"] = value.lower()

    def gvb_orb_max_iter(self, value="show"):
        '''
Name: GVB_ORB_MAX_ITER
Type: INTEGER

Description: Controls the number of orbital iterations allowed in GVB-CC calculations. Some jobs, particularly unrestricted PP jobs can require 500-1000 iterations.
Recommendation: : Default is typically adequate, but some jobs, particularly UPP jobs, can  require 500-1000 iterations if converged tightly.    '''
        if value == "":
            if "GVB_ORB_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_MAX_ITER"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_ORB_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_MAX_ITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_ORB_MAX_ITER"] = value.lower()

    def gvb_orb_scale(self, value="show"):
        '''
Name: GVB_ORB_SCALE
Type: INTEGER

Description: Scales the default orbital step size by n/1000.
Recommendation: : Default is usually fine, but for some stretched geometries it can help with convergence to use smaller values.    '''
        if value == "":
            if "GVB_ORB_SCALE" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_SCALE"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_ORB_SCALE" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_SCALE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_ORB_SCALE"] = value.lower()

    def gvb_restart(self, value="show"):
        '''
Name: GVB_RESTART
Type: STRING
Default: FALSE

Options:
    'FALSE'......................... FALSE

Description: Restart a job from previously-converged GVB-CC orbitals.
Recommendation: : Useful when trying to converge to the same GVB solution at slightly different geometries, for example.    '''
        if value == "":
            if "GVB_RESTART" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_RESTART"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_RESTART" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_RESTART"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_RESTART"] = value.lower()

    def gvb_unrestricted(self, value="show"):
        '''
Name: GVB_UNRESTRICTED
Type: STRING
Default: same value as UNRESTRICTED

Options:
    'same value as UNRESTRICTED'.... same value as UNRESTRICTED

Description: Controls restricted versus unrestricted PP jobs. Usually handled  automatically.
Recommendation: : Set this variable explicitly only to do a UPP job from an RHF or ROHF initial guess.    '''
        if value == "":
            if "GVB_UNRESTRICTED" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_UNRESTRICTED"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_UNRESTRICTED" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_UNRESTRICTED"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_UNRESTRICTED"] = value.lower()

    def gvb_print(self, value="show"):
        '''
Name: GVB_PRINT
Type: INTEGER

Description: Controls the amount of information printed during a GVB-CC job.
Recommendation: : Should never need to go above 0 or 1.    '''
        if value == "":
            if "GVB_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_PRINT"] = value.lower()

    def nvo_lin_convergence(self, value="show"):
        '''
Name: NVO_LIN_CONVERGENCE
Type: INTEGER

Description: Target error factor in the preconditioned conjugate gradient solver of the single-excitation amplitude equations.
Recommendation: : Solution of the single-excitation amplitude equations is considered converged if the maximum residual is less than 10-n multiplied by the current DIIS error. For the ARS correction, n is automatically set to 1 since the locally-projected DIIS error is normally several orders of magnitude smaller than the full DIIS error.    '''
        if value == "":
            if "NVO_LIN_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_LIN_CONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_LIN_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_LIN_CONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_LIN_CONVERGENCE"] = value.lower()

    def nvo_lin_max_ite(self, value="show"):
        '''
Name: NVO_LIN_MAX_ITE
Type: INTEGER

Description: Maximum number of iterations in the preconditioned conjugate gradient solver of the single-excitation amplitude equations.
    '''
        if value == "":
            if "NVO_LIN_MAX_ITE" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_LIN_MAX_ITE"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_LIN_MAX_ITE" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_LIN_MAX_ITE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_LIN_MAX_ITE"] = value.lower()

    def nvo_truncate_precond(self, value="show"):
        '''
Name: NVO_TRUNCATE_PRECOND
Type: INTEGER

Description: Specifies which atomic blocks of the Fock matrix are used to construct the preconditioner. This variable is used only if NVO_TRUNCATE_DIST is set to -2.
Recommendation: : Use default. Increasing n improves convergence of the PCG algorithm but overall may slow down the calculations.    '''
        if value == "":
            if "NVO_TRUNCATE_PRECOND" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_TRUNCATE_PRECOND"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_TRUNCATE_PRECOND" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_TRUNCATE_PRECOND"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_TRUNCATE_PRECOND"] = value.lower()

    def nvo_uvv_maxpwr(self, value="show"):
        '''
Name: NVO_UVV_MAXPWR
Type: INTEGER

Description: Controls convergence of the Taylor series when calculating the Uvv block from the single-excitation amplitudes. If the series is not converged at the nth term, more expensive direct inversion is used to calculate the Uvv block.
    '''
        if value == "":
            if "NVO_UVV_MAXPWR" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_UVV_MAXPWR"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_UVV_MAXPWR" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_UVV_MAXPWR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_UVV_MAXPWR"] = value.lower()

    def nvo_uvv_precision(self, value="show"):
        '''
Name: NVO_UVV_PRECISION
Type: INTEGER

Description: Controls convergence of the Taylor series when calculating the Uvv block from the single-excitation amplitudes. Series is considered converged when the maximum element of the term is less than 10-n.
Recommendation: : NVO_UVV_PRECISION must be the same as or larger than THRESH.    '''
        if value == "":
            if "NVO_UVV_PRECISION" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_UVV_PRECISION"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_UVV_PRECISION" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_UVV_PRECISION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_UVV_PRECISION"] = value.lower()

    def print_dist_matrix(self, value="show"):
        '''
Name: PRINT_DIST_MATRIX
Type: INTEGER

Description: Controls the printing of the inter-atomic distance matrix
Recommendation: : Use default unless distances are required for large systems    '''
        if value == "":
            if "PRINT_DIST_MATRIX" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_DIST_MATRIX"]
                print("Keyword removed.")
        elif value == "show":
            if "PRINT_DIST_MATRIX" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_DIST_MATRIX"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PRINT_DIST_MATRIX"] = value.lower()

    def rc_r0(self, value="show"):
        '''
Name: RC_R0
Type: INTEGER

Description: Determines the parameter in the Gaussian weight function used to smooth the  density at the nuclei.
Recommendation: : We recommend value of 250 for a typical spit valence basis.  For basis sets  with increased flexibility in the nuclear vicinity the smaller values of r_0  also yield adequate spin density.    '''
        if value == "":
            if "RC_R0" in self.dict_of_keywords:
                del self.dict_of_keywords["RC_R0"]
                print("Keyword removed.")
        elif value == "show":
            if "RC_R0" in self.dict_of_keywords:
                return self.dict_of_keywords["RC_R0"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RC_R0"] = value.lower()

    def svp_cavity_conv(self, value="show"):
        '''
Name: SVP_CAVITY_CONV
Type: INTEGER

Description: Determines the convergence value of the iterative isodensity cavity procedure.
Recommendation: : The default value unless convergence problems arise.    '''
        if value == "":
            if "SVP_CAVITY_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_CAVITY_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP_CAVITY_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_CAVITY_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP_CAVITY_CONV"] = value.lower()

    def rpath_max_cycles(self, value="show"):
        '''
Name: RPATH_MAX_CYCLES
Type: INTEGER
Default: 20
Options: Range from 1 to 500

Description: Specifies the maximum number of points to find on the reaction path.
Recommendation: : Use more points if the minimum is desired, but not reached using the default.    '''
        if value == "":
            if "RPATH_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_MAX_CYCLES"] = value.lower()

    def rpath_tol_displacement(self, value="show"):
        '''
Name: RPATH_TOL_DISPLACEMENT
Type: INTEGER
Factor: 0.0001
Default: 50 [=0.0050]
Options: Range from 1 [=0.0001] to 2000 [=0.2000]

Description: Specifies the convergence threshold (in a.u.) for the step. If a step size is chosen by the algorithm that is smaller than this, the path is deemed to have reached the minimum.
    '''
        if value == "":
            if "RPATH_TOL_DISPLACEMENT" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_TOL_DISPLACEMENT" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"] = value.lower()

    def rpath_print(self, value="show"):
        '''
Name: RPATH_PRINT
Type: INTEGER
Default: 2
Options: Range from 1 to 5

Description: Controls the level of output for a reaction coordinate calculation.
Recommendation: : Use default, little additional information is printed at higher levels. Most of the output arises from the multiple single point calculations that are performed along the reaction pathway.    '''
        if value == "":
            if "RPATH_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_PRINT"] = value.lower()

    def scf_convergence(self, value="show"):
        '''
Name: SCF_CONVERGENCE
Type: INTEGER
Default: 5
Options: Range from 0 to 12

Description: SCF is considered converged when the wavefunction error is less that 10-SCF_CONVERGENCE. Adjust the value of THRESH at the same time. Note that in Q-Chem 3.0 the DIIS error is measured by the maximum error rather than the RMS error as in previous versions.
Recommendation: : Tighter criteria for geometry optimization and vibration analysis. Larger values provide more significant figures, at greater computational cost.    '''
        if value == "":
            if "SCF_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_CONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_CONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_CONVERGENCE"] = value.lower()

    def aimd_steps(self, value="show"):
        '''
Name: AIMD_STEPS
Type: INTEGER
Default: 0
Options: Range from 0 to 500

Description: Specifies the requested number of molecular dynamics steps.
    '''
        if value == "":
            if "AIMD_STEPS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_STEPS"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_STEPS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_STEPS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_STEPS"] = value.lower()

    def ao2mo_disk(self, value="show"):
        '''
Name: AO2MO_DISK
Type: INTEGER
Default: 2000
Options: Range from 0 to 8000

Description: Sets the amount of disk space (in megabytes) available for MP2 calculations. 

Recommendation: : Should be set as large as possible, as discussed in the manual.    '''
        if value == "":
            if "AO2MO_DISK" in self.dict_of_keywords:
                del self.dict_of_keywords["AO2MO_DISK"]
                print("Keyword removed.")
        elif value == "show":
            if "AO2MO_DISK" in self.dict_of_keywords:
                return self.dict_of_keywords["AO2MO_DISK"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AO2MO_DISK"] = value.lower()

    def cc_canonize_final(self, value="show"):
        '''
Name: CC_CANONIZE_FINAL
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether to semi-canonicalize orbitals at the end of the ground state calculation.
Recommendation: : Should not normally have to be altered.    '''
        if value == "":
            if "CC_CANONIZE_FINAL" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE_FINAL"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CANONIZE_FINAL" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE_FINAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CANONIZE_FINAL"] = value.lower()

    def cc_canonize(self, value="show"):
        '''
Name: CC_CANONIZE
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether to semi-canonicalize orbitals at the start of the calculation (i.e. Fock matrix is diagonalized in each orbital subspace)
Recommendation: : Should not normally have to be altered.    '''
        if value == "":
            if "CC_CANONIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CANONIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CANONIZE"] = value.lower()

    def cc_convergence(self, value="show"):
        '''
Name: CC_CONVERGENCE
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: Overall convergence criterion for the coupled-cluster codes. This is designed to ensure at least n significant digits in the calculated energy, and automatically sets the other convergence-related variables (CC_E_CONV, CC_T_CONV, CC_THETA_CONV, CC_THETA_GRAD_CONV, CC_Z_CONV) [10-n].
    '''
        if value == "":
            if "CC_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CONVERGENCE"] = value.lower()

    def cc_dconvergence(self, value="show"):
        '''
Name: CC_DCONVERGENCE
Type: INTEGER
Default: 5
Options: Range from 0 to 12

Description: Convergence criterion for the RMS residuals of excited state vectors
Recommendation: : Use default. Should normally be set to the same value as CC_DTHRESHOLD.    '''
        if value == "":
            if "CC_DCONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DCONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DCONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DCONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DCONVERGENCE"] = value.lower()

    def cdft(self, value="show"):
        '''
Name: CDFT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Initiates a constrained DFT calculation
Recommendation: : Set to TRUE if a Constrained DFT calculation is desired.    '''
        if value == "":
            if "CDFT" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT"]
                print("Keyword removed.")
        elif value == "show":
            if "CDFT" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CDFT"] = value.lower()

    def cd_algorithm(self, value="show"):
        '''
Name: CD_ALGORITHM
Type: STRING
Default: Program determined

Options:
    'Program determined'............ Program determined
    'Direct'........................ Direct
    'Semi-direct'................... Semi-direct
    'Local-occupied'................ Local-occupied

Description: Determines the algorithm for MP2 integral transformations.
Recommendation: : Semi-direct is usually most efficient, and will normally be chosen by default.    '''
        if value == "":
            if "CD_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["CD_ALGORITHM"]
                print("Keyword removed.")
        elif value == "show":
            if "CD_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["CD_ALGORITHM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CD_ALGORITHM"] = value.lower()

    def cdft_thresh(self, value="show"):
        '''
Name: CDFT_THRESH
Type: INTEGER
Default: 5
Options: Range from 0 to 12

Description: Determines how tightly the constraint must be satisfied.

Recommendation: : Use default unless problems occur.    '''
        if value == "":
            if "CDFT_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "CDFT_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CDFT_THRESH"] = value.lower()

    def cdft_postdiis(self, value="show"):
        '''
Name: CDFT_POSTDIIS
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether the connstraint is enforced after DIIS extrapolation.
Recommentation: Use default unless convergence problems arise, in which case it may be benificial to turn this option off.  If selected, energies should be variational after the first iteration.
    '''
        if value == "":
            if "CDFT_POSTDIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_POSTDIIS"]
                print("Keyword removed.")
        elif value == "show":
            if "CDFT_POSTDIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_POSTDIIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CDFT_POSTDIIS"] = value.lower()

    def cdft_prediis(self, value="show"):
        '''
Name: CDFT_PREDIIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls wheter the constraint is enforced before DIIS extrapolation.

Recommendation: : Use default unless problems arise, in which case it might be beneficial to turn this option on.    '''
        if value == "":
            if "CDFT_PREDIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_PREDIIS"]
                print("Keyword removed.")
        elif value == "show":
            if "CDFT_PREDIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_PREDIIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CDFT_PREDIIS"] = value.lower()

    def cfmm_order(self, value="show"):
        '''
Name: CFMM_ORDER
Type: INTEGER
Default: 15
Options: Range from 5 to 30

Description: Controls the order of the multipole expansions in CFMM calculation.
    '''
        if value == "":
            if "CFMM_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["CFMM_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "CFMM_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["CFMM_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CFMM_ORDER"] = value.lower()

    def chemsol_nn(self, value="show"):
        '''
Name: CHEMSOL_NN
Type: INTEGER
Default: 5
Options: Range from 1 to 20

Description: Sets the number of grids used to calculate the average hydration free energy.
    '''
        if value == "":
            if "CHEMSOL_NN" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_NN"]
                print("Keyword removed.")
        elif value == "show":
            if "CHEMSOL_NN" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_NN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHEMSOL_NN"] = value.lower()

    def cis_convergence(self, value="show"):
        '''
Name: CIS_CONVERGENCE
Type: INTEGER
Default: 6
Options: Range from 0 to 12

Description: CIS is considered converged when error is less than 10-CIS_CONVERGENCE
    '''
        if value == "":
            if "CIS_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_CONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_CONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_CONVERGENCE"] = value.lower()

    def cis_guess_disk(self, value="show"):
        '''
Name: CIS_GUESS_DISK
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Read the CIS guess from disk (previous calculation)
Recommendation: : Requires a guess from previous calculation.    '''
        if value == "":
            if "CIS_GUESS_DISK" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_GUESS_DISK"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_GUESS_DISK" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_GUESS_DISK"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_GUESS_DISK"] = value.lower()

    def cis_relaxed_density(self, value="show"):
        '''
Name: CIS_RELAXED_DENSITY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Determines whether or not to use the relaxed CIS density for attachment/detachment density analysis
    '''
        if value == "":
            if "CIS_RELAXED_DENSITY" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RELAXED_DENSITY"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RELAXED_DENSITY" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RELAXED_DENSITY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RELAXED_DENSITY"] = value.lower()

    def cis_singlets(self, value="show"):
        '''
Name: CIS_SINGLETS
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Solve for singlet excited states in RCIS calculations (ignored for UCIS)
    '''
        if value == "":
            if "CIS_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_SINGLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_SINGLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_SINGLETS"] = value.lower()

    def cis_triplets(self, value="show"):
        '''
Name: CIS_TRIPLETS
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Solve for triplet excited states in RCIS calculations (ignored for UCIS)
    '''
        if value == "":
            if "CIS_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_TRIPLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_TRIPLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_TRIPLETS"] = value.lower()

    def core_character(self, value="show"):
        '''
Name: CORE_CHARACTER
Type: INTEGER
Default: 0
Options: Range from 0 to 4

Description: Selects how the core orbitals are determined in the frozen-core approximation.
Recommendation: : Use default, unless performing calculations on molecules with heavy elements.    '''
        if value == "":
            if "CORE_CHARACTER" in self.dict_of_keywords:
                del self.dict_of_keywords["CORE_CHARACTER"]
                print("Keyword removed.")
        elif value == "show":
            if "CORE_CHARACTER" in self.dict_of_keywords:
                return self.dict_of_keywords["CORE_CHARACTER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CORE_CHARACTER"] = value.lower()

    def cpscf_nseg(self, value="show"):
        '''
Name: CPSCF_NSEG
Type: INTEGER
Default: 0
Options: Range from 0 to 20

Description: Controls the number of segments used to calculate the CPSCF equations.
Recommendation: : Use default unless too much memory is requested.  Increasing this option reduces memory requirements.    '''
        if value == "":
            if "CPSCF_NSEG" in self.dict_of_keywords:
                del self.dict_of_keywords["CPSCF_NSEG"]
                print("Keyword removed.")
        elif value == "show":
            if "CPSCF_NSEG" in self.dict_of_keywords:
                return self.dict_of_keywords["CPSCF_NSEG"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CPSCF_NSEG"] = value.lower()

    def deuterate(self, value="show"):
        '''
Name: DEUTERATE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Requests that all hydrogen atoms be replaces with deuterium.
Recommendation: : Replacing hydrogen atoms reduces the fastest vibrational frequencies by a factor of 1.4, which allow for a larger fictitious mass and time step in ELMD calculations. There is no reason to replace hydrogens in BOMD calculations.    '''
        if value == "":
            if "DEUTERATE" in self.dict_of_keywords:
                del self.dict_of_keywords["DEUTERATE"]
                print("Keyword removed.")
        elif value == "show":
            if "DEUTERATE" in self.dict_of_keywords:
                return self.dict_of_keywords["DEUTERATE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DEUTERATE"] = value.lower()

    def dma_midpoints(self, value="show"):
        '''
Name: DMA_MIDPOINTS
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies whether to include bond midpoints in the DMA expansion.
    '''
        if value == "":
            if "DMA_MIDPOINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["DMA_MIDPOINTS"]
                print("Keyword removed.")
        elif value == "show":
            if "DMA_MIDPOINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["DMA_MIDPOINTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DMA_MIDPOINTS"] = value.lower()

    def dual_basis_energy(self, value="show"):
        '''
Name: DUAL_BASIS_ENERGY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activates dual-basis SCF (HF or DFT) energy correction. 
Recommendation: : Use Dual-Basis to capture large-basis effects at smaller basis cost. Particularly useful with RI-MP2, in which HF often dominates. Use only proper subsets for small-basis calculation.    '''
        if value == "":
            if "DUAL_BASIS_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["DUAL_BASIS_ENERGY"]
                print("Keyword removed.")
        elif value == "show":
            if "DUAL_BASIS_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["DUAL_BASIS_ENERGY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DUAL_BASIS_ENERGY"] = value.lower()

    def eda_bsse(self, value="show"):
        '''
Name: EDA_BSSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Calculates the BSSE correction when performing the energy decomposition analysis.
Recommendation: : Set to TRUE unless a very large basis set is used.    '''
        if value == "":
            if "EDA_BSSE" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_BSSE"]
                print("Keyword removed.")
        elif value == "show":
            if "EDA_BSSE" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_BSSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EDA_BSSE"] = value.lower()

    def eda_print_covp(self, value="show"):
        '''
Name: EDA_PRINT_COVP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Replace the final MOs with the CVOP orbitals in the end of the run.
Recommendation: : Set to TRUE to print COVP orbitals instead of conventional MOs.    '''
        if value == "":
            if "EDA_PRINT_COVP" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_PRINT_COVP"]
                print("Keyword removed.")
        elif value == "show":
            if "EDA_PRINT_COVP" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_PRINT_COVP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EDA_PRINT_COVP"] = value.lower()

    def epao_iterate(self, value="show"):
        '''
Name: EPAO_ITERATE
Type: INTEGER
Factor: 10
Default: 0 [=0]
Options: Range from 0 [=0] to 100 [=1000]

Description: Controls iterations for EPAO calculations (see PAO_METHOD).
Recommendation: : Use default. For molecules that are not too large, one can test the sensitivity of the results to the type of minimal functions by the use of optimized EPAOs in which case a value of 500 is reasonable.    '''
        if value == "":
            if "EPAO_ITERATE" in self.dict_of_keywords:
                del self.dict_of_keywords["EPAO_ITERATE"]
                print("Keyword removed.")
        elif value == "show":
            if "EPAO_ITERATE" in self.dict_of_keywords:
                return self.dict_of_keywords["EPAO_ITERATE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EPAO_ITERATE"] = value.lower()

    def fast_xc(self, value="show"):
        '''
Name: FAST_XC
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls direct variable thresholds to accelerate the calculation of exchange and correlation (XC) in DFT.
Recommendation: : This option improves the speed of a DFT calculation, but may occasionally cause the SCF calculation to diverge.    '''
        if value == "":
            if "FAST_XC" in self.dict_of_keywords:
                del self.dict_of_keywords["FAST_XC"]
                print("Keyword removed.")
        elif value == "show":
            if "FAST_XC" in self.dict_of_keywords:
                return self.dict_of_keywords["FAST_XC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FAST_XC"] = value.lower()

    def frgm_lpcorr(self, value="show"):
        '''
Name: FRGM_LPCORR
Type: STRING
Default: None

Options:
    'None'.......................... None
    'ARS'........................... ARS
    'RS'............................ RS
    'EXACT_SCF'..................... EXACT_SCF
    'ARS_EXACT_SCF'................. ARS_EXACT_SCF
    'RS_EXACT_SCF'.................. RS_EXACT_SCF

Description: Specifies a correction method performed after the locally-projected equations are converged.
Recommendation: : For large basis sets use ARS, use RS if ARS fails.    '''
        if value == "":
            if "FRGM_LPCORR" in self.dict_of_keywords:
                del self.dict_of_keywords["FRGM_LPCORR"]
                print("Keyword removed.")
        elif value == "show":
            if "FRGM_LPCORR" in self.dict_of_keywords:
                return self.dict_of_keywords["FRGM_LPCORR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FRGM_LPCORR"] = value.lower()

    def frgm_method(self, value="show"):
        '''
Name: FRGM_METHOD
Type: STRING
Default: None

Options:
    'None'.......................... None
    'STOLL'......................... STOLL
    'GIA'........................... GIA
    'NOSCF_RS'...................... NOSCF_RS
    'NOSCF_ARS'..................... NOSCF_ARS
    'NOSCF_DRS'..................... NOSCF_DRS
    'NOSCF_RS_FOCK'................. NOSCF_RS_FOCK

Description: Specifies the locally-projected method.
Recommendation: : STOLL and GIA - variational optimization of the ALMOs. NOSCF options are for computationally fast corrections of the FRAGMO initial guess.     '''
        if value == "":
            if "FRGM_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["FRGM_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "FRGM_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["FRGM_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FRGM_METHOD"] = value.lower()

    def ftc_class_thresh_mult(self, value="show"):
        '''
Name: FTC_CLASS_THRESH_MULT
Type: INTEGER
Default: 5
Options: Range from 1 to 9

Description: Together with FTC_CLASS_THRESH_ORDER, determines the cutoff threshold for included a shell-pair in the dd class, i.e. the class that is expanded in terms of plane waves. 
    '''
        if value == "":
            if "FTC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_CLASS_THRESH_MULT"]
                print("Keyword removed.")
        elif value == "show":
            if "FTC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_CLASS_THRESH_MULT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FTC_CLASS_THRESH_MULT"] = value.lower()

    def ftc_class_thresh_order(self, value="show"):
        '''
Name: FTC_CLASS_THRESH_ORDER
Type: INTEGER
Default: 5
Options: Range from 1 to 9

Description: Together with FTC_CLASS_THRESH_MULT, determines the cutoff threshold for included a shell-pair in the dd class, i.e. the class that is expanded in terms of plane waves.
    '''
        if value == "":
            if "FTC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "FTC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"] = value.lower()

    def qui_geom_opt_fallback(self, value="show"):
        '''
Name: QUI_GEOM_OPT_FALLBACK
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Sets whether or not to fall back to cartesian coordinates if the optimization in internal or Z-matrix coordinates fails.
    '''
        if value == "":
            if "QUI_GEOM_OPT_FALLBACK" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_GEOM_OPT_FALLBACK" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"] = value.lower()

    def geom_opt_linear_angle(self, value="show"):
        '''
Name: GEOM_OPT_LINEAR_ANGLE
Type: INTEGER
Default: 165
Options: Range from 150 to 180

Description: Threshold for near linear bond angles (in degrees).
Recommendation: : Use default.    '''
        if value == "":
            if "GEOM_OPT_LINEAR_ANGLE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_LINEAR_ANGLE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"] = value.lower()

    def geom_opt_max_cycles(self, value="show"):
        '''
Name: GEOM_OPT_MAX_CYCLES
Type: INTEGER
Default: 50
Options: Range from 1 to 200

Description: Maximum number of optimization cycles.
Recommendation: : The default should be sufficient for most cases. Increase if the initial guess geometry is poor, or for systems with shallow potential wells.    '''
        if value == "":
            if "GEOM_OPT_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"] = value.lower()

    def geom_opt_mode(self, value="show"):
        '''
Name: GEOM_OPT_MODE
Type: INTEGER
Default: 0
Options: Range from 0 to 200

Description: Determines which Hessian mode is followed during a transition state search.
    '''
        if value == "":
            if "GEOM_OPT_MODE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MODE"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_MODE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MODE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_MODE"] = value.lower()

    def geom_opt_print(self, value="show"):
        '''
Name: GEOM_OPT_PRINT
Type: INTEGER
Default: 3
Options: Range from 0 to 7

Description: Controls the amount of optimization print output.
Recommendation: : Use the default.    '''
        if value == "":
            if "GEOM_OPT_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_PRINT"] = value.lower()

    def geom_print(self, value="show"):
        '''
Name: GEOM_PRINT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the printing of additional geometric information at each step.
Recommendation: : Use if you want to be able to quickly examine geometric parameters at the beginning and end of optimizations. Only prints in the beginning of single point energy calculations.    '''
        if value == "":
            if "GEOM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_PRINT"] = value.lower()

    def gvb_amp_scale(self, value="show"):
        '''
Name: GVB_AMP_SCALE
Type: INTEGER
Factor: 0.001
Default: 1000 [=1.000]
Options: Range from 1 [=0.001] to 1000 [=1.000]

Description: Scales the default orbital amplitude iteration step size for IP/RCC. PP amplitude equations are solved analytically, so this parameter does not affect PP.
Recommendation: : Default is usually fine, but in some highly-correlated systems it can help with convergence to use smaller values.    '''
        if value == "":
            if "GVB_AMP_SCALE" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_AMP_SCALE"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_AMP_SCALE" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_AMP_SCALE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_AMP_SCALE"] = value.lower()

    def gvb_guess_mix(self, value="show"):
        '''
Name: GVB_GUESS_MIX
Type: INTEGER
Default: 0
Options: Range from 0 to 100

Description: Similar to SCF_GUESS_MIX, it breaks alpha-beta symmetry for UPP by mixing the alpha HOMO and LUMO orbitals according to the user-defined fraction of LUMO to add the HOMO. 100 corresponds to a 1:1 ratio of HOMO and LUMO in the mixed orbitals.
Recommendation: : 25 often works well to break symmetry without overly impeding convergence.    '''
        if value == "":
            if "GVB_GUESS_MIX" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_GUESS_MIX"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_GUESS_MIX" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_GUESS_MIX"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_GUESS_MIX"] = value.lower()

    def gvb_local(self, value="show"):
        '''
Name: GVB_LOCAL
Type: STRING
Default: Pipek-Mezey

Options:
    'Boys localized'................ Boys localized
    'Pipek-Mezey'................... Pipek-Mezey

Description: Sets the localization scheme used in the initial guess wave function.
Recommendation: : Different initial guesses can sometimes lead to different solutions. It can be helpful to try both to ensure the global minimum has been found.    '''
        if value == "":
            if "GVB_LOCAL" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_LOCAL"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_LOCAL" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_LOCAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_LOCAL"] = value.lower()

    def gvb_n_pairs(self, value="show"):
        '''
Name: GVB_N_PAIRS
Type: INTEGER
Default: 0
Options: Range from 0 to 100

Description: Alternative to CC_REST_OCC and CC_REST_VIR for setting active space size in GVB and valence coupled cluster methods.
Recommendation: : Use default unless one wants to study a special active space. When using small active spaces, it is important to ensure that the proper orbitals are incorporated in the active space. If not, use the $reordermo feature to adjust the SCF orbitals appropriately.    '''
        if value == "":
            if "GVB_N_PAIRS" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_N_PAIRS"]
                print("Keyword removed.")
        elif value == "show":
            if "GVB_N_PAIRS" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_N_PAIRS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GVB_N_PAIRS"] = value.lower()

    def incdft(self, value="show"):
        '''
Name: INCDFT
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Toggles the use of the IncDFT procedure for DFT energy calculations.
Recommendation: : Turning this option on can lead to faster SCF calculations, particularly towards the end of the SCF. Please note that for some systems use of this option may lead to convergence problems.    '''
        if value == "":
            if "INCDFT" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT"]
                print("Keyword removed.")
        elif value == "show":
            if "INCDFT" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCDFT"] = value.lower()

    def incdft_dendiff_thresh(self, value="show"):
        '''
Name: INCDFT_DENDIFF_THRESH
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: Sets the threshold for screening density matrix values in the IncDFT procedure. 
Recommendation: : If the default value causes convergence problems, set this value higher to tighten the threshold.    '''
        if value == "":
            if "INCDFT_DENDIFF_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_DENDIFF_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "INCDFT_DENDIFF_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_DENDIFF_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCDFT_DENDIFF_THRESH"] = value.lower()

    def incdft_dendiff_varthresh(self, value="show"):
        '''
Name: INCDFT_DENDIFF_VARTHRESH
Type: INTEGER
Default: 0
Options: Range from 0 to 12

Description: Sets the lower bound for the variable threshold for screening density matrix values in the IncDFT procedure. The threshold will begin at this value and then vary depending on the error in the current SCF iteration until the value specified by INCDFT_DENDIFF_THRESH is reached. This means this value must be set lower than INCDFT_DENDIFF_THRESH.
Recommendation: : If the default value causes convergence problems, set this value higher to tighten accuracy. If this fails, set to 0 and use a static threshold.    '''
        if value == "":
            if "INCDFT_DENDIFF_VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "INCDFT_DENDIFF_VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"] = value.lower()

    def incdft_griddiff_thresh(self, value="show"):
        '''
Name: INCDFT_GRIDDIFF_THRESH
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: Sets the threshold for screening functional values in the IncDFT procedure
Recommendation: : If the default value causes convergence problems, set this value higher to tighten the threshold.    '''
        if value == "":
            if "INCDFT_GRIDDIFF_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "INCDFT_GRIDDIFF_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"] = value.lower()

    def incfock(self, value="show"):
        '''
Name: INCFOCK
Type: INTEGER
Default: 1
Options: Range from 0 to 100

Description: Iteration number after which the incremental Fock matrix algorithm is initiated.  Setting this to 0 turns INCFOCK off.
Recommendation: : May be necessary to allow several iterations before switching on INCFOCK.    '''
        if value == "":
            if "INCFOCK" in self.dict_of_keywords:
                del self.dict_of_keywords["INCFOCK"]
                print("Keyword removed.")
        elif value == "show":
            if "INCFOCK" in self.dict_of_keywords:
                return self.dict_of_keywords["INCFOCK"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCFOCK"] = value.lower()

    def integrals_buffer(self, value="show"):
        '''
Name: INTEGRALS_BUFFER
Type: INTEGER
Default: 15
Options: Range from 1 to 128

Description: Controls the size of in-core integral storage buffer (in megabytes).
Recommendation: : Use the default, or consult your systems administrator for hardware limits.    '''
        if value == "":
            if "INTEGRALS_BUFFER" in self.dict_of_keywords:
                del self.dict_of_keywords["INTEGRALS_BUFFER"]
                print("Keyword removed.")
        elif value == "show":
            if "INTEGRALS_BUFFER" in self.dict_of_keywords:
                return self.dict_of_keywords["INTEGRALS_BUFFER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTEGRALS_BUFFER"] = value.lower()

    def lin_k(self, value="show"):
        '''
Name: LIN_K
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether linear scaling evaluation of exact exchange (LinK) is used.
Recommendation: : Use for HF and hybrid DFT calculations with large numbers of atoms.    '''
        if value == "":
            if "LIN_K" in self.dict_of_keywords:
                del self.dict_of_keywords["LIN_K"]
                print("Keyword removed.")
        elif value == "show":
            if "LIN_K" in self.dict_of_keywords:
                return self.dict_of_keywords["LIN_K"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["LIN_K"] = value.lower()

    def max_sub_file_num(self, value="show"):
        '''
Name: MAX_SUB_FILE_NUM
Type: INTEGER
Default: 16
Options: Range from 1 to 64

Description: Sets the maximum number of sub files allowed.
Recommendation: : Leave as default, or adjust according to your system limits.    '''
        if value == "":
            if "MAX_SUB_FILE_NUM" in self.dict_of_keywords:
                del self.dict_of_keywords["MAX_SUB_FILE_NUM"]
                print("Keyword removed.")
        elif value == "show":
            if "MAX_SUB_FILE_NUM" in self.dict_of_keywords:
                return self.dict_of_keywords["MAX_SUB_FILE_NUM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MAX_SUB_FILE_NUM"] = value.lower()

    def qui_frozen_core(self, value="show"):
        '''
Name: QUI_FROZEN_CORE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_FROZEN_CORE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_FROZEN_CORE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_FROZEN_CORE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_FROZEN_CORE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_FROZEN_CORE"] = value.lower()

    def xc_smart_grid(self, value="show"):
        '''
Name: XC_SMART_GRID
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Uses SG-0 (where available) for early SCF cycles, and switches to the (larger) grid specified by XC_GRID (which defaults to SG-1, if not otherwise specified) for final cycles of the SCF.
Recommendation: : The use of the smart grid can save some time on initial SCF cycles.    '''
        if value == "":
            if "XC_SMART_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["XC_SMART_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "XC_SMART_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["XC_SMART_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XC_SMART_GRID"] = value.lower()

    def xopt_seam_only(self, value="show"):
        '''
Name: XOPT_SEAM_ONLY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Orders an intersection seam search only, no minimization is to perform.
Recommendation: : In systems with a large number of degrees of freedom it might be useful to locate the seam first setting this option to TRUE and use that geometry as a starting point for the minimization.    '''
        if value == "":
            if "XOPT_SEAM_ONLY" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_SEAM_ONLY"]
                print("Keyword removed.")
        elif value == "show":
            if "XOPT_SEAM_ONLY" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_SEAM_ONLY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XOPT_SEAM_ONLY"] = value.lower()

    def xcis(self, value="show"):
        '''
Name: XCIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Do an XCIS calculation in addition to a CIS calculation
    '''
        if value == "":
            if "XCIS" in self.dict_of_keywords:
                del self.dict_of_keywords["XCIS"]
                print("Keyword removed.")
        elif value == "show":
            if "XCIS" in self.dict_of_keywords:
                return self.dict_of_keywords["XCIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XCIS"] = value.lower()

    def wavefunction_analysis(self, value="show"):
        '''
Name: WAVEFUNCTION_ANALYSIS
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the running of the default wavefunction analysis tasks.
    '''
        if value == "":
            if "WAVEFUNCTION_ANALYSIS" in self.dict_of_keywords:
                del self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"]
                print("Keyword removed.")
        elif value == "show":
            if "WAVEFUNCTION_ANALYSIS" in self.dict_of_keywords:
                return self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"] = value.lower()

    def vibman_print(self, value="show"):
        '''
Name: VIBMAN_PRINT
Type: INTEGER
Default: 1
Options: Range from 1 to 7

Description: Controls level of extra print out for vibrational analysis.
Recommendation: : Use default.    '''
        if value == "":
            if "VIBMAN_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["VIBMAN_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "VIBMAN_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["VIBMAN_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["VIBMAN_PRINT"] = value.lower()

    def vci(self, value="show"):
        '''
Name: VCI
Type: INTEGER
Default: 0
Options: Range from 0 to 10

Description: Specifies the number of quanta involved in the VCI calculation.
Recommendation: : The availability depends on the memory of the machine.  For example, a machine with 1.5 GB memory and for molecules with fewer than 4 atoms, VCI(10) can be carried out, for molecule containing fewer than 5 atoms, VCI(6) can be carried out, for molecule containing fewer than 6 atoms, VCI(5) can be carried out. For molecules containing fewer than 50 atoms, VCI(2) is available. VCI(1) and VCI(3) usually overestimated the true energy while VCI(4) usually gives an answer close to the converged energy.    '''
        if value == "":
            if "VCI" in self.dict_of_keywords:
                del self.dict_of_keywords["VCI"]
                print("Keyword removed.")
        elif value == "show":
            if "VCI" in self.dict_of_keywords:
                return self.dict_of_keywords["VCI"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["VCI"] = value.lower()

    def stability_analysis(self, value="show"):
        '''
Name: STABILITY_ANALYSIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Performs stability analysis for a HF or DFT solution.
Recommendation: : Set to TRUE when a HF or DFT solution is suspected to be unstable.    '''
        if value == "":
            if "STABILITY_ANALYSIS" in self.dict_of_keywords:
                del self.dict_of_keywords["STABILITY_ANALYSIS"]
                print("Keyword removed.")
        elif value == "show":
            if "STABILITY_ANALYSIS" in self.dict_of_keywords:
                return self.dict_of_keywords["STABILITY_ANALYSIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["STABILITY_ANALYSIS"] = value.lower()

    def scf_print_frgm(self, value="show"):
        '''
Name: SCF_PRINT_FRGM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the output of Q-Chem jobs on isolated fragments.
Recommendation: : Use TRUE if details about isolated fragments are important.    '''
        if value == "":
            if "SCF_PRINT_FRGM" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_PRINT_FRGM"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_PRINT_FRGM" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_PRINT_FRGM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_PRINT_FRGM"] = value.lower()

    def chemsol_print(self, value="show"):
        '''
Name: CHEMSOL_PRINT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Increases the amount of ChemSol output.
    '''
        if value == "":
            if "CHEMSOL_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "CHEMSOL_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHEMSOL_PRINT"] = value.lower()

    def qui_charge(self, value="show"):
        '''
Name: QUI_CHARGE
Type: INTEGER
Default: 0
Options: Range from -100 to 100

Description: Sets the total charge of the system.
    '''
        if value == "":
            if "QUI_CHARGE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_CHARGE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_CHARGE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_CHARGE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_CHARGE"] = value.lower()

    def unrestricted(self, value="show"):
        '''
Name: UNRESTRICTED
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of restricted or unrestricted orbitals.
Recommendation: : Use default unless ROHF is desired. Note that for unrestricted calculations on systems with an even number of electrons it is usually necessary to break alpha-beta symmetry in the initial guess, by using SCF_GUESS_MIX or providing $occupied information.    '''
        if value == "":
            if "UNRESTRICTED" in self.dict_of_keywords:
                del self.dict_of_keywords["UNRESTRICTED"]
                print("Keyword removed.")
        elif value == "show":
            if "UNRESTRICTED" in self.dict_of_keywords:
                return self.dict_of_keywords["UNRESTRICTED"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["UNRESTRICTED"] = value.lower()

    def write_wfn(self, value="show"):
        '''
Name: WRITE_WFN
Type: STRING
Default: qchem

Options:
    'qchem'......................... qchem

Description: Specifies whether or not a wfn file is created, which is suitable for use with AIMPAC. Note that the output to this file is currently limited to f orbitals, which is the highest angular momentum implemented in AIMPAC.
    '''
        if value == "":
            if "WRITE_WFN" in self.dict_of_keywords:
                del self.dict_of_keywords["WRITE_WFN"]
                print("Keyword removed.")
        elif value == "show":
            if "WRITE_WFN" in self.dict_of_keywords:
                return self.dict_of_keywords["WRITE_WFN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["WRITE_WFN"] = value.lower()

    def aimd_moments(self, value="show"):
        '''
Name: AIMD_MOMENTS
Type: INTEGER
Default: 0
Options: Range from 0 to 20

Description: Specifies the order of multipole moments that are output at each time step.  Setting this to 0 disables the printing of moments.
    '''
        if value == "":
            if "AIMD_MOMENTS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_MOMENTS"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_MOMENTS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_MOMENTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_MOMENTS"] = value.lower()

    def n_frozen_core(self, value="show"):
        '''
Name: N_FROZEN_CORE
Type: INTEGER
Default: 0
Options: Range from 0 to 100

Description: Sets the number of frozen core orbitals in a post-Hartree-Fock calculation.
Recommendation: : While the default is not to freeze orbitals, MP2 calculations are more efficient with frozen core orbitals. Use FC if possible.    '''
        if value == "":
            if "N_FROZEN_CORE" in self.dict_of_keywords:
                del self.dict_of_keywords["N_FROZEN_CORE"]
                print("Keyword removed.")
        elif value == "show":
            if "N_FROZEN_CORE" in self.dict_of_keywords:
                return self.dict_of_keywords["N_FROZEN_CORE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["N_FROZEN_CORE"] = value.lower()

    def n_frozen_virtual(self, value="show"):
        '''
Name: N_FROZEN_VIRTUAL
Type: INTEGER
Default: 0
Options: Range from 0 to 500

Description: Sets the number of frozen virtual orbitals in a post-Hartree-Fock calculation.
    '''
        if value == "":
            if "N_FROZEN_VIRTUAL" in self.dict_of_keywords:
                del self.dict_of_keywords["N_FROZEN_VIRTUAL"]
                print("Keyword removed.")
        elif value == "show":
            if "N_FROZEN_VIRTUAL" in self.dict_of_keywords:
                return self.dict_of_keywords["N_FROZEN_VIRTUAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["N_FROZEN_VIRTUAL"] = value.lower()

    def chemsol(self, value="show"):
        '''
Name: CHEMSOL
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of ChemSol in Q-Chem.
    '''
        if value == "":
            if "CHEMSOL" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL"]
                print("Keyword removed.")
        elif value == "show":
            if "CHEMSOL" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHEMSOL"] = value.lower()

    def ftc(self, value="show"):
        '''
Name: FTC
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the overall use of the FTC.
Recommendation: : Use FTC when bigger and/or diffuse basis sets are used.     '''
        if value == "":
            if "FTC" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC"]
                print("Keyword removed.")
        elif value == "show":
            if "FTC" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FTC"] = value.lower()

    def qui_cfmm(self, value="show"):
        '''
Name: QUI_CFMM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_CFMM" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_CFMM"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_CFMM" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_CFMM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_CFMM"] = value.lower()

    def qui_largemol_none(self, value="show"):
        '''
Name: QUI_LARGEMOL_NONE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_LARGEMOL_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_LARGEMOL_NONE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_LARGEMOL_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_LARGEMOL_NONE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_LARGEMOL_NONE"] = value.lower()

    def qui_solvent_onsager(self, value="show"):
        '''
Name: QUI_SOLVENT_ONSAGER
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_SOLVENT_ONSAGER" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_ONSAGER"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_ONSAGER" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_ONSAGER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_SOLVENT_ONSAGER"] = value.lower()

    def qui_plots_points(self, value="show"):
        '''
Name: QUI_PLOTS_POINTS
Type: INTEGER
Default: 1
Options: Range from 1 to 10000

Description: 
    '''
        if value == "":
            if "QUI_PLOTS_POINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_PLOTS_POINTS"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_PLOTS_POINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_PLOTS_POINTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_PLOTS_POINTS"] = value.lower()

    def ssg(self, value="show"):
        '''
Name: SSG
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the calculation of the SSG wavefunction.
Recommendation: : See also the UNRESTRICTED and DIIS_SUBSPACE_SIZE $rem variables.    '''
        if value == "":
            if "SSG" in self.dict_of_keywords:
                del self.dict_of_keywords["SSG"]
                print("Keyword removed.")
        elif value == "show":
            if "SSG" in self.dict_of_keywords:
                return self.dict_of_keywords["SSG"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SSG"] = value.lower()

    def qui_eom_method(self, value="show"):
        '''
Name: QUI_EOM_METHOD
Type: STRING
Default: Spin Flip

Options:
    'Spin Flip'..................... Spin Flip
    'EA (Diffuse Orbital)'.......... EA (Diffuse Orbital)
    'IP (Diffuse Orbital)'.......... IP (Diffuse Orbital)
    'IP (Proper)'................... IP (Proper)
    'DIP'........................... DIP

Description: Specifies the type of EOM calculation to perform
    '''
        if value == "":
            if "QUI_EOM_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_METHOD"] = value.lower()

    def geom_opt_tol_displacement(self, value="show"):
        '''
Name: GEOM_OPT_TOL_DISPLACEMENT
Type: INTEGER
Default: 1200
Options: Range from 0 to 5000

Description: Convergence on maximum atomic displacement (in micro-angstroms).
Recommendation: : Use the default. To converge the gradient and either one of the energy or displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_DISPLACEMENT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_TOL_DISPLACEMENT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"] = value.lower()

    def geom_opt_tol_energy(self, value="show"):
        '''
Name: GEOM_OPT_TOL_ENERGY
Type: INTEGER
Default: 100
Options: Range from 1 to 500

Description: Convergence on energy change of successive optimization cycles (x 10-8).
Recommendation: : Use the default. To converge the gradient and either one of the energy or displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_TOL_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"] = value.lower()

    def geom_opt_tol_gradient(self, value="show"):
        '''
Name: GEOM_OPT_TOL_GRADIENT
Type: INTEGER
Default: 300
Options: Range from 1 to 1000

Description: Convergence on maximum gradient component.
Recommendation: : Use the default. To converge the gradient and either one of the energy and displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_GRADIENT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_TOL_GRADIENT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"] = value.lower()

    def skip_cis_rpa(self, value="show"):
        '''
Name: SKIP_CIS_RPA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Skips the solution of the CIS, RPA, TDA or TDDFT equations for wavefunction analysis.

Recommendation: : Set to true to speed up the generation of plot data if the same calculation has ben run previously.    '''
        if value == "":
            if "SKIP_CIS_RPA" in self.dict_of_keywords:
                del self.dict_of_keywords["SKIP_CIS_RPA"]
                print("Keyword removed.")
        elif value == "show":
            if "SKIP_CIS_RPA" in self.dict_of_keywords:
                return self.dict_of_keywords["SKIP_CIS_RPA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SKIP_CIS_RPA"] = value.lower()

    def qui_coordinates(self, value="show"):
        '''
Name: QUI_COORDINATES
Type: STRING
Default: Cartesian

Options:
    'Cartesian'..................... Cartesian
    'Z-matrix'...................... Z-matrix
    'Z-matrix (compact)'............ Z-matrix (compact)

Description: Controls the format of the geometry in the output file.
    '''
        if value == "":
            if "QUI_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_COORDINATES"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_COORDINATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_COORDINATES"] = value.lower()

    def pseudo_canonical(self, value="show"):
        '''
Name: PSEUDO_CANONICAL
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: When SCF_ALGORITHM = DM, this controls the way the initial step, and steps after subspace resets are taken.
Recommendation: : The default is usually more efficient, but choosing TRUE sometimes avoids problems with orbital reordering.    '''
        if value == "":
            if "PSEUDO_CANONICAL" in self.dict_of_keywords:
                del self.dict_of_keywords["PSEUDO_CANONICAL"]
                print("Keyword removed.")
        elif value == "show":
            if "PSEUDO_CANONICAL" in self.dict_of_keywords:
                return self.dict_of_keywords["PSEUDO_CANONICAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PSEUDO_CANONICAL"] = value.lower()

    def purecart(self, value="show"):
        '''
Name: PURECART
Type: STRING
Default: 1111

Options:
    '1111'.......................... All pure
    '2222'.......................... All cartesian
    '2111'.......................... Pure d, f and g
    '2211'.......................... Pure d and f
    '2221'.......................... Pure d only

Description: Controls the use of pure (spherical harmonic) or Cartesian angular forms.
    '''
        if value == "":
            if "PURECART" in self.dict_of_keywords:
                del self.dict_of_keywords["PURECART"]
                print("Keyword removed.")
        elif value == "show":
            if "PURECART" in self.dict_of_keywords:
                return self.dict_of_keywords["PURECART"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PURECART"] = value.lower()

    def moprop(self, value="show"):
        '''
Name: MOPROP
Type: STRING
Default: 0

Options:
    '0'............................. None
    '1'............................. NMR Shielding Tensors
    '2'............................. Static Polarizability
    '100'........................... Dynamic Polarizability
    '101'........................... Hyperpolarizability
    '102'........................... Hyperpolarizability (Read)
    '103'........................... Hyperpolarizabilty (Wigner)
    '104'........................... Hyperpolarizability (Wigner+Read)

Description: Specifies the job for mopropman.  Note that for hyperpolarizabilities, the Wigner option uses the (2n+1) rule, and the Read option takes first order results from disk.
    '''
        if value == "":
            if "MOPROP" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP"] = value.lower()

    def ecp(self, value="show"):
        '''
Name: ECP
Type: STRING
Default: None

Options:
    'None'.......................... None

    'CRENBL'........................ CRENBL
    'CRENBS'........................ CRENBS
    'HWMB'.......................... HWMB
    'HWVDZ'......................... HWVDZ
    'LACVP'......................... LACVP
    'LANL2DZ'....................... LANL2DZ
    'SBKJC'......................... SBKJC
    'SRLC'.......................... SRLC
    'SRSC'.......................... SRSC

    'gen'........................... User-defined

Description: Defines the effective core potential and associated basis set to be used
Recommendation: : Pseudopotentials are recommended for first row transition metals and heavier elements. Consult the reviews for more details.    '''
        if value == "":
            if "ECP" in self.dict_of_keywords:
                del self.dict_of_keywords["ECP"]
                print("Keyword removed.")
        elif value == "show":
            if "ECP" in self.dict_of_keywords:
                return self.dict_of_keywords["ECP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ECP"] = value.lower()

    def aimd_fictitious_mass(self, value="show"):
        '''
Name: AIMD_FICTITIOUS_MASS
Type: INTEGER
Default: 100
Options: Range from 1 to 500

Description: Specifies the value of the fictious electronic mass is atomic units.  This quantity has dimensions (energy)x(time)2. 

Recommendation: : Values in the range 50-200 a.u. have been employed.    '''
        if value == "":
            if "AIMD_FICTITIOUS_MASS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FICTITIOUS_MASS"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_FICTITIOUS_MASS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FICTITIOUS_MASS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_FICTITIOUS_MASS"] = value.lower()

    def anhar(self, value="show"):
        '''
Name: ANHAR
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether or not to perform various nuclear vibrational theory (TOSH, VPT2, VCI) calculations to obtain vibrational anharmonic frequencies. 

Recommendation: : Since this calculation involves the third and fourth derivatives at the minimum of the potential energy surface, it is recommended that the geometry optimization tolerances (displacement, gradient and energy) be set tighter.  Note that VPT2 calculations may fail if the system involves accidental degenerate resonances. See the VCI $rem variable for more details about increasing the accuracy of anharmonic calculations.    '''
        if value == "":
            if "ANHAR" in self.dict_of_keywords:
                del self.dict_of_keywords["ANHAR"]
                print("Keyword removed.")
        elif value == "show":
            if "ANHAR" in self.dict_of_keywords:
                return self.dict_of_keywords["ANHAR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ANHAR"] = value.lower()

    def aimd_initial_velocities(self, value="show"):
        '''
Name: AIMD_INITIAL_VELOCITIES
Type: STRING
Default: 0

Options:
    '0'............................. Read
    'ZPE'........................... ZPE
    'Random'........................ Random
    'Thermal'....................... Thermal

Description: Specifies the method for selecting initial nuclear velocities. 
Recommendation: : This variable need only be specified in the event that velocities are not specified explicitly in a $velocity  section.    '''
        if value == "":
            if "AIMD_INITIAL_VELOCITIES" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_INITIAL_VELOCITIES" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"] = value.lower()

    def aimd_temperature(self, value="show"):
        '''
Name: AIMD_TEMPERATURE
Type: INTEGER
Default: 0
Options: Range from 0 to 1000

Description: Specifies a temperature (in Kelvin) for Maxwell-Boltzmann velocity sampling.
Recommendation: : This variable is only useful in conjunction with AIMD velocit initialisation set to Thermal. Note that the simulations are run at constant energy, rather than constant temperature, so the mean nuclear kinetic energy will fluctuate in the course of the simulation.    '''
        if value == "":
            if "AIMD_TEMP" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_TEMP"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_TEMP" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_TEMP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_TEMP"] = value.lower()

    def anharmonic(self, value="show"):
        '''
Name: ANHARMONIC
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether or not to perform various nuclear vibrational theory (TOSH, VPT2, VCI) calculations to obtain vibrational anharmonic frequencies. 

Recommendation: : Since this calculation involves the third and fourth derivatives at the minimum of the potential energy surface, it is recommended that the geometry optimization tolerances (displacement, gradient and energy) be set tighter.  Note that VPT2 calculations may fail if the system involves accidental degenerate resonances. See the VCI $rem variable for more details about increasing the accuracy of anharmonic calculations.    '''
        if value == "":
            if "ANHARMONIC" in self.dict_of_keywords:
                del self.dict_of_keywords["ANHARMONIC"]
                print("Keyword removed.")
        elif value == "show":
            if "ANHARMONIC" in self.dict_of_keywords:
                return self.dict_of_keywords["ANHARMONIC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ANHARMONIC"] = value.lower()

    def basis_linear_dependence_thresh(self, value="show"):
        '''
Name: BASIS_LINEAR_DEPENDENCE_THRESH
Type: INTEGER
Default: 6
Options: Range from 1 to 10

Description: Sets the threshold for determining linear dependence in the basis set.  The threshold is set to 10-n.
Recommendation: : Set to 5 or smaller if you have a poorly behaved SCF and you suspect linear dependence in you basis set. Lower values (larger thresholds) may affect the accuracy of the calculation.    '''
        if value == "":
            if "BASIS_LINEAR_DEPENDENCE_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS_LINEAR_DEPENDENCE_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "BASIS_LINEAR_DEPENDENCE_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS_LINEAR_DEPENDENCE_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "BASIS_LINEAR_DEPENDENCE_THRESH"] = value.lower()

    def cc_amplitude_response(self, value="show"):
        '''
Name: CC_AMPLITUDE_RESPONSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If set to TRUE, adds amplitude response terms to one-particle and two-particle CCSD density matrices before calculation of properties. CC_PROP must be set to TRUE. 
Recommendation: : The cost is always about the cost of an analytic gradient calculation, independent of whether or not the two-particle properties are requested. Besides, adding amplitude response terms without orbital response will unlikely improve the quality of the properties. However, it can be used for debugging purposes.    '''
        if value == "":
            if "CC_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"] = value.lower()

    def cc_properties(self, value="show"):
        '''
Name: CC_PROPERTIES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments < X2>, < Y2>, and < Z2> of electron density, and the total < R2> = < X2> +< Y2> +< Z2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.
Recommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PROPERTIES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PROPERTIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PROPERTIES"] = value.lower()

    def cc_two_particle_properties(self, value="show"):
        '''
Name: CC_TWO_PARTICLE_PROPERTIES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Request for calculation of non-relaxed two-particle CCSD properties. The two-particle properties currently include. The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation is inferior compared to the cost of 2>. The variable CC_PROPERTIES must be also set to TRUE.
Recommendation: : The two-particle properties are extremely computationally expensive, since they require calculation and use of the two-particle density matrix (the cost is approximately the same as the cost of an analytic gradient calculation). Do not request the two-particle properties unless you really need them.    '''
        if value == "":
            if "CC_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"] = value.lower()

    def cc_block_tensor_buffer_size(self, value="show"):
        '''
Name: CC_BLOCK_TENSOR_BUFFER_SIZE
Type: INTEGER
Default: 1000
Options: Range from 0 to 8000

Description: Specifies the maximum size, in Mb, of the buffers for in-core storage of block-tensors.
Recommendation: : Larger values can give better I/O performance and are recommended for systems with large memory (add to your .qchemrc file)    '''
        if value == "":
            if "CC_BLOCK_TENSOR_BUFFER_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_BLOCK_TENSOR_BUFFER_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"] = value.lower()

    def cc_diis_maximum_overlap(self, value="show"):
        '''
Name: CC_DIIS_MAXIMUM_OVERLAP
Type: INTEGER
Factor: 0.01
Default: 100 [=1.00]
Options: Range from 1 [=0.01] to 100 [=1.00]

Description: 
    '''
        if value == "":
            if "CC_DIIS_MAXIMUM_OVERLAP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS_MAXIMUM_OVERLAP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"] = value.lower()

    def cc_diis12_switch(self, value="show"):
        '''
Name: CC_DIIS12_SWITCH
Type: INTEGER
Default: 5
Options: Range from 0 to 12

Description: When to switch from DIIS 2 to DIIS 1 procedure, or when DIIS 2 procedure is required to generate DIIS guesses less frequently. Total value of DIIS error vector must be less than 10-n, where n is the value of this option.
    '''
        if value == "":
            if "CC_DIIS12_SWITCH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS12_SWITCH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS12_SWITCH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS12_SWITCH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS12_SWITCH"] = value.lower()

    def cc_diis_extrapolation_frequency(self, value="show"):
        '''
Name: CC_DIIS_EXTRAPOLATION_FREQUENCY
Type: INTEGER
Default: 2
Options: Range from 0 to 100

Description: DIIS extrapolation will be attempted every n iterations. However, DIIS2 will be attempted every iteration while total error vector exceeds CC_DIIS12_SWITCH. DIIS1 cannot generate guesses more frequently than every 2 iterations. 
    '''
        if value == "":
            if "CC_DIIS_EXTRAPOLATION_FREQUENCY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_EXTRAPOLATION_FREQUENCY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS_EXTRAPOLATION_FREQUENCY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_EXTRAPOLATION_FREQUENCY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "CC_DIIS_EXTRAPOLATION_FREQUENCY"] = value.lower()

    def cc_diis_minimum_overlap(self, value="show"):
        '''
Name: CC_DIIS_MINIMUM_OVERLAP
Type: INTEGER
Default: 11
Options: Range from 0 to 12

Description: The DIIS procedure will be halted when the square root of smallest element of the error overlap matrix is less than 10-n, where n is the value of this option. Small values of the B matrix mean it will become near-singular, making the DIIS equations difficult to solve.
    '''
        if value == "":
            if "CC_DIIS_MINIMUM_OVERLAP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS_MINIMUM_OVERLAP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"] = value.lower()

    def cc_diis_size(self, value="show"):
        '''
Name: CC_DIIS_SIZE
Type: INTEGER
Default: 7
Options: Range from 1 to 50

Description: Specifies the maximum size of the DIIS space.
Recommendation: : Larger values involve larger amounts of disk storage.    '''
        if value == "":
            if "CC_DIIS_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS_SIZE"] = value.lower()

    def cc_diis_start(self, value="show"):
        '''
Name: CC_DIIS_START
Type: INTEGER
Default: 3
Options: Range from 1 to 500

Description: Iteration number when DIIS is turned on. Set to a large number to disable DIIS.
Recommendation: : Occasionally DIIS can cause optimized orbital coupled-cluster calculations to diverge through large orbital changes. If this is seen, DIIS should be disabled.    '''
        if value == "":
            if "CC_DIIS_START" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_START"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS_START" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_START"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS_START"] = value.lower()

    def cc_dmaxiter(self, value="show"):
        '''
Name: CC_DMAXITER
Type: INTEGER
Default: 30
Options: Range from 1 to 100

Description: Maximum number of iteration allowed for Davidson diagonalization procedure. 
Recommendation: : Default is usually sufficient    '''
        if value == "":
            if "CC_DMAXITER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DMAXITER"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DMAXITER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DMAXITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DMAXITER"] = value.lower()

    def cc_do_disconnected(self, value="show"):
        '''
Name: CC_DO_DISCONNECTED
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Determines whether disconnected terms included in the EOM-OD equations
Recommendation: : Inclusion of disconnected terms has very small effects and is not necessary.    '''
        if value == "":
            if "CC_DO_DISCONNECTED" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DISCONNECTED"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DO_DISCONNECTED" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DISCONNECTED"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DO_DISCONNECTED"] = value.lower()

    def cc_do_cisdt(self, value="show"):
        '''
Name: CC_DO_CISDT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the calculation of full CISDT
    '''
        if value == "":
            if "CC_DO_CISDT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_CISDT"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DO_CISDT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_CISDT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DO_CISDT"] = value.lower()

    def cc_do_dyson_ee(self, value="show"):
        '''
Name: CC_DO_DYSON_EE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether excited state Dyson orbitals will be calculated for EOM-IP/EA-CCSD calculations.
Recommendation: : none    '''
        if value == "":
            if "CC_DO_DYSON_EE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DYSON_EE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DO_DYSON_EE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DYSON_EE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DO_DYSON_EE"] = value.lower()

    def cc_do_dyson(self, value="show"):
        '''
Name: CC_DO_DYSON
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether ground state Dyson orbitals will be calculated for EOM-IP/EA-CCSD calculations.
Recommendation: : none    '''
        if value == "":
            if "CC_DO_DYSON" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DYSON"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DO_DYSON" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DYSON"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DO_DYSON"] = value.lower()

    def cc_ea(self, value="show"):
        '''
Name: CC_EA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If TRUE, calculates EOM-EA-CCSD excitation energies and properties using the diffuse orbital trick. A very diffuse orbital is added to the basis set, excitations from which correspond to electron attachment.
    '''
        if value == "":
            if "CC_EA" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EA"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EA" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EA"] = value.lower()

    def cc_ip(self, value="show"):
        '''
Name: CC_IP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If TRUE, calculates EOM-IP-CCSD excitation energies and properties using the diffuse orbital trick. A very diffuse orbital is added to the basis set, excitations to which correspond to ionization. CC_IP and CC_IP_PROPER keywords are mutually exclusive. 
    '''
        if value == "":
            if "CC_IP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_IP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_IP"] = value.lower()

    def cc_ip_filter(self, value="show"):
        '''
Name: CC_IP_FILTER
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If TRUE, filters the EOM-IP-CCSD amplitudes obtained using the diffuse orbital implementation, to be used in conjunction with CC_IP keyword. Helps with convergence.
    '''
        if value == "":
            if "CC_IP_FILTER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP_FILTER"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_IP_FILTER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP_FILTER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_IP_FILTER"] = value.lower()

    def cc_ip_proper(self, value="show"):
        '''
Name: CC_IP_PROPER
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If TRUE, calculates proper EOM-IP-CCSD excitation energies and properties. This implementation does not use the diffuse orbital and is the recommended method of doing EOM-IP-CCSD calculations. CC_IP and CC_IP_PROPER keywords are mutually exclusive.
    '''
        if value == "":
            if "CC_IP_PROPER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP_PROPER"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_IP_PROPER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP_PROPER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_IP_PROPER"] = value.lower()

    def cc_mp2no_grad(self, value="show"):
        '''
Name: CC_MP2NO_GRAD
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If CC_MP2NO_GUESS is TRUE, what kind of one-particle density matrix is used to make the guess orbitals? 
Recommendation: : The two definitions give generally similar performance.    '''
        if value == "":
            if "CC_MP2NO_GRAD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MP2NO_GRAD"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_MP2NO_GRAD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MP2NO_GRAD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_MP2NO_GRAD"] = value.lower()

    def cc_mp2no_guess(self, value="show"):
        '''
Name: CC_MP2NO_GUESS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Will guess orbitals be natural orbitals of the MP1 wavefunction? Alternatively, it is possible to use an effective one-particle density matrix to define the natural orbitals. 
    '''
        if value == "":
            if "CC_MP2NO_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MP2NO_GUESS"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_MP2NO_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MP2NO_GUESS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_MP2NO_GUESS"] = value.lower()

    def cc_preconv_doubles(self, value="show"):
        '''
Name: CC_PRECONV_DOUBLES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: When TRUE, doubly-excited vectors are converged prior to a full excited states calculation. 
Recommendation: : Occasionally necessary to ensure a doubly excited state is found.    '''
        if value == "":
            if "CC_PRECONV_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONV_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONV_DOUBLES"] = value.lower()

    def cc_preconv_singles(self, value="show"):
        '''
Name: CC_PRECONV_SINGLES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: When TRUE, singly-excited vectors are converged prior to a full excited states calculation. 
    '''
        if value == "":
            if "CC_PRECONV_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_SINGLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONV_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_SINGLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONV_SINGLES"] = value.lower()

    def cc_prop(self, value="show"):
        '''
Name: CC_PROP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments < X2>, < Y2>, and < Z2> of electron density, and the total < R2> = < X2> +< Y2> +< Z2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.
Recommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PROP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PROP"] = value.lower()

    def cc_restart(self, value="show"):
        '''
Name: CC_RESTART
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Allows an optimized orbital coupled cluster calculation to begin with an initial guess for the orbital transformation matrix U other than the unit vector. The scratch file from a previous run must be available for the U matrix to be read successfully. 
Recommendation: : Useful for restarting a job that did not converge, if files were saved.    '''
        if value == "":
            if "CC_RESTART" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTART"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_RESTART" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTART"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_RESTART"] = value.lower()

    def cc_restart_no_scf(self, value="show"):
        '''
Name: CC_RESTART_NO_SCF
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Should an optimized orbital coupled cluster calculation begin with optimized orbitals from a previous calculation? When TRUE, molecular orbitals are initially orthogonalized, and CC_PRECONV_T2Z and CC_CANONIZE are set to TRUE while other guess options are set to FALSE.
    '''
        if value == "":
            if "CC_RESTART_NO_SCF" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTART_NO_SCF"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_RESTART_NO_SCF" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTART_NO_SCF"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_RESTART_NO_SCF"] = value.lower()

    def cc_sd_3(self, value="show"):
        '''
Name: CC_SD_3
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: This keyword initializes calculation of non-iterative triples corrections (fT) and (dT) for EE or SF after the EOM-CCSD is done.
    '''
        if value == "":
            if "CC_SD_3" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SD_3"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_SD_3" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SD_3"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_SD_3"] = value.lower()

    def cc_spin_flip(self, value="show"):
        '''
Name: CC_SPIN_FLIP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether do perform a standard excited state calculation, or a spin-flip calculation. Spin multiplicity should be set to 3 for systems with an even number of electrons, and 4 for systems with an odd number of electrons.
    '''
        if value == "":
            if "CC_SPIN_FLIP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SPIN_FLIP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_SPIN_FLIP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SPIN_FLIP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_SPIN_FLIP"] = value.lower()

    def cc_symmetry(self, value="show"):
        '''
Name: CC_SYMMETRY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of symmetry in coupled-cluster calculations
Recommendation: : It is automatically turned off for any finite difference calculations, e.g. second derivatives.    '''
        if value == "":
            if "CC_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SYMMETRY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SYMMETRY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_SYMMETRY"] = value.lower()

    def cc_eom_amplitude_response(self, value="show"):
        '''
Name: CC_EOM_AMPLITUDE_RESPONSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If set to TRUE, adds amplitude response terms to one-particle and two-particle EOM-CCSD density matrices before calculation of properties. CC_EXSTATES_PROP must be set to TRUE.
Recommendation: : The cost is always about the cost of an analytic gradient calculation for each state, independent of whether or not the two-particle properties are requested. Besides, adding amplitude response terms without orbital response will unlikely improve the quality of the properties. However, it can be used for debugging purposes.    '''
        if value == "":
            if "CC_EOM_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"] = value.lower()

    def cc_eom_properties(self, value="show"):
        '''
Name: CC_EOM_PROPERTIES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the non-relaxed (expectation value) one-particle EOM-CCSD target state properties will be calculated. The properties currently include permanent dipole moment, the second moments 2>, 2>, and 2> of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.
Recommendation: : Additional equations (EOM-CCSD equations for the left eigenvectors) need to be solved for properties, approximately doubling the cost of calculation for each irrep. Sometimes the equations for left and right eigenvectors converge to different sets of target states. In this case, the simultaneous iterations of left and right vectors will diverge, and the properties for several or all the target states may be incorrect! The problem can be solved by varying the number of requested states, specified with CC_NLOWSPIN and CC_NHIGHSPIN, or the number of guess vectors (CC_NGUESS_SINGLES). The cost of the one-particle properties calculation itself is low. The one-particle density of an EOM-CCSD target state can be analyzed with NBO package by specifying the state with CC_REFSYM and CC_STATE_DERIV and requesting NBO=TRUE and CC_EXSTATES_PROP=TRUE.    '''
        if value == "":
            if "CC_EOM_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROPERTIES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROPERTIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_PROPERTIES"] = value.lower()

    def cc_eom_full_response(self, value="show"):
        '''
Name: CC_EOM_FULL_RESPONSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If set to TRUE, adds both amplitude and orbital response terms to one- and two-particle EOM-CCSD density matrices before calculation of the properties. CC_EXSTATES_PROP must be set to TRUE. If both CC_EOM_AMPL_RESP=TRUE and CC_EOM_FULL_RESP=TRUE, the CC_EOM_AMPL_RESP=TRUE will be ignored.
Recommendation: : The cost for the full response properties calculation is about the same as the cost of the analytic gradient for each state. Adding full response terms improves quality of calculated properties, but usually it is a small but expensive correction. Use it only if you really need accurate properties.    '''
        if value == "":
            if "CC_EOM_FULL_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_FULL_RESPONSE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_FULL_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_FULL_RESPONSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_FULL_RESPONSE"] = value.lower()

    def cc_canonize_frequency(self, value="show"):
        '''
Name: CC_CANONIZE_FREQUENCY
Type: INTEGER
Default: 50
Options: Range from 1 to 100

Description: The orbitals will be semi-canonicalized every n theta resets. The thetas (orbital rotation angles) are reset every CC_RESET_THETA iterations. The counting of iterations differs for active space (VOD, VQCCD) calculations, where the orbitals are always canonicalized at the first theta-reset.
Recommendation: : Smaller values can be tried in cases that do not converge.    '''
        if value == "":
            if "CC_CANONIZE_FREQUENCY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE_FREQUENCY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CANONIZE_FREQUENCY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE_FREQUENCY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CANONIZE_FREQUENCY"] = value.lower()

    def cc_eom_transition_properties(self, value="show"):
        '''
Name: CC_EOM_TRANSITION_PROPERTIES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_REFSYM and CC_STATE_DERIV must be specified for this state.
Recommendation: : Additional equations (for the left EOM-CCSD eigenvectors plus lambda CCSD equations in case if transition properties between the CCSD reference and EOM-CCSD target states are requested) need to be solved for transition properties, approximately doubling the computational cost. The cost of the transition properties calculation itself is low.    '''
        if value == "":
            if "CC_EOM_TRANSITION_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TRANSITION_PROPERTIES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_TRANSITION_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TRANSITION_PROPERTIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "CC_EOM_TRANSITION_PROPERTIES"] = value.lower()

    def cc_convergence_energy(self, value="show"):
        '''
Name: CC_CONVERGENCE_ENERGY
Type: INTEGER
Default: 10
Options: Range from 0 to 12

Description: Convergence desired on the change in total energy, between iterations.
    '''
        if value == "":
            if "CC_CONVERGENCE_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_ENERGY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CONVERGENCE_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_ENERGY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CONVERGENCE_ENERGY"] = value.lower()

    def cc_convergence_zvector(self, value="show"):
        '''
Name: CC_CONVERGENCE_ZVECTOR
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: Convergence criterion on the RMS difference between successive doubles Z-vector amplitudes [10-n].
Recommendation: : Use Default    '''
        if value == "":
            if "CC_CONVERGENCE_ZVECTOR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CONVERGENCE_ZVECTOR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"] = value.lower()

    def cc_convergence_amplitudes(self, value="show"):
        '''
Name: CC_CONVERGENCE_AMPLITUDES
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: Convergence criterion on the RMS difference between successive sets of coupled-cluster doubles amplitudes [10-n]
Recommendation: : Use default    '''
        if value == "":
            if "CC_CONVERGENCE_AMPLITUDES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_CONVERGENCE_AMPLITUDES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"] = value.lower()

    def cc_full_response(self, value="show"):
        '''
Name: CC_FULL_RESPONSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: If set to TRUE, adds both amplitude and orbital response terms to one- and two-particle CCSD density matrices before calculation of the properties. CC_PROP must be set to TRUE. If both CC_AMPL_RESP=TRUE and CC_FULL_RESP=TRUE, the CC_AMPL_RESP=TRUE will be ignored.
Recommendation: : The cost for the full response properties calculation is about the same as the cost of the analytic gradient. Adding full response terms improves quality of calculated properties, but usually it is a small but expensive correction. Use it only if you really need accurate properties.    '''
        if value == "":
            if "CC_FULL_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_FULL_RESPONSE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_FULL_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_FULL_RESPONSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_FULL_RESPONSE"] = value.lower()

    def cc_hessian_thresh(self, value="show"):
        '''
Name: CC_HESSIAN_THRESH
Type: INTEGER
Factor: 0.001
Default: 10 [=0.010]
Options: Range from 1 [=0.001] to 1000 [=1.000]

Description: Minimum alloed value for the orbital Hessian.  Smaller values are replaced with this constant.
    '''
        if value == "":
            if "CC_HESSIAN_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_HESSIAN_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_HESSIAN_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_HESSIAN_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_HESSIAN_THRESH"] = value.lower()

    def cc_preconverge_doubles(self, value="show"):
        '''
Name: CC_PRECONVERGE_DOUBLES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: When TRUE, doubly-excited vectors are converged prior to a full excited states calculation. 
Recommendation: : Occasionally necessary to ensure a doubly excited state is found.    '''
        if value == "":
            if "CC_PRECONVERGE_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"] = value.lower()

    def cc_preconverge_sd(self, value="show"):
        '''
Name: CC_PRECONVERGE_SD
Type: INTEGER

Description: Solves the EOM-CCSD equations, prints energies, then uses EOM-CCSD vectors as initial vectors in EOM-CC(2,3). Very convenient for calculations using energy additivity schemes.
Recommendation: : Turning this option on is recommended    '''
        if value == "":
            if "CC_PRECONVERGE_SD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_SD"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_SD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_SD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_SD"] = value.lower()

    def cc_preconverge_singles(self, value="show"):
        '''
Name: CC_PRECONVERGE_SINGLES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: When TRUE, singly-excited vectors are converged prior to a full excited states calculation. 
    '''
        if value == "":
            if "CC_PRECONVERGE_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_SINGLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_SINGLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_SINGLES"] = value.lower()

    def qui_solvent_none(self, value="show"):
        '''
Name: QUI_SOLVENT_NONE
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Checking this disables all solvent models.
    '''
        if value == "":
            if "QUI_SOLVENT_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_NONE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_NONE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_SOLVENT_NONE"] = value.lower()

    def chemsol_efield(self, value="show"):
        '''
Name: CHEMSOL_EFIELD
Type: STRING
Default: 1

Options:
    '1'............................. Exact
    '0'............................. Mulliken

Description: Determines how the solute charge distribution is approximated in evaluating the electrostatic field of the solute.  Either the exact solute charge distribution is used, or the charge distribution is approximated by Mulliken atomic charges. 
Recommentation: Mulliken charges are faster, but less rigorous.
    '''
        if value == "":
            if "CHEMSOL_EFIELD" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_EFIELD"]
                print("Keyword removed.")
        elif value == "show":
            if "CHEMSOL_EFIELD" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_EFIELD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHEMSOL_EFIELD"] = value.lower()

    def cis_state_derivative(self, value="show"):
        '''
Name: CIS_STATE_DERIVATIVE
Type: INTEGER
Default: 0
Options: Range from 0 to 200

Description: Sets which CIS state to use for excited state optimizations and vibrational analysis.
Recommendation: : Check to see that the states do no change order during an optimization.    '''
        if value == "":
            if "CIS_STATE_DERIVATIVE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_STATE_DERIVATIVE"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_STATE_DERIVATIVE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_STATE_DERIVATIVE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_STATE_DERIVATIVE"] = value.lower()

    def rpa(self, value="show"):
        '''
Name: RPA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Do an RPA calculation in addition to a CIS calculation
    '''
        if value == "":
            if "RPA" in self.dict_of_keywords:
                del self.dict_of_keywords["RPA"]
                print("Keyword removed.")
        elif value == "show":
            if "RPA" in self.dict_of_keywords:
                return self.dict_of_keywords["RPA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPA"] = value.lower()

    def cis_ras_cutoff_occupied(self, value="show"):
        '''
Name: CIS_RAS_CUTOFF_OCCUPIED
Type: INTEGER
Factor: 0.01
Default: 50 [=0.50]
Options: Range from 0 [=0.00] to 200 [=2.00]

Description: Specifies the occupied orbital cutoff
    '''
        if value == "":
            if "CIS_RAS_CUTOFF_OCCUPIED" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS_CUTOFF_OCCUPIED" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"] = value.lower()

    def cis_ras_cutoff_virtual(self, value="show"):
        '''
Name: CIS_RAS_CUTOFF_VIRTUAL
Type: INTEGER
Factor: 0.01
Default: 0 [=0]
Options: Range from 0 [=0] to 100 [=1.00]

Description: Specifies the virtual orbital cutoff.
    '''
        if value == "":
            if "CIS_RAS_CUTOFF_VIRTUAL" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS_CUTOFF_VIRTUAL" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"] = value.lower()

    def cis_ras(self, value="show"):
        '''
Name: CIS_RAS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether reduced single excitation space is used
    '''
        if value == "":
            if "CIS_RAS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS"] = value.lower()

    def cis_ras_n_solute_atoms(self, value="show"):
        '''
Name: CIS_RAS_N_SOLUTE_ATOMS
Type: INTEGER

Description: Specifies number of atoms or orbitals in solute
    '''
        if value == "":
            if "CIS_RAS_N_SOLUTE_ATOMS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS_N_SOLUTE_ATOMS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"] = value.lower()

    def cis_ras_print(self, value="show"):
        '''
Name: CIS_RAS_PRINT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether or not to print additional output.
    '''
        if value == "":
            if "CIS_RAS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS_PRINT"] = value.lower()

    def cis_ras_type(self, value="show"):
        '''
Name: CIS_RAS_TYPE
Type: STRING
Default: 1

Options:
    '1'............................. Localized Orbitals
    '2'............................. User-defined

Description: Controls how reduced subspace is specified
    '''
        if value == "":
            if "CIS_RAS_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_TYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_RAS_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_TYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_RAS_TYPE"] = value.lower()

    def diis_error_metric(self, value="show"):
        '''
Name: DIIS_ERROR_METRIC
Type: STRING
Default: false

Options:
    'false'......................... Maximum
    'true'.......................... RMS

Description: Changes the DIIS convergence metric from the maximum to the RMS error.
Recommendation: : Use default, the maximum error provides a more reliable criterion.    '''
        if value == "":
            if "DIIS_ERROR_METRIC" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_ERROR_METRIC"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_ERROR_METRIC" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_ERROR_METRIC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_ERROR_METRIC"] = value.lower()

    def dma(self, value="show"):
        '''
Name: DMA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies whether to perform Distributed Multipole Analysis.
    '''
        if value == "":
            if "DMA" in self.dict_of_keywords:
                del self.dict_of_keywords["DMA"]
                print("Keyword removed.")
        elif value == "show":
            if "DMA" in self.dict_of_keywords:
                return self.dict_of_keywords["DMA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DMA"] = value.lower()

    def raman(self, value="show"):
        '''
Name: RAMAN
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls calculation of Raman intensities. Only relevant for a frequency calculation.
    '''
        if value == "":
            if "RAMAN" in self.dict_of_keywords:
                del self.dict_of_keywords["RAMAN"]
                print("Keyword removed.")
        elif value == "show":
            if "RAMAN" in self.dict_of_keywords:
                return self.dict_of_keywords["RAMAN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RAMAN"] = value.lower()

    def dscf_convergence_level_1(self, value="show"):
        '''
Name: DSCF_CONVERGENCE_LEVEL_1
Type: INTEGER
Default: 4
Options: Range from 0 to 10

Description: Sets the convergence criterion for the level-1 iterations. This preconditions the density for the level-2 calculation, and does not include any two-electron integrals. 
Recommendation: : The criterion for level-1 convergence must be less than or equal to the level-2 criterion, otherwise the D-CPSCF will not converge.    '''
        if value == "":
            if "DSCF_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"]
                print("Keyword removed.")
        elif value == "show":
            if "DSCF_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"] = value.lower()

    def dscf_convergence_level_2(self, value="show"):
        '''
Name: DSCF_CONVERGENCE_LEVEL_2
Type: INTEGER
Default: 4
Options: Range from 0 to 10

Description: Sets the convergence criterion for the level-2 iterations.
    '''
        if value == "":
            if "DSCF_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"]
                print("Keyword removed.")
        elif value == "show":
            if "DSCF_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"] = value.lower()

    def dscf_diis_subspace(self, value="show"):
        '''
Name: DSCF_DIIS_SUBSPACE
Type: INTEGER
Default: 11
Options: Range from 0 to 25

Description: Specifies the number of matrices to use in the DIIS extrapolation in the D-CPSCF.
Recommendation: : Use the default.    '''
        if value == "":
            if "DSCF_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_DIIS_SUBSPACE"]
                print("Keyword removed.")
        elif value == "show":
            if "DSCF_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_DIIS_SUBSPACE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSCF_DIIS_SUBSPACE"] = value.lower()

    def dcpscf_pertnum(self, value="show"):
        '''
Name: DCPSCF_PERTNUM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies whether to do the perturbations all together.
    '''
        if value == "":
            if "DCPSCF_PERTNUM" in self.dict_of_keywords:
                del self.dict_of_keywords["DCPSCF_PERTNUM"]
                print("Keyword removed.")
        elif value == "show":
            if "DCPSCF_PERTNUM" in self.dict_of_keywords:
                return self.dict_of_keywords["DCPSCF_PERTNUM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DCPSCF_PERTNUM"] = value.lower()

    def dcpscf_perturbations(self, value="show"):
        '''
Name: DCPSCF_PERTURBATIONS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies whether to do the perturbations all together.
    '''
        if value == "":
            if "DCPSCF_PERTURBATIONS" in self.dict_of_keywords:
                del self.dict_of_keywords["DCPSCF_PERTURBATIONS"]
                print("Keyword removed.")
        elif value == "show":
            if "DCPSCF_PERTURBATIONS" in self.dict_of_keywords:
                return self.dict_of_keywords["DCPSCF_PERTURBATIONS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DCPSCF_PERTURBATIONS"] = value.lower()

    def fd_derivative_type(self, value="show"):
        '''
Name: FD_DERIVATIVE_TYPE
Type: STRING
Default: 0

Options:
    '0'............................. Energies Only
    '1'............................. Gradients Only
    '2'............................. Hessians Only
    '3 '............................ Enegies, Gradients and Hessians

Description: Controls what types of gradient information are used to compute higher derivatives. The default uses a combination of energy, gradient and Hessian information, which makes the force field calculation faster. 
Recommendation: : When the molecule is larger than benzene with small basis set, using only Hessian information may be faster. Note that this option will be set lower if analytic derivatives of the requested order are not available.     '''
        if value == "":
            if "FD_DERIVATIVE_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["FD_DERIVATIVE_TYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "FD_DERIVATIVE_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["FD_DERIVATIVE_TYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FD_DERIVATIVE_TYPE"] = value.lower()

    def fd_step_size(self, value="show"):
        '''
Name: FD_STEP_SIZE
Type: INTEGER
Factor: 0.0001
Default: 10 [=0.0010]
Options: Range from 1 [=0.0001 ] to 100 [=0.0100]

Description: Displacement used for calculating derivatives by finite difference.
Recommendation: : Use default, unless on a very flat potential, in which case a larger value may be required.    '''
        if value == "":
            if "FD_STEP_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["FD_STEP_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "FD_STEP_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["FD_STEP_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FD_STEP_SIZE"] = value.lower()

    def analytic_derivative_order(self, value="show"):
        '''
Name: ANALYTIC_DERIVATIVE_ORDER
Type: INTEGER
Default: 2
Options: Range from 0 to 2

Description: Controls the order of derivatives that are evaluated analytically. The user is not normally required to specify a value, unless numerical derivatives are desired. The derivatives will be evaluated numerically if this option is set lower than the type of job requires.
Recommendation: : Usually set to the maximum possible for efficiency. Note that this option will be set lower if analytic derivatives of the requested order are not available.    '''
        if value == "":
            if "ANALYTIC_DERIVATIVE_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "ANALYTIC_DERIVATIVE_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"] = value.lower()

    def pao_algorithm(self, value="show"):
        '''
Name: PAO_ALGORITHM
Type: STRING
Default: 0

Options:
    '0'............................. Efficient
    '1'............................. Conservative

Description: Algorithm used to optimize polarized atomic orbitals (see PAO_METHOD)
    '''
        if value == "":
            if "PAO_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["PAO_ALGORITHM"]
                print("Keyword removed.")
        elif value == "show":
            if "PAO_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["PAO_ALGORITHM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PAO_ALGORITHM"] = value.lower()

    def pao_method(self, value="show"):
        '''
Name: PAO_METHOD
Type: STRING
Default: EPAO 

Options:
    'EPAO '......................... EPAO 
    'PAO'........................... PAO

Description: Controls evaluation of polarized atomic orbitals (PAOs).
    '''
        if value == "":
            if "PAO_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["PAO_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "PAO_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["PAO_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PAO_METHOD"] = value.lower()

    def epao_weights(self, value="show"):
        '''
Name: EPAO_WEIGHTS
Type: STRING
Default: 115

Options:
    '115'........................... 1st & 2nd order
    '15'............................ 1st order only

Description: Controls algorithm and weights for EPAO calculations (see PAO_METHOD).
Recommendation: : Use default, unless convergence failure is encountered.    '''
        if value == "":
            if "EPAO_WEIGHTS" in self.dict_of_keywords:
                del self.dict_of_keywords["EPAO_WEIGHTS"]
                print("Keyword removed.")
        elif value == "show":
            if "EPAO_WEIGHTS" in self.dict_of_keywords:
                return self.dict_of_keywords["EPAO_WEIGHTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EPAO_WEIGHTS"] = value.lower()

    def aimd_fock_extrapolation_order(self, value="show"):
        '''
Name: AIMD_FOCK_EXTRAPOLATION_ORDER
Type: INTEGER
Default: 0
Options: Range from 0 to 20

Description: Specifies the polynomial order N for Fock matrix extrapolation.
    '''
        if value == "":
            if "AIMD_FOCK_EXTRAPOLATION_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_FOCK_EXTRAPOLATION_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "AIMD_FOCK_EXTRAPOLATION_ORDER"] = value.lower()

    def ftc_fast(self, value="show"):
        '''
Name: FTC_FAST
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether or not the operator is evaluated on a large grid and stored in memory to speed up the calculation.
Recommendation: : Use the default if possible.  Turning this option off conserves some memory, but causes a slow down in speed.    '''
        if value == "":
            if "FTC_FAST" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_FAST"]
                print("Keyword removed.")
        elif value == "show":
            if "FTC_FAST" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_FAST"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["FTC_FAST"] = value.lower()

    def cis_max_cycles(self, value="show"):
        '''
Name: CIS_MAX_CYCLES
Type: INTEGER
Default: 30
Options: Range from 0 to 100

Description: Maximum number of CIS iterative cycles allowed
Recommendation: : Default is usually sufficient.    '''
        if value == "":
            if "CIS_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_MAX_CYCLES"] = value.lower()

    def dscf_max_cycles_level_2(self, value="show"):
        '''
Name: DSCF_MAX_CYCLES_LEVEL_2
Type: INTEGER
Default: 30
Options: Range from 1 to 500

Description: Sets the maximum number of level-2 iterations.
Recommendation: : Use default.    '''
        if value == "":
            if "DSCF_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"]
                print("Keyword removed.")
        elif value == "show":
            if "DSCF_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"] = value.lower()

    def dscf_max_cycles_level_1(self, value="show"):
        '''
Name: DSCF_MAX_CYCLES_LEVEL_1
Type: INTEGER
Factor: 10
Default: 10 [=100]
Options: Range from 0 [=0] to 50 [=500]

Description: Sets the maximum number of level-1 iterations.
Recommendation: : Use default.    '''
        if value == "":
            if "DSCF_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"]
                print("Keyword removed.")
        elif value == "show":
            if "DSCF_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"] = value.lower()

    def geom_opt_diis_subspace(self, value="show"):
        '''
Name: GEOM_OPT_DIIS_SUBSPACE
Type: INTEGER
Default: 0
Options: Range from -1 to 50

Description: Controls maximum size of subspace for GDIIS. 0 turns off GDIIS and -1 causes the program to select the default size.
    '''
        if value == "":
            if "GEOM_OPT_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"] = value.lower()

    def geom_opt_symmetry(self, value="show"):
        '''
Name: GEOM_OPT_SYMMETRY
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of point-group symmetry in the optimization.
    '''
        if value == "":
            if "GEOM_OPT_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_SYMMETRY"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_SYMMETRY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_SYMMETRY"] = value.lower()

    def geom_opt_coordinates(self, value="show"):
        '''
Name: GEOM_OPT_COORDINATES
Type: STRING
Default: 0

Options:
    '0'............................. Cartesian
    '1'............................. Internal
    '2'............................. Z-matrix

Description: Controls the type of optimization coordinates.
Recommendation: : Use the default; delocalized internals are more efficient.    '''
        if value == "":
            if "GEOM_OPT_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_COORDINATES"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_COORDINATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_COORDINATES"] = value.lower()

    def geom_opt_max_step_size(self, value="show"):
        '''
Name: GEOM_OPT_MAX_STEP_SIZE
Type: INTEGER
Factor: 0.001
Default: 300 [=0.300]
Options: Range from 1 [=0.001] to 999 [=0.999]

Description: Maximum allowed step size in the geometry optimization.
    '''
        if value == "":
            if "GEOM_OPT_MAX_STEP_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_MAX_STEP_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"] = value.lower()

    def geom_opt_hessian_update(self, value="show"):
        '''
Name: GEOM_OPT_HESSIAN_UPDATE
Type: STRING
Default: -1

Options:
    '-1'............................ Default
    '0'............................. No Update
    '1'............................. Murtagh-Sargent
    '2'............................. Powell
    '3'............................. Powell/Murtagh-Sargent
    '4'............................. BFGS
    '5'............................. BFGS w/ safeguards

Description: Controls the Hessian update algorithm.  The default dpends on the type of job.
    '''
        if value == "":
            if "GEOM_OPT_HESSIAN_UPDATE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_HESSIAN_UPDATE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"] = value.lower()

    def cfmm_grain(self, value="show"):
        '''
Name: CFMM_GRAIN
Type: STRING
Default: -1

Options:
    '-1'............................ Automatic
    '1'............................. Off
    '8'............................. 8
    '9'............................. 9
    '10'............................ 10
    '11'............................ 11
    '12'............................ 12
    '13'............................ 13
    '14'............................ 14
    '15'............................ 15
    '16'............................ 16

Description: Controls the number of lowest-level boxes in one dimension for CFMM.
Recommendation: : This is an expert option; either use the default, or use a value of 1 if CFMM is not desired.    '''
        if value == "":
            if "CFMM_GRAIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CFMM_GRAIN"]
                print("Keyword removed.")
        elif value == "show":
            if "CFMM_GRAIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CFMM_GRAIN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CFMM_GRAIN"] = value.lower()

    def plots_property(self, value="show"):
        '''
Name: PLOTS_PROPERTY
Type: STRING
Default: None

Options:
    'None'.......................... None
    '0'............................. ESP Only
    '1'............................. ESP & EFIELD
    '2'............................. EFIELD Only

Description: Triggers the calculation of the electrostatic potential and/or the electric field at the points given in the file ESPGrid. 
Recommendation: : Must use this option when IGDESP is specified.    '''
        if value == "":
            if "PLOTS_PROPERTY" in self.dict_of_keywords:
                del self.dict_of_keywords["PLOTS_PROPERTY"]
                print("Keyword removed.")
        elif value == "show":
            if "PLOTS_PROPERTY" in self.dict_of_keywords:
                return self.dict_of_keywords["PLOTS_PROPERTY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PLOTS_PROPERTY"] = value.lower()

    def intracule_conserve_memory(self, value="show"):
        '''
Name: INTRACULE_CONSERVE_MEMORY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Reduce memory required in the evaluation of W(u,v). 
Recommendation: : The low memory option is slower, use default unless memory is limited.    '''
        if value == "":
            if "INTRACULE_CONSERVE_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_CONSERVE_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"] = value.lower()

    def intracule(self, value="show"):
        '''
Name: INTRACULE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether intracule properties are calculated.  Setting this option causes the data in $intracule to be activated.
    '''
        if value == "":
            if "INTRACULE" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE"] = value.lower()

    def intracule_grid(self, value="show"):
        '''
Name: INTRACULE_GRID
Type: STRING
Default: 194

Options:
    '6'............................. 6
    '18'............................ 18
    '26'............................ 26
    '38'............................ 38
    '50'............................ 50
    '74'............................ 74
    '86'............................ 86
    '110'........................... 110
    '146'........................... 146
    '170'........................... 170
    '194'........................... 194
    '230'........................... 230
    '266'........................... 266
    '302'........................... 302
    '350'........................... 350
    '434'........................... 434
    '590'........................... 590
    '770'........................... 770
    '974'........................... 974
    '1202'.......................... 1202
    '1454'.......................... 1454
    '1730'.......................... 1730
    '2030'.......................... 2030
    '2354'.......................... 2354
    '2702'.......................... 2702
    '3074'.......................... 3074
    '3470'.......................... 3470
    '3890'.......................... 3890
    '4334'.......................... 4334
    '4802'.......................... 4802
    '5294'.......................... 5294

Description: Specify angular Lebedev grid for Wigner intracule calculations.
Recommendation: : Larger grids if high accuracy required.    '''
        if value == "":
            if "INTRACULE_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE_GRID"] = value.lower()

    def intracule_wigner_series_limit(self, value="show"):
        '''
Name: INTRACULE_WIGNER_SERIES_LIMIT
Type: INTEGER
Default: 10
Options: Range from 1 to 100

Description: Sets summation limit for Wigner integrals.
Recommendation: : Increase n for greater accuracy.    '''
        if value == "":
            if "INTRACULE_WIGNER_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_WIGNER_SERIES_LIMIT"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_WIGNER_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_WIGNER_SERIES_LIMIT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "INTRACULE_WIGNER_SERIES_LIMIT"] = value.lower()

    def intracule_j_series_limit(self, value="show"):
        '''
Name: INTRACULE_J_SERIES_LIMIT
Type: INTEGER
Default: 40
Options: Range from 1 to 100

Description: Sets summation limit for series expansion evaluation of j_n(x).
Recommendation: : Lower values speed up the calculation, but may affect accuracy.    '''
        if value == "":
            if "INTRACULE_J_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_J_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"] = value.lower()

    def intracule_i_series_limit(self, value="show"):
        '''
Name: INTRACULE_I_SERIES_LIMIT
Type: INTEGER
Default: 40
Options: Range from 1 to 100

Description: Sets summation limit for series expansion evaluation of i_n(x).
Recommendation: : Lower values speed up the calculation, but may affect accuracy.    '''
        if value == "":
            if "INTRACULE_I_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_I_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"] = value.lower()

    def intracule_method(self, value="show"):
        '''
Name: INTRACULE_METHOD
Type: STRING
Default: 0

Options:
    '0'............................. Series
    '1'............................. Quadrature

Description: Use Lebedev quadrature to evaluate Wigner integrals.
    '''
        if value == "":
            if "INTRACULE_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "INTRACULE_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INTRACULE_METHOD"] = value.lower()

    def rca_max_cycles(self, value="show"):
        '''
Name: RCA_MAX_CYCLES
Type: INTEGER
Default: 50
Options: Range from 0 to 100

Description: The maximum number of RCA iterations before switching to DIIS when the SCF algorithm is RCA_DIIS
    '''
        if value == "":
            if "RCA_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "RCA_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RCA_MAX_CYCLES"] = value.lower()

    def scf_max_cycles(self, value="show"):
        '''
Name: SCF_MAX_CYCLES
Type: INTEGER
Default: 50
Options: Range from 0 to 200

Description: Controls the maximum number of SCF iterations permitted.
Recommendation: : Increase for slowly converging systems such as those containing transition metals.    '''
        if value == "":
            if "SCF_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_MAX_CYCLES"] = value.lower()

    def mm_charges(self, value="show"):
        '''
Name: MM_CHARGES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Requests the calculation of multipole-derived charges (MDCs).
Recommendation: : Set to TRUE if MDCs or the traceless form of the multipole moments are desired. The calculation does not take long.    '''
        if value == "":
            if "MM_CHARGES" in self.dict_of_keywords:
                del self.dict_of_keywords["MM_CHARGES"]
                print("Keyword removed.")
        elif value == "show":
            if "MM_CHARGES" in self.dict_of_keywords:
                return self.dict_of_keywords["MM_CHARGES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MM_CHARGES"] = value.lower()

    def molden_format(self, value="show"):
        '''
Name: MOLDEN_FORMAT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Requests a $molden-formatted input file containing information from a Q-Chem job.
    '''
        if value == "":
            if "MOLDEN_FORMAT" in self.dict_of_keywords:
                del self.dict_of_keywords["MOLDEN_FORMAT"]
                print("Keyword removed.")
        elif value == "show":
            if "MOLDEN_FORMAT" in self.dict_of_keywords:
                return self.dict_of_keywords["MOLDEN_FORMAT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOLDEN_FORMAT"] = value.lower()

    def mom_print(self, value="show"):
        '''
Name: MOM_PRINT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Switches printing on within the MOM procedure.
    '''
        if value == "":
            if "MOM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["MOM_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "MOM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["MOM_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOM_PRINT"] = value.lower()

    def moprop_convergence_level_1(self, value="show"):
        '''
Name: MOPROP_CONVERGENCE_LEVEL_1
Type: INTEGER
Default: 6
Options: Range from 0 to 12

Description: Sets the convergence criteria for CPSCF and 1st order TDSCF.
    '''
        if value == "":
            if "MOPROP_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"] = value.lower()

    def moprop_convergence_level_2(self, value="show"):
        '''
Name: MOPROP_CONVERGENCE_LEVEL_2
Type: INTEGER
Default: 6
Options: Range from 0 to 12

Description: Sets the convergence criterium for second-order TDSCF.
    '''
        if value == "":
            if "MOPROP_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"] = value.lower()

    def moprop_diis(self, value="show"):
        '''
Name: MOPROP_DIIS
Type: INTEGER
Factor: 5
Default: 1 [=5]
Options: Range from 0 [=0] to 1 [=5]

Description: Controls the use of Pulays DIIS.
    '''
        if value == "":
            if "MOPROP_DIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_DIIS"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_DIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_DIIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_DIIS"] = value.lower()

    def moprop_diis_subspace(self, value="show"):
        '''
Name: MOPROP_DIIS_SUBSPACE
Type: INTEGER
Default: 20
Options: Range from 0 to 50

Description: Specified the DIIS subspace dimension.
    '''
        if value == "":
            if "MOPROP_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"] = value.lower()

    def moprop_max_cycles_level_2(self, value="show"):
        '''
Name: MOPROP_MAX_CYCLES_LEVEL_2
Type: INTEGER
Default: 50
Options: Range from 1 to 500

Description: The maximal number of iterations for second-order TDSCF.
Recommendation: : Use default.    '''
        if value == "":
            if "MOPROP_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"] = value.lower()

    def moprop_max_cycles_level_1(self, value="show"):
        '''
Name: MOPROP_MAX_CYCLES_LEVEL_1
Type: INTEGER
Default: 50
Options: Range from 1 to 500

Description: The maximal number of iterations for CPSCF and first-order TDSCF.
Recommendation: : Use default.    '''
        if value == "":
            if "MOPROP_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"] = value.lower()

    def moprop_perturbations(self, value="show"):
        '''
Name: MOPROP_PERTURBATIONS
Type: INTEGER
Default: 0
Options: Range from 0 to 20

Description: Set the number of perturbed densities that will to be treated together.
Recommendation: : Use default    '''
        if value == "":
            if "MOPROP_PERTURBATIONS" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_PERTURBATIONS"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_PERTURBATIONS" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_PERTURBATIONS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_PERTURBATIONS"] = value.lower()

    def nbo(self, value="show"):
        '''
Name: NBO
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of the NBO package.
    '''
        if value == "":
            if "NBO" in self.dict_of_keywords:
                del self.dict_of_keywords["NBO"]
                print("Keyword removed.")
        elif value == "show":
            if "NBO" in self.dict_of_keywords:
                return self.dict_of_keywords["NBO"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NBO"] = value.lower()

    def qmmm_charges(self, value="show"):
        '''
Name: QMMM_CHARGES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the printing of QM charges to file.
Recommendation: : Use default unless running calculations with $charmm where charges on the QM region need to be saved.    '''
        if value == "":
            if "QMMM_CHARGES" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_CHARGES"]
                print("Keyword removed.")
        elif value == "show":
            if "QMMM_CHARGES" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_CHARGES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QMMM_CHARGES"] = value.lower()

    def qmmm_print(self, value="show"):
        '''
Name: QMMM_PRINT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the amount of output printed from a QM/MM job.
Recommendation: : Use default unless running calculations with $charmm.    '''
        if value == "":
            if "QMMM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "QMMM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QMMM_PRINT"] = value.lower()

    def multipole_order(self, value="show"):
        '''
Name: MULTIPOLE_ORDER
Type: INTEGER
Default: 4
Options: Range from 0 to 20

Description: Determines the order to which the multipole expansion of the solute charge density is carried out.
Recommendation: : Use default unless higher (or lower) precision is desired.    '''
        if value == "":
            if "MULTIPOLE_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["MULTIPOLE_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "MULTIPOLE_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["MULTIPOLE_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MULTIPOLE_ORDER"] = value.lower()

    def plots_grid(self, value="show"):
        '''
Name: PLOTS_GRID
Type: STRING
Default: None

Options:
    'None'.......................... None
    '-1'............................ User-defined
    '0'............................. Nuclear positions
    '1'............................. Read from file

Description: Controls evaluation of the electrostatic potential on a grid of points. If enabled, the output is in an ACSII file, plot.esp, in the format x, y, z, esp for each point.
    '''
        if value == "":
            if "PLOTS_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["PLOTS_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "PLOTS_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["PLOTS_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PLOTS_GRID"] = value.lower()

    def mulliken(self, value="show"):
        '''
Name: MULLIKEN
Type: STRING
Default: 0

Options:
    '0'............................. None
    '1'............................. Atomic
    '2'............................. Shell

Description: 
    '''
        if value == "":
            if "MULLIKEN" in self.dict_of_keywords:
                del self.dict_of_keywords["MULLIKEN"]
                print("Keyword removed.")
        elif value == "show":
            if "MULLIKEN" in self.dict_of_keywords:
                return self.dict_of_keywords["MULLIKEN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MULLIKEN"] = value.lower()

    def core_character_print(self, value="show"):
        '''
Name: CORE_CHARACTER_PRINT
Type: STRING
Default: 0

Options:
    '0'............................. None
    '1'............................. Occupied MOs
    '2'............................. MOs and AOs

Description: Determines the print level for the CORE_CHARACTER option.
Recommendation: : Use default, unless you are uncertain about what the core character is.    '''
        if value == "":
            if "CORE_CHARACTER_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CORE_CHARACTER_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "CORE_CHARACTER_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CORE_CHARACTER_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CORE_CHARACTER_PRINT"] = value.lower()

    def print_distance_matrix(self, value="show"):
        '''
Name: PRINT_DISTANCE_MATRIX
Type: INTEGER

Description: Controls the printing of the inter-atomic distance matrix
Recommendation: : Use default unless distances are required for large systems    '''
        if value == "":
            if "PRINT_DISTANCE_MATRIX" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_DISTANCE_MATRIX"]
                print("Keyword removed.")
        elif value == "show":
            if "PRINT_DISTANCE_MATRIX" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_DISTANCE_MATRIX"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PRINT_DISTANCE_MATRIX"] = value.lower()

    def qui_print_orbitals(self, value="show"):
        '''
Name: QUI_PRINT_ORBITALS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_PRINT_ORBITALS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_PRINT_ORBITALS"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_PRINT_ORBITALS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_PRINT_ORBITALS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_PRINT_ORBITALS"] = value.lower()

    def qmmm(self, value="show"):
        '''
Name: QMMM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Turns on the Q-Chem/CHARMM interface.
Recommendation: : Use default unless running calculations with $charmm.    '''
        if value == "":
            if "QMMM" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM"]
                print("Keyword removed.")
        elif value == "show":
            if "QMMM" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QMMM"] = value.lower()

    def rca_print(self, value="show"):
        '''
Name: RCA_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 3

Description: Controls the amount of output from a RCA SCF optimization.
    '''
        if value == "":
            if "RCA_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "RCA_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RCA_PRINT"] = value.lower()

    def rpath_coordinates(self, value="show"):
        '''
Name: RPATH_COORDINATES
Type: STRING
Default: 1

Options:
    '1'............................. Mass-weighted
    '2'............................. Z-matrix

Description: Determines which coordinate system to use in the IRC search.
Recommendation: : Mass weighted coordinates are usually more effective.    '''
        if value == "":
            if "RPATH_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_COORDINATES"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_COORDINATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_COORDINATES"] = value.lower()

    def rpath_direction(self, value="show"):
        '''
Name: RPATH_DIRECTION
Type: STRING
Default: 1

Options:
    '1'............................. Positive
    '-1'............................ Negative

Description: Determines the direction of the eigen mode to follow. This will not usually be known prior to the Hessian diagonalization and thereforr both directions will have to be considered.
    '''
        if value == "":
            if "RPATH_DIRECTION" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_DIRECTION"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_DIRECTION" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_DIRECTION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_DIRECTION"] = value.lower()

    def rpath_max_stepsize(self, value="show"):
        '''
Name: RPATH_MAX_STEPSIZE
Type: INTEGER
Factor: 0.001
Default: 150 [=0.150]
Options: Range from 1 [=0.001] to 500 [=0.500]

Description: Specifies the maximum step size to be taken (in a.u.).
    '''
        if value == "":
            if "RPATH_MAX_STEPSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_MAX_STEPSIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "RPATH_MAX_STEPSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_MAX_STEPSIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RPATH_MAX_STEPSIZE"] = value.lower()

    def geom_opt_scf_guess_always(self, value="show"):
        '''
Name: GEOM_OPT_SCF_GUESS_ALWAYS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Switch to force the regeneration of a new initial guess for each series of SCF iterations (for use in geometry optimization).
Recommendation: : Use default unless SCF convergence issues arise    '''
        if value == "":
            if "GEOM_OPT_SCF_GUESS_ALWAYS" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_SCF_GUESS_ALWAYS" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"] = value.lower()

    def scf_guess_print(self, value="show"):
        '''
Name: SCF_GUESS_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 2

Description: Controls printing of guess MOs, Fock and density matrices.
    '''
        if value == "":
            if "SCF_GUESS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_GUESS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_GUESS_PRINT"] = value.lower()

    def scf_print(self, value="show"):
        '''
Name: SCF_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 3

Description: Controls level of output from SCF procedure to Q-Chem output file.
Recommendation: : Proceed with care; can result in extremely large output files at level 2 or higher. These levels are primarily for program debugging.    '''
        if value == "":
            if "SCF_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_PRINT"] = value.lower()

    def solute_radius(self, value="show"):
        '''
Name: SOLUTE_RADIUS
Type: INTEGER
Factor: 0.0001
Default: 0 [=0.0000]
Options: Range from 0 [=0.0000] to 999998 [=99.9999]

Description: Sets the Onsager solvent model cavity radius.
Recommendation: : Use equation (\ref{eq1000}).    '''
        if value == "":
            if "SOLUTE_RADIUS" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLUTE_RADIUS"]
                print("Keyword removed.")
        elif value == "show":
            if "SOLUTE_RADIUS" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLUTE_RADIUS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SOLUTE_RADIUS"] = value.lower()

    def solvent_dielectric(self, value="show"):
        '''
Name: SOLVENT_DIELECTRIC
Type: INTEGER
Factor: 0.0001
Default: 0 [=0.0000]
Options: Range from 0 [=0.0000] to 999998 [=99.9999]

Description: Sets the dielectric constant of the Onsager solvent continuum.
Recommendation: : As per required solvent.    '''
        if value == "":
            if "SOLVENT_DIELECTRIC" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLVENT_DIELECTRIC"]
                print("Keyword removed.")
        elif value == "show":
            if "SOLVENT_DIELECTRIC" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLVENT_DIELECTRIC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SOLVENT_DIELECTRIC"] = value.lower()

    def rca_switch_thresh(self, value="show"):
        '''
Name: RCA_SWITCH_THRESH
Type: INTEGER
Default: 3
Options: Range from 0 to 12

Description: The threshold for switching between RCA and DIIS when the SCF algorithm is set to RCA_DIIS
    '''
        if value == "":
            if "RCA_SWITCH_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_SWITCH_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "RCA_SWITCH_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_SWITCH_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["RCA_SWITCH_THRESH"] = value.lower()

    def qui_angular_grid(self, value="show"):
        '''
Name: QUI_ANGULAR_GRID
Type: STRING
Default: SG-1

Options:
    'SG-0'.......................... SG-0
    'SG-1'.......................... SG-1

    '6'............................. 6
    '18'............................ 18
    '26'............................ 26
    '38'............................ 38
    '50'............................ 50
    '74'............................ 74
    '86'............................ 86
    '110'........................... 110
    '146'........................... 146
    '170'........................... 170
    '194'........................... 194
    '230'........................... 230
    '266'........................... 266
    '302'........................... 302
    '350'........................... 350
    '434'........................... 434
    '590'........................... 590
    '770'........................... 770
    '974'........................... 974
    '1202'.......................... 1202
    '1454'.......................... 1454
    '1730'.......................... 1730
    '2030'.......................... 2030
    '2354'.......................... 2354
    '2702'.......................... 2702
    '3074'.......................... 3074
    '3470'.......................... 3470
    '3890'.......................... 3890
    '4334'.......................... 4334
    '4802'.......................... 4802
    '5294'.......................... 5294

Description: Specifies the quadrature grid to be used for evaluating the exchange-correlation component of the energy.  Either a standard grid should be selected, or a Lebedev grid with the corresponding number of points.

Recommendation: : Use the default unless convergence difficulties arise.  Larger grids are required for calculations involving derivatives and excited states.    '''
        if value == "":
            if "QUI_ANGULAR_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_ANGULAR_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_ANGULAR_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_ANGULAR_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_ANGULAR_GRID"] = value.lower()

    def qui_radial_grid(self, value="show"):
        '''
Name: QUI_RADIAL_GRID
Type: INTEGER
Default: 50
Options: Range from 1 to 200

Description: Specifies the number of radial point for the exchange-correlation quadrature.
    '''
        if value == "":
            if "QUI_RADIAL_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_RADIAL_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_RADIAL_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_RADIAL_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_RADIAL_GRID"] = value.lower()

    def isotopes(self, value="show"):
        '''
Name: ISOTOPES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies if non-default masses are to be used in the frequency calculation.  If this option is selected the $isotopes section is read.
    '''
        if value == "":
            if "ISOTOPES" in self.dict_of_keywords:
                del self.dict_of_keywords["ISOTOPES"]
                print("Keyword removed.")
        elif value == "show":
            if "ISOTOPES" in self.dict_of_keywords:
                return self.dict_of_keywords["ISOTOPES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ISOTOPES"] = value.lower()

    def memory_static(self, value="show"):
        '''
Name: MEMORY_STATIC
Type: INTEGER
Default: 64
Options: Range from 1 to 512

Description: Sets the memory (in megabytes) for individual program modules.
Recommendation: : For direct and semi-direct MP2 calculations, this must exceed OVN + requirements for AO integral evaluation (32-160 Mb).    '''
        if value == "":
            if "MEMORY_STATIC" in self.dict_of_keywords:
                del self.dict_of_keywords["MEMORY_STATIC"]
                print("Keyword removed.")
        elif value == "show":
            if "MEMORY_STATIC" in self.dict_of_keywords:
                return self.dict_of_keywords["MEMORY_STATIC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MEMORY_STATIC"] = value.lower()

    def memory_total(self, value="show"):
        '''
Name: MEMORY_TOTAL
Type: INTEGER
Factor: 10
Default: 200 [=2000]
Options: Range from 12 [=128] to 800 [=8000]

Description: Sets the total memory available to Q-Chem, in megabytes.
Recommendation: : Use default, or set to the physical memory of your machine.    '''
        if value == "":
            if "MEMORY_TOTAL" in self.dict_of_keywords:
                del self.dict_of_keywords["MEMORY_TOTAL"]
                print("Keyword removed.")
        elif value == "show":
            if "MEMORY_TOTAL" in self.dict_of_keywords:
                return self.dict_of_keywords["MEMORY_TOTAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MEMORY_TOTAL"] = value.lower()

    def aimd_method(self, value="show"):
        '''
Name: AIMD_METHOD
Type: STRING
Default: BOMD

Options:
    'BOMD'.......................... BOMD
    'Curvy'......................... Curvy

Description: Selects an ab initio molecular dynamics algorithm.
Recommendation: : Born-oppenheimer MD (BOMD) yields exact classical molecular dynamics, provided that the energy is tolerably conserved. Curvy-steps extended Lagrangian MD (Curvy) is an approximation to exact classical dynamics whose validity should be tested for the properties of interest.     '''
        if value == "":
            if "AIMD_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_METHOD"] = value.lower()

    def cc_dthreshold(self, value="show"):
        '''
Name: CC_DTHRESHOLD
Type: INTEGER
Factor: 0.00001
Default: 1 [=0.00001]
Options: Range from 0 [=0.000001] to 99999 [=1.00000]

Description: Specifies threshold for including a new expansion vector in the iterative Davidson diagonalization. Their norm must be above this threshold. 
Recommendation: : Use default unless converge problems are encountered. Should normally be set to the same values as CC_DCONVERGENCE, if convergence problems arise try setting to a value less than CC_DCONVERGENCE.    '''
        if value == "":
            if "CC_DTHRESHOLD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DTHRESHOLD"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DTHRESHOLD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DTHRESHOLD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DTHRESHOLD"] = value.lower()

    def cc_state_derivative(self, value="show"):
        '''
Name: CC_STATE_DERIVATIVE
Type: INTEGER

Description: Selects which EOM or CIS(D) state is to be considered for optimization or property calculations.
    '''
        if value == "":
            if "CC_STATE_DERIVATIVE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_DERIVATIVE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_STATE_DERIVATIVE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_DERIVATIVE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_STATE_DERIVATIVE"] = value.lower()

    def gui(self, value="show"):
        '''
Name: GUI
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the output of auxiliary information for third party packages.
Recommendation: : Use default unless the additional information is required. Please note that any existing Test.FChk file will be overwritten.    '''
        if value == "":
            if "GUI" in self.dict_of_keywords:
                del self.dict_of_keywords["GUI"]
                print("Keyword removed.")
        elif value == "show":
            if "GUI" in self.dict_of_keywords:
                return self.dict_of_keywords["GUI"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GUI"] = value.lower()

    def cc_eom_two_particle_properties(self, value="show"):
        '''
Name: CC_EOM_TWO_PARTICLE_PROPERTIES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Request for calculation of non-relaxed two-particle EOM-CCSD target state properties.  The two-part properties currently include<>. The one-particle properties will also be calculated since the additional cost of these is small in comparison.  The vairable CC_EXSTATES_PROP must also be set.

Recommendation:  Two-particle properties are extremely computationally expensive since the require calculation and use of the two-particle density matrix (the cost being about the same as the cost of an analytic gradient calculation for each state.    '''
        if value == "":
            if "CC_EOM_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TWO_PARTICLE_PROPERTIES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TWO_PARTICLE_PROPERTIES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "CC_EOM_TWO_PARTICLE_PROPERTIES"] = value.lower()

    def qui_section_opt(self, value="show"):
        '''
Name: QUI_SECTION_OPT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Adds the $opt section for specifying constraints in the geometry optimization
    '''
        if value == "":
            if "QUI_SECTION_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SECTION_OPT"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SECTION_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SECTION_OPT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_SECTION_OPT"] = value.lower()

    def dft_d_a(self, value="show"):
        '''
Name: DFT_D_A
Type: INTEGER
Factor: 0.01
Default: 600 [=6.00]
Options: Range from 1 [=0.01] to 10000 [=100.00]

Description: Controls the strength of the dispersion corrections in the Chai-Head-Gordon scheme.   The default value should be apprpriate for most systems.

    '''
        if value == "":
            if "DFT_D_A" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D_A"]
                print("Keyword removed.")
        elif value == "show":
            if "DFT_D_A" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D_A"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFT_D_A"] = value.lower()

    def scf_guess_mix(self, value="show"):
        '''
Name: SCF_GUESS_MIX
Type: INTEGER
Factor: 10
Default: 0 [=0]
Options: Range from 0 [=0] to 10 [=100]

Description: Controls mixing of LUMO and HOMO to break symmetry in the initial guess. For unrestricted jobs, the mixing is performed only for the alpha orbitals.
Recommendation: : When performing unrestricted calculations on molecules with an even number of electrons, it is often necessary to break alpha-beta symmetry in the initial guess with this option, or by specifying input for $occupied.    '''
        if value == "":
            if "SCF_GUESS_MIX" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS_MIX"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_GUESS_MIX" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS_MIX"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_GUESS_MIX"] = value.lower()

    def cc_print(self, value="show"):
        '''
Name: CC_PRINT
Type: INTEGER

Description: Controls the output from post-MP2 coupled-cluster module of QChem
Recommendation: : Increase if you need more output and don't like trees    '''
        if value == "":
            if "CC_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRINT"] = value.lower()

    def moprop_save_last_gpx(self, value="show"):
        '''
Name: MOPROP_SAVE_LAST_GPX
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Save last G[P]x when calculating dynamic polarizabilities in order to call mopropman in a second run with MOPROP = 102.
    '''
        if value == "":
            if "MOPROP_SAVE_LAST_GPX" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"]
                print("Keyword removed.")
        elif value == "show":
            if "MOPROP_SAVE_LAST_GPX" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"] = value.lower()

    def qui_title(self, value="show"):
        '''
Name: QUI_TITLE
Type: STRING
Default:  

Options:
    ' '.............................  

Description: Sets the lable for this section of the input file.
    '''
        if value == "":
            if "QUI_TITLE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_TITLE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_TITLE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_TITLE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_TITLE"] = value.lower()

    def smx_solvent(self, value="show"):
        '''
Name: SMX_SOLVENT
Type: STRING
Default: water

Options:
    '111trichloroethane'............ 111trichloroethane
    '112trichloroethane'............ 112trichloroethane
    '11dichloroethane'.............. 11dichloroethane
    '124trimethylbenzene'........... 124trimethylbenzene
    '14dioxane'..................... 14dioxane
    '1bromo2methylpropane'.......... 1bromo2methylpropane
    '1bromopentane'................. 1bromopentane
    '1bromopropane'................. 1bromopropane
    '1butanol'...................... 1butanol
    '1chloropentane'................ 1chloropentane
    '1chloropropane'................ 1chloropropane
    '1decanol'...................... 1decanol
    '1fluorooctane'................. 1fluorooctane
    '1heptanol'..................... 1heptanol
    '1hexanol'...................... 1hexanol
    '1hexene'....................... 1hexene
    '1hexyne'....................... 1hexyne
    '1iodobutane'................... 1iodobutane
    '1iodopentene'.................. 1iodopentene
    '1iodopropane'.................. 1iodopropane
    '1nitropropane'................. 1nitropropane
    '1nonanol'...................... 1nonanol
    '1octanol'...................... 1octanol
    '1pentanol'..................... 1pentanol
    '1pentene'...................... 1pentene
    '1pentyne'...................... 1pentyne
    '1propanol'..................... 1propanol
    '222trifluoroethanol'........... 222trifluoroethanol
    '224trimethylpentane'........... 224trimethylpentane
    '24dimethylpentane'............. 24dimethylpentane
    '24dimethylpyridine'............ 24dimethylpyridine
    '26dimethylpyridine'............ 26dimethylpyridine
    '2bromopropane'................. 2bromopropane
    '2chlorobutane'................. 2chlorobutane
    '2heptanone'.................... 2heptanone
    '2hexanone'..................... 2hexanone
    '2methylpentane'................ 2methylpentane
    '2methylpyridine'............... 2methylpyridine
    '2nitropropane'................. 2nitropropane
    '2octanone'..................... 2octanone
    '2pentanone'.................... 2pentanone
    '2propanol'..................... 2propanol
    '2propen1ol'.................... 2propen1ol
    '3methylpyridine'............... 3methylpyridine
    '3pentanone'.................... 3pentanone
    '4heptanone'.................... 4heptanone
    '4methyl2pentanone'............. 4methyl2pentanone
    '4methylpyridine'............... 4methylpyridine
    '5nonanone'..................... 5nonanone
    'aceticacid'.................... aceticacid
    'acetone'....................... acetone
    'acetonitrile'.................. acetonitrile
    'aniline'....................... aniline
    'anisole'....................... anisole
    'benzaldehyde'.................. benzaldehyde
    'benzene'....................... benzene
    'benzonitrile'.................. benzonitrile
    'benzylalcohol'................. benzylalcohol
    'bromobenzene'.................. bromobenzene
    'bromoethane'................... bromoethane
    'bromooctane'................... bromooctane
    'butanal'....................... butanal
    'butanoicacid'.................. butanoicacid
    'butanone'...................... butanone
    'butanonitrile'................. butanonitrile
    'butylethanoate'................ butylethanoate
    'butylamine'.................... butylamine
    'butylbenzene'.................. butylbenzene
    'carbondisulfide'............... carbondisulfide
    'carbontet'..................... carbontet
    'chlorobenzene'................. chlorobenzene
    'chlorotoluene'................. chlorotoluene
    'cis12dimethylcyclohexane'...... cis12dimethylcyclohexane
    'decalin'....................... decalin
    'cyclohexane'................... cyclohexane
    'cyclohexanone'................. cyclohexanone
    'cyclopentane'.................. cyclopentane
    'cyclopentanol'................. cyclopentanol
    'cyclopentanone'................ cyclopentanone
    'decane'........................ decane
    'dibromomethane'................ dibromomethane
    'dibutylether'.................. dibutylether
    'dichloromethane'............... dichloromethane
    'diethylether'.................. diethylether
    'diethylsulfide'................ diethylsulfide
    'diethylamine'.................. diethylamine
    'diiodomethane'................. diiodomethane
    'dimethyldisulfide'............. dimethyldisulfide
    'dimethylacetamide'............. dimethylacetamide
    'dimethylformamide'............. dimethylformamide
    'dimethylpyridine'.............. dimethylpyridine
    'dmso'.......................... dmso
    'dipropylamine'................. dipropylamine
    'dodecane'...................... dodecane
    'E12dichloroethene'............. E12dichloroethene
    'E2pentene'..................... E2pentene
    'ethanethiol'................... ethanethiol
    'ethanol'....................... ethanol
    'ethylethanoate'................ ethylethanoate
    'ethylmethanoate'............... ethylmethanoate
    'ethylphenylether'.............. ethylphenylether
    'ethylbenzene'.................. ethylbenzene
    'ethyleneglycol'................ ethyleneglycol
    'fluorobenzene'................. fluorobenzene
    'formamide'..................... formamide
    'formicacid'.................... formicacid
    'hexadecyliodide'............... hexadecyliodide
    'hexanoic'...................... hexanoic
    'iodobenzene'................... iodobenzene
    'iodoethane'.................... iodoethane
    'iodomethane'................... iodomethane
    'isobutanol'.................... isobutanol
    'isopropylether'................ isopropylether
    'isopropylbenzene'.............. isopropylbenzene
    'isopropyltoluene'.............. isopropyltoluene
    'mcresol'....................... mcresol
    'mesitylene'.................... mesitylene
    'methanol'...................... methanol
    'methylbenzoate'................ methylbenzoate
    'methylethanoate'............... methylethanoate
    'methylmethanoate'.............. methylmethanoate
    'methylphenylketone'............ methylphenylketone
    'methylpropanoate'.............. methylpropanoate
    'methylbutanoate'............... methylbutanoate
    'methylcyclohexane'............. methylcyclohexane
    'methylformamide'............... methylformamide
    'methylformamide'............... methylformamide
    'heptane'....................... heptane
    'hexadecane'.................... hexadecane
    'hexane'........................ hexane
    'nitrobenzene'.................. nitrobenzene
    'nitroethane'................... nitroethane
    'nitromethane'.................. nitromethane
    'methylaniline'................. methylaniline
    'nonane'........................ nonane
    'octane'........................ octane
    'pentane'....................... pentane
    'ochlorotoluene'................ ochlorotoluene
    'ocresol'....................... ocresol
    'odichlorobenzene'.............. odichlorobenzene
    'onitrotoluene'................. onitrotoluene
    'oxylene'....................... oxylene
    'pentadecane'................... pentadecane
    'pentanal'...................... pentanal
    'pentanoicacid'................. pentanoicacid
    'pentylethanoate'............... pentylethanoate
    'pentylamine'................... pentylamine
    'perfluorobenzene'.............. perfluorobenzene
    'phenyletherphenylether'........ phenyletherphenylether
    'propanal'...................... propanal
    'propanoicacid'................. propanoicacid
    'propanonitrile'................ propanonitrile
    'propylethanoate'............... propylethanoate
    'propylamine'................... propylamine
    'pxylene'....................... pxylene
    'pyridine'...................... pyridine
    'pyrrolidine'................... pyrrolidine
    'secbutanol'.................... secbutanol
    'tbutanol'...................... tbutanol
    'tbutylbenzene'................. tbutylbenzene
    'tetrachloroethene'............. tetrachloroethene
    'thf'........................... thf
    'tetrahyrothiophenedioxide'..... tetrahyrothiophenedioxide
    'tetralin'...................... tetralin
    'thiophene'..................... thiophene
    'thiophenol'.................... thiophenol
    'toluene'....................... toluene
    'transdecalin'.................. transdecalin
    'tribromomethane'............... tribromomethane
    'tributylphosphate'............. tributylphosphate
    'trichloroethene'............... trichloroethene
    'trichloromethane'.............. trichloromethane
    'triethylamine'................. triethylamine
    'undecane'...................... undecane
    'water'......................... water
    'Z12dichloroethene'............. Z12dichloroethene

Description: Specifies which solvent to use in the SM8 model.
    '''
        if value == "":
            if "SMX_SOLVENT" in self.dict_of_keywords:
                del self.dict_of_keywords["SMX_SOLVENT"]
                print("Keyword removed.")
        elif value == "show":
            if "SMX_SOLVENT" in self.dict_of_keywords:
                return self.dict_of_keywords["SMX_SOLVENT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SMX_SOLVENT"] = value.lower()

    def smx_solvation(self, value="show"):
        '''
Name: SMX_SOLVATION
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Sets whether or not to use the SM8 solvation model.
    '''
        if value == "":
            if "SMX_SOLVATION" in self.dict_of_keywords:
                del self.dict_of_keywords["SMX_SOLVATION"]
                print("Keyword removed.")
        elif value == "show":
            if "SMX_SOLVATION" in self.dict_of_keywords:
                return self.dict_of_keywords["SMX_SOLVATION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SMX_SOLVATION"] = value.lower()

    def link_atom_projection(self, value="show"):
        '''
Name: LINK_ATOM_PROJECTION
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether to perform a link-atom projection, which is necessary in a full QM/MM hessian evaluation on a system with link atoms.
    '''
        if value == "":
            if "LINK_ATOM_PROJECTION" in self.dict_of_keywords:
                del self.dict_of_keywords["LINK_ATOM_PROJECTION"]
                print("Keyword removed.")
        elif value == "show":
            if "LINK_ATOM_PROJECTION" in self.dict_of_keywords:
                return self.dict_of_keywords["LINK_ATOM_PROJECTION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["LINK_ATOM_PROJECTION"] = value.lower()

    def qmmm_full_hessian(self, value="show"):
        '''
Name: QMMM_FULL_HESSIAN
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QMMM_FULL_HESSIAN" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_FULL_HESSIAN"]
                print("Keyword removed.")
        elif value == "show":
            if "QMMM_FULL_HESSIAN" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_FULL_HESSIAN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QMMM_FULL_HESSIAN"] = value.lower()

    def gaussian_blur(self, value="show"):
        '''
Name: GAUSSIAN_BLUR
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Enables the use of Gaussian-delocalized external charges in a QM/MM calculation.  If set to FALSE, then regular point charges are used.
    '''
        if value == "":
            if "GAUSSIAN_BLUR" in self.dict_of_keywords:
                del self.dict_of_keywords["GAUSSIAN_BLUR"]
                print("Keyword removed.")
        elif value == "show":
            if "GAUSSIAN_BLUR" in self.dict_of_keywords:
                return self.dict_of_keywords["GAUSSIAN_BLUR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GAUSSIAN_BLUR"] = value.lower()

    def hess_proj_trm(self, value="show"):
        '''
Name: HESS_PROJ_TRM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether or not to project out the rotational and translational degrees of freedom in a frequency calculation.
    '''
        if value == "":
            if "HESS_PROJ_TRM" in self.dict_of_keywords:
                del self.dict_of_keywords["HESS_PROJ_TRM"]
                print("Keyword removed.")
        elif value == "show":
            if "HESS_PROJ_TRM" in self.dict_of_keywords:
                return self.dict_of_keywords["HESS_PROJ_TRM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["HESS_PROJ_TRM"] = value.lower()

    def geom_opt_iproj(self, value="show"):
        '''
Name: GEOM_OPT_IPROJ
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Allows the molecule to reorient during a geometry optimization.  Turn this option off if using external charges.
    '''
        if value == "":
            if "GEOM_OPT_IPROJ" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_IPROJ"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_IPROJ" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_IPROJ"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_IPROJ"] = value.lower()

    def qui_qchem_executable(self, value="show"):
        '''
Name: QUI_QCHEM_EXECUTABLE
Type: STRING
Default: qchem_s.exe

Options:
    'qcprog.exe'.................... qcprog.exe
    'qchem_s.exe'................... qchem_s.exe

Description: 
    '''
        if value == "":
            if "QUI_QCHEM_EXECUTABLE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_QCHEM_EXECUTABLE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"] = value.lower()

    def qui_avogadro_visualize_file(self, value="show"):
        '''
Name: QUI_AVOGADRO_VISUALIZE_FILE
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Determines which file is passed as an argument to Avogadro.  If true, then the .out file is passed, if false the .Fchk file is passed.
    '''
        if value == "":
            if "QUI_AVOGADRO_VISUALIZE_FILE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_AVOGADRO_VISUALIZE_FILE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"] = value.lower()

    def qui_windows_directory(self, value="show"):
        '''
Name: QUI_WINDOWS_DIRECTORY
Type: STRING
Default: /Windows/System32

Options:
    '/Windows/System32'............. /Windows/System32

Description: Sets the directory for tskill or taskkill
    '''
        if value == "":
            if "QUI_WINDOWS_DIRECTORY" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_WINDOWS_DIRECTORY" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"] = value.lower()

    def qui_windows_kill_command(self, value="show"):
        '''
Name: QUI_WINDOWS_KILL_COMMAND
Type: STRING
Default: tskill qchem_s

Options:
    'tskill qchem_s'................ tskill qchem_s
    'taskkill /IM qchem_s.exe /F'... taskkill /IM qchem_s.exe /F

Description: The command required to kill qchem jobs, only used on Windows
    '''
        if value == "":
            if "QUI_WINDOWS_KILL_COMMAND" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_WINDOWS_KILL_COMMAND" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"] = value.lower()

    def pdb_print(self, value="show"):
        '''
Name: PDB_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 2

Description: Prints final coordinates at the end of the output file using the PDB format.
    '''
        if value == "":
            if "PDB_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["PDB_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "PDB_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["PDB_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PDB_PRINT"] = value.lower()

    def aaaa(self, value="show"):
        '''
Name: AAAA
Type: STRING
Default: New option

Options:
    '0'............................. 0
    'New option'.................... New option
    '1'............................. 1
    '1'............................. 1
    '1'............................. 1

Description: 
    '''
        if value == "":
            if "AAAA" in self.dict_of_keywords:
                del self.dict_of_keywords["AAAA"]
                print("Keyword removed.")
        elif value == "show":
            if "AAAA" in self.dict_of_keywords:
                return self.dict_of_keywords["AAAA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AAAA"] = value.lower()

    def threads(self, value="show"):
        '''
Name: THREADS
Type: INTEGER
Default: 1
Options: Range from 1 to 1024

Description: Number of threads in shared memory parallel calculations.
    '''
        if value == "":
            if "THREADS" in self.dict_of_keywords:
                del self.dict_of_keywords["THREADS"]
                print("Keyword removed.")
        elif value == "show":
            if "THREADS" in self.dict_of_keywords:
                return self.dict_of_keywords["THREADS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["THREADS"] = value.lower()

    def cc_max_iter(self, value="show"):
        '''
Name: CC_MAX_ITER
Type: INTEGER
Default: 200
Options: Range from 1 to 1000

Description: Maximum number of iterations to optimize the coupled-cluster energy. 
    '''
        if value == "":
            if "CC_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MAX_ITER"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MAX_ITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_MAX_ITER"] = value.lower()

    def cc_memory(self, value="show"):
        '''
Name: CC_MEMORY
Type: INTEGER
Factor: 10
Default: 150 [=1500]
Options: Range from 19 [=192] to 100000 [=1000000]

Description: Specifies the maximum size, in Mb, of the buffers for in-core storage of block-tensors in CCMAN and CCMAN2.
    '''
        if value == "":
            if "CC_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MEMORY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MEMORY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_MEMORY"] = value.lower()

    def mem_total(self, value="show"):
        '''
Name: MEM_TOTAL
Type: INTEGER
Factor: 10
Default: 200 [=2000]
Options: Range from 20 [=200] to 100000 [=1000000]

Description: Sets the total memory available to Q-Chem, in megabytes.
    '''
        if value == "":
            if "MEM_TOTAL" in self.dict_of_keywords:
                del self.dict_of_keywords["MEM_TOTAL"]
                print("Keyword removed.")
        elif value == "show":
            if "MEM_TOTAL" in self.dict_of_keywords:
                return self.dict_of_keywords["MEM_TOTAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MEM_TOTAL"] = value.lower()

    def mem_static(self, value="show"):
        '''
Name: MEM_STATIC
Type: INTEGER
Factor: 10
Default: 24 [=240]
Options: Range from 3 [=32] to 100000 [=1000000]

Description: Sets the memory (in megabytes) for individual fortran program modules.
Recommendation: : For direct and semi-direct MP2 calculations, this must exceed OVN + requirements for AO integral evaluation (32-160 Mb).    '''
        if value == "":
            if "MEM_STATIC" in self.dict_of_keywords:
                del self.dict_of_keywords["MEM_STATIC"]
                print("Keyword removed.")
        elif value == "show":
            if "MEM_STATIC" in self.dict_of_keywords:
                return self.dict_of_keywords["MEM_STATIC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MEM_STATIC"] = value.lower()

    def cc_incl_core_corr(self, value="show"):
        '''
Name: CC_INCL_CORE_CORR
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether to include the correlation contribution from frozen core orbitals in non iterative (2) corrections, such as OD(2) and CCSD(2).
Recommendation: : Use default unless no core-valence or core correlation is desired (e.g., for comparison with other methods or because the basis used cannot describe core correlation).    '''
        if value == "":
            if "CC_INCL_CORE_CORR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_INCL_CORE_CORR"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_INCL_CORE_CORR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_INCL_CORE_CORR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_INCL_CORE_CORR"] = value.lower()

    def cc_ref_prop(self, value="show"):
        '''
Name: CC_REF_PROP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments , , and  of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.
Recommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_REF_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REF_PROP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REF_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REF_PROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REF_PROP"] = value.lower()

    def cc_nguess_doubles(self, value="show"):
        '''
Name: CC_NGUESS_DOUBLES
Type: INTEGER

Description: Specifies number of excited state guess vectors which are double excitations. 
Recommendation: : This should be set to the expected number of doubly excited states (see also EOM_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "CC_NGUESS_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_NGUESS_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_NGUESS_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_NGUESS_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_NGUESS_DOUBLES"] = value.lower()

    def ccman2(self, value="show"):
        '''
Name: CCMAN2
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "CCMAN2" in self.dict_of_keywords:
                del self.dict_of_keywords["CCMAN2"]
                print("Keyword removed.")
        elif value == "show":
            if "CCMAN2" in self.dict_of_keywords:
                return self.dict_of_keywords["CCMAN2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CCMAN2"] = value.lower()

    def cc_t_conv(self, value="show"):
        '''
Name: CC_T_CONV
Type: INTEGER
Default: 8
Options: Range from 0 to 12

Description: 
    '''
        if value == "":
            if "CC_T_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_T_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_T_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_T_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_T_CONV"] = value.lower()

    def cc_ref_prop_te(self, value="show"):
        '''
Name: CC_REF_PROP_TE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Request for calculation of non-relaxed two-particle CCSD properties. The two-particle properties currently include . The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation is
inferior compared to the cost of . The variable CC_REF_PROP must be also set to TRUE.
    '''
        if value == "":
            if "CC_REF_PROP_TE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REF_PROP_TE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REF_PROP_TE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REF_PROP_TE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REF_PROP_TE"] = value.lower()

    def eom_corr(self, value="show"):
        '''
Name: EOM_CORR
Type: STRING
Default: CIS

Options:
    'CIS'........................... CIS
    'CIS(D)'........................ CIS(D)
    'SDT'........................... SDT
    'DT'............................ DT
    'SD'............................ SD
    'SD(dT)'........................ SD(dT)
    'SD(fT)'........................ SD(fT)
    'SD(sT)'........................ SD(sT)

Description: Specifies the correlation level
    '''
        if value == "":
            if "EOM_CORR" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_CORR"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_CORR" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_CORR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_CORR"] = value.lower()

    def eom_davidson_max_iter(self, value="show"):
        '''
Name: EOM_DAVIDSON_MAX_ITER
Type: INTEGER
Default: 30
Options: Range from 0 to 100

Description: Maximum number of iteration allowed for Davidson diagonalization procedure.
n User-defined number of iterations
    '''
        if value == "":
            if "EOM_DAVIDSON_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_DAVIDSON_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"] = value.lower()

    def eom_ngues_doubles(self, value="show"):
        '''
Name: EOM_NGUES_DOUBLES
Type: INTEGER
Default: 0
Options: Range from 0 to 1000

Description: Specifies number of excited state guess vectors which are double excitations. 
Options: n Include n guess vectors that are double excitations
Recommendation: : This should be set to the expected number of doubly excited states (see also CC_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "EOM_NGUES_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUES_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_NGUES_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUES_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_NGUES_DOUBLES"] = value.lower()

    def eom_nguess_doubles(self, value="show"):
        '''
Name: EOM_NGUESS_DOUBLES
Type: INTEGER
Default: 0
Options: Range from 0 to 1000

Description: Specifies number of excited state guess vectors which are double excitations. 
Options: n Include n guess vectors that are double excitations
Recommendation: : This should be set to the expected number of doubly excited states (see also CC_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "EOM_NGUESS_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUESS_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_NGUESS_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUESS_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_NGUESS_DOUBLES"] = value.lower()

    def eom_nguess_singles(self, value="show"):
        '''
Name: EOM_NGUESS_SINGLES
Type: INTEGER
Default: 0
Options: Range from 0 to 1000

Description: Specifies number of excited state guess vectors which are single excitations. 
Options: n Include n guess vectors that are single excitations
Recommendation: : Should be greater or equal than the number of excited states requested.    '''
        if value == "":
            if "EOM_NGUESS_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUESS_SINGLES"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_NGUESS_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUESS_SINGLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_NGUESS_SINGLES"] = value.lower()

    def cc_eom_trans_prop(self, value="show"):
        '''
Name: CC_EOM_TRANS_PROP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_REFSYM and CC_STATE_DERIV must be specified for this state.
Recommendation: : Additional equations (for the left EOM-CCSD eigenvectors plus lambda CCSD equations in case if transition properties between the CCSD reference and EOM-CCSD target states are requested) need to be solved for transition properties, approximately doubling the computational cost. The cost of the transition properties calculation itself is low.    '''
        if value == "":
            if "CC_EOM_TRANS_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TRANS_PROP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_TRANS_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TRANS_PROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_TRANS_PROP"] = value.lower()

    def cc_eom_prop(self, value="show"):
        '''
Name: CC_EOM_PROP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the non-relaxed (expectation value) one-particle EOM-CCSD target state properties will be calculated. The properties currently include permanent dipole moment, the second moments 2>, 2>, and 2> of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.
Recommendation: : Additional equations (EOM-CCSD equations for the left eigenvectors) need to be solved for properties, approximately doubling the cost of calculation for each irrep. Sometimes the equations for left and right eigenvectors converge to different sets of target states. In this case, the simultaneous iterations of left and right vectors will diverge, and the properties for several or all the target states may be incorrect! The problem can be solved by varying the number of requested states, specified with CC_NLOWSPIN and CC_NHIGHSPIN, or the number of guess vectors (CC_NGUESS_SINGLES). The cost of the one-particle properties calculation itself is low. The one-particle density of an EOM-CCSD target state can be analyzed with NBO package by specifying the state with CC_REFSYM and CC_STATE_DERIV and requesting NBO=TRUE and CC_EXSTATES_PROP=TRUE.    '''
        if value == "":
            if "CC_EOM_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_PROP"] = value.lower()

    def cc_eom_prop_te(self, value="show"):
        '''
Name: CC_EOM_PROP_TE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Request for calculation of non-relaxed two-particle EOM-CC properties. The two-particle properties currently include . The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation
is inferior compared to the cost of . The variable CC_EOM_PROP must be also set to TRUE. Alternatively, CC_CALC_SSQ can be used to request  calculation.
    '''
        if value == "":
            if "CC_EOM_PROP_TE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROP_TE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_PROP_TE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROP_TE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_PROP_TE"] = value.lower()

    def cc_fullresponse(self, value="show"):
        '''
Name: CC_FULLRESPONSE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Fully relaxed properties (including orbital relaxation terms) will be computed. The variable CC EOM PROP must be also set to TRUE.
Recommendation: : Not available for non-UHF/RHF references. Only available for EOM/CI methods for which analytic gradients are available.    '''
        if value == "":
            if "CC_FULLRESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_FULLRESPONSE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_FULLRESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_FULLRESPONSE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_FULLRESPONSE"] = value.lower()

    def eom_davidson_threshold(self, value="show"):
        '''
Name: EOM_DAVIDSON_THRESHOLD
Type: STRING
Default: 00105

Options:
    '00105'......................... 00105

Description: Species threshold for including a new expansion vector in the iterative Davidson diagonalization. Their norm must be above this threshold.
abcde Integer code is mapped to abc x 10^-de
    '''
        if value == "":
            if "EOM_DAVIDSON_THRESHOLD" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_DAVIDSON_THRESHOLD" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"] = value.lower()

    def cc_diis(self, value="show"):
        '''
Name: CC_DIIS
Type: INTEGER
Default: 0
Options: Range from 0 to 2

Description: Specify the version of Pulay's Direct Inversion of the Iterative Subspace (DIIS) convergence accelerator to be used in the coupled{cluster code.
0 Activates procedure 2 initially, and procedure 1 when gradients are smaller
than DIIS12 SWITCH.
1 Uses error vectors dened as dierences between parameter vectors from
successive iterations. Most ecient near convergence.
2 Error vectors are dened as gradients scaled by square root of the
approximate diagonal Hessian. Most ecient far from convergence.
Recommendation: : DIIS1 can be more stable. If DIIS problems are encountered in the early stages of a calculation (when gradients are large) try DIIS 1.    '''
        if value == "":
            if "CC_DIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DIIS"] = value.lower()

    def cc_dov_thresh(self, value="show"):
        '''
Name: CC_DOV_THRESH
Type: STRING
Default: 2502

Options:
    '2502'.......................... 2502

Description: Specifies the minimum allowed values for the coupled-cluster energy denominators.  Smaller values are replaced by this constant during the early iterations only, so the final results are unaffected, but initial convergence is improved when the guess is poor.OPTIONS: abcde Integer code is mapped to abc x 10^-de
RECOMMENDATION:

Recommendation: : Increase to 0.5 or 0.75 for non-convergent coupled-cluster calculations.    '''
        if value == "":
            if "CC_DOV_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DOV_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DOV_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DOV_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DOV_THRESH"] = value.lower()

    def cc_nguess_singles(self, value="show"):
        '''
Name: CC_NGUESS_SINGLES
Type: INTEGER
Default: 0
Options: Range from 0 to 200

Description: Specifies number of excited state guess vectors which are single excitations. 
Recommendation: : Should be greater or equal than the number of excited states requested.    '''
        if value == "":
            if "CC_NGUESS_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_NGUESS_SINGLES"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_NGUESS_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_NGUESS_SINGLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_NGUESS_SINGLES"] = value.lower()

    def cc_trans_prop(self, value="show"):
        '''
Name: CC_TRANS_PROP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_STATE_TO_OPT must be specified for this state.
    '''
        if value == "":
            if "CC_TRANS_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TRANS_PROP"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_TRANS_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TRANS_PROP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_TRANS_PROP"] = value.lower()

    def qui_multiplicity(self, value="show"):
        '''
Name: QUI_MULTIPLICITY
Type: INTEGER
Factor: 2
Default: 0 [=1]
Options: Range from 0 [=0] to 10 [=20]

Description: Sets the multiplicity of the system
    '''
        if value == "":
            if "QUI_MULTIPLICITY" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_MULTIPLICITY"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_MULTIPLICITY" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_MULTIPLICITY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_MULTIPLICITY"] = value.lower()

    def efp(self, value="show"):
        '''
Name: EFP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: The keyword should be present if excited state calculation is requested.
    '''
        if value == "":
            if "EFP" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP"]
                print("Keyword removed.")
        elif value == "show":
            if "EFP" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EFP"] = value.lower()

    def efp_fragments_only(self, value="show"):
        '''
Name: EFP_FRAGMENTS_ONLY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Set to true if there is no QM part to the calculation.
    '''
        if value == "":
            if "EFP_FRAGMENTS_ONLY" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP_FRAGMENTS_ONLY"]
                print("Keyword removed.")
        elif value == "show":
            if "EFP_FRAGMENTS_ONLY" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP_FRAGMENTS_ONLY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EFP_FRAGMENTS_ONLY"] = value.lower()

    def efp_input(self, value="show"):
        '''
Name: EFP_INPUT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: True indicates the new format without a dummy atom in the $molecule section.  False indicates the old format which requires a dummy atom (e.g. He) in the $molecule section for an EFP-only calculation.
    '''
        if value == "":
            if "EFP_INPUT" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP_INPUT"]
                print("Keyword removed.")
        elif value == "show":
            if "EFP_INPUT" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP_INPUT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EFP_INPUT"] = value.lower()

    def chemsol_read_vdw(self, value="show"):
        '''
Name: CHEMSOL_READ_VDW
Type: STRING
Default: false

Options:
    'false'......................... Default
    'true'.......................... User-defined

Description: Controls the input of user-defined atomic radii for a ChemSol calculation.
    '''
        if value == "":
            if "CHEMSOL_READ_VDW" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_READ_VDW"]
                print("Keyword removed.")
        elif value == "show":
            if "CHEMSOL_READ_VDW" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_READ_VDW"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHEMSOL_READ_VDW"] = value.lower()

    def pcm_print(self, value="show"):
        '''
Name: PCM_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 5

Description: Controls the print level during PCM calculations.
    '''
        if value == "":
            if "PCM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["PCM_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "PCM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["PCM_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PCM_PRINT"] = value.lower()

    def qui_solvent_cosmo(self, value="show"):
        '''
Name: QUI_SOLVENT_COSMO
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_SOLVENT_COSMO" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_COSMO"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_COSMO" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_COSMO"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_SOLVENT_COSMO"] = value.lower()

    def qui_solvent_pcm(self, value="show"):
        '''
Name: QUI_SOLVENT_PCM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Use an apparent surface charge polarizable continuum solvent model.
    '''
        if value == "":
            if "QUI_SOLVENT_PCM" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_PCM"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_PCM" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_PCM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_SOLVENT_PCM"] = value.lower()

    def solvent_method(self, value="show"):
        '''
Name: SOLVENT_METHOD
Type: STRING
Default: SCRF

Options:
    'SCRF'.......................... SCRF
    'PCM'........................... PCM
    'COSMO'......................... COSMO

Description: Sets the preferred solvent model.
    '''
        if value == "":
            if "SOLVENT_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLVENT_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "SOLVENT_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLVENT_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SOLVENT_METHOD"] = value.lower()

    def sol_order(self, value="show"):
        '''
Name: SOL_ORDER
Type: INTEGER
Default: 15
Options: Range from 1 to 25

Description: Determines the order to which the multipole expansion of the solute charge density is carried out.
    '''
        if value == "":
            if "SOL_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["SOL_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "SOL_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["SOL_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SOL_ORDER"] = value.lower()

    def svp(self, value="show"):
        '''
Name: SVP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Sets whether to perform the isodensity solvation procedure.
    '''
        if value == "":
            if "SVP" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP"] = value.lower()

    def svp_charge_conv(self, value="show"):
        '''
Name: SVP_CHARGE_CONV
Type: INTEGER

Description: Determines the convergence value for the charges on the cavity. When the change in charges fall below this value, if the electron density is converged, then the calculation is considered converged.
Recommendation: : The default value unless convergence problems arise.    '''
        if value == "":
            if "SVP_CHARGE_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_CHARGE_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP_CHARGE_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_CHARGE_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP_CHARGE_CONV"] = value.lower()

    def svp_guess(self, value="show"):
        '''
Name: SVP_GUESS
Type: STRING
Default: 0

Options:
    '0'............................. No Guess
    '1'............................. Read
    '2'............................. Specify

Description: Specifies how and if the solvation module will use a given guess for the charges and cavity points.
Recommendation: : It is helpful to also set SCF_GUESS to READ when using a guess from a previous Q-Chem run.     '''
        if value == "":
            if "SVP_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_GUESS"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_GUESS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP_GUESS"] = value.lower()

    def svp_memory(self, value="show"):
        '''
Name: SVP_MEMORY
Type: INTEGER
Default: 125
Options: Range from 32 to 2048

Description: Specifies the amount of memory for use by the solvation module.
Recommendation: :     '''
        if value == "":
            if "SVP_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_MEMORY"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_MEMORY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP_MEMORY"] = value.lower()

    def svp_path(self, value="show"):
        '''
Name: SVP_PATH
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Specifies whether to run a gas phase computation prior to performing the solvation procedure.
Recommendation: : Running the gas-phase calculation provides a good guess to start the solvation stage and provides a more complete set of solvated properties.    '''
        if value == "":
            if "SVP_PATH" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_PATH"]
                print("Keyword removed.")
        elif value == "show":
            if "SVP_PATH" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_PATH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SVP_PATH"] = value.lower()

    def symmetry_decomposition(self, value="show"):
        '''
Name: SYMMETRY_DECOMPOSITION
Type: INTEGER

Description: Determines symmetry decompositions to calculate.
    '''
        if value == "":
            if "SYMMETRY_DECOMPOSITION" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_DECOMPOSITION"]
                print("Keyword removed.")
        elif value == "show":
            if "SYMMETRY_DECOMPOSITION" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_DECOMPOSITION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SYMMETRY_DECOMPOSITION"] = value.lower()

    def qui_solvent_dielectric_cosmo(self, value="show"):
        '''
Name: QUI_SOLVENT_DIELECTRIC_COSMO
Type: INTEGER
Factor: 0.0001
Default: 0 [=0.0000]
Options: Range from 0 [=0.0000] to 999998 [=99.9999]

Description: Sets the dielectric constant of the solvent
    '''
        if value == "":
            if "QUI_SOLVENT_DIELECTRIC_COSMO" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_COSMO"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_DIELECTRIC_COSMO" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_COSMO"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "QUI_SOLVENT_DIELECTRIC_COSMO"] = value.lower()

    def qui_solvent_dielectric_onsager(self, value="show"):
        '''
Name: QUI_SOLVENT_DIELECTRIC_ONSAGER
Type: INTEGER
Factor: 0.0001
Default: 0 [=0.0000]
Options: Range from 0 [=0.0000] to 999998 [=99.9999]

Description: Sets the dielectric constant of the solvent
    '''
        if value == "":
            if "QUI_SOLVENT_DIELECTRIC_ONSAGER" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_ONSAGER"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SOLVENT_DIELECTRIC_ONSAGER" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_ONSAGER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "QUI_SOLVENT_DIELECTRIC_ONSAGER"] = value.lower()

    def cholesky_tol(self, value="show"):
        '''
Name: CHOLESKY_TOL
Type: INTEGER
Default: 3
Options: Range from 0 to 16

Description: Tolerance for the Cholesky decomposition of two-electron integrals.

Recommendation: :
2 - qualitative calculations
3 - appropriate for most cases
4 - quantatative (error < 10-6 Eh)    '''
        if value == "":
            if "CHOLESKY_TOL" in self.dict_of_keywords:
                del self.dict_of_keywords["CHOLESKY_TOL"]
                print("Keyword removed.")
        elif value == "show":
            if "CHOLESKY_TOL" in self.dict_of_keywords:
                return self.dict_of_keywords["CHOLESKY_TOL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHOLESKY_TOL"] = value.lower()

    def qui_integral_decomposition_none(self, value="show"):
        '''
Name: QUI_INTEGRAL_DECOMPOSITION_NONE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: No integral decomposition
    '''
        if value == "":
            if "QUI_INTEGRAL_DECOMPOSITION_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_INTEGRAL_DECOMPOSITION_NONE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_INTEGRAL_DECOMPOSITION_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_INTEGRAL_DECOMPOSITION_NONE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "QUI_INTEGRAL_DECOMPOSITION_NONE"] = value.lower()

    def direct_ri(self, value="show"):
        '''
Name: DIRECT_RI
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
Recommendation: :
By default, all integrals are used in decomposed format allowing significant reduction of memory use.  If all integrals are transformed back (TRUE option) no memory reduction is achieved and decompostion error is introduced.  However, the integral transformation is performed significantly faster and conventional CC/EOM algorithms are used.    '''
        if value == "":
            if "DIRECT_RI" in self.dict_of_keywords:
                del self.dict_of_keywords["DIRECT_RI"]
                print("Keyword removed.")
        elif value == "show":
            if "DIRECT_RI" in self.dict_of_keywords:
                return self.dict_of_keywords["DIRECT_RI"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIRECT_RI"] = value.lower()

    def print_general_basis(self, value="show"):
        '''
Name: PRINT_GENERAL_BASIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls print out of built in basis sets in input format
Recommendation: : Useful for modification of standard basis sets.    '''
        if value == "":
            if "PRINT_GENERAL_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_GENERAL_BASIS"]
                print("Keyword removed.")
        elif value == "show":
            if "PRINT_GENERAL_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_GENERAL_BASIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PRINT_GENERAL_BASIS"] = value.lower()

    def print_orbitals(self, value="show"):
        '''
Name: PRINT_ORBITALS
Type: STRING
Default: 0

Options:
    '0'............................. 0
    '99'............................ 99
    '0'............................. 0
    '1'............................. 1

Description: Prints orbitals coefficients with atom labels in the analysis part of the output
Recommendation: : Usually not required as the orbitals are more easily examined visually via the formatted checkpoint file.    '''
        if value == "":
            if "PRINT_ORBITALS" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_ORBITALS"]
                print("Keyword removed.")
        elif value == "show":
            if "PRINT_ORBITALS" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_ORBITALS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["PRINT_ORBITALS"] = value.lower()

    def scf_algorithm(self, value="show"):
        '''
Name: SCF_ALGORITHM
Type: STRING
Default: DIIS

Options:
    'DIIS'.......................... DIIS
    'DM'............................ DM
    'GDM'........................... GDM
    'RCA'........................... RCA
    'ROOTHAAN'...................... ROOTHAAN
    'DIIS_DM'....................... DIIS_DM
    'DIIS_GDM'...................... DIIS_GDM
    'RCA_DIIS'...................... RCA_DIIS

Description: Selects the algorithm to use for converging the SCF.
Recommendation: : Use DIIS unless performing a restricted open-shell calculation, in which case GDM is recommended. If DIIS fails to find a reasonable approximate solution in the initial iterations, RCA_DIIS is the recommended fallback option. If DIIS approaches the correct solution but fails to finally converge, DIIS_GDM is the recommended fallback.     '''
        if value == "":
            if "SCF_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_ALGORITHM"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_ALGORITHM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_ALGORITHM"] = value.lower()

    def thresh(self, value="show"):
        '''
Name: THRESH
Type: INTEGER
Default: 8
Options: Range from 0 to 14

Description: Cutoff for neglect of two electron integrals. 10-THRESH (THRESH ? 14).
Recommendation: : Should be at least three greater than the SCF convergence setting. Increase for more significant figures, at greater computational cost.    '''
        if value == "":
            if "THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["THRESH"] = value.lower()

    def basis_projection_type(self, value="show"):
        '''
Name: BASIS_PROJECTION_TYPE
Type: STRING
Default: FOPPROJECTION

Options:
    'FOPPROJECTION'................. Fock Matrix
    'OVPROJECTION'.................. MO Overlap

Description: Determines which method to use when projecting the density matrix for the basis set projection guess.
    '''
        if value == "":
            if "BASIS_PROJECTION_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS_PROJECTION_TYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "BASIS_PROJECTION_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS_PROJECTION_TYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["BASIS_PROJECTION_TYPE"] = value.lower()

    def diis_max_cycles(self, value="show"):
        '''
Name: DIIS_MAX_CYCLES
Type: INTEGER
Default: 50
Options: Range from 1 to 100

Description: The maximum number of DIIS iterations before switching to (geometric) direct minimization when SCF_ALGORITHM is DIIS_GDM or DIIS_DM. See also THRESH_DIIS_SWITCH. 
    '''
        if value == "":
            if "DIIS_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_MAX_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_MAX_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_MAX_CYCLES"] = value.lower()

    def diis_separate_errvec(self, value="show"):
        '''
Name: DIIS_SEPARATE_ERRVEC
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the optimization of DIIS error vectors in unrestricted caclulations.  When using DIIS in QChem a convenient optimization for unrestricted calculations is to sume the alpha and beta error vectors into a single vector which is used for extrapolation.  This is often extremely effective, but in some pathological systems with symmetry breaking, can lead to false solutions being detected, where the alpha and beta components of the error vector cancel exactly giving a zero DIIS error.  

Recommendation: : While an extremely uncommon 
occurrence, if cancelation is suspected, set to TRUE to check.    '''
        if value == "":
            if "DIIS_SEPARATE_ERRVEC" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_SEPARATE_ERRVEC"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_SEPARATE_ERRVEC" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_SEPARATE_ERRVEC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_SEPARATE_ERRVEC"] = value.lower()

    def diis_subspace_size(self, value="show"):
        '''
Name: DIIS_SUBSPACE_SIZE
Type: INTEGER
Default: 15
Options: Range from 1 to 50

Description: Controls the size of the DIIS and/or RCA subspace during the SCF.
    '''
        if value == "":
            if "DIIS_SUBSPACE_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_SUBSPACE_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_SUBSPACE_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_SUBSPACE_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_SUBSPACE_SIZE"] = value.lower()

    def diis_switch_thresh(self, value="show"):
        '''
Name: DIIS_SWITCH_THRESH
Type: INTEGER

Description: The threshold for switching between DIIS extrapolation and direct minimization of the SCF energy is 10-THRESH_DIIS_SWITCH when SCF_ALGORITHM is DIIS_GDM or DIIS_DM. See also MAX_DIIS_MAX_CYCLES
    '''
        if value == "":
            if "DIIS_SWITCH_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_SWITCH_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_SWITCH_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_SWITCH_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_SWITCH_THRESH"] = value.lower()

    def direct_scf(self, value="show"):
        '''
Name: DIRECT_SCF
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls direct SCF, the default is determined by the program.
Recommendation: : Use default; direct SCF switches off in-core integrals.    '''
        if value == "":
            if "DIRECT_SCF" in self.dict_of_keywords:
                del self.dict_of_keywords["DIRECT_SCF"]
                print("Keyword removed.")
        elif value == "show":
            if "DIRECT_SCF" in self.dict_of_keywords:
                return self.dict_of_keywords["DIRECT_SCF"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIRECT_SCF"] = value.lower()

    def meteco(self, value="show"):
        '''
Name: METECO
Type: STRING
Default: 1

Options:
    '1'............................. Machine precision
    '2'............................. Integral Thresh

Description: Sets the threshold criteria for discarding shell-pairs.
Recommendation: : Use default.    '''
        if value == "":
            if "METECO" in self.dict_of_keywords:
                del self.dict_of_keywords["METECO"]
                print("Keyword removed.")
        elif value == "show":
            if "METECO" in self.dict_of_keywords:
                return self.dict_of_keywords["METECO"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["METECO"] = value.lower()

    def symmetry_ignore(self, value="show"):
        '''
Name: SYMMETRY_IGNORE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether or not Q-Chem determines the point group of the molecule and reorients the molecule to the standard orientation.
Recommendation: : Use default unless you do not want the molecule to be reoriented. Note that symmetry usage is disabled for RIMP2 jobs.    '''
        if value == "":
            if "SYMMETRY_IGNORE" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_IGNORE"]
                print("Keyword removed.")
        elif value == "show":
            if "SYMMETRY_IGNORE" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_IGNORE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SYMMETRY_IGNORE"] = value.lower()

    def symmetry_integral(self, value="show"):
        '''
Name: SYMMETRY_INTEGRAL
Type: LOGICAL
Default: TRUE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the efficiency through the use of point group symmetry for calculating integrals.
Recommendation: : Use default unless benchmarking. Note that symmetry usage is disabled for RIMP2, FFT and QM/MM jobs.    '''
        if value == "":
            if "SYMMETRY_INTEGRAL" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_INTEGRAL"]
                print("Keyword removed.")
        elif value == "show":
            if "SYMMETRY_INTEGRAL" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_INTEGRAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SYMMETRY_INTEGRAL"] = value.lower()

    def symmetry_tolerance(self, value="show"):
        '''
Name: SYMMETRY_TOLERANCE
Type: INTEGER
Default: 5
Options: Range from 0 to 10

Description: Controls the tolerance for determining point group symmetry. Differences in atom locations less than 10-SYM_TOL are treated as zero.
Recommendation: : Use the default unless the molecule has high symmetry which is not being correctly identified. Note that relaxing this tolerance too much may introduce errors into the calculation.    '''
        if value == "":
            if "SYMMETRY_TOLERANCE" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_TOLERANCE"]
                print("Keyword removed.")
        elif value == "show":
            if "SYMMETRY_TOLERANCE" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_TOLERANCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SYMMETRY_TOLERANCE"] = value.lower()

    def dftvdw_alpha1(self, value="show"):
        '''
Name: DFTVDW_ALPHA1
Type: INTEGER
Factor: 10
Default: 8 [=83]
Options: Range from 1 [=10] to 100 [=1000]

Description: Parameter in XDM calculations with higher-order terms
    '''
        if value == "":
            if "DFTVDW_ALPHA1" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_ALPHA1"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_ALPHA1" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_ALPHA1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_ALPHA1"] = value.lower()

    def dftvdw_alpha2(self, value="show"):
        '''
Name: DFTVDW_ALPHA2
Type: INTEGER
Factor: 10
Default: 15 [=155]
Options: Range from 1 [=10] to 100 [=1000]

Description: Parameter in XDM calculations with higher-order terms.
    '''
        if value == "":
            if "DFTVDW_ALPHA2" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_ALPHA2"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_ALPHA2" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_ALPHA2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_ALPHA2"] = value.lower()

    def dftvdw_jobnumber(self, value="show"):
        '''
Name: DFTVDW_JOBNUMBER
Type: STRING
Default: 0

Options:
    '0'............................. None
    '1'............................. Post SCF
    '2'............................. Full SCF

Description: Basic vdW job control
    '''
        if value == "":
            if "DFTVDW_JOBNUMBER" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_JOBNUMBER"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_JOBNUMBER" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_JOBNUMBER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_JOBNUMBER"] = value.lower()

    def dftvdw_kai(self, value="show"):
        '''
Name: DFTVDW_KAI
Type: INTEGER
Factor: 10
Default: 80 [=800]
Options: Range from 1 [=10] to 100 [=1000]

Description: Damping factor K for C6 only damping functions
    '''
        if value == "":
            if "DFTVDW_KAI" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_KAI"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_KAI" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_KAI"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_KAI"] = value.lower()

    def dftvdw_mol1natoms(self, value="show"):
        '''
Name: DFTVDW_MOL1NATOMS
Type: INTEGER
Default: 0
Options: Range from 0 to 9999

Description: The number of atoms in the first monomer in a dimer calculation.
    '''
        if value == "":
            if "DFTVDW_MOL1NATOMS" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_MOL1NATOMS"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_MOL1NATOMS" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_MOL1NATOMS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_MOL1NATOMS"] = value.lower()

    def dftvdw_print(self, value="show"):
        '''
Name: DFTVDW_PRINT
Type: STRING
Default: 1

Options:
    '0'............................. None
    '1'............................. Minimal
    '2'............................. Debug

Description: Controls the amount of output from the VDW code.
    '''
        if value == "":
            if "DFTVDW_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_PRINT"] = value.lower()

    def lrc_dft(self, value="show"):
        '''
Name: LRC_DFT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the application of long-range-corrected DFT

Recommendation: :  Long-range correction is available for any combination of HF, B88 and PBE exchange.    '''
        if value == "":
            if "LRC_DFT" in self.dict_of_keywords:
                del self.dict_of_keywords["LRC_DFT"]
                print("Keyword removed.")
        elif value == "show":
            if "LRC_DFT" in self.dict_of_keywords:
                return self.dict_of_keywords["LRC_DFT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["LRC_DFT"] = value.lower()

    def omega(self, value="show"):
        '''
Name: OMEGA
Type: INTEGER
Factor: 0.001
Default: 200 [=0.200]
Options: Range from 1 [=0.001] to 9999 [=9.999]

Description: Controls the degree of attenuation of the Coulomb operator.
    '''
        if value == "":
            if "OMEGA" in self.dict_of_keywords:
                del self.dict_of_keywords["OMEGA"]
                print("Keyword removed.")
        elif value == "show":
            if "OMEGA" in self.dict_of_keywords:
                return self.dict_of_keywords["OMEGA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["OMEGA"] = value.lower()

    def dftvdw_method(self, value="show"):
        '''
Name: DFTVDW_METHOD
Type: STRING
Default: 1

Options:
    '1'............................. C6 Term Only
    '2'............................. Include C8 & C10

Description: Selects the damping function used in XDM
    '''
        if value == "":
            if "DFTVDW_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["DFTVDW_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "DFTVDW_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["DFTVDW_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFTVDW_METHOD"] = value.lower()

    def dft_d(self, value="show"):
        '''
Name: DFT_D
Type: STRING
Default: 0

Options:
    '0'............................. None
    'EMPIRICAL_GRIMME'.............. Grimme
    'EMPIRICAL_GRIMME3'............. DFT-D3
    'EMPIRICAL_CHG'................. Chai Head-Gordon

Description: Specifies what dispersion correction to use within a DFT calculation.
    '''
        if value == "":
            if "DFT_D" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D"]
                print("Keyword removed.")
        elif value == "show":
            if "DFT_D" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFT_D"] = value.lower()

    def dft_d3_s6(self, value="show"):
        '''
Name: DFT_D3_S6
Type: INTEGER
Factor: 0.001
Default: 1000 [=1.000]
Options: Range from 0 [=0] to 9999 [=9.999]

Description: Controls the strength of dispersion corrections, s6, in Grimme?s DFT-D3 method
    '''
        if value == "":
            if "DFT_D3_S6" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D3_S6"]
                print("Keyword removed.")
        elif value == "show":
            if "DFT_D3_S6" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D3_S6"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFT_D3_S6"] = value.lower()

    def dft_d3_s8(self, value="show"):
        '''
Name: DFT_D3_S8
Type: INTEGER
Factor: 0.001
Default: 1000 [=1.000]
Options: Range from 0 [=0] to 9999 [=9.999]

Description: Controls the strength of dispersion corrections, s8 , in Grimme?s DFT-D3 method.
    '''
        if value == "":
            if "DFT_D3_S8" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D3_S8"]
                print("Keyword removed.")
        elif value == "show":
            if "DFT_D3_S8" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D3_S8"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFT_D3_S8"] = value.lower()

    def dft_d3_sr6(self, value="show"):
        '''
Name: DFT_D3_SR6
Type: INTEGER
Factor: 0.001
Default: 1000 [=1.000]
Options: Range from 0 [=0] to 9999 [=9.999]

Description: Controls the strength of dispersion corrections, sr6 , in the Grimme?s DFT-D3 method.
    '''
        if value == "":
            if "DFT_D3_SR6" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D3_SR6"]
                print("Keyword removed.")
        elif value == "show":
            if "DFT_D3_SR6" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D3_SR6"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFT_D3_SR6"] = value.lower()

    def nl_correlation(self, value="show"):
        '''
Name: NL_CORRELATION
Type: STRING
Default: None

Options:
    'None'.......................... None
    'vdW-DF-04'..................... vdW-DF-04
    'vdW-DF-10'..................... vdW-DF-10
    'VV09'.......................... VV09
    'VV10'.......................... VV10

Description: Speci?es a non-local correlation functional that includes non-empirical dispersion.

Recommendation: :  Do not forget to add the LSDA correlation (PW92 is recommended) when using vdW-DF-04, vdW-DF-10, or VV09. VV10 should be used with PBE correlation. Choose exchange functionals carefully: HF, rPW86, revPBE, and some of the LRC exchange functionals are among the recommended choices.    '''
        if value == "":
            if "NL_CORRELATION" in self.dict_of_keywords:
                del self.dict_of_keywords["NL_CORRELATION"]
                print("Keyword removed.")
        elif value == "show":
            if "NL_CORRELATION" in self.dict_of_keywords:
                return self.dict_of_keywords["NL_CORRELATION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NL_CORRELATION"] = value.lower()

    def nl_grid(self, value="show"):
        '''
Name: NL_GRID
Type: STRING
Default: 0

Options:
    '0'............................. SG-0
    '1'............................. SG-1
    '2'............................. El-cheeso

    '35,110'........................ 35,110
    '50,194'........................ 50,194
    '75,302'........................ 75,302
    '99,590'........................ 99,590

Description: Speci?es the grid to use for non-local correlation.

Recommendation: : Use default unless computational cost becomes prohibitive, in which case SG-0 may be 
used.  XC_GRID should generally be ?ner than NL_GRID.    '''
        if value == "":
            if "NL_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["NL_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "NL_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["NL_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NL_GRID"] = value.lower()

    def nl_vv_b(self, value="show"):
        '''
Name: NL_VV_B
Type: INTEGER
Factor: 0.01
Default: 590 [=5.90]
Options: Range from 1 [=0.01] to 9999 [=99.99]

Description: Sets the parameter b in VV10. This parameter controls the short range behavior of the nonlocal correlation energy.

Recommendation: : The optimal value depends strongly on the exchange functional used. b = 5.9 is recommended for rPW86.    '''
        if value == "":
            if "NL_VV_B" in self.dict_of_keywords:
                del self.dict_of_keywords["NL_VV_B"]
                print("Keyword removed.")
        elif value == "show":
            if "NL_VV_B" in self.dict_of_keywords:
                return self.dict_of_keywords["NL_VV_B"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NL_VV_B"] = value.lower()

    def nl_vv_c(self, value="show"):
        '''
Name: NL_VV_C
Type: INTEGER
Factor: 0.00001
Default: 88 [=0.00089]
Options: Range from 1 [=0.00001] to 99999 [=1.00000]

Description: Sets the parameter C in VV09 and VV10.  This parameter is fitted to asymptotic van der Waals C6 coefficients.

Recommendation: : C = 0.0093 is recommended when a semilocal exchange functional is used. C = 0.0089 
is recommended when a long-range corrected (LRC) hybrid functional is used.    '''
        if value == "":
            if "NL_VV_C" in self.dict_of_keywords:
                del self.dict_of_keywords["NL_VV_C"]
                print("Keyword removed.")
        elif value == "show":
            if "NL_VV_C" in self.dict_of_keywords:
                return self.dict_of_keywords["NL_VV_C"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NL_VV_C"] = value.lower()

    def local_interp_order(self, value="show"):
        '''
Name: LOCAL_INTERP_ORDER
Type: INTEGER
Default: 6
Options: Range from 1 to 10

Description: Controls the order of the B-spline.
    '''
        if value == "":
            if "LOCAL_INTERP_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["LOCAL_INTERP_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "LOCAL_INTERP_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["LOCAL_INTERP_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["LOCAL_INTERP_ORDER"] = value.lower()

    def mrxc(self, value="show"):
        '''
Name: MRXC
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the use of multi-resolution XC.

Recommendation: : MRXC is very ef?cient for medium and large molecules, especially when medium and 
large basis sets are used.    '''
        if value == "":
            if "MRXC" in self.dict_of_keywords:
                del self.dict_of_keywords["MRXC"]
                print("Keyword removed.")
        elif value == "show":
            if "MRXC" in self.dict_of_keywords:
                return self.dict_of_keywords["MRXC"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MRXC"] = value.lower()

    def mrxc_class_thresh_mult(self, value="show"):
        '''
Name: MRXC_CLASS_THRESH_MULT
Type: INTEGER
Default: 1
Options: Range from 1 to 9

Description: Controls the prefactor of the smoothness precision.
    '''
        if value == "":
            if "MRXC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                del self.dict_of_keywords["MRXC_CLASS_THRESH_MULT"]
                print("Keyword removed.")
        elif value == "show":
            if "MRXC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                return self.dict_of_keywords["MRXC_CLASS_THRESH_MULT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MRXC_CLASS_THRESH_MULT"] = value.lower()

    def mrxc_class_thresh_order(self, value="show"):
        '''
Name: MRXC_CLASS_THRESH_ORDER
Type: INTEGER
Default: 6
Options: Range from 1 to 9

Description: Controls the exponent of the smoothness precision.
    '''
        if value == "":
            if "MRXC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["MRXC_CLASS_THRESH_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "MRXC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["MRXC_CLASS_THRESH_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MRXC_CLASS_THRESH_ORDER"] = value.lower()

    def exchange(self, value="show"):
        '''
Name: EXCHANGE
Type: STRING
Default: HF

Options:
    'HF'............................ HF

    'Slater'........................ Slater

    'B86'........................... B86
    'B88'........................... B88
    'B97-D'......................... B97-D
    'GG99'.......................... GG99
    'Gill'.......................... Gill96
    'PBE'........................... PBE
    'PW86'.......................... PW86
    'PW91'.......................... PW91
    'revPBE'........................ revPBE
    'rPW86'......................... rPW86

    'BOP'........................... BOP
    'EDF1'.......................... EDF1
    'PBEOP'......................... PBEOP
    'SOGGA'......................... SOGGA
    'SOGGA11'....................... SOGGA11

    'B3LYP'......................... B3LYP
    'B3LYP5'........................ B3LYP5
    'B3P86'......................... B3P86
    'B3PW91'........................ B3PW91
    'B5050LYP'...................... B5050LYP
    'B97'........................... B97
    'B97-1'......................... B97-1
    'B97-2'......................... B97-2
    'BHHLYP'........................ BHHLYP
    'BR89B94h'...................... BR89B94h
    'EDF2'.......................... EDF2
    'HCTH'.......................... HCTH
    'HCTH-120'...................... HCTH-120
    'HCTH-147'...................... HCTH-147
    'HCTH-407'...................... HCTH-407
    'mPW1B95'....................... mPW1B95
    'mPW1LYP'....................... mPW1LYP
    'mPW1PBE'....................... mPW1PBE
    'mPW1PW91'...................... mPW1PW91
    'mPWB1K'........................ mPWB1K
    'PBE1PBE'....................... PBE0
    'PBE50'......................... PBE50
    'SOGGA11X'...................... SOGGA11X
    'X3LYP'......................... X3LYP

    'CAM-B3LYP'..................... CAM-B3LYP
    'LRC-wPBEPBE'................... LRC-wPBEPBE
    'LRC-wPBEhPBE'.................. LRC-wPBEhPBE
    'muB88'......................... muB88
    'muPBE'......................... muPBE
    'omegaB97'...................... omegaB97
    'omegaB97X'..................... omegaB97X
    'omegaB97X-D'................... omegaB97X-D
    'omegaB97X-2(LP)'............... omegaB97X-2(LP)
    'omegaB97X-2(TQZ)'.............. omegaB97X-2(TQZ)
    'wPBE'.......................... wPBE

    'B1B95'......................... B1B95
    'BMK'........................... BMK
    'BR89'.......................... BR89
    'M05'........................... M05
    'M052X'......................... M052X
    'M06'........................... M06
    'M06L'.......................... M06L
    'M06HF'......................... M06HF
    'M062X'......................... M062X
    'M08HX'......................... M08HX
    'M08SO'......................... M08SO
    'M11'........................... M11
    'M11L'.......................... M11L
    'TPSS'.......................... TPSS
    'TPSSH'......................... TPSSH
    'VSXC'.......................... VSXC

    'B05'........................... B05
    'B3tLAP'........................ B3tLAP
    'BM05'.......................... BM05
    'MCY2'.......................... MCY2
    'PSTS'.......................... PSTS

    'LXYGJOS'....................... LXYGJOS
    'XYG3'.......................... XYG3
    'XYGJOS'........................ XYGJOS

    'gen'........................... User-defined

Description: Specifies the exchange level of theory.
Recommendation: : Consult the literature and reviews for guidence    '''
        if value == "":
            if "EXCHANGE" in self.dict_of_keywords:
                del self.dict_of_keywords["EXCHANGE"]
                print("Keyword removed.")
        elif value == "show":
            if "EXCHANGE" in self.dict_of_keywords:
                return self.dict_of_keywords["EXCHANGE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EXCHANGE"] = value.lower()

    def dfpt_exchange(self, value="show"):
        '''
Name: DFPT_EXCHANGE
Type: STRING
Default: B3LYP

Options:
    'B3LYP'......................... B3LYP
    'B3LYP5'........................ B3LYP5
    'B3P86'......................... B3P86
    'B3PW91'........................ B3PW91
    'B5050LYP'...................... B5050LYP
    'B97'........................... B97
    'B98-1'......................... B98-1
    'B97-2'......................... B97-2
    'BHHLYP'........................ BHHLYP
    'BR89B96h'...................... BR89B96h
    'EFD2'.......................... EFD2
    'HCTH'.......................... HCTH
    'HCTH-120'...................... HCTH-120
    'HCTH-147'...................... HCTH-147
    'HCTH-407'...................... HCTH-407
    'mPW1B95'....................... mPW1B95
    'mPW1LYP'....................... mPW1LYP
    'mPW1PBE'....................... mPW1PBE
    'mPW1PW91'...................... mPW1PW91
    'mPWB1K'........................ mPWB1K
    'PBE0'.......................... PBE0
    'SOGGA11X'...................... SOGGA11X
    'X3LYP'......................... X3LYP

    'B1B95'......................... B1B95
    'BMK'........................... BMK
    'BR89'.......................... BR89
    'M05'........................... M05
    'M052X'......................... M052X
    'M06'........................... M06
    'M06L'.......................... M06L
    'M06HF'......................... M06HF
    'M062X'......................... M062X
    'M08HX'......................... M08HX
    'M08SO'......................... M08SO
    'M11'........................... M11
    'M11L'.......................... M11L
    'TPSS'.......................... TPSS
    'TPSSH'......................... TPSSH

    'B05'........................... B05
    'B3tLAP'........................ B3tLAP
    'BM05'.......................... BM05
    'MCY2'.......................... MCY2
    'PSTS'.......................... PSTS

    'LXYGJOS'....................... LXYGJOS
    'XYG3'.......................... XYG3
    'XYGJOS'........................ XYGJOS

    'gen'........................... User-defined

Description: Specifies the secondary functional in a density functional perturbation theory (DFPT) calculation.

    '''
        if value == "":
            if "DFPT_EXCHANGE" in self.dict_of_keywords:
                del self.dict_of_keywords["DFPT_EXCHANGE"]
                print("Keyword removed.")
        elif value == "show":
            if "DFPT_EXCHANGE" in self.dict_of_keywords:
                return self.dict_of_keywords["DFPT_EXCHANGE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFPT_EXCHANGE"] = value.lower()

    def dfpt_xc_grid(self, value="show"):
        '''
Name: DFPT_XC_GRID
Type: STRING
Default: 1

Options:
    '1'............................. SG-1
    '75,302'........................ 75,302
    '99,590'........................ 99,590

Description: Specifies the secondary grid in a density functional perturbation theory (DFPT) calculation.

Recommendation: : This should be at least as large as the XC_GRID    '''
        if value == "":
            if "DFPT_XC_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["DFPT_XC_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "DFPT_XC_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["DFPT_XC_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DFPT_XC_GRID"] = value.lower()

    def hfpt(self, value="show"):
        '''
Name: HFPT
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: This is set implicitly by setting HFPT_BASIS
    '''
        if value == "":
            if "HFPT" in self.dict_of_keywords:
                del self.dict_of_keywords["HFPT"]
                print("Keyword removed.")
        elif value == "show":
            if "HFPT" in self.dict_of_keywords:
                return self.dict_of_keywords["HFPT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["HFPT"] = value.lower()

    def xc_grid(self, value="show"):
        '''
Name: XC_GRID
Type: STRING
Default: 1

Options:
    '0'............................. SG-0
    '1'............................. SG-1
    '2'............................. El-cheeso

    '35,110'........................ 35,110
    '50,194'........................ 50,194
    '75,302'........................ 75,302
    '99,590'........................ 99,590

Description: Specifies the quadrature grid to be used for evaluating the exchange-correlation component of the energy.  Either a standard grid should be selected, or a Lebedev grid with the corresponding number of points.

Recommendation: : Use the default unless convergence difficulties arise.  Larger grids are required for calculations involving derivatives and excited states.    '''
        if value == "":
            if "XC_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["XC_GRID"]
                print("Keyword removed.")
        elif value == "show":
            if "XC_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["XC_GRID"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XC_GRID"] = value.lower()

    def auxiliary_basis(self, value="show"):
        '''
Name: AUXILIARY_BASIS
Type: STRING
Default: None

Options:
    'None'.......................... None
    'RIMP2-VDZ'..................... RIMP2-VDZ
    'RIMP2-TZVPP'................... RIMP2-TZVPP
    'RIMP2-cc-PVDZ'................. RIMP2-cc-PVDZ
    'RIMP2-cc-PVTZ'................. RIMP2-cc-PVTZ
    'RIMP2-cc-PVQZ'................. RIMP2-cc-PVQZ
    'RIMP2-aug-cc-PVDZ'............. RIMP2-aug-cc-PVDZ
    'RIMP2-aug-cc-PVTZ'............. RIMP2-aug-cc-PVTZ
    'RIMP2-aug-cc-PVQZ'............. RIMP2-aug-cc-PVQZ

Description: Specifies the auxiliary basis to be used in a RI-MP2 calculation
    '''
        if value == "":
            if "AUXILIARY_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["AUXILIARY_BASIS"]
                print("Keyword removed.")
        elif value == "show":
            if "AUXILIARY_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["AUXILIARY_BASIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AUXILIARY_BASIS"] = value.lower()

    def basis(self, value="show"):
        '''
Name: BASIS
Type: STRING
Default: 6-31G

Options:
    'STO-3G'........................ STO-3G
    'STO-6G'........................ STO-6G

    '3-21G'......................... 3-21G
    '4-31G'......................... 4-31G
    '6-31G'......................... 6-31G
    '6-31G*'........................ 6-31G*
    '6-31+G*'....................... 6-31+G*
    '6-31G**'....................... 6-31G**
    '6-31++G**'..................... 6-31++G**
    '6-311G'........................ 6-311G
    '6-311G*'....................... 6-311G*
    '6-311+G*'...................... 6-311+G*
    '6-311G**'...................... 6-311G**
    '6-311++G**'.................... 6-311++G**
    '6-311++G(3df,3pd)'............. 6-311++G(3df,3pd)

    'pc-0'.......................... pc-0
    'pc-1'.......................... pc-1
    'pc-2'.......................... pc-2
    'pc-3'.......................... pc-3
    'pc-4'.......................... pc-4
    'pcJ-0'......................... pcJ-0
    'pcJ-1'......................... pcJ-1
    'pcJ-2'......................... pcJ-2
    'pcJ-3'......................... pcJ-3
    'pcJ-4'......................... pcJ-4
    'pcS-0'......................... pcS-0
    'pcS-1'......................... pcS-1
    'pcS-2'......................... pcS-2
    'pcS-3'......................... pcS-3
    'pcS-4'......................... pcS-4

    'cc-pVDZ'....................... cc-pVDZ
    'cc-pVTZ'....................... cc-pVTZ
    'cc-pVQZ'....................... cc-pVQZ
    'cc-pcVDZ'...................... cc-pcVDZ
    'cc-pcVTZ'...................... cc-pcVTZ
    'cc-pcVQZ'...................... cc-pcVQZ
    'aug-cc-pVDZ'................... aug-cc-pVDZ
    'aug-cc-pVTZ'................... aug-cc-pVTZ
    'aug-cc-pVQZ'................... aug-cc-pVQZ
    'aug-cc-pcVDZ'.................. aug-cc-pcVDZ
    'aug-cc-pcVTZ'.................. aug-cc-pcVTZ
    'aug-cc-pcVQZ'.................. aug-cc-pcVQZ

    'G3Large'....................... G3Large
    'G3MP2Large'.................... G3MP2Large

    'CRENBL'........................ CRENBL
    'CRENBS'........................ CRENBS
    'HWMB'.......................... HWMB
    'HWVDZ'......................... HWVDZ
    'LACVP'......................... LACVP
    'LANL2DZ'....................... LANL2DZ
    'SBKJC'......................... SBKJC
    'SRLC'.......................... SRLC
    'SRSC'.......................... SRSC

    'gen'........................... User-defined
    'Mixed'......................... Mixed

Description: Specifies the basis sets to be used.
Recommendation: : Consult literature and reviews to aid your selection.    '''
        if value == "":
            if "BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS"]
                print("Keyword removed.")
        elif value == "show":
            if "BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["BASIS"] = value.lower()

    def chelpg(self, value="show"):
        '''
Name: CHELPG
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the calculation of CHELPG charges. 

Recommendation: : For large molecules there is some overhead associated with computing CHELPG charges, especially if the number of grid points is large.    '''
        if value == "":
            if "CHELPG" in self.dict_of_keywords:
                del self.dict_of_keywords["CHELPG"]
                print("Keyword removed.")
        elif value == "show":
            if "CHELPG" in self.dict_of_keywords:
                return self.dict_of_keywords["CHELPG"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CHELPG"] = value.lower()

    def cis_moments(self, value="show"):
        '''
Name: CIS_MOMENTS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls the calculation of excited-date (CIS or TDDFT) multipole moments.
    '''
        if value == "":
            if "CIS_MOMENTS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_MOMENTS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_MOMENTS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_MOMENTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_MOMENTS"] = value.lower()

    def cis_ampl_anal(self, value="show"):
        '''
Name: CIS_AMPL_ANAL
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Perform additional analysis of CIS and TDDFT excitation amplitudes, including generation of natural transition orbitals, excited-state multipole moments, and Mulliken analysis of the excited state densities and particle/hole density matrices.
    '''
        if value == "":
            if "CIS_AMPL_ANAL" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_AMPL_ANAL"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_AMPL_ANAL" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_AMPL_ANAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_AMPL_ANAL"] = value.lower()

    def spin_flip(self, value="show"):
        '''
Name: SPIN_FLIP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects whether to perform a spin-?ip calculation.  Spin multiplicity should be set to 3 for systems with an even number of electrons, and 4 for systems with an odd number of electrons.
    '''
        if value == "":
            if "SPIN_FLIP" in self.dict_of_keywords:
                del self.dict_of_keywords["SPIN_FLIP"]
                print("Keyword removed.")
        elif value == "show":
            if "SPIN_FLIP" in self.dict_of_keywords:
                return self.dict_of_keywords["SPIN_FLIP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SPIN_FLIP"] = value.lower()

    def spin_flip_xcis(self, value="show"):
        '''
Name: SPIN_FLIP_XCIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Do a SF-XCIS
    '''
        if value == "":
            if "SPIN_FLIP_XCIS" in self.dict_of_keywords:
                del self.dict_of_keywords["SPIN_FLIP_XCIS"]
                print("Keyword removed.")
        elif value == "show":
            if "SPIN_FLIP_XCIS" in self.dict_of_keywords:
                return self.dict_of_keywords["SPIN_FLIP_XCIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SPIN_FLIP_XCIS"] = value.lower()

    def cis_mulliken(self, value="show"):
        '''
Name: CIS_MULLIKEN
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls Mulliken and Loewdin population analyses for excited-state particle and hole density matrices. 

    '''
        if value == "":
            if "CIS_MULLIKEN" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_MULLIKEN"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_MULLIKEN" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_MULLIKEN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_MULLIKEN"] = value.lower()

    def jobtype(self, value="show"):
        '''
Name: JOBTYPE
Type: STRING
Default: SP

Options:
    'SP'............................ Energy
    'Force'......................... Forces
    'Optimization'.................. Geometry
    'PES_Scan'...................... PES Scan
    'RPath'......................... Reaction Path
    'TS'............................ Transition State
    'Frequency'..................... Frequencies
    'NMR'........................... Chemical Shifts
    'ISSC'.......................... Spin-Spin Couplings
    'AIMD'.......................... Ab Initio MD
    'SP'............................ Properties

Description: Specifies the type of calculation to run. 
    '''
        if value == "":
            if "JOBTYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["JOBTYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "JOBTYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["JOBTYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["JOBTYPE"] = value.lower()

    def hfpt_basis(self, value="show"):
        '''
Name: HFPT_BASIS
Type: STRING
Default: 6-31++G**

Options:
    '6-31++G**'..................... 6-31++G**
    '6-311+G*'...................... 6-311+G*
    '6-311G**'...................... 6-311G**
    '6-311++G**'.................... 6-311++G**
    '6-311++G(3df,3pd)'............. 6-311++G(3df,3pd)

    'cc-pVTZ'....................... cc-pVTZ
    'cc-pVQZ'....................... cc-pVQZ
    'cc-pcVTZ'...................... cc-pcVTZ
    'cc-pcVQZ'...................... cc-pcVQZ
    'aug-cc-pVTZ'................... aug-cc-pVTZ
    'aug-cc-pVQZ'................... aug-cc-pVQZ
    'aug-cc-pcVTZ'.................. aug-cc-pcVTZ
    'aug-cc-pcVQZ'.................. aug-cc-pcVQZ

    'G3Large'....................... G3Large
    'G3MP2Large'.................... G3MP2Large

    'gen'........................... User-defined

Description: Specifies the secondary basis in a HFPC/DFPC calculation
    '''
        if value == "":
            if "HFPT_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["HFPT_BASIS"]
                print("Keyword removed.")
        elif value == "show":
            if "HFPT_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["HFPT_BASIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["HFPT_BASIS"] = value.lower()

    def hirshfeld(self, value="show"):
        '''
Name: HIRSHFELD
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Compute Hirshfeld populations
    '''
        if value == "":
            if "HIRSHFELD" in self.dict_of_keywords:
                del self.dict_of_keywords["HIRSHFELD"]
                print("Keyword removed.")
        elif value == "show":
            if "HIRSHFELD" in self.dict_of_keywords:
                return self.dict_of_keywords["HIRSHFELD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["HIRSHFELD"] = value.lower()

    def qui_use_case(self, value="show"):
        '''
Name: QUI_USE_CASE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Selects the atenuated coulomb operator (CASE approximation).
    '''
        if value == "":
            if "QUI_USE_CASE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_USE_CASE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_USE_CASE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_USE_CASE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_USE_CASE"] = value.lower()

    def qui_cis_guess(self, value="show"):
        '''
Name: QUI_CIS_GUESS
Type: STRING
Default: New Guess

Options:
    'New Guess'..................... New Guess
    'Read Singlets'................. Read Singlets
    'Read Triplets'................. Read Triplets
    'Read Singlets & Triplets'...... Read Singlets & Triplets

Description: Determines what guess should be used for a CIS calculation.  Note that the read options require a guess from a previous calculation.
    '''
        if value == "":
            if "QUI_CIS_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_CIS_GUESS"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_CIS_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_CIS_GUESS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_CIS_GUESS"] = value.lower()

    def sts_mom(self, value="show"):
        '''
Name: STS_MOM
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Compute transition moments between excited states in the CIS and TD-DFT calculations.
    '''
        if value == "":
            if "STS_MOM" in self.dict_of_keywords:
                del self.dict_of_keywords["STS_MOM"]
                print("Keyword removed.")
        elif value == "show":
            if "STS_MOM" in self.dict_of_keywords:
                return self.dict_of_keywords["STS_MOM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["STS_MOM"] = value.lower()

    def max_cis_subspace(self, value="show"):
        '''
Name: MAX_CIS_SUBSPACE
Type: INTEGER
Default: 0
Options: Range from 0 to 9999

Description: Maximum number of subspace vectors allowed in the CIS iterations:

Recommendation: : The default is usually appropriate, unless a large number of states are requested for a large molecule. The total memory required to store the subspace vectors is bounded above by 2nOV , where O and V represent the number of occupied and virtual orbitals, respectively. n can be reduced to save memory, at the cost of a larger number of CIS iterations. Convergence may be impaired if n is not much larger than CIS_N_ROOTS.    '''
        if value == "":
            if "MAX_CIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["MAX_CIS_SUBSPACE"]
                print("Keyword removed.")
        elif value == "show":
            if "MAX_CIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["MAX_CIS_SUBSPACE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MAX_CIS_SUBSPACE"] = value.lower()

    def cis_dynamic_memory(self, value="show"):
        '''
Name: CIS_DYNAMIC_MEMORY
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether to use static or dynamic memory in CIS and TDDFT calculations.

Recommendation: : The default control requires static memory (MEM_STATIC) to hold a temporary array whose minimum size is OV ? CIS_N_ROOTS. For a large calculation, one has to specify a large value for MEM_STATIC, which is not recommended. Therefore, it is recommended to use dynamic memory for large calculations. 
    '''
        if value == "":
            if "CIS_DYNAMIC_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_DYNAMIC_MEMORY"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_DYNAMIC_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_DYNAMIC_MEMORY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_DYNAMIC_MEMORY"] = value.lower()

    def roks(self, value="show"):
        '''
Name: ROKS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Controls whether a restricted open-shell Kohn-Sham calculation will be performed.
    '''
        if value == "":
            if "ROKS" in self.dict_of_keywords:
                del self.dict_of_keywords["ROKS"]
                print("Keyword removed.")
        elif value == "show":
            if "ROKS" in self.dict_of_keywords:
                return self.dict_of_keywords["ROKS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ROKS"] = value.lower()

    def roks_level_shift(self, value="show"):
        '''
Name: ROKS_LEVEL_SHIFT
Type: INTEGER
Factor: 0.01
Default: 0 [=0.00]
Options: Range from 0 [=0.00] to 999 [=9.99]

Description: Introduce a level-shift (in Hartree) to aid convergence.
    '''
        if value == "":
            if "ROKS_LEVEL_SHIFT" in self.dict_of_keywords:
                del self.dict_of_keywords["ROKS_LEVEL_SHIFT"]
                print("Keyword removed.")
        elif value == "show":
            if "ROKS_LEVEL_SHIFT" in self.dict_of_keywords:
                return self.dict_of_keywords["ROKS_LEVEL_SHIFT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ROKS_LEVEL_SHIFT"] = value.lower()

    def cc_dmaxvectors(self, value="show"):
        '''
Name: CC_DMAXVECTORS
Type: INTEGER
Default: 60
Options: Range from 1 to 500

Description: Specifies maximum number of vectors in the subspace for the Davidson diagonalization. 
Recommendation: : Larger values increase disk storage but accelerate and stabilize convergence.    '''
        if value == "":
            if "CC_DMAXVECTORS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DMAXVECTORS"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_DMAXVECTORS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DMAXVECTORS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_DMAXVECTORS"] = value.lower()

    def eda_covp(self, value="show"):
        '''
Name: EDA_COVP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Perform COVP analysis when evaluating the RS or ARS charge-transfer correction. COVP analysis is currently implemented only for systems of two fragments.
Recommendation: : Set to TRUE to perform COVP analysis of the CT term in an EDA or SCF MI(RS) job.    '''
        if value == "":
            if "EDA_COVP" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_COVP"]
                print("Keyword removed.")
        elif value == "show":
            if "EDA_COVP" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_COVP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EDA_COVP"] = value.lower()

    def eom_davidson_convergence(self, value="show"):
        '''
Name: EOM_DAVIDSON_CONVERGENCE
Type: INTEGER
Default: 5
Options: Range from 0 to 12

Description: Convergence criterion for the RMS residuals of excited state vectors.
n Corresponding to 10^-n convergence criterion
Use default. Should normally be set to the same value as
EOM DAVIDSON THRESHOLD.
    '''
        if value == "":
            if "EOM_DAVIDSON_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_DAVIDSON_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"] = value.lower()

    def eom_davidson_maxvectors(self, value="show"):
        '''
Name: EOM_DAVIDSON_MAXVECTORS
Type: INTEGER
Default: 60
Options: Range from 0 to 200

Description: Species maximum number of vectors in the subspace for the Davidson diagonalization.
n Up to n vectors per root before the subspace is reset
RECOMMENDATION: 
Larger values increase disk storage but accelerate and stabilize convergence.
    '''
        if value == "":
            if "EOM_DAVIDSON_MAXVECTORS" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_MAXVECTORS"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_DAVIDSON_MAXVECTORS" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_MAXVECTORS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_DAVIDSON_MAXVECTORS"] = value.lower()

    def eom_ref_prop_te(self, value="show"):
        '''
Name: EOM_REF_PROP_TE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Request for calculation of non-relaxed two-particle EOM-CC properties. The two-particle properties currently include 
?S2?. The one-particle properties also will be calculated, since their cost is much less than the two-particle properties. The variable CC_EOM_PROP must be also set to TRUE. Alternatively, CC_CALC_SSQ can be used to request ?S2? calculation. 

Recommendation: :  The two-particle properties are computationally expensive since they require calculation 
and use of the two-particle density matrix (the cost is approximately the same as the cost of an analytic gradient calculation). Do not request the two-particle properties unless you really need them.     '''
        if value == "":
            if "EOM_REF_PROP_TE" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_REF_PROP_TE"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_REF_PROP_TE" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_REF_PROP_TE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_REF_PROP_TE"] = value.lower()

    def qui_eom_ea(self, value="show"):
        '''
Name: QUI_EOM_EA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_EOM_EA" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_EA"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_EA" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_EA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_EA"] = value.lower()

    def qui_eom_ip(self, value="show"):
        '''
Name: QUI_EOM_IP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_EOM_IP" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_IP"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_IP" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_IP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_IP"] = value.lower()

    def qui_eom_sf(self, value="show"):
        '''
Name: QUI_EOM_SF
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_EOM_SF" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_SF"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_SF" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_SF"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_SF"] = value.lower()

    def cc_theta_stepsize(self, value="show"):
        '''
Name: CC_THETA_STEPSIZE
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Scale factor for the orbital rotation step size. The optimal rotation steps should be approximately equal to the gradient vector.
Recommendation: : Try a smaller value in cases of poor convergence and very large orbital gradients. For example, a value of 01001 translates to 0.1    '''
        if value == "":
            if "CC_THETA_STEPSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_STEPSIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_THETA_STEPSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_STEPSIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_THETA_STEPSIZE"] = value.lower()

    def qui_eom_ee(self, value="show"):
        '''
Name: QUI_EOM_EE
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_EOM_EE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_EE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_EE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_EE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_EE"] = value.lower()

    def qui_eom_dip(self, value="show"):
        '''
Name: QUI_EOM_DIP
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "QUI_EOM_DIP" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_DIP"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_DIP" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_DIP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_DIP"] = value.lower()

    def aimd_fock_extrapolation_points(self, value="show"):
        '''
Name: AIMD_FOCK_EXTRAPOLATION_POINTS
Type: INTEGER
Default: 0
Options: Range from 0 to 20

Description: Specifies the number M of old Fock matrices that are retained for use in extrapolation. 
Recommendation: : Higher-order extrapolations with more saved Fock matrices are faster and conserve energy better than low-order extrapolations, up to a point. In many cases, the scheme (N = 6, M = 12), in conjunction with SCF_CONVERGENCE = 6, is found to provide about a 50% savings in computational cost while still conserving energy.    '''
        if value == "":
            if "AIMD_FOCK_EXTRAPOLATION_POINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_POINTS"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_FOCK_EXTRAPOLATION_POINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_POINTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "AIMD_FOCK_EXTRAPOLATION_POINTS"] = value.lower()

    def aimd_time_step(self, value="show"):
        '''
Name: AIMD_TIME_STEP
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Specifies the molecular dynamics time step, in atomic units (1 a.u. = 0.0242 fs).
Recommendation: : Smaller time steps lead to better energy conservation; too large a time step may cause the job to fail entirely. Make the time step as large as possible, consistent with tolerable energy conservation.    '''
        if value == "":
            if "AIMD_TIME_STEP" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_TIME_STEP"]
                print("Keyword removed.")
        elif value == "show":
            if "AIMD_TIME_STEP" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_TIME_STEP"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIMD_TIME_STEP"] = value.lower()

    def cc_eom_state_to_opt(self, value="show"):
        '''
Name: CC_EOM_STATE_TO_OPT
Type: undefined
Default: [0,0]

Options:
    '[0,0]'......................... [0,0]

Description: Specifies which state to optimize.
[i,j] optimize the jth state of the ith irrep.
    '''
        if value == "":
            if "CC_EOM_STATE_TO_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_STATE_TO_OPT"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_EOM_STATE_TO_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_STATE_TO_OPT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_EOM_STATE_TO_OPT"] = value.lower()

    def cc_high_spin(self, value="show"):
        '''
Name: CC_HIGH_SPIN
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of high-spin excited state roots to find.  Works only for singlet reference states and triplet excited states.  The default is to not look for any excited states.  Setting this to [i,j,k...] looks for i excited states in the first irrep, j states in the second irrep and so on.
    '''
        if value == "":
            if "CC_HIGH_SPIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_HIGH_SPIN"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_HIGH_SPIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_HIGH_SPIN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_HIGH_SPIN"] = value.lower()

    def cc_iterate_ov(self, value="show"):
        '''
Name: CC_ITERATE_OV
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: In active space calculations, use a mixed iteration procedure if the value is greater than 0. Then, if the RMS orbital gradient is larger than the value of CC_THETA_GRAD_THRESH, micro-iterations will be performed to converge the occupied-virtual mixing angles for the current active space. The maximum number of such iterations is given by this option. 
Recommendation: : Can be useful for non-convergent active space calculations.    '''
        if value == "":
            if "CC_ITERATE_OV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ITERATE_OV"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_ITERATE_OV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ITERATE_OV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_ITERATE_OV"] = value.lower()

    def cc_low_spin(self, value="show"):
        '''
Name: CC_LOW_SPIN
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of low-spin excited state roots to find.  In the cas of closed-shell reference states, excited singlet states will be found.  For any other reference state all states (e.g. both singlet and triplet) will be calculated.  The default is to not look for any excited states.  Setting this to [i,j,k...] looks for i excited states in the first irrep, j states in the second irrep and so on.
    '''
        if value == "":
            if "CC_LOW_SPIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_LOW_SPIN"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_LOW_SPIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_LOW_SPIN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_LOW_SPIN"] = value.lower()

    def cc_preconverge_fz(self, value="show"):
        '''
Name: CC_PRECONVERGE_FZ
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: In active space methods, whether to preconverge other wavefunction variables for fixed initial guess of active space.
Recommendation: :     '''
        if value == "":
            if "CC_PRECONVERGE_FZ" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_FZ"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_FZ" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_FZ"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_FZ"] = value.lower()

    def cc_preconverge_t2z(self, value="show"):
        '''
Name: CC_PRECONVERGE_T2Z
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Whether to pre-converge the cluster amplitudes before beginning orbital optimization in optimized orbital cluster methods.
Recommendation: : Experiment with this option in cases of convergence failure.    '''
        if value == "":
            if "CC_PRECONVERGE_T2Z" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_T2Z"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_T2Z" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_T2Z"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_T2Z"] = value.lower()

    def cc_preconverge_t2z_each(self, value="show"):
        '''
Name: CC_PRECONVERGE_T2Z_EACH
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Whether to pre-converge the cluster amplitudes before each change of the orbitals in optimized orbital coupled-cluster methods. The maximum number of iterations in this pre-convergence procedure is given by the value of this parameter.
Recommendation: : A very slow last resort option for jobs that do not converge.    '''
        if value == "":
            if "CC_PRECONVERGE_T2Z_EACH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONVERGE_T2Z_EACH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"] = value.lower()

    def cc_preconv_t2z(self, value="show"):
        '''
Name: CC_PRECONV_T2Z
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Whether to pre-converge the cluster amplitudes before beginning orbital optimization in optimized orbital cluster methods.
Recommendation: : Experiment with this option in cases of convergence failure.    '''
        if value == "":
            if "CC_PRECONV_T2Z" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_T2Z"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONV_T2Z" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_T2Z"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONV_T2Z"] = value.lower()

    def cc_preconv_t2z_each(self, value="show"):
        '''
Name: CC_PRECONV_T2Z_EACH
Type: INTEGER
Default: 0
Options: Range from 0 to 0

Description: Whether to pre-converge the cluster amplitudes before each change of the orbitals in optimized orbital coupled-cluster methods. The maximum number of iterations in this pre-convergence procedure is given by the value of this parameter.
Recommendation: : A very slow last resort option for jobs that do not converge.    '''
        if value == "":
            if "CC_PRECONV_T2Z_EACH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_T2Z_EACH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_PRECONV_T2Z_EACH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_T2Z_EACH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_PRECONV_T2Z_EACH"] = value.lower()

    def cc_qccd_theta_switch(self, value="show"):
        '''
Name: CC_QCCD_THETA_SWITCH
Type: INTEGER
Default: 0
Options: Range from 2 to 0

Description: QCCD calculations switch from OD to QCCD when the rotation gradient is below this threshold [10-n]
    '''
        if value == "":
            if "CC_QCCD_THETA_SWITCH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_QCCD_THETA_SWITCH"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_QCCD_THETA_SWITCH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_QCCD_THETA_SWITCH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_QCCD_THETA_SWITCH"] = value.lower()

    def cc_reference_symmetry(self, value="show"):
        '''
Name: CC_REFERENCE_SYMMETRY
Type: INTEGER
Default: -1
Options: Range from -1 to 0

Description: Together with CC_STATE_DERIV, selects which EOM state is to be considered for optimization or property calculations. When transition properties are requested, the transition properties will be calculated between this state and all other EOM states.
    '''
        if value == "":
            if "CC_REFERENCE_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REFERENCE_SYMMETRY"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REFERENCE_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REFERENCE_SYMMETRY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REFERENCE_SYMMETRY"] = value.lower()

    def cc_refsym(self, value="show"):
        '''
Name: CC_REFSYM
Type: INTEGER
Default: -1
Options: Range from -1 to 0

Description: Together with CC_STATE_DERIV, selects which EOM state is to be considered for optimization or property calculations. When transition properties are requested, the transition properties will be calculated between this state and all other EOM states.
    '''
        if value == "":
            if "CC_REFSYM" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REFSYM"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_REFSYM" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REFSYM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_REFSYM"] = value.lower()

    def cc_spin_flip_ms(self, value="show"):
        '''
Name: CC_SPIN_FLIP_MS
Type: INTEGER
Default: 0
Options: Range from 0 to 2

Description: This option is only used in EOM-SF using quintet references and including triple excitations. By default, SF flips the spin of one ? electron. One can ask to flip the spins of two ? electrons by specifying CC_SPIN_FLIP_MS = 1
Recommendation: : This option can be useful when starting from quintet references - though this is not typical for EOM-SF.    '''
        if value == "":
            if "CC_SPIN_FLIP_MS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SPIN_FLIP_MS"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_SPIN_FLIP_MS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SPIN_FLIP_MS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_SPIN_FLIP_MS"] = value.lower()

    def cc_state_deriv(self, value="show"):
        '''
Name: CC_STATE_DERIV
Type: INTEGER

Description: Selects which EOM or CIS(D) state is to be considered for optimization or property calculations.
    '''
        if value == "":
            if "CC_STATE_DERIV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_DERIV"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_STATE_DERIV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_DERIV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_STATE_DERIV"] = value.lower()

    def cc_state_to_opt(self, value="show"):
        '''
Name: CC_STATE_TO_OPT
Type: undefined
Default: [0,0]

Options:
    '[0,0]'......................... [0,0]

Description: Species which state to optimize.
[i,j] optimize the jth state of the ith irrep.
    '''
        if value == "":
            if "CC_STATE_TO_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_TO_OPT"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_STATE_TO_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_TO_OPT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_STATE_TO_OPT"] = value.lower()

    def cc_theta_conv(self, value="show"):
        '''
Name: CC_THETA_CONV
Type: INTEGER
Default: 0
Options: Range from 5 to 0

Description: Convergence criterion on the RMS difference between successive sets of orbital rotation angles [10-n].
Recommendation: : Use default    '''
        if value == "":
            if "CC_THETA_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_THETA_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_THETA_CONV"] = value.lower()

    def cc_theta_grad_conv(self, value="show"):
        '''
Name: CC_THETA_GRAD_CONV
Type: INTEGER
Default: 0
Options: Range from 7 to 0

Description: Convergence desired on the RMS gradient of the energy with respect to orbital rotation angles [10-n]. 
Recommendation: : Use default    '''
        if value == "":
            if "CC_THETA_GRAD_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_GRAD_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "CC_THETA_GRAD_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_GRAD_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CC_THETA_GRAD_CONV"] = value.lower()

    def cis_guess_disk_type(self, value="show"):
        '''
Name: CIS_GUESS_DISK_TYPE
Type: INTEGER
Default: -1
Options: Range from -1 to 2

Description: Determines the type of guesses to be read from disk
Recommendation: : Must be specified if a CIS guess in to be read from disk.    '''
        if value == "":
            if "CIS_GUESS_DISK_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_GUESS_DISK_TYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_GUESS_DISK_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_GUESS_DISK_TYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_GUESS_DISK_TYPE"] = value.lower()

    def cis_n_roots(self, value="show"):
        '''
Name: CIS_N_ROOTS
Type: INTEGER
Default: 0
Options: Range from 0 to 200

Description: Sets the number of CI-Singles (CIS) excited state roots to find
    '''
        if value == "":
            if "CIS_N_ROOTS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_N_ROOTS"]
                print("Keyword removed.")
        elif value == "show":
            if "CIS_N_ROOTS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_N_ROOTS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CIS_N_ROOTS"] = value.lower()

    def diis_print(self, value="show"):
        '''
Name: DIIS_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 4

Description: Controls the output from DIIS SCF optimization:
1: Chosen method and DIIS coefficients
2: Level 1 + print changes in multipole moments
3: Level 2 + multipole moments
4: Level 3 + extrapolated Fock matrices
    '''
        if value == "":
            if "DIIS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "DIIS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIIS_PRINT"] = value.lower()

    def dip_singlets(self, value="show"):
        '''
Name: DIP_SINGLETS
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of singlet DIP roots to find. Works only for closed-shell references.
[i, j, k...] Find i DIP singlet states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "DIP_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["DIP_SINGLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "DIP_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["DIP_SINGLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIP_SINGLETS"] = value.lower()

    def dip_states(self, value="show"):
        '''
Name: DIP_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of DIP roots to find. For closed-shell reference, defaults into EOM_DIP_SINGLETS. For open-shell references, speccies all low-lying states.
[i, j, k...] Find i DIP states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "DIP_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["DIP_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "DIP_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["DIP_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIP_STATES"] = value.lower()

    def dip_triplets(self, value="show"):
        '''
Name: DIP_TRIPLETS
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of triplet DIP roots to find. Works only for closed-shell references.
[i, j, k...] Find i DIP triplet states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "DIP_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["DIP_TRIPLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "DIP_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["DIP_TRIPLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DIP_TRIPLETS"] = value.lower()

    def dsf_states(self, value="show"):
        '''
Name: DSF_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of double spin-?ip target states roots to ?nd.
[i, j, k . . .] Find i SF states in the ?rst irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "DSF_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["DSF_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "DSF_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["DSF_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["DSF_STATES"] = value.lower()

    def ee_singlets(self, value="show"):
        '''
Name: EE_SINGLETS
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "EE_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EE_SINGLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "EE_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EE_SINGLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EE_SINGLETS"] = value.lower()

    def ee_states(self, value="show"):
        '''
Name: EE_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "EE_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EE_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "EE_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EE_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EE_STATES"] = value.lower()

    def ee_triplets(self, value="show"):
        '''
Name: EE_TRIPLETS
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "EE_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EE_TRIPLETS"]
                print("Keyword removed.")
        elif value == "show":
            if "EE_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EE_TRIPLETS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EE_TRIPLETS"] = value.lower()

    def eom_ea_alpha(self, value="show"):
        '''
Name: EOM_EA_ALPHA
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of attached target states derived by attaching alpha electron (Ms = 1/2).
[i, j, k...] Find i EA states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_EA_ALPHA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_ALPHA"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_EA_ALPHA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_ALPHA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_EA_ALPHA"] = value.lower()

    def eom_ea_beta(self, value="show"):
        '''
Name: EOM_EA_BETA
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of attached target states derived by attaching beta electron (Ms = -1/2).
[i, j, k...] Find i EA states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_EA_BETA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_BETA"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_EA_BETA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_BETA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_EA_BETA"] = value.lower()

    def eom_ea_states(self, value="show"):
        '''
Name: EOM_EA_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of attached target states roots to find. By default, alpha electron will be attached (see EOM_EA_ALPHA).
[i; j; k...] Find i EA states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_EA_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_EA_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_EA_STATES"] = value.lower()

    def eom_ip_alpha(self, value="show"):
        '''
Name: EOM_IP_ALPHA
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of ionized target states derived by removing alpha electron (Ms = -1/2).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_IP_ALPHA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_ALPHA"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_IP_ALPHA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_ALPHA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_IP_ALPHA"] = value.lower()

    def eom_ip_beta(self, value="show"):
        '''
Name: EOM_IP_BETA
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of ionized target states derived by removing beta electron (Ms = 1/2, default for EOM-IP).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_IP_BETA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_BETA"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_IP_BETA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_BETA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_IP_BETA"] = value.lower()

    def eom_ip_states(self, value="show"):
        '''
Name: EOM_IP_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of ionized target states roots to find. By default, B electron will be removed (see EOM_IP_BETA).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "EOM_IP_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "EOM_IP_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["EOM_IP_STATES"] = value.lower()

    def incdft_griddiff_varthresh(self, value="show"):
        '''
Name: INCDFT_GRIDDIFF_VARTHRESH
Type: INTEGER
Default: 0
Options: Range from 0 to 12

Description: Sets the lower bound for the variable threshold for screening the functional values in the IncDFT procedure. The threshold will begin at this value and then vary depending on the error in the current SCF iteration until the value specified by INCDFT_GRIDDIFF_THRESH is reached. This means that this value must be set lower than INCDFT_GRIDDIFF_THRESH.
Recommendation: : If the default value causes convergence problems, set this value higher to tighten accuracy. If this fails, set to 0 and use a static threshold.    '''
        if value == "":
            if "INCDFT_GRIDDIFF_VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "INCDFT_GRIDDIFF_VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"] = value.lower()

    def mom_start(self, value="show"):
        '''
Name: MOM_START
Type: INTEGER
Default: 0
Options: Range from 0 to 999

Description: Determines when MOM is switched on to stabilize DIIS iterations.
Recommendation: : Set to 1 if preservation of initial orbitals is desired. If MOM is to be used to aid convergence, an SCF without MOM should be run to determine when the SCF starts oscillating. MOM should be set to start just before the oscillations.    '''
        if value == "":
            if "MOM_START" in self.dict_of_keywords:
                del self.dict_of_keywords["MOM_START"]
                print("Keyword removed.")
        elif value == "show":
            if "MOM_START" in self.dict_of_keywords:
                return self.dict_of_keywords["MOM_START"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MOM_START"] = value.lower()

    def nvo_method(self, value="show"):
        '''
Name: NVO_METHOD
Type: INTEGER
Default: 0
Options: Range from 9 to 0

Description: Sets method to be used to converge solution of the single-excitation amplitude equations.
Recommendation: : Experimental option. Use default.    '''
        if value == "":
            if "NVO_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_METHOD"] = value.lower()

    def nvo_truncate_dist(self, value="show"):
        '''
Name: NVO_TRUNCATE_DIST
Type: INTEGER
Default: -1
Options: Range from -2 to 999

Description: Specifies which atomic blocks of the Fock matrix are used to construct the preconditioner.
Recommendation: : This option does not affect the final result. However, it affects the rate of the PCG algorithm convergence. For small systems use default.    '''
        if value == "":
            if "NVO_TRUNCATE_DIST" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_TRUNCATE_DIST"]
                print("Keyword removed.")
        elif value == "show":
            if "NVO_TRUNCATE_DIST" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_TRUNCATE_DIST"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NVO_TRUNCATE_DIST"] = value.lower()

    def qui_eom_states1(self, value="show"):
        '''
Name: QUI_EOM_STATES1
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "QUI_EOM_STATES1" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_STATES1"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_STATES1" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_STATES1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_STATES1"] = value.lower()

    def qui_eom_states2(self, value="show"):
        '''
Name: QUI_EOM_STATES2
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "QUI_EOM_STATES2" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_STATES2"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_EOM_STATES2" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_STATES2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_EOM_STATES2"] = value.lower()

    def scf_final_print(self, value="show"):
        '''
Name: SCF_FINAL_PRINT
Type: INTEGER
Default: 0
Options: Range from 0 to 3

Description: Controls level of output from SCF procedure to Q-Chem output file at the end of the SCF.
Recommendation: : The break-down of energies is often useful (level 1).    '''
        if value == "":
            if "SCF_FINAL_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_FINAL_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_FINAL_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_FINAL_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_FINAL_PRINT"] = value.lower()

    def scf_guess(self, value="show"):
        '''
Name: SCF_GUESS
Type: STRING
Default: SAD

Options:
    'SAD'........................... SAD
    'CORE'.......................... CORE
    'GWH'........................... GWH
    'READ'.......................... READ

Description: Specifies the initial guess procedure to use for the SCF.
Recommendation: : SAD guess for standard basis sets. For general basis sets, it is best to use the BASIS2 $rem. Alternatively, try the GWH or core Hamiltonian guess. For ROHF it can be useful to READ guesses from an SCF calculation on the corresponding cation or anion. Note that because the density is made spherical, this may favor an undesired state for atomic systems, especially transition metals.    '''
        if value == "":
            if "SCF_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS"]
                print("Keyword removed.")
        elif value == "show":
            if "SCF_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SCF_GUESS"] = value.lower()

    def sf_states(self, value="show"):
        '''
Name: SF_STATES
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Sets the number of spin-?ip target states roots to ?nd.
[i, j, k . . .] Find i SF states in the ?rst irrep, j states in the second irrep etc.
    '''
        if value == "":
            if "SF_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["SF_STATES"]
                print("Keyword removed.")
        elif value == "show":
            if "SF_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["SF_STATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SF_STATES"] = value.lower()

    def varthresh(self, value="show"):
        '''
Name: VARTHRESH
Type: INTEGER
Default: 0
Options: Range from 0 to 12

Description: Controls the temporary integral cut-off threshold. The variable threshold is set to 10-VARTHRESH? DIIS_error
Recommendation: : 3 has been found to be a practical level, and can slightly speed up SCF evaluation.    '''
        if value == "":
            if "VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["VARTHRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["VARTHRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["VARTHRESH"] = value.lower()

    def xopt_state_1(self, value="show"):
        '''
Name: XOPT_STATE_1
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Specify two electronic states the intersection of which will be searched.
[spin, irrep, state]
spin = 0 Addresses states with low spin, see also EOM EE SINGLETS.
spin = 1 Addresses states with high spin, see also EOM EE TRIPLETS.
irrep Species the irreducible representation to which the state belongs, for C2v point group symmetry 
irrep = 1 for A1, irrep = 2 for A2,
irrep = 3 for B1, irrep = 4 for B2.
state Species the state number within the irreducible
representation, state = 1 means the lowest excited
state, state = 2 is the second excited state, etc.
0, 0, -1 Ground state.
    '''
        if value == "":
            if "XOPT_STATE_1" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_STATE_1"]
                print("Keyword removed.")
        elif value == "show":
            if "XOPT_STATE_1" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_STATE_1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XOPT_STATE_1"] = value.lower()

    def xopt_state_2(self, value="show"):
        '''
Name: XOPT_STATE_2
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: Specify two electronic states the intersection of which will be searched.
[spin, irrep, state]
spin = 0 Addresses states with low spin, see also EOM EE SINGLETS.
spin = 1 Addresses states with high spin, see also EOM EE TRIPLETS.
irrep Species the irreducible representation to which the state belongs, for C2v point group symmetry 
irrep = 1 for A1, irrep = 2 for A2,
irrep = 3 for B1, irrep = 4 for B2.
state Species the state number within the irreducible
representation, state = 1 means the lowest excited
state, state = 2 is the second excited state, etc.
0, 0, -1 Ground state.
    '''
        if value == "":
            if "XOPT_STATE_2" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_STATE_2"]
                print("Keyword removed.")
        elif value == "show":
            if "XOPT_STATE_2" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_STATE_2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XOPT_STATE_2"] = value.lower()

    def qui_section_swap_occupied_virtual(self, value="show"):
        '''
Name: QUI_SECTION_SWAP_OCCUPIED_VIRTUAL
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Change the occupancies of the guess orbitals (not compatible with the SAD guess)
    '''
        if value == "":
            if "QUI_SECTION_SWAP_OCCUPIED_VIRTUAL" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SECTION_SWAP_OCCUPIED_VIRTUAL"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_SECTION_SWAP_OCCUPIED_VIRTUAL" in self.dict_of_keywords:
                return self.dict_of_keywords[
                    "QUI_SECTION_SWAP_OCCUPIED_VIRTUAL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords[
                "QUI_SECTION_SWAP_OCCUPIED_VIRTUAL"] = value.lower()

    def adc_davidson_maxiter(self, value="show"):
        '''
Name: ADC_DAVIDSON_MAXITER
Type: INTEGER
Default: 60
Options: Range from 1 to 500

Description: Maximum number of iterations to determine the eigenstates in an ADC calculation using the Davidson algorithm.
    '''
        if value == "":
            if "ADC_DAVIDSON_MAXITER" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DAVIDSON_MAXITER"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DAVIDSON_MAXITER" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DAVIDSON_MAXITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DAVIDSON_MAXITER"] = value.lower()

    def adc_diis_maxiter(self, value="show"):
        '''
Name: ADC_DIIS_MAXITER
Type: STRING
Default: 

Options:
    ''.............................. 

Description: Maximum number of iterations to determine the eigenstates in an ADC calculation using the DIIS algorithm.
    '''
        if value == "":
            if "ADC_DIIS_MAXITER" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DIIS_MAXITER"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DIIS_MAXITER" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DIIS_MAXITER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DIIS_MAXITER"] = value.lower()

    def adc_davidson_maxsubspace(self, value="show"):
        '''
Name: ADC_DAVIDSON_MAXSUBSPACE
Type: INTEGER
Default: 40
Options: Range from 1 to 500

Description: Maximum number of ADC amplitudes in the subspace for the Davidson diagonalization. 
Recommendation: : Larger values increase disk storage but accelerate and stabilize convergence.
    '''
        if value == "":
            if "ADC_DAVIDSON_MAXSUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DAVIDSON_MAXSUBSPACE"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DAVIDSON_MAXSUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DAVIDSON_MAXSUBSPACE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DAVIDSON_MAXSUBSPACE"] = value.lower()

    def adc_c_c(self, value="show"):
        '''
Name: ADC_C_C
Type: INTEGER
Default: 1000
Options: Range from 0 to 9999

Description: Scaling factor cC for spin-opposite scaled ADC calculations. The parameter value is devided by 1000 to obtain the proper scaling factor. 

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_C_C" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_C_C"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_C_C" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_C_C"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_C_C"] = value.lower()

    def adc_c_t(self, value="show"):
        '''
Name: ADC_C_T
Type: INTEGER
Default: 1300
Options: Range from 0 to 9999

Description: Scaling factor cT  for MP(2) T amplitudes in spin-opposite scaled ADC calculations. The parameter value is devided by 1000 to obtain the proper scaling factor. 

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_C_T" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_C_T"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_C_T" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_C_T"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_C_T"] = value.lower()

    def adc_c_x(self, value="show"):
        '''
Name: ADC_C_X
Type: INTEGER
Default: 900
Options: Range from 0 to 9999

Description: Scaling factor cX  in spin-opposite scaled ADC calculations. The parameter value is devided by 1000 to obtain the proper scaling factor. 

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_C_X" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_C_X"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_C_X" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_C_X"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_C_X"] = value.lower()

    def adc_davidson_conv(self, value="show"):
        '''
Name: ADC_DAVIDSON_CONV
Type: INTEGER
Default: 6
Options: Range from 0 to 14

Description: Convergence criterion on the RMS difference between successive sets of ADC amplitudes during the Davidson diagonalization [10-n].
Recommendation: : Use default.    '''
        if value == "":
            if "ADC_DAVIDSON_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DAVIDSON_CONV"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DAVIDSON_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DAVIDSON_CONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DAVIDSON_CONV"] = value.lower()

    def adc_davidson_thresh(self, value="show"):
        '''
Name: ADC_DAVIDSON_THRESH
Type: INTEGER
Default: 10
Options: Range from 0 to 14

Description: Numerical threshold to suppress noise in the ADC amplitudes during the Davidson diagonalization [10-n].
Recommendation: : Use default.    '''
        if value == "":
            if "ADC_DAVIDSON_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DAVIDSON_THRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DAVIDSON_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DAVIDSON_THRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DAVIDSON_THRESH"] = value.lower()

    def adc_diis_econv(self, value="show"):
        '''
Name: ADC_DIIS_ECONV
Type: INTEGER
Default: 6
Options: Range from 0 to 14

Description: Convergence criterion on the RMS difference between successive sets of ADC amplitudes during the Davidson diagonalization [10-n].
Recommendation: : Use default.    '''
        if value == "":
            if "ADC_DIIS_ECONV" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DIIS_ECONV"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DIIS_ECONV" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DIIS_ECONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DIIS_ECONV"] = value.lower()

    def adc_diis_rconv(self, value="show"):
        '''
Name: ADC_DIIS_RCONV
Type: INTEGER
Default: 6
Options: Range from 0 to 14

Description: Convergence criterion on the RMS difference between ADC amplitudes in successive steps of the DIIS solver [10-n].
Recommendation: : Use default.    '''
        if value == "":
            if "ADC_DIIS_RCONV" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DIIS_RCONV"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DIIS_RCONV" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DIIS_RCONV"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DIIS_RCONV"] = value.lower()

    def adc_diis_size(self, value="show"):
        '''
Name: ADC_DIIS_SIZE
Type: INTEGER
Default: 5
Options: Range from 1 to 50

Description: Specifies the maximum size of the DIIS space in an ADC(2) calculation.
Recommendation: : Larger values involve larger amounts of disk storage.    '''
        if value == "":
            if "ADC_DIIS_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DIIS_SIZE"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DIIS_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DIIS_SIZE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DIIS_SIZE"] = value.lower()

    def adc_diis_start(self, value="show"):
        '''
Name: ADC_DIIS_START
Type: INTEGER
Default: 1
Options: Range from 0 to 200

Description: Iteration number when the DIIS steps are turned on in an ADC(2)  calculation. Set to a large number to disable DIIS and use the Jacobi algorithm. 

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_DIIS_START" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DIIS_START"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DIIS_START" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DIIS_START"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DIIS_START"] = value.lower()

    def adc_do_diis(self, value="show"):
        '''
Name: ADC_DO_DIIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Compute the ADC(2) eigenvalues and eigenvectors by solving non-linear equation systems using the DIIS algorithm. 

Recommendation: : Use only with extreme care!    '''
        if value == "":
            if "ADC_DO_DIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_DO_DIIS"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_DO_DIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_DO_DIIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_DO_DIIS"] = value.lower()

    def adc_ecorr(self, value="show"):
        '''
Name: ADC_ECORR
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the computation of higher-order energy corrections for ADC excitation energies. [EXPERIMENTAL]
    '''
        if value == "":
            if "ADC_ECORR" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_ECORR"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_ECORR" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_ECORR"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_ECORR"] = value.lower()

    def adc_extended(self, value="show"):
        '''
Name: ADC_EXTENDED
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Use the extended second order ADC variant ADC(2)-x. 

Recommendation: : This keyword is deprecated. Use METHOD instead.    '''
        if value == "":
            if "ADC_EXTENDED" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_EXTENDED"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_EXTENDED" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_EXTENDED"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_EXTENDED"] = value.lower()

    def adc_cvs(self, value="show"):
        '''
Name: ADC_CVS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the core-valence separation approximation in ADC calculations to compute core-excited states. In addition the parameter CC_REST_OCC has to be set to define the core orbitals.

Recommendation: : This keyword is deprecated. Use METHOD instead.    '''
        if value == "":
            if "ADC_CVS" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_CVS"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_CVS" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_CVS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_CVS"] = value.lower()

    def adc_order(self, value="show"):
        '''
Name: ADC_ORDER
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate an ADC calculation and set the order of ADC to use.

Recommendation: : This keyword is deprecated. Use METHOD instead.    '''
        if value == "":
            if "ADC_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_ORDER"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_ORDER"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_ORDER"] = value.lower()

    def adc_nguess_doubles(self, value="show"):
        '''
Name: ADC_NGUESS_DOUBLES
Type: INTEGER
Default: 0
Options: Range from 0 to 500

Description: Set the number of guesses from the double excitation manifold.

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_NGUESS_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_NGUESS_DOUBLES"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_NGUESS_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_NGUESS_DOUBLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_NGUESS_DOUBLES"] = value.lower()

    def adc_nguess_singles(self, value="show"):
        '''
Name: ADC_NGUESS_SINGLES
Type: INTEGER
Default: 0
Options: Range from 0 to 500

Description: Set the number of guesses from the single excitation manifold.

Recommendation: : Use default ( = number of states to requested).    '''
        if value == "":
            if "ADC_NGUESS_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_NGUESS_SINGLES"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_NGUESS_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_NGUESS_SINGLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_NGUESS_SINGLES"] = value.lower()

    def adc_prop_es2es(self, value="show"):
        '''
Name: ADC_PROP_ES2ES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the calculation of state-to-state transition properties in an ADC calculation.
    '''
        if value == "":
            if "ADC_PROP_ES2ES" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_PROP_ES2ES"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_PROP_ES2ES" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_PROP_ES2ES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_PROP_ES2ES"] = value.lower()

    def adc_prop_tpa(self, value="show"):
        '''
Name: ADC_PROP_TPA
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the calculation of two-photon absorption cross-sections for the excited states in an ADC calculation.
    '''
        if value == "":
            if "ADC_PROP_TPA" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_PROP_TPA"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_PROP_TPA" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_PROP_TPA"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_PROP_TPA"] = value.lower()

    def adc_prop_es(self, value="show"):
        '''
Name: ADC_PROP_ES
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the calculation of excited state properties in an ADC calculation (one-particle transition properties from the ground state are activated by default).
    '''
        if value == "":
            if "ADC_PROP_ES" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_PROP_ES"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_PROP_ES" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_PROP_ES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_PROP_ES"] = value.lower()

    def adc_sos(self, value="show"):
        '''
Name: ADC_SOS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Activate the spin-opposite scaled variant of the second order ADC methods.

Recommendation: : This keyword is deprecated. Use METHOD instead.    '''
        if value == "":
            if "ADC_SOS" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_SOS"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_SOS" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_SOS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_SOS"] = value.lower()

    def adc_print(self, value="show"):
        '''
Name: ADC_PRINT
Type: INTEGER
Default: 2
Options: Range from 0 to 6

Description: Set the print level in the ADC part of the calculation.

Recommendation: : Use default.    '''
        if value == "":
            if "ADC_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["ADC_PRINT"]
                print("Keyword removed.")
        elif value == "show":
            if "ADC_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["ADC_PRINT"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["ADC_PRINT"] = value.lower()

    def state_analysis(self, value="show"):
        '''
Name: STATE_ANALYSIS
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: Performs certain excited state analyses for CIS/TD-DFT, ADC, and CC excited states. [EXPERIMENTAL]
    '''
        if value == "":
            if "STATE_ANALYSIS" in self.dict_of_keywords:
                del self.dict_of_keywords["STATE_ANALYSIS"]
                print("Keyword removed.")
        elif value == "show":
            if "STATE_ANALYSIS" in self.dict_of_keywords:
                return self.dict_of_keywords["STATE_ANALYSIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["STATE_ANALYSIS"] = value.lower()

    def qui_adc_states1(self, value="show"):
        '''
Name: QUI_ADC_STATES1
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "QUI_ADC_STATES1" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_ADC_STATES1"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_ADC_STATES1" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_ADC_STATES1"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_ADC_STATES1"] = value.lower()

    def qui_adc_states2(self, value="show"):
        '''
Name: QUI_ADC_STATES2
Type: undefined
Default: [0]

Options:
    '[0]'........................... [0]

Description: 
    '''
        if value == "":
            if "QUI_ADC_STATES2" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_ADC_STATES2"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_ADC_STATES2" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_ADC_STATES2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_ADC_STATES2"] = value.lower()

    def qui_adc_core(self, value="show"):
        '''
Name: QUI_ADC_CORE
Type: INTEGER
Default: 1
Options: Range from 1 to 500

Description: Set the number of core orbitals in an CVS-ADC calculation.
    '''
        if value == "":
            if "QUI_ADC_CORE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_ADC_CORE"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_ADC_CORE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_ADC_CORE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_ADC_CORE"] = value.lower()

    def geom_opt_hessian(self, value="show"):
        '''
Name: GEOM_OPT_HESSIAN
Type: STRING
Default: READ

Options:
    '0'............................. None
    'READ'.......................... READ
    'Diagonal'...................... Diagonal

Description: Determines the initial Hessian status.
Recommendation: : An accurate initial Hessian will improve the performance of the optimizer, but is expensive to compute.    '''
        if value == "":
            if "GEOM_OPT_HESSIAN" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_HESSIAN"]
                print("Keyword removed.")
        elif value == "show":
            if "GEOM_OPT_HESSIAN" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_HESSIAN"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["GEOM_OPT_HESSIAN"] = value.lower()

    def mp2v(self, value="show"):
        '''
Name: MP2V
Type: LOGICAL
Default: FALSE

Options:
    'FALSE'......................... FALSE
    'TRUE'.......................... TRUE

Description: 
    '''
        if value == "":
            if "MP2V" in self.dict_of_keywords:
                del self.dict_of_keywords["MP2V"]
                print("Keyword removed.")
        elif value == "show":
            if "MP2V" in self.dict_of_keywords:
                return self.dict_of_keywords["MP2V"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MP2V"] = value.lower()

    def correlation(self, value="show"):
        '''
Name: CORRELATION
Type: STRING
Default: None

Options:
    'None'.......................... None

    'B94'........................... B94
    'B94hyb'........................ B94hyb
    'B95'........................... B95
    'LYP'........................... LYP
    'LYP(EDF1)'..................... LYP(EDF1)
    'P86'........................... P86
    'PBE'........................... PBE
    'PK09'.......................... PK09
    'PK09'.......................... PK09
    'PW91'.......................... PW91
    'PW92'.......................... PW92
    'PZ81'.......................... PZ81
    'TPSS'.......................... TPSS
    'VWN'........................... VWN
    'Wigner'........................ Wigner
    '(PBE)OP'....................... (PBE)OP
    '(B88)OP'....................... (B88)OP

    'MP2'........................... MP2
    'CCMP2'......................... CCMP2
    'MP4SDQ'........................ MP4SDQ
    'ZAPT2'......................... ZAPT2
    'Local_MP2'..................... Local_MP2
    'RI-MP2'........................ RI-MP2
    'SOSMP2'........................ SOSMP2
    'ATTMP2'........................ ATT-MP2
    'MOSMP2'........................ MOSMP2
    'RILMP2'........................ RILMP2

    'MP3'........................... MP3
    'MP4'........................... MP4
    'CCD'........................... CCD
    'CCD(2)'........................ CCD(2)
    'CCSD'.......................... CCSD
    'CCSD(T)'....................... CCSD(T)
    'CCSD(2)'....................... CCSD(2)
    'CCSD(fT)'...................... CCSD(fT)
    'CCSD(dT)'...................... CCSD(dT)
    'QCCD'.......................... QCCD
    'QCCD(T)'....................... QCCD(T)
    'QCCD(2)'....................... QCCD(2)
    'QCISD'......................... QCISD
    'QCISD(T)'...................... QCISD(T)
    'OD'............................ OD
    'OD(T)'......................... OD(T)
    'OD(2)'......................... OD(2)
    'VOD'........................... VOD
    'VOD(2)'........................ VOD(2)
    'VQCCD'......................... VQCCD
    'VQCCD(T)'...................... VQCCD(T)
    'VQCCD(2)'...................... VQCCD(2)

    'PP'............................ PP
    'CCVB'.......................... CCVB
    'GVB_IP'........................ GVB_IP
    'GVB_SIP'....................... GVB_SIP
    'GVB_DIP'....................... GVB_DIP
    'OP'............................ OP
    'NP'............................ NP
    '2P'............................ 2P

Description: Specifies the correlation level of theory, either DFT or wavefunction-based.
Recommendation: : Consult the literature and reviews for guidence    '''
        if value == "":
            if "CORRELATION" in self.dict_of_keywords:
                del self.dict_of_keywords["CORRELATION"]
                print("Keyword removed.")
        elif value == "show":
            if "CORRELATION" in self.dict_of_keywords:
                return self.dict_of_keywords["CORRELATION"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["CORRELATION"] = value.lower()

    def qui_primary_basis(self, value="show"):
        '''
Name: QUI_PRIMARY_BASIS
Type: STRING
Default: 6-31G

Options:
    'STO-3G'........................ STO-3G
    'STO-6G'........................ STO-6G

    '3-21G'......................... 3-21G
    '4-31G'......................... 4-31G
    '6-31G'......................... 6-31G
    '6-31G*'........................ 6-31G*
    '6-31+G*'....................... 6-31+G*
    '6-31G**'....................... 6-31G**
    '6-31++G**'..................... 6-31++G**
    '6-311G'........................ 6-311G
    '6-311G*'....................... 6-311G*
    '6-311+G*'...................... 6-311+G*
    '6-311G**'...................... 6-311G**
    '6-311++G**'.................... 6-311++G**
    '6-311++G(3df,3pd)'............. 6-311++G(3df,3pd)

    'pc-0'.......................... pc-0
    'pc-1'.......................... pc-1
    'pc-2'.......................... pc-2
    'pc-3'.......................... pc-3
    'pc-4'.......................... pc-4
    'pcJ-0'......................... pcJ-0
    'pcJ-1'......................... pcJ-1
    'pcJ-2'......................... pcJ-2
    'pcJ-3'......................... pcJ-3
    'pcJ-4'......................... pcJ-4
    'pcS-0'......................... pcS-0
    'pcS-1'......................... pcS-1
    'pcS-2'......................... pcS-2
    'pcS-3'......................... pcS-3
    'pcS-4'......................... pcS-4

    'cc-pVDZ'....................... cc-pVDZ
    'cc-pVTZ'....................... cc-pVTZ
    'cc-pVQZ'....................... cc-pVQZ
    'cc-pcVDZ'...................... cc-pcVDZ
    'cc-pcVTZ'...................... cc-pcVTZ
    'cc-pcVQZ'...................... cc-pcVQZ
    'aug-cc-pVDZ'................... aug-cc-pVDZ
    'aug-cc-pVTZ'................... aug-cc-pVTZ
    'aug-cc-pVQZ'................... aug-cc-pVQZ
    'aug-cc-pcVDZ'.................. aug-cc-pcVDZ
    'aug-cc-pcVTZ'.................. aug-cc-pcVTZ
    'aug-cc-pcVQZ'.................. aug-cc-pcVQZ

    'G3Large'....................... G3Large
    'G3MP2Large'.................... G3MP2Large

    'CRENBL'........................ CRENBL
    'CRENBS'........................ CRENBS
    'HWMB'.......................... HWMB
    'HWVDZ'......................... HWVDZ
    'LACVP'......................... LACVP
    'LANL2DZ'....................... LANL2DZ
    'SBKJC'......................... SBKJC
    'SRLC'.......................... SRLC
    'SRSC'.......................... SRSC

    'gen'........................... User-defined
    'Mixed'......................... Mixed

Description: Specifies the basis sets to be used for the inital SCF to determine the occupied orbitals.
Recommendation: : The primary basis should be smaller than the target basis.    '''
        if value == "":
            if "QUI_PRIMARY_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_PRIMARY_BASIS"]
                print("Keyword removed.")
        elif value == "show":
            if "QUI_PRIMARY_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_PRIMARY_BASIS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["QUI_PRIMARY_BASIS"] = value.lower()

    def basis2_save(self, value="show"):
        '''
Name: BASIS2_SAVE
Type: STRING
Default: None

Options:
    'None'.......................... None

    'STO-3G'........................ STO-3G
    'STO-6G'........................ STO-6G
    '3-21G'......................... 3-21G
    '4-21G'......................... 4-21G
    '6-31G'......................... 6-31G
    '6-31G(d)'...................... 6-31G(d)
    'cc-pVDZ'....................... cc-pVDZ

    'r64G'.......................... r64G
    '6-31G'......................... 6-31G
    '6-31G*'........................ 6-31G*
    '6-311G*'....................... 6-311G*
    '6-311+G*'...................... 6-311+G*
    'rcc-pVTZ'...................... rcc-pVTZ
    'rcc-pVQZ'...................... rcc-pVQZ
    'racc-pVDZ'..................... racc-pVDZ
    'racc-pVTZ'..................... racc-pVTZ
    'racc-pVQZ'..................... racc-pVQZ

Description: Selects either a small basis set to use in basis set projection for the initial guess, or a subset basis for dual basis set calculations.
    '''
        if value == "":
            if "BASIS2_SAVE" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS2_SAVE"]
                print("Keyword removed.")
        elif value == "show":
            if "BASIS2_SAVE" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS2_SAVE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["BASIS2_SAVE"] = value.lower()

    def basis2(self, value="show"):
        '''
Name: BASIS2
Type: STRING
Default: None

Options:
    'None'.......................... None

    'STO-3G'........................ STO-3G
    'STO-6G'........................ STO-6G

    '3-21G'......................... 3-21G
    '4-31G'......................... 4-31G
    '6-31G'......................... 6-31G
    '6-31G*'........................ 6-31G*
    '6-31+G*'....................... 6-31+G*
    '6-31G**'....................... 6-31G**
    '6-31++G**'..................... 6-31++G**
    '6-311G'........................ 6-311G
    '6-311G*'....................... 6-311G*
    '6-311+G*'...................... 6-311+G*
    '6-311G**'...................... 6-311G**
    '6-311++G**'.................... 6-311++G**
    '6-311++G(3df,3pd)'............. 6-311++G(3df,3pd)

    'rcc-pVTZ'...................... rcc-pVTZ
    'rcc-pVQZ'...................... rcc-pVQZ
    'racc-pVDZ'..................... racc-pVDZ
    'racc-pVTZ'..................... racc-pVTZ
    'racc-pVQZ'..................... racc-pVQZ

    'pc-0'.......................... pc-0
    'pc-1'.......................... pc-1
    'pc-2'.......................... pc-2
    'pc-3'.......................... pc-3
    'pc-4'.......................... pc-4
    'pcJ-0'......................... pcJ-0
    'pcJ-1'......................... pcJ-1
    'pcJ-2'......................... pcJ-2
    'pcJ-3'......................... pcJ-3
    'pcJ-4'......................... pcJ-4
    'pcS-0'......................... pcS-0
    'pcS-1'......................... pcS-1
    'pcS-2'......................... pcS-2
    'pcS-3'......................... pcS-3
    'pcS-4'......................... pcS-4

    'cc-pVDZ'....................... cc-pVDZ
    'cc-pVTZ'....................... cc-pVTZ
    'cc-pVQZ'....................... cc-pVQZ
    'cc-pcVDZ'...................... cc-pcVDZ
    'cc-pcVTZ'...................... cc-pcVTZ
    'cc-pcVQZ'...................... cc-pcVQZ
    'aug-cc-pVDZ'................... aug-cc-pVDZ
    'aug-cc-pVTZ'................... aug-cc-pVTZ
    'aug-cc-pVQZ'................... aug-cc-pVQZ
    'aug-cc-pcVDZ'.................. aug-cc-pcVDZ
    'aug-cc-pcVTZ'.................. aug-cc-pcVTZ
    'aug-cc-pcVQZ'.................. aug-cc-pcVQZ

    'G3Large'....................... G3Large
    'G3MP2Large'.................... G3MP2Large

    'CRENBL'........................ CRENBL
    'CRENBS'........................ CRENBS
    'HWMB'.......................... HWMB
    'HWVDZ'......................... HWVDZ
    'LACVP'......................... LACVP
    'LANL2DZ'....................... LANL2DZ
    'SBKJC'......................... SBKJC
    'SRLC'.......................... SRLC
    'SRSC'.......................... SRSC

    'gen'........................... User-defined
    'Mixed'......................... Mixed

Description: Specifies the basis sets to be used.
Recommendation: : Consult literature and reviews to aid your selection.    '''
        if value == "":
            if "BASIS2" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS2"]
                print("Keyword removed.")
        elif value == "show":
            if "BASIS2" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS2"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["BASIS2"] = value.lower()

    def method(self, value="show"):
        '''
Name: METHOD
Type: STRING
Default: HF

Options:
    'Custom'........................ Custom
    'HF'............................ HF

    'B3LYP'......................... B3LYP
    'M06-2X'........................ M06-2X
    'wB97XD'........................ Omega-B97X-D
    'wB97XV'........................ Omega-B97X-V
    'BLYP'.......................... BLYP
    'CAMB3LYP'...................... CAM-B3LYP
    'EDF1'.......................... EDF1
    'B97-D'......................... B97-D
    'EDF2'.......................... EDF2
    'LDA'........................... LDA
    'M06'........................... M06
    'PBE'........................... PBE
    'PBE50'......................... PBE50

    'MP2'........................... MP2
    'MP2'........................... MP2[V]
    'RI-MP2'........................ RI-MP2
    'ATT-MP2'....................... ATT-MP2
    'SOS-MP2'....................... SOS-MP2
    'CCSD'.......................... CCSD
    'CCSD(T)'....................... CCSD(T)

    'TD-DFT'........................ TD-DFT
    'CIS'........................... CIS
    'CIS(D)'........................ CIS(D)
    'RI-CIS(D)'..................... RI-CIS(D)
    'SOS-CIS(D)'.................... SOS-CIS(D)
    'SOS-CIS(D0)'................... SOS-CIS(D0)

    'EOM-CCSD'...................... EOM-CCSD

    'ADC(1)'........................ ADC(1)
    'ADC(2)'........................ ADC(2)
    'ADC(2)-x'...................... ADC(2)-x
    'ADC(3)'........................ ADC(3)
    'CVS-ADC(1)'.................... CVS-ADC(1)
    'CVS-ADC(2)'.................... CVS-ADC(2)
    'CVS-ADC(2)-x'.................. CVS-ADC(2)-x
    'SOS-ADC(2)'.................... SOS-ADC(2)
    'SOS-ADC(2)-x'.................. SOS-ADC(2)-x

Description: The level of theory used in the calculation.
    '''
        if value == "":
            if "METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["METHOD"]
                print("Keyword removed.")
        elif value == "show":
            if "METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["METHOD"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["METHOD"] = value.lower()

    # --------------- End of computer generated keyword list ------------------

    # --------------  REM keywords added by hand (obsolete?) -----------------

    def aifdem(self, value="show"):
        if value == "":
            if "AIFDEM" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM"] = value.lower()

    def aifdem_ntothresh(self, value="show"):
        if value == "":
            if "AIFDEM_NTOTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM_NTOTHRESH"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM_NTOTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM_NTOTHRESH"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM_NTOTHRESH"] = value.lower()

    def aifdem_embed_range(self, value="show"):
        if value == "":
            if "AIFDEM_EMBED_RANGE" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM_EMBED_RANGE"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM_EMBED_RANGE" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM_EMBED_RANGE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM_EMBED_RANGE"] = value.lower()

    def aifdem_ctstates(self, value="show"):
        if value == "":
            if "AIFDEM_CTSTATES" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM_CTSTATES"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM_CTSTATES" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM_CTSTATES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM_CTSTATES"] = value.lower()

    def xpol(self, value="show"):
        if value == "":
            if "XPOL" in self.dict_of_keywords:
                del self.dict_of_keywords["XPOL"]
                print("Keyword removed.")
        elif value == "show":
            if "XPOL" in self.dict_of_keywords:
                return self.dict_of_keywords["XPOL"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XPOL"] = value.lower()

    def xpol_noscf(self, value="show"):
        if value == "":
            if "XPOL_NOSCF" in self.dict_of_keywords:
                del self.dict_of_keywords["XPOL_NOSCF"]
                print("Keyword removed.")
        elif value == "show":
            if "XPOL_NOSCF" in self.dict_of_keywords:
                return self.dict_of_keywords["XPOL_NOSCF"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XPOL_NOSCF"] = value.lower()

    def xpol_charge_type(self, value="show"):
        if value == "":
            if "XPOL_CHARGE_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["XPOL_CHARGE_TYPE"]
                print("Keyword removed.")
        elif value == "show":
            if "XPOL_CHARGE_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["XPOL_CHARGE_TYPE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["XPOL_CHARGE_TYPE"] = value.lower()

    def nto_pairs(self, value="show"):
        if value == "":
            if "NTO_PAIRS" in self.dict_of_keywords:
                del self.dict_of_keywords["NTO_PAIRS"]
                print("Keyword removed.")
        elif value == "show":
            if "NTO_PAIRS" in self.dict_of_keywords:
                return self.dict_of_keywords["NTO_PAIRS"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["NTO_PAIRS"] = value.lower()

    def max_scf_cycles(self, value="show"):
        if value == "":
            if "MAX_SCF_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["MAX_SCF_CYCLES"]
                print("Keyword removed.")
        elif value == "show":
            if "MAX_SCF_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["MAX_SCF_CYCLES"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["MAX_SCF_CYCLES"] = value.lower()

    def symmetry(self, value="show"):
        if value == "":
            if "SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY"]
                print("Keyword removed.")
        elif value == "show":
            if "SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["SYMMETRY"] = value.lower()

    def aifdem_frgm_read(self, value="show"):
        if value == "":
            if "AIFDEM_FRGM_READ" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM_FRGM_READ"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM_FRGM_READ" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM_FRGM_READ"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM_FRGM_READ"] = value.lower()

    def aifdem_frgm_write(self, value="show"):
        if value == "":
            if "AIFDEM_FRGM_WRITE" in self.dict_of_keywords:
                del self.dict_of_keywords["AIFDEM_FRGM_WRITE"]
                print("Keyword removed.")
        elif value == "show":
            if "AIFDEM_FRGM_WRITE" in self.dict_of_keywords:
                return self.dict_of_keywords["AIFDEM_FRGM_WRITE"]
            else:
                print("Value not set.")
        else:
            self.dict_of_keywords["AIFDEM_FRGM_WRITE"] = value.lower()

    # ------------------------ End of manual keyword list ----------------------

    def add(self, keyword, value):
        '''\nFor rem values without documenation herein, please add keyword and value manually'''
        self.dict_of_keywords[keyword.upper()] = value.lower()

    def remove(self, keyword):
        del self.dict_of_keywords[keyword.upper()]

    def clear(self):
        '''Removes all keywords from array.'''
        self.dict_of_keywords.clear()

    def __str__(self):
        str_ret = "$rem\n"
        for key, value in self.dict_of_keywords.items():
            str_ret += key.upper() + (
                    rem_array.__tabstop - len(key)) * " " + value + "\n"
        str_ret += "$end\n"
        return str_ret

    def info(self):
        print("Type: rem array")
        print("Keywords: " + str(len(self.dict_of_keywords)))


######################### REM_FRGM FRAGMENT ##############################

class rem_frgm_array(rem_array):
    __tabstop = 30

    def __init__(self, rem_init=""):
        self.dict_of_keywords = {}
        rem_init = rem_init.splitlines()
        if len(rem_init) != 0:
            for i in rem_init:
                i = i.split(" ")
                if len(i) == 0:
                    i = i.split("=")
                if i[0].startswith("$"):
                    continue
                self.add(i[0], i[1])

    def __str__(self):
        str_ret = "$rem_frgm\n"
        for key, value in self.dict_of_keywords.items():
            str_ret += key.upper() + (rem_frgm_array.__tabstop - len(
                key)) * " " + value + "\n"
        str_ret += "$end\n"
        return str_ret

    def info(self):
        print("Type: rem_frgm array")
        print("Keywords: " + str(len(self.dict_of_keywords)))
