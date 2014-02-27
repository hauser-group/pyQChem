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
import constants
import running_scripts

############################# RUNDATA  ###############################

class _rundata(object):
    def __init__(self):
        self.name=''
        self.loc53=''
        self.qchem=''
        self.nt=1
        self.np=1
        self.timestamp=False

    def __str__(self):
        ret_str = "Submission status summary:\n" + 26*"-" + "\n\n"
        if self.name!='':
            ret_str +=  "Filename is " + self.name + "\n"
        else:
            ret_str += "No filename provided, will use timestamp instead\n"
        if self.loc53!='':
            ret_str += "53.0 is stored at \'" + self.loc53 + "\'\n"
        if self.qchem!='':
            ret_str += "Q-Chem version is " + self.qchem + "\n"
        if self.nt>0:
            ret_str += "Using " + str(self.nt) + " threads\n"
        if self.np>0:
            ret_str += "Using " + str(self.np) + " cores\n"
        return ret_str

    def info(self):
        print self

########################### MULTIFILE  ##############################

class multifile(object):

    def __init__(self, jobs=[]):
        self.list_of_jobs=[]
        self.list_of_content=[]
        for k in jobs:
            self.add(k)

    def add(self,new_job):
        ''' Adds an inputfile to your batch object.'''
        if type(new_job) == type(inputfile()):
            self.list_of_jobs.append(new_job)
            self.list_of_content.append(new_job._jtype)
        else:
            print "Can only add inputfiles."

    def remove(self,position=0): #if not specified delete last
        ''' Removes an inputfile from your batch object. If no other specified the last is removed.'''
        del self.list_of_content[position-1] 
        del self.list_of_jobs[position-1] 

    def __str__(self):
        if self.list_of_jobs==[]:
            ret_str =  "empty"
        else:
            ret_str = self.list_of_jobs[0].__str__()
        if len(self.list_of_jobs)>1:
            for k in self.list_of_jobs[1:]:
                ret_str += "@@@\n\n" + k.__str__()
        return ret_str

    def write(self,filename):
        '''Writes the batch jobfile to disk.'''
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

    def run(self,name='',loc53='',qchem='',nt=1,np=1,timestamp=False):
        '''Makes Q-Chem process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        loc53 ...... 53.0 file location
        nt ......... number of threads
        np ......... number of processors.

        If nothing specified, pyQChem will fall back on information in the corresponding runinfo object.'''
        if name == '':
            name = self.runinfo.name
            loc53 = self.runinfo.loc53
            qchem = self.runinfo.qchem
            nt = self.runinfo.nt
            np = self.runinfo.np
            timestamp = self.runinfo.timestamp
        running_scripts._run(self,name,loc53,qchem,nt,np,timestamp)

########################### INPUTFILE  ##############################

class inputfile(object):
    
    def __init__(self, arrays=[]):
        self.list_of_arrays=[]
        self.list_of_content=[]
        self.runinfo=_rundata()
        self.__jtype="undef"
        for k in arrays:
            self.add(k)
            
    def add(self,new_array):
        ''' Adds an array to your inputfile object.'''
        if type(new_array) == type(rem_array()):
            self.rem = new_array
            if "rem" in self.list_of_content:
                index = self.list_of_content.index("rem")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("rem")
                self.list_of_arrays.append(new_array)
            self._jtype = new_array.jobtype() #rem variable "jobtype" defines type
        
        elif type(new_array) == type(mol_array()):
            self.molecule = new_array
            if "molecule" in self.list_of_content:
                index = self.list_of_content.index("molecule")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("molecule")
                self.list_of_arrays.append(new_array)
        
        elif type(new_array) == type(comment_array()):
            self.list_of_content.append("comment")
            self.list_of_arrays.append(new_array)

        elif type(new_array) == type(basis_array()):
            self.basis = new_array
            self.list_of_content.append("basis")
            self.list_of_arrays.append(new_array)

        elif type(new_array) == type(ecp_array()):
            self.ecp = new_array
            self.list_of_content.append("ecp")
            self.list_of_arrays.append(new_array)

        elif type(new_array) == type(_unsupported_array()):
            self.list_of_content.append(str(new_array.type))
            self.list_of_arrays.append(new_array)

        else:
            print "Array type unknown."

    def remove(self,position=0): #if not specified delete last
        ''' Removes an array from your inputfile object. If no other specified the last is removed.'''
        del self.list_of_content[position-1] 
        del self.list_of_arrays[position-1] 
             
    def __str__(self):
        ret_str = ""
        for k in self.list_of_arrays:
            ret_str += k.__str__() + "\n" 
        return ret_str
 
    def write(self,filename):
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

    def info(self):
        '''A quick overview of your inputfile.''' # Health check could be put here
        if "rem" and "molecule" in self.list_of_content:
            status = "valid"
        else:
            status = "invalid"

        print "Type: inputfile"
        print "Status: " + status

    def run(self,name='',loc53='',qchem='',nt=1,np=1,timestamp=False):
        '''Makes Q-Chem process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        loc53 ...... 53.0 file location
        nt ......... number of threads
        np ......... number of processors.

        If nothing specified, pyQChem will fall back on information in the corresponding runinfo object.'''
        if name == '':
            name = self.runinfo.name
            loc53 = self.runinfo.loc53
            qchem = self.runinfo.qchem
            nt = self.runinfo.nt
            np = self.runinfo.np
            timestamp = self.runinfo.timestamp
        running_scripts._run(self,name,loc53,qchem,nt,np,timestamp)
    
######################## INPUT FRAGMENTS ############################

class _array(object):

    def __init__(self,content="undef"):
        self.content = content

    def __str__(self):
        ret_str = "$undef\n"
        ret_str += str(self.content) + "\n"
        ret_str += "$end\n"
        return ret_str
    
    def write(self,filename):
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

##################### UNSUPPORTED FRAGMENT ##########################

class _unsupported_array(_array):

    def __init__(self,arraytype="undef"):
        self.content = []
        self.type = arraytype

    def add_line(self,line):
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
        
    def add_atom(self,line,position=0):  #if not specified add at the end
        '''Adds an atom to your Z-Matrix.'''
        if position ==0:
            self.__lines.append(line)
        else:
            self.__lines.insert(position-1,line)
        self.__Natoms += 1

    def remove_atom(self,position=0): #if not specified delete last
        '''Removes an atom from your Z-Matrix. Takes the last if no other specified.'''
        del self.__lines[position-1] 
        self.__Natoms -= 1

    def variable(self,key="variable",value="show"):
        '''Adds, changes or removes variable definitions.'''
        if value == "" and key in self.__variables:
            del self.__variables[key]
        elif value == "show":
            return self.__variables[key]
        else:
            self.__variables[key]=value        
        
    def __str__(self):
        ret_str = ""
        for k in self.__lines:
            ret_str += k + "\n"
        ret_str += "\n"
        for key,value in self.__variables.iteritems():
            ret_str += key + " "*(zmat.__tabstop-len(key)) + value + "\n"
        return ret_str
    

####################### CARTESIAN FRAGMENT ##########################

class cartesian(_array):
    def __init__(self,title="",atom_list=[]):
        import copy
        self.__title = title
        self.__Natoms = 0
        self.xyzs=[]
        self.com=_np.array([0.0,0.0,0.0])
        self.centroid=_np.array([0.0,0.0,0.0])
        self.list_of_atoms = copy.deepcopy(atom_list)
        if atom_list!=[]:
            self.__Natoms=len(atom_list)
            for i in xrange(self.__Natoms):
                x=self.list_of_atoms[i][1]
                y=self.list_of_atoms[i][2]
                z=self.list_of_atoms[i][3]
                self.xyzs.append(_np.array([float(x),float(y),float(z)]))
            self.xyzs=_np.array(self.xyzs)
            self.__center_of_mass()
    def fix(self):
        """This fixes any odd errors resulting from modifying the number of atoms"""
        self.__Natoms=len(self.list_of_atoms)
        if self.__Natoms==0:
            return
        self.xyzs=[]
        for i in xrange(self.__Natoms):
            x=self.list_of_atoms[i][1]
            y=self.list_of_atoms[i][2]
            z=self.list_of_atoms[i][3]
            self.xyzs.append(_np.array([float(x),float(y),float(z)]))
        self.xyzs=_np.array(self.xyzs)
        self.__center_of_mass()
        
    def __center_of_mass(self):
        """This computes the centroid and center of mass using standard atomic masses"""
        #print self.xyzs, self.__Natoms
        self.com=_np.array([0.0,0.0,0.0])
        self.centroid=_np.array([0.0,0.0,0.0])
        if len(self.xyzs)==0:
            return   
        total_mass=0.0
        self.centroid=sum(self.xyzs)/len(self.xyzs)
        wts=[constants.dict_of_atomic_masses[self.list_of_atoms[i][0].replace("@","")]  for i in xrange(self.__Natoms)]
        for i,atom in enumerate(self.xyzs):
            wt=wts[i]
            total_mass=total_mass+wt
            self.com=self.com+atom*wt
        self.centroid=_np.array([i/self.__Natoms for i in self.centroid])
        self.com=_np.array([i/total_mass for i in self.com])
    
    def title(self,title="show"):
        if title == "show":
            return self.__title
        else:
            self.__title=title
        
    def add_atom(self,name="H",x="0",y="0",z="0"):
        self.list_of_atoms.append([name,x,y,z])
        self.fix()
        self.__center_of_mass()
    
    def remove_atom(self,position=0):
        del self.list_of_atoms[position-1]  # First atom is atom 1
        self.fix()
        self.__center_of_mass()
    
    def ghost(self):
        atoms=[]
        for i in xrange(self.__Natoms):
            atoms.append(['@'+self.list_of_atoms[i][0],self.list_of_atoms[i][1],self.list_of_atoms[i][2],self.list_of_atoms[i][3]])
        return atoms

    def atoms(self):
        for i, k in enumerate(self.list_of_atoms):
            print str(i+1) + ":\t" +  k[0] + "\t" + k[1] + "\t" + k[2] + "\t" + k[3]

    def print_centroid(self):
        print str(self.centroid[0])+'\t'+str(self.centroid[1])+'\t'+str(self.centroid[2])
        
    def print_center_of_mass(self):
        print str(self.com[0])+'\t'+str(self.com[1])+'\t'+str(self.com[2])
        
    def move(self,dir,amt=1.0):
        dir=_np.array(dir)
        for i in xrange(self.__Natoms):
            self.xyzs[i]=self.xyzs[i]+dir*amt
            self.list_of_atoms[i][1]=str(self.xyzs[i][0])
            self.list_of_atoms[i][2]=str(self.xyzs[i][1])
            self.list_of_atoms[i][3]=str(self.xyzs[i][2])
        self.__center_of_mass()

    def __str__(self):
        str_ret = str(self.__Natoms) + "\n" + self.__title + "\n"
        for k in self.list_of_atoms:
            str_ret += k[0] + "    " + k[1] + "    " + k[2] + "    " + k[3] + "\n"
        return str_ret
    def __add__(self,other):
        if type(other)==type([]):  #let's move the atoms
            self.move(other,1.0)
        if type(other)==type(self.com):  #let's move the atoms using a numpy array
            self.move(other,1.0)
        if type(other)==type(self):                      #merge two cartesians
            atoms=self.list_of_atoms+other.list_of_atoms  
            return cartesian(atom_list=atoms)
    def __radd__(self,other):  #reverse of above
        if type(other)==type([]):
            self.move(other)
        if type(other)==type(self):
            atoms=self.list_of_atoms+other.list_of_atoms
            return cartesian(atom_list=atoms)
            

####################### TINKER FRAGMENT ##########################

class tinker(cartesian):

    def __init__(self,title=""):
        self.__title = title
        self.__Natoms = 0
        self.list_of_atoms = []
        self.dict_of_types = {}  # dictionary of atom types

    def title(self,title="show"):
        if title == "show":
            return self.__title
        else:
            self.__title=title
            
    def remove_atom(self,position):
        del self.list_of_atoms[position-1]  # First atom is atom 1
        self.__Natoms -= 1


    def add_atom(self,name="H",x="0",y="0",z="0",atomtype="0",con1="0",con2="0",con3="0",con4="0"):
        if not self.dict_of_types.has_key(name):
            self.dict_of_types[name]=atomtype
        self.list_of_atoms.append([name,x,y,z,atomtype,con1,con2,con3,con4])
        self.__Natoms += 1

    def change_type(self,atomtype="none",value="show"):
        '''Changes the definition number of atom "atomtype" to "value".'''
        if value == "" and self.dict_of_types.has_key(atomtype):
            del self.dict_of_types[atomtype]
        elif value == "show":
            return self.dict_of_types
        else:
            self.dict_of_types[atomtype]=value 
            # atomtype definition has changed, list_of_atoms needs to be updated:
            Nchanges = 0
            for k in self.list_of_atoms:
                if k[0]==atomtype:
                    k[4] = value
                    Nchanges += 1
            print "Atomtype definition has changed, " + str(Nchanges) + " atoms updated in list_of_atoms."
    
    def atoms(self):
        for i, k in enumerate(self.list_of_atoms):
            print str(i+1) + ":\t" +  k[0] + "    " + k[1] + "    " + k[2] + \
            "    " + k[3] + "    " + k[4] + "    " + k[5] + "    " + k[6] + "    " + k[7] + "    " + k[8]
        
    def __str__(self):
        str_ret = str(self.__Natoms) + "\t" + self.__title + "\n"
        for i,k in enumerate(self.list_of_atoms):
            str_ret += str(i+1) + "\t" + k[0] + "    " + k[1] + "    " + k[2] + \
            "    " + k[3] + "    " + k[4] + "    " + k[5] + "    " + k[6] + "    " + \
            k[7] + "    " + k[8] + "\n"
        return str_ret
        
######################### MOL FRAGMENT ##############################

class mol_array(_array):
            
    def __init__(self,geometry=""):
        if geometry == "":
            geometry = cartesian()
        self.content = {"CHARGE":"0","MULTIPLICITY":"1","GEOMETRY":geometry}
     
    def charge(self,value="show"):
        '''Total charge of the molecule.'''
        if value == "show":
            return self.content["CHARGE"]
        else:
            self.content["CHARGE"]=value
             
    def multiplicity(self,value="show"):
        '''Spin multiplicity of the molecule.'''
        if value == "show":
            return self.content["MULTIPLICITY"]
        else:
            self.content["MULTIPLICITY"]=value
        
    def geometry(self,value="show"):
        '''Reads xyz, txyz or zmat coordinate array.'''
        if value == "show":
            return self.content["GEOMETRY"]
        elif value == "read":
            self.content["GEOMETRY"]=value
        else:
            if type(value)==type(cartesian()) or type(value)==type(zmat()) or type(value)==type(tinker()):
                self.content["GEOMETRY"]=value
            else:
                print "Only cartesian, tinker or zmat arrays can be added here."
       
    def clear(self):
        self.content = {"CHARGE":"0","MULTIPLICITY":"1","GEOMETRY":""}
        
    def __str__(self):
        if self.content["GEOMETRY"]=="read":
            str_ret = "$molecule\nread\n$end\n"
        else:
            str_ret = "$molecule\n" + self.content["CHARGE"] + " " + self.content["MULTIPLICITY"] + "\n"
            if type(self.content["GEOMETRY"])==type(cartesian()):
                for k in (self.content["GEOMETRY"]).list_of_atoms:
                    str_ret += k[0] + "    " + k[1] + "    " + k[2] + "    " + k[3] + "\n"
            elif type(self.content["GEOMETRY"])==type(zmat()):
                str_ret += (self.content["GEOMETRY"]).__str__() 
            elif type(self.content["GEOMETRY"])==type(tinker()):
                for k in (self.content["GEOMETRY"]).list_of_atoms:
                    str_ret +=  k[0] + "    " + k[1] + "   " + k[2] + \
                    "    " + k[3] + "    " + k[4] + "    " + k[5] + \
                    "    " + k[6] + "    " + k[7] + "    " + k[8] + "\n"
            str_ret += "$end\n"
        return str_ret
         
    def info(self):
        switch = 0
        if type(self.content["GEOMETRY"])==type(cartesian()):
            coor_type = "cartesian coordinates"
            switch = 1
        elif type(self.content["GEOMETRY"])==type(zmat()):
            coor_type = "Z-Matrix"
            switch = 1
        elif type(self.content["GEOMETRY"])==type(tinker()):
            coor_type = "Tinker"
            switch = 1
        else:
            coor_type = "empty"
        print "Type: molecule array, " + coor_type
        if switch == 1:
            print "Number of atoms: " + str(len((self.content["GEOMETRY"]).list_of_atoms))

######################### BASIS FRAGMENT ############################

class basis_array(_array):

    def __init__(self):
        self.dict_of_atoms = {}
    
    def add(self,atom,line):
        if self.dict_of_atoms.has_key(atom):
            self.dict_of_atoms[atom].append(line)
        else:
            self.dict_of_atoms[atom]=list()
            self.dict_of_atoms[atom].append(line)

    def remove(self,atom):
        if self.dict_of_atoms.has_key(atom):
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

    def __init__(self,content=""):
        self.content = content

    def __str__(self):
        ret_str = "$comment\n" 
        lines = filter(bool,self.content.split("\n")) #remove empty list entries
        for line in lines:
            ret_str += line.strip() + "\n" 
        ret_str += "$end\n"
        return ret_str

######################### REM FRAGMENT ##############################

class rem_array(_array):
    
    __tabstop = 30
        
    def __init__(self):
        self.dict_of_keywords = {"JOBTYPE":"sp","EXCHANGE":"hf"}
    
    # -------------- Computer-generated List of REM keywords  -----------------

    def cc_dip(self,value="show"):
        '''\nName: CC_DIP\nType: INTEGER\nDefault: 0\nOptions: 0:1\nDescription: Initializes a EOM-DIP-CCSD calculation\n    '''
        if value == "":
            if "CC_DIP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIP"]=value.lower()


    def cc_do_triples(self,value="show"):
        '''\nName: CC_DO_TRIPLES\nType: INTEGER\nDefault: 0\nOptions: 0:1\nDescription: This keyword initializes a EOM-CC(2,3) calculation. If {CC_IP_PROPER} is set then EOM-IP-CC(2,3) states are calculated.\n    '''
        if value == "":
            if "CC_DO_TRIPLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_TRIPLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DO_TRIPLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_TRIPLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DO_TRIPLES"]=value.lower()


    def cc_iterate_on(self,value="show"):
        '''\nName: CC_ITERATE_ON\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: In active space calculations, use a mixed iteration procedure if the value is greater than 0.  Then if the RMS orbital gradient is larger than the value of CC_THETA_GRAD_THRESH, micro-iterations will be performed to converge the occupied-virtual mixing angles for the current active-space. The maximum number of space iterations is given by this option.\nRecommendation: : Can be useful for non-convergent active space calculations    '''
        if value == "":
            if "CC_ITERATE_ON" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ITERATE_ON"]
                print "Keyword removed."
        elif value == "show":
            if "CC_ITERATE_ON" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ITERATE_ON"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_ITERATE_ON"]=value.lower()


    def cc_iterate_ov(self,value="show"):
        '''\nName: CC_ITERATE_OV\nType: INTEGER\nDefault: 0\nOptions: 0   \nDescription: In active space calculations, use a mixed iteration procedure if the value is greater than 0. Then, if the RMS orbital gradient is larger than the value of CC_THETA_GRAD_THRESH, micro-iterations will be performed to converge the occupied-virtual mixing angles for the current active space. The maximum number of such iterations is given by this option. \nRecommendation: : Can be useful for non-convergent active space calculations.    '''
        if value == "":
            if "CC_ITERATE_OV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ITERATE_OV"]
                print "Keyword removed."
        elif value == "show":
            if "CC_ITERATE_OV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ITERATE_OV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_ITERATE_OV"]=value.lower()


    def cc_orbs_per_block(self,value="show"):
        '''\nName: CC_ORBS_PER_BLOCK\nType: INTEGER\nDefault: 0\nOptions: 16\nDescription: Specifies target (and maximum) size of blocks in orbital space.\n    '''
        if value == "":
            if "CC_ORBS_PER_BLOCK" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_ORBS_PER_BLOCK"]
                print "Keyword removed."
        elif value == "show":
            if "CC_ORBS_PER_BLOCK" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_ORBS_PER_BLOCK"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_ORBS_PER_BLOCK"]=value.lower()


    def cc_preconv_sd(self,value="show"):
        '''\nName: CC_PRECONV_SD\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Solves the EOM-CCSD equations, prints energies, then uses EOM-CCSD vectors  as initial vectors in EOM-CC(2,3). Very convenient for calculations using energy additivity schemes.\nRecommendation: : Turning this option on is recommended    '''
        if value == "":
            if "CC_PRECONV_SD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_SD"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONV_SD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_SD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONV_SD"]=value.lower()


    def cc_preconv_t2z_each(self,value="show"):
        '''\nName: CC_PRECONV_T2Z_EACH\nType: INTEGER\nDefault: 1\nOptions: 0:0   \nDescription: Whether to pre-converge the cluster amplitudes before each change of the  orbitals in optimized orbital coupled-cluster methods. The maximum number of  iterations in this pre-convergence procedure is given by the value of this  parameter.\nRecommendation: : A very slow last resort option for jobs that do not converge.    '''
        if value == "":
            if "CC_PRECONV_T2Z_EACH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_T2Z_EACH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONV_T2Z_EACH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_T2Z_EACH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONV_T2Z_EACH"]=value.lower()


    def cc_preconv_t2z(self,value="show"):
        '''\nName: CC_PRECONV_T2Z\nType: INTEGER\nDefault: 1\nOptions: 0:0  \nDescription: Whether to pre-converge the cluster amplitudes before beginning orbital  optimization in optimized orbital cluster methods.\nRecommendation: : Experiment with this option in cases of convergence failure.    '''
        if value == "":
            if "CC_PRECONV_T2Z" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_T2Z"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONV_T2Z" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_T2Z"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONV_T2Z"]=value.lower()


    def cc_qccd_theta_switch(self,value="show"):
        '''\nName: CC_QCCD_THETA_SWITCH\nType: INTEGER\nDefault: 0\nOptions: 2   \nDescription: QCCD calculations switch from OD to QCCD when the rotation gradient is below this threshold [10-n]\n    '''
        if value == "":
            if "CC_QCCD_THETA_SWITCH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_QCCD_THETA_SWITCH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_QCCD_THETA_SWITCH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_QCCD_THETA_SWITCH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_QCCD_THETA_SWITCH"]=value.lower()


    def cc_refsym(self,value="show"):
        '''\nName: CC_REFSYM\nType: INTEGER\nDefault: 0\nOptions: -1 \nDescription: Together with CC_STATE_DERIV, selects which EOM  state is to be considered  for optimization or property calculations. When transition properties are requested, the transition properties will be calculated between this state and all other EOM states.\n    '''
        if value == "":
            if "CC_REFSYM" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REFSYM"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REFSYM" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REFSYM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REFSYM"]=value.lower()


    def cc_reset_theta(self,value="show"):
        '''\nName: CC_RESET_THETA\nType: INTEGER\nDefault: 0\nOptions: 15\nDescription: The reference MO coefficient matrix is reset every n iterations to help   overcome problems associated with the theta metric as theta becomes large. \n    '''
        if value == "":
            if "CC_RESET_THETA" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESET_THETA"]
                print "Keyword removed."
        elif value == "show":
            if "CC_RESET_THETA" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESET_THETA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_RESET_THETA"]=value.lower()


    def cc_restr_ampl(self,value="show"):
        '''\nName: CC_RESTR_AMPL\nType: INTEGER\nDefault: 1\nOptions: 0:1\nDescription: Controls the restriction on amplitudes is there are restricted orbitals\n    '''
        if value == "":
            if "CC_RESTR_AMPL" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTR_AMPL"]
                print "Keyword removed."
        elif value == "show":
            if "CC_RESTR_AMPL" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTR_AMPL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_RESTR_AMPL"]=value.lower()


    def cc_restr_triples(self,value="show"):
        '''\nName: CC_RESTR_TRIPLES\nType: INTEGER\nDefault: 0\nOptions: 0:1\nDescription: Controls which space the triples correction is computed in\n    '''
        if value == "":
            if "CC_RESTR_TRIPLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTR_TRIPLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_RESTR_TRIPLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTR_TRIPLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_RESTR_TRIPLES"]=value.lower()


    def cc_rest_occ(self,value="show"):
        '''\nName: CC_REST_OCC\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Sets the number of restricted occupied orbitals including frozen occupied  orbitals.\n    '''
        if value == "":
            if "CC_REST_OCC" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REST_OCC"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REST_OCC" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REST_OCC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REST_OCC"]=value.lower()


    def cc_rest_vir(self,value="show"):
        '''\nName: CC_REST_VIR\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Sets the number of restricted virtual orbitals including frozen virtual  orbitals.\n    '''
        if value == "":
            if "CC_REST_VIR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REST_VIR"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REST_VIR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REST_VIR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REST_VIR"]=value.lower()


    def cc_spin_flip_ms(self,value="show"):
        '''\nName: CC_SPIN_FLIP_MS\nType: INTEGER\nDefault: 0\nOptions: 0 :2\nDescription: This option is only used in EOM-SF using quintet references  and including triple excitations. By default, SF flips the spin of one  electron.  One can ask to flip the spins of two   electrons by specifying CC_SPIN_FLIP_MS = 1\nRecommendation: : This option can be useful when starting from quintet references - though this is not typical for EOM-SF.    '''
        if value == "":
            if "CC_SPIN_FLIP_MS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SPIN_FLIP_MS"]
                print "Keyword removed."
        elif value == "show":
            if "CC_SPIN_FLIP_MS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SPIN_FLIP_MS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_SPIN_FLIP_MS"]=value.lower()


    def cc_state_deriv(self,value="show"):
        '''\nName: CC_STATE_DERIV\nType: INTEGER\nDefault: 0\nOptions: -1 \nDescription: Selects which EOM or CIS(D) state is to be considered  for optimization or property calculations.\n    '''
        if value == "":
            if "CC_STATE_DERIV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_DERIV"]
                print "Keyword removed."
        elif value == "show":
            if "CC_STATE_DERIV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_DERIV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_STATE_DERIV"]=value.lower()


    def cc_theta_conv(self,value="show"):
        '''\nName: CC_THETA_CONV\nType: INTEGER\nDefault: 0\nOptions: 5 \nDescription: Convergence criterion on the RMS difference between successive sets of  orbital rotation angles [10-n].\nRecommendation: : Use default    '''
        if value == "":
            if "CC_THETA_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "CC_THETA_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_THETA_CONV"]=value.lower()


    def cc_theta_grad_conv(self,value="show"):
        '''\nName: CC_THETA_GRAD_CONV\nType: INTEGER\nDefault: 0\nOptions: 7   \nDescription: Convergence desired on the RMS gradient of the energy with respect to   orbital rotation angles [10-n]. \nRecommendation: : Use default    '''
        if value == "":
            if "CC_THETA_GRAD_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_GRAD_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "CC_THETA_GRAD_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_GRAD_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_THETA_GRAD_CONV"]=value.lower()


    def cc_theta_grad_thresh(self,value="show"):
        '''\nName: CC_THETA_GRAD_THRESH\nType: INTEGER\nDefault: 0\nOptions: 2\nDescription: RMS orbital gradient threshold [10-n] above which mixed iterations  are performed in active space calculations if CC_ITERATE_OV is  TRUE.\nRecommendation: : Can be made smaller if convergence difficulties are encountered.    '''
        if value == "":
            if "CC_THETA_GRAD_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_GRAD_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_THETA_GRAD_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_GRAD_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_THETA_GRAD_THRESH"]=value.lower()


    def cc_tmpbuffsize(self,value="show"):
        '''\nName: CC_TMPBUFFSIZE\nType: INTEGER\nDefault: 0\nOptions: 3\% of \remvar{MEM\_TOTAL}\nDescription: Maximum size, in Mb, of additional buffers for temporary arrays used to work with individual blocks or matrices.\nRecommendation: : Should not be smaller than the size of the largest possible block.    '''
        if value == "":
            if "CC_TMPBUFFSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TMPBUFFSIZE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_TMPBUFFSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TMPBUFFSIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_TMPBUFFSIZE"]=value.lower()


    def gvb_orb_conv(self,value="show"):
        '''\nName: GVB_ORB_CONV\nType: INTEGER\nDefault: 0\nOptions: 5\nDescription: The GVB-CC wave function is considered converged when the root-mean-square  orbital gradient and orbital step sizes are less than  10-GVB_ORB_CONV. Adjust THRESH simultaneously.\nRecommendation: : Use 6 for PP(2) jobs or geometry optimizations. Tighter convergence (i.e. 7 or higher) cannot always be reliably achieved.    '''
        if value == "":
            if "GVB_ORB_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_ORB_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_ORB_CONV"]=value.lower()


    def gvb_orb_max_iter(self,value="show"):
        '''\nName: GVB_ORB_MAX_ITER\nType: INTEGER\nDefault: 0\nOptions: 256\nDescription: Controls the number of orbital iterations allowed in GVB-CC calculations. Some jobs, particularly unrestricted PP jobs can require 500-1000 iterations.\nRecommendation: : Default is typically adequate, but some jobs, particularly UPP jobs, can  require 500-1000 iterations if converged tightly.    '''
        if value == "":
            if "GVB_ORB_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_MAX_ITER"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_ORB_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_MAX_ITER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_ORB_MAX_ITER"]=value.lower()


    def gvb_orb_scale(self,value="show"):
        '''\nName: GVB_ORB_SCALE\nType: INTEGER\nDefault: 0\nOptions: 1000  \nDescription: Scales the default orbital step size by n/1000.\nRecommendation: : Default is usually fine, but for some stretched geometries it can help with convergence to use smaller values.    '''
        if value == "":
            if "GVB_ORB_SCALE" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_ORB_SCALE"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_ORB_SCALE" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_ORB_SCALE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_ORB_SCALE"]=value.lower()


    def gvb_restart(self,value="show"):
        '''\nName: GVB_RESTART\nType: STRING\nDefault: 0\nOptions: FALSE\nDescription: Restart a job from previously-converged GVB-CC orbitals.\nRecommendation: : Useful when trying to converge to the same GVB solution at slightly different geometries, for example.    '''
        if value == "":
            if "GVB_RESTART" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_RESTART"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_RESTART" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_RESTART"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_RESTART"]=value.lower()


    def gvb_unrestricted(self,value="show"):
        '''\nName: GVB_UNRESTRICTED\nType: STRING\nDefault: 0\nOptions: same value as UNRESTRICTED\nDescription: Controls restricted versus unrestricted PP jobs. Usually handled  automatically.\nRecommendation: : Set this variable explicitly only to do a UPP job from an RHF or ROHF initial guess.    '''
        if value == "":
            if "GVB_UNRESTRICTED" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_UNRESTRICTED"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_UNRESTRICTED" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_UNRESTRICTED"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_UNRESTRICTED"]=value.lower()


    def gvb_print(self,value="show"):
        '''\nName: GVB_PRINT\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Controls the amount of information printed during a GVB-CC job.\nRecommendation: : Should never need to go above 0 or 1.    '''
        if value == "":
            if "GVB_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_PRINT"]=value.lower()


    def mom_start(self,value="show"):
        '''\nName: MOM_START\nType: INTEGER\nDefault: 0\nOptions: 0 (FALSE)\nDescription: Determines when MOM is switched on to stabilize DIIS iterations.\nRecommendation: : Set to 1 if preservation of  initial orbitals is desired.  If MOM is to be  used to aid convergence, an SCF without MOM should be run to determine when  the SCF starts oscillating.  MOM should be set to start just before the  oscillations.    '''
        if value == "":
            if "MOM_START" in self.dict_of_keywords:
                del self.dict_of_keywords["MOM_START"]
                print "Keyword removed."
        elif value == "show":
            if "MOM_START" in self.dict_of_keywords:
                return self.dict_of_keywords["MOM_START"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOM_START"]=value.lower()


    def nvo_lin_convergence(self,value="show"):
        '''\nName: NVO_LIN_CONVERGENCE\nType: INTEGER\nDefault: 0\nOptions: 3\nDescription: Target error factor in the preconditioned conjugate gradient solver of the single-excitation amplitude equations.\nRecommendation: : Solution of the single-excitation amplitude equations is considered converged if the maximum residual is less than 10-n multiplied by the current DIIS error. For the ARS correction, n is automatically set to 1 since the locally-projected DIIS error is normally several orders of magnitude smaller than the full DIIS error.    '''
        if value == "":
            if "NVO_LIN_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_LIN_CONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_LIN_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_LIN_CONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_LIN_CONVERGENCE"]=value.lower()


    def nvo_lin_max_ite(self,value="show"):
        '''\nName: NVO_LIN_MAX_ITE\nType: INTEGER\nDefault: 0\nOptions: 30\nDescription: Maximum number of iterations in the preconditioned conjugate gradient solver of the single-excitation amplitude equations.\n    '''
        if value == "":
            if "NVO_LIN_MAX_ITE" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_LIN_MAX_ITE"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_LIN_MAX_ITE" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_LIN_MAX_ITE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_LIN_MAX_ITE"]=value.lower()


    def nvo_method(self,value="show"):
        '''\nName: NVO_METHOD\nType: INTEGER\nDefault: 0\nOptions: 9\nDescription: Sets method to be used to converge solution of the single-excitation amplitude equations.\nRecommendation: : Experimental option. Use default.    '''
        if value == "":
            if "NVO_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_METHOD"]=value.lower()


    def nvo_truncate_dist(self,value="show"):
        '''\nName: NVO_TRUNCATE_DIST\nType: INTEGER\nDefault: 0\nOptions: -1:0:1:2\nDescription: Specifies which atomic blocks of the Fock matrix are used to construct the preconditioner.\nRecommendation: : This option does not affect the final result. However, it affects the rate of the PCG algorithm convergence. For small systems use default.    '''
        if value == "":
            if "NVO_TRUNCATE_DIST" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_TRUNCATE_DIST"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_TRUNCATE_DIST" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_TRUNCATE_DIST"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_TRUNCATE_DIST"]=value.lower()


    def nvo_truncate_precond(self,value="show"):
        '''\nName: NVO_TRUNCATE_PRECOND\nType: INTEGER\nDefault: 0\nOptions: 2\nDescription: Specifies which atomic blocks of the Fock matrix are used to construct the preconditioner. This variable is used only if NVO_TRUNCATE_DIST is set to -2.\nRecommendation: : Use default. Increasing n improves convergence of the PCG algorithm but overall may slow down the calculations.    '''
        if value == "":
            if "NVO_TRUNCATE_PRECOND" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_TRUNCATE_PRECOND"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_TRUNCATE_PRECOND" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_TRUNCATE_PRECOND"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_TRUNCATE_PRECOND"]=value.lower()


    def nvo_uvv_maxpwr(self,value="show"):
        '''\nName: NVO_UVV_MAXPWR\nType: INTEGER\nDefault: 0\nOptions: 10\nDescription: Controls convergence of the Taylor series when calculating the Uvv block from the single-excitation amplitudes. If the series is not converged at the nth term, more expensive direct inversion is used to calculate the Uvv block.\n    '''
        if value == "":
            if "NVO_UVV_MAXPWR" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_UVV_MAXPWR"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_UVV_MAXPWR" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_UVV_MAXPWR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_UVV_MAXPWR"]=value.lower()


    def nvo_uvv_precision(self,value="show"):
        '''\nName: NVO_UVV_PRECISION\nType: INTEGER\nDefault: 0\nOptions: 11\nDescription: Controls convergence of the Taylor series when calculating the Uvv block from the single-excitation amplitudes. Series is considered converged when the maximum element of the term is less than 10-n.\nRecommendation: : NVO_UVV_PRECISION must be the same as or larger than THRESH.    '''
        if value == "":
            if "NVO_UVV_PRECISION" in self.dict_of_keywords:
                del self.dict_of_keywords["NVO_UVV_PRECISION"]
                print "Keyword removed."
        elif value == "show":
            if "NVO_UVV_PRECISION" in self.dict_of_keywords:
                return self.dict_of_keywords["NVO_UVV_PRECISION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NVO_UVV_PRECISION"]=value.lower()


    def print_dist_matrix(self,value="show"):
        '''\nName: PRINT_DIST_MATRIX\nType: INTEGER\nDefault: 1\nOptions: 0:15\nDescription: Controls the printing of the inter-atomic distance matrix\nRecommendation: : Use default unless distances are required for large systems    '''
        if value == "":
            if "PRINT_DIST_MATRIX" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_DIST_MATRIX"]
                print "Keyword removed."
        elif value == "show":
            if "PRINT_DIST_MATRIX" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_DIST_MATRIX"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PRINT_DIST_MATRIX"]=value.lower()


    def rc_r0(self,value="show"):
        '''\nName: RC_R0\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Determines the parameter in the Gaussian weight function used to smooth the  density at the nuclei.\nRecommendation: : We recommend value of 250 for a typical spit valence basis.  For basis sets  with increased flexibility in the nuclear vicinity the smaller values of r_0  also yield adequate spin density.    '''
        if value == "":
            if "RC_R0" in self.dict_of_keywords:
                del self.dict_of_keywords["RC_R0"]
                print "Keyword removed."
        elif value == "show":
            if "RC_R0" in self.dict_of_keywords:
                return self.dict_of_keywords["RC_R0"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RC_R0"]=value.lower()


    def svp_cavity_conv(self,value="show"):
        '''\nName: SVP_CAVITY_CONV\nType: INTEGER\nDefault: 0\nOptions: 10\nDescription: Determines the convergence value of the iterative isodensity cavity procedure.\nRecommendation: : The default value unless convergence problems arise.    '''
        if value == "":
            if "SVP_CAVITY_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_CAVITY_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "SVP_CAVITY_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_CAVITY_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP_CAVITY_CONV"]=value.lower()


    def rpath_max_cycles(self,value="show"):
        '''\nName: RPATH_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 1:500:20:1\nDescription: Specifies the maximum number of points to find on the reaction path.\nRecommendation: : Use more points if the minimum is desired, but not reached using the default.    '''
        if value == "":
            if "RPATH_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_MAX_CYCLES"]=value.lower()


    def rpath_tol_displacement(self,value="show"):
        '''\nName: RPATH_TOL_DISPLACEMENT\nType: INTEGER\nDefault: 2\nOptions: 0.0001:0.2000:0.0050:0.0001\nDescription: Specifies the convergence threshold (in a.u.) for the step. If a step size is chosen by the algorithm that is smaller than this, the path is deemed to have reached the minimum.\n    '''
        if value == "":
            if "RPATH_TOL_DISPLACEMENT" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_TOL_DISPLACEMENT" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_TOL_DISPLACEMENT"]=value.lower()


    def rpath_print(self,value="show"):
        '''\nName: RPATH_PRINT\nType: INTEGER\nDefault: 2\nOptions: 1:5:2:1\nDescription: Controls the level of output for a reaction coordinate calculation.\nRecommendation: : Use default, little additional information is printed at higher levels. Most of the output arises from the multiple single point calculations that are performed along the reaction pathway.    '''
        if value == "":
            if "RPATH_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_PRINT"]=value.lower()


    def scf_convergence(self,value="show"):
        '''\nName: SCF_CONVERGENCE\nType: INTEGER\nDefault: 2\nOptions: 0:12:5:1\nDescription: SCF is considered converged when the wavefunction error is less that 10-SCF_CONVERGENCE. Adjust the value of THRESH at the same time. Note that in Q-Chem 3.0 the DIIS error is measured by the maximum error rather than the RMS error as in previous versions.\nRecommendation: : Tighter criteria for geometry optimization and vibration analysis. Larger values provide more significant figures, at greater computational cost.    '''
        if value == "":
            if "SCF_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_CONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_CONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_CONVERGENCE"]=value.lower()


    def scf_guess(self,value="show"):
        '''\nName: SCF_GUESS\nType: STRING\nDefault: 0\nOptions: SAD  :CORE:GWH:READ\nDescription: Specifies the initial guess procedure to use for the SCF.\nRecommendation: : SAD guess for standard basis sets. For general basis sets, it is best to use the BASIS2 $rem. Alternatively, try the GWH or core Hamiltonian guess. For ROHF it can be useful to READ guesses from an SCF calculation on the corresponding cation or anion. Note that because the density is made spherical, this may favor an undesired state for atomic systems, especially transition metals.    '''
        if value == "":
            if "SCF_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_GUESS"]=value.lower()


    def scf_algorithm(self,value="show"):
        '''\nName: SCF_ALGORITHM\nType: STRING\nDefault: 0\nOptions: DIIS      :DM:GDM:RCA:ROOTHAAN:DIIS_DM:DIIS_GDM:RCA_DIIS\nDescription: Selects the algorithm to use for converging the SCF.\nRecommendation: : Use DIIS unless performing a restricted open-shell calculation, in which case GDM is recommended. If DIIS fails to find a reasonable approximate solution in the initial iterations, RCA_DIIS is the recommended fallback option. If DIIS approaches the correct solution but fails to finally converge, DIIS_GDM is the recommended fallback.     '''
        if value == "":
            if "SCF_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_ALGORITHM"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_ALGORITHM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_ALGORITHM"]=value.lower()


    def diis_print(self,value="show"):
        '''\nName: DIIS_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:4 :0:1\nDescription: Controls the output from DIIS SCF optimization:
1: Chosen method and DIIS coefficients
2: Level 1 + print changes in multipole moments
3: Level 2 + multipole moments
4: Level 3 + extrapolated Fock matrices\n    '''
        if value == "":
            if "DIIS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "DIIS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIIS_PRINT"]=value.lower()


    def varthresh(self,value="show"):
        '''\nName: VARTHRESH\nType: INTEGER\nDefault: 2\nOptions: 0   :12:0:1\nDescription: Controls the temporary integral cut-off threshold. The variable threshold is set to 10-VARTHRESH? DIIS_error\nRecommendation: : 3 has been found to be a practical level, and can slightly speed up SCF evaluation.    '''
        if value == "":
            if "VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["VARTHRESH"]
                print "Keyword removed."
        elif value == "show":
            if "VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["VARTHRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["VARTHRESH"]=value.lower()


    def aimd_steps(self,value="show"):
        '''\nName: AIMD_STEPS\nType: INTEGER\nDefault: 2\nOptions: 0:500:0:1\nDescription: Specifies the requested number of molecular dynamics steps.\n    '''
        if value == "":
            if "AIMD_STEPS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_STEPS"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_STEPS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_STEPS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_STEPS"]=value.lower()


    def ao2mo_disk(self,value="show"):
        '''\nName: AO2MO_DISK\nType: INTEGER\nDefault: 2\nOptions: 0:8000:2000:1\nDescription: Sets the amount of disk space (in megabytes) available for MP2 calculations. 
\nRecommendation: : Should be set as large as possible, as discussed in the manual.    '''
        if value == "":
            if "AO2MO_DISK" in self.dict_of_keywords:
                del self.dict_of_keywords["AO2MO_DISK"]
                print "Keyword removed."
        elif value == "show":
            if "AO2MO_DISK" in self.dict_of_keywords:
                return self.dict_of_keywords["AO2MO_DISK"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AO2MO_DISK"]=value.lower()


    def cc_canonize_final(self,value="show"):
        '''\nName: CC_CANONIZE_FINAL\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether to semi-canonicalize orbitals at the end of the ground state calculation.\nRecommendation: : Should not normally have to be altered.    '''
        if value == "":
            if "CC_CANONIZE_FINAL" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE_FINAL"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CANONIZE_FINAL" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE_FINAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CANONIZE_FINAL"]=value.lower()


    def cc_canonize(self,value="show"):
        '''\nName: CC_CANONIZE\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Whether to semi-canonicalize orbitals at the start of the calculation (i.e. Fock matrix is diagonalized in each orbital subspace)\nRecommendation: : Should not normally have to be altered.    '''
        if value == "":
            if "CC_CANONIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CANONIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CANONIZE"]=value.lower()


    def cc_convergence(self,value="show"):
        '''\nName: CC_CONVERGENCE\nType: INTEGER\nDefault: 2\nOptions: 0:12:8:1\nDescription: Overall convergence criterion for the coupled-cluster codes. This is designed to ensure at least n significant digits in the calculated energy, and automatically sets the other convergence-related variables (CC_E_CONV, CC_T_CONV, CC_THETA_CONV, CC_THETA_GRAD_CONV, CC_Z_CONV) [10-n].\n    '''
        if value == "":
            if "CC_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CONVERGENCE"]=value.lower()


    def cc_dconvergence(self,value="show"):
        '''\nName: CC_DCONVERGENCE\nType: INTEGER\nDefault: 2\nOptions: 0:12:5:1\nDescription: Convergence criterion for the RMS residuals of excited state vectors\nRecommendation: : Use default. Should normally be set to the same value as CC_DTHRESHOLD.    '''
        if value == "":
            if "CC_DCONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DCONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DCONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DCONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DCONVERGENCE"]=value.lower()


    def cdft(self,value="show"):
        '''\nName: CDFT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Initiates a constrained DFT calculation\nRecommendation: : Set to TRUE if a Constrained DFT calculation is desired.    '''
        if value == "":
            if "CDFT" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT"]
                print "Keyword removed."
        elif value == "show":
            if "CDFT" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CDFT"]=value.lower()


    def cd_algorithm(self,value="show"):
        '''\nName: CD_ALGORITHM\nType: STRING\nDefault: 0\nOptions: Program determined:Direct:Semi-direct:Local-occupied\nDescription: Determines the algorithm for MP2 integral transformations.\nRecommendation: : Semi-direct is usually most efficient, and will normally be chosen by default.    '''
        if value == "":
            if "CD_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["CD_ALGORITHM"]
                print "Keyword removed."
        elif value == "show":
            if "CD_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["CD_ALGORITHM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CD_ALGORITHM"]=value.lower()


    def cdft_thresh(self,value="show"):
        '''\nName: CDFT_THRESH\nType: INTEGER\nDefault: 2\nOptions: 0:12:5:1\nDescription: Determines how tightly the constraint must be satisfied.
\nRecommendation: : Use default unless problems occur.    '''
        if value == "":
            if "CDFT_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "CDFT_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CDFT_THRESH"]=value.lower()


    def cdft_postdiis(self,value="show"):
        '''\nName: CDFT_POSTDIIS\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls whether the connstraint is enforced after DIIS extrapolation.
Recommentation: Use default unless convergence problems arise, in which case it may be benificial to turn this option off.  If selected, energies should be variational after the first iteration.\n    '''
        if value == "":
            if "CDFT_POSTDIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_POSTDIIS"]
                print "Keyword removed."
        elif value == "show":
            if "CDFT_POSTDIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_POSTDIIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CDFT_POSTDIIS"]=value.lower()


    def cdft_prediis(self,value="show"):
        '''\nName: CDFT_PREDIIS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls wheter the constraint is enforced before DIIS extrapolation.
\nRecommendation: : Use default unless problems arise, in which case it might be beneficial to turn this option on.    '''
        if value == "":
            if "CDFT_PREDIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CDFT_PREDIIS"]
                print "Keyword removed."
        elif value == "show":
            if "CDFT_PREDIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CDFT_PREDIIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CDFT_PREDIIS"]=value.lower()


    def cfmm_order(self,value="show"):
        '''\nName: CFMM_ORDER\nType: INTEGER\nDefault: 2\nOptions: 5  :30:15:1\nDescription: Controls the order of the multipole expansions in CFMM calculation.\n    '''
        if value == "":
            if "CFMM_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["CFMM_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "CFMM_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["CFMM_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CFMM_ORDER"]=value.lower()


    def chemsol_nn(self,value="show"):
        '''\nName: CHEMSOL_NN\nType: INTEGER\nDefault: 2\nOptions: 1:20:5:1\nDescription: Sets the number of grids used to calculate the average hydration free energy.\n    '''
        if value == "":
            if "CHEMSOL_NN" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_NN"]
                print "Keyword removed."
        elif value == "show":
            if "CHEMSOL_NN" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_NN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHEMSOL_NN"]=value.lower()


    def cis_convergence(self,value="show"):
        '''\nName: CIS_CONVERGENCE\nType: INTEGER\nDefault: 2\nOptions: 0   :12:6:1\nDescription: CIS is considered converged when error is less than 10-CIS_CONVERGENCE\n    '''
        if value == "":
            if "CIS_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_CONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_CONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_CONVERGENCE"]=value.lower()


    def cis_guess_disk(self,value="show"):
        '''\nName: CIS_GUESS_DISK\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Read the CIS guess from disk (previous calculation)\nRecommendation: : Requires a guess from previous calculation.    '''
        if value == "":
            if "CIS_GUESS_DISK" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_GUESS_DISK"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_GUESS_DISK" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_GUESS_DISK"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_GUESS_DISK"]=value.lower()


    def cis_n_roots(self,value="show"):
        '''\nName: CIS_N_ROOTS\nType: INTEGER\nDefault: 2\nOptions: 0   :200:0:1\nDescription: Sets the number of CI-Singles (CIS) excited state roots to find\n    '''
        if value == "":
            if "CIS_N_ROOTS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_N_ROOTS"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_N_ROOTS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_N_ROOTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_N_ROOTS"]=value.lower()


    def cis_relaxed_density(self,value="show"):
        '''\nName: CIS_RELAXED_DENSITY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Determines whether or not to use the relaxed CIS density for attachment/detachment density analysis\n    '''
        if value == "":
            if "CIS_RELAXED_DENSITY" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RELAXED_DENSITY"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RELAXED_DENSITY" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RELAXED_DENSITY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RELAXED_DENSITY"]=value.lower()


    def cis_singlets(self,value="show"):
        '''\nName: CIS_SINGLETS\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Solve for singlet excited states in RCIS calculations (ignored for UCIS)\n    '''
        if value == "":
            if "CIS_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_SINGLETS"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_SINGLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_SINGLETS"]=value.lower()


    def cis_triplets(self,value="show"):
        '''\nName: CIS_TRIPLETS\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Solve for triplet excited states in RCIS calculations (ignored for UCIS)\n    '''
        if value == "":
            if "CIS_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_TRIPLETS"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_TRIPLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_TRIPLETS"]=value.lower()


    def core_character(self,value="show"):
        '''\nName: CORE_CHARACTER\nType: INTEGER\nDefault: 2\nOptions: 0:4:0:1\nDescription: Selects how the core orbitals are determined in the frozen-core approximation.\nRecommendation: : Use default, unless performing calculations on molecules with heavy elements.    '''
        if value == "":
            if "CORE_CHARACTER" in self.dict_of_keywords:
                del self.dict_of_keywords["CORE_CHARACTER"]
                print "Keyword removed."
        elif value == "show":
            if "CORE_CHARACTER" in self.dict_of_keywords:
                return self.dict_of_keywords["CORE_CHARACTER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CORE_CHARACTER"]=value.lower()


    def cpscf_nseg(self,value="show"):
        '''\nName: CPSCF_NSEG\nType: INTEGER\nDefault: 2\nOptions: 0:20:0:1\nDescription: Controls the number of segments used to calculate the CPSCF equations.\nRecommendation: : Use default unless too much memory is requested.  Increasing this option reduces memory requirements.    '''
        if value == "":
            if "CPSCF_NSEG" in self.dict_of_keywords:
                del self.dict_of_keywords["CPSCF_NSEG"]
                print "Keyword removed."
        elif value == "show":
            if "CPSCF_NSEG" in self.dict_of_keywords:
                return self.dict_of_keywords["CPSCF_NSEG"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CPSCF_NSEG"]=value.lower()


    def deuterate(self,value="show"):
        '''\nName: DEUTERATE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Requests that all hydrogen atoms be replaces with deuterium.\nRecommendation: : Replacing hydrogen atoms reduces the fastest vibrational frequencies by a factor of 1.4, which allow for a larger fictitious mass and time step in ELMD calculations. There is no reason to replace hydrogens in BOMD calculations.    '''
        if value == "":
            if "DEUTERATE" in self.dict_of_keywords:
                del self.dict_of_keywords["DEUTERATE"]
                print "Keyword removed."
        elif value == "show":
            if "DEUTERATE" in self.dict_of_keywords:
                return self.dict_of_keywords["DEUTERATE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DEUTERATE"]=value.lower()


    def dma_midpoints(self,value="show"):
        '''\nName: DMA_MIDPOINTS\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Specifies whether to include bond midpoints in the DMA expansion.\n    '''
        if value == "":
            if "DMA_MIDPOINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["DMA_MIDPOINTS"]
                print "Keyword removed."
        elif value == "show":
            if "DMA_MIDPOINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["DMA_MIDPOINTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DMA_MIDPOINTS"]=value.lower()


    def direct_scf(self,value="show"):
        '''\nName: DIRECT_SCF\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls direct SCF.\nRecommendation: : Use default; direct SCF switches off in-core integrals.    '''
        if value == "":
            if "DIRECT_SCF" in self.dict_of_keywords:
                del self.dict_of_keywords["DIRECT_SCF"]
                print "Keyword removed."
        elif value == "show":
            if "DIRECT_SCF" in self.dict_of_keywords:
                return self.dict_of_keywords["DIRECT_SCF"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIRECT_SCF"]=value.lower()


    def dual_basis_energy(self,value="show"):
        '''\nName: DUAL_BASIS_ENERGY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Activates dual-basis SCF (HF or DFT) energy correction. \nRecommendation: : Use Dual-Basis to capture large-basis effects at smaller basis cost. Particularly useful with RI-MP2, in which HF often dominates. Use only proper subsets for small-basis calculation.    '''
        if value == "":
            if "DUAL_BASIS_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["DUAL_BASIS_ENERGY"]
                print "Keyword removed."
        elif value == "show":
            if "DUAL_BASIS_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["DUAL_BASIS_ENERGY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DUAL_BASIS_ENERGY"]=value.lower()


    def eda_bsse(self,value="show"):
        '''\nName: EDA_BSSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Calculates the BSSE correction when performing the energy decomposition analysis.\nRecommendation: : Set to TRUE unless a very large basis set is used.    '''
        if value == "":
            if "EDA_BSSE" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_BSSE"]
                print "Keyword removed."
        elif value == "show":
            if "EDA_BSSE" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_BSSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EDA_BSSE"]=value.lower()


    def eda_covp(self,value="show"):
        '''\nName: EDA_COVP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Perform COVP analysis when evaluating the RS or ARS charge-transfer correction. COVP analysis is currently implemented only for systems of two fragments.\nRecommendation: : Set to TRUE to perform COVP analysis of the CT term in an EDA or SCF MI(RS) job.    '''
        if value == "":
            if "EDA_COVP" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_COVP"]
                print "Keyword removed."
        elif value == "show":
            if "EDA_COVP" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_COVP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EDA_COVP"]=value.lower()


    def eda_print_covp(self,value="show"):
        '''\nName: EDA_PRINT_COVP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Replace the final MOs with the CVOP orbitals in the end of the run.\nRecommendation: : Set to TRUE to print COVP orbitals instead of conventional MOs.    '''
        if value == "":
            if "EDA_PRINT_COVP" in self.dict_of_keywords:
                del self.dict_of_keywords["EDA_PRINT_COVP"]
                print "Keyword removed."
        elif value == "show":
            if "EDA_PRINT_COVP" in self.dict_of_keywords:
                return self.dict_of_keywords["EDA_PRINT_COVP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EDA_PRINT_COVP"]=value.lower()


    def epao_iterate(self,value="show"):
        '''\nName: EPAO_ITERATE\nType: INTEGER\nDefault: 2\nOptions: 0:1000:0:10\nDescription: Controls iterations for EPAO calculations (see PAO_METHOD).\nRecommendation: : Use default. For molecules that are not too large, one can test the sensitivity of the results to the type of minimal functions by the use of optimized EPAOs in which case a value of 500 is reasonable.    '''
        if value == "":
            if "EPAO_ITERATE" in self.dict_of_keywords:
                del self.dict_of_keywords["EPAO_ITERATE"]
                print "Keyword removed."
        elif value == "show":
            if "EPAO_ITERATE" in self.dict_of_keywords:
                return self.dict_of_keywords["EPAO_ITERATE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EPAO_ITERATE"]=value.lower()


    def fast_xc(self,value="show"):
        '''\nName: FAST_XC\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls direct variable thresholds to accelerate the calculation of exchange and correlation (XC) in DFT.\nRecommendation: : This option improves the speed of a DFT calculation, but may occasionally cause the SCF calculation to diverge.    '''
        if value == "":
            if "FAST_XC" in self.dict_of_keywords:
                del self.dict_of_keywords["FAST_XC"]
                print "Keyword removed."
        elif value == "show":
            if "FAST_XC" in self.dict_of_keywords:
                return self.dict_of_keywords["FAST_XC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FAST_XC"]=value.lower()


    def frgm_lpcorr(self,value="show"):
        '''\nName: FRGM_LPCORR\nType: STRING\nDefault: 0\nOptions: None:ARS:RS:EXACT_SCF:ARS_EXACT_SCF:RS_EXACT_SCF\nDescription: Specifies a correction method performed after the locally-projected equations are converged.\nRecommendation: : For large basis sets use ARS, use RS if ARS fails.    '''
        if value == "":
            if "FRGM_LPCORR" in self.dict_of_keywords:
                del self.dict_of_keywords["FRGM_LPCORR"]
                print "Keyword removed."
        elif value == "show":
            if "FRGM_LPCORR" in self.dict_of_keywords:
                return self.dict_of_keywords["FRGM_LPCORR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FRGM_LPCORR"]=value.lower()


    def frgm_method(self,value="show"):
        '''\nName: FRGM_METHOD\nType: STRING\nDefault: 0\nOptions: None:STOLL:GIA:NOSCF_RS:NOSCF_ARS:NOSCF_DRS:NOSCF_RS_FOCK\nDescription: Specifies the locally-projected method.\nRecommendation: : STOLL and GIA - variational optimization of the ALMOs. NOSCF options are for computationally fast corrections of the FRAGMO initial guess.     '''
        if value == "":
            if "FRGM_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["FRGM_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "FRGM_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["FRGM_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FRGM_METHOD"]=value.lower()


    def ftc_class_thresh_mult(self,value="show"):
        '''\nName: FTC_CLASS_THRESH_MULT\nType: INTEGER\nDefault: 2\nOptions: 1 :9:5:1\nDescription: Together with FTC_CLASS_THRESH_ORDER, determines the cutoff threshold for included a shell-pair in the dd class, i.e. the class that is expanded in terms of plane waves. \n    '''
        if value == "":
            if "FTC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_CLASS_THRESH_MULT"]
                print "Keyword removed."
        elif value == "show":
            if "FTC_CLASS_THRESH_MULT" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_CLASS_THRESH_MULT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FTC_CLASS_THRESH_MULT"]=value.lower()


    def ftc_class_thresh_order(self,value="show"):
        '''\nName: FTC_CLASS_THRESH_ORDER\nType: INTEGER\nDefault: 2\nOptions: 1:9:5:1\nDescription: Together with FTC_CLASS_THRESH_MULT, determines the cutoff threshold for included a shell-pair in the dd class, i.e. the class that is expanded in terms of plane waves.\n    '''
        if value == "":
            if "FTC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "FTC_CLASS_THRESH_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FTC_CLASS_THRESH_ORDER"]=value.lower()


    def qui_geom_opt_fallback(self,value="show"):
        '''\nName: QUI_GEOM_OPT_FALLBACK\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Sets whether or not to fall back to cartesian coordinates if the optimization in internal or Z-matrix coordinates fails.\n    '''
        if value == "":
            if "QUI_GEOM_OPT_FALLBACK" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_GEOM_OPT_FALLBACK" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_GEOM_OPT_FALLBACK"]=value.lower()


    def geom_opt_hessian(self,value="show"):
        '''\nName: GEOM_OPT_HESSIAN\nType: STRING\nDefault: 0\nOptions: Diagonal:READ\nDescription: Determines the initial Hessian status.\nRecommendation: : An accurate initial Hessian will improve the performance of the optimizer, but is expensive to compute.    '''
        if value == "":
            if "GEOM_OPT_HESSIAN" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_HESSIAN"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_HESSIAN" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_HESSIAN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_HESSIAN"]=value.lower()


    def geom_opt_linear_angle(self,value="show"):
        '''\nName: GEOM_OPT_LINEAR_ANGLE\nType: INTEGER\nDefault: 2\nOptions: 150:180:165:1\nDescription: Threshold for near linear bond angles (in degrees).\nRecommendation: : Use default.    '''
        if value == "":
            if "GEOM_OPT_LINEAR_ANGLE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_LINEAR_ANGLE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_LINEAR_ANGLE"]=value.lower()


    def geom_opt_max_cycles(self,value="show"):
        '''\nName: GEOM_OPT_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 1:200:50:1\nDescription: Maximum number of optimization cycles.\nRecommendation: : The default should be sufficient for most cases. Increase if the initial guess geometry is poor, or for systems with shallow potential wells.    '''
        if value == "":
            if "GEOM_OPT_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_MAX_CYCLES"]=value.lower()


    def geom_opt_mode(self,value="show"):
        '''\nName: GEOM_OPT_MODE\nType: INTEGER\nDefault: 2\nOptions: 0:200:0:1\nDescription: Determines which Hessian mode is followed during a transition state search.\n    '''
        if value == "":
            if "GEOM_OPT_MODE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MODE"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_MODE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MODE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_MODE"]=value.lower()


    def geom_opt_print(self,value="show"):
        '''\nName: GEOM_OPT_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:7:3:1\nDescription: Controls the amount of optimization print output.\nRecommendation: : Use the default.    '''
        if value == "":
            if "GEOM_OPT_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_PRINT"]=value.lower()


    def geom_print(self,value="show"):
        '''\nName: GEOM_PRINT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the printing of additional geometric information at each step.\nRecommendation: : Use if you want to be able to quickly examine geometric parameters at the beginning and end of optimizations. Only prints in the beginning of single point energy calculations.    '''
        if value == "":
            if "GEOM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_PRINT"]=value.lower()


    def gvb_amp_scale(self,value="show"):
        '''\nName: GVB_AMP_SCALE\nType: INTEGER\nDefault: 2\nOptions: 0.001:1.000:1.000:0.001\nDescription: Scales the default orbital amplitude iteration step size for IP/RCC. PP amplitude equations are solved analytically, so this parameter does not affect PP.\nRecommendation: : Default is usually fine, but in some highly-correlated systems it can help with convergence to use smaller values.    '''
        if value == "":
            if "GVB_AMP_SCALE" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_AMP_SCALE"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_AMP_SCALE" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_AMP_SCALE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_AMP_SCALE"]=value.lower()


    def gvb_guess_mix(self,value="show"):
        '''\nName: GVB_GUESS_MIX\nType: INTEGER\nDefault: 2\nOptions: 0:100:0:1\nDescription: Similar to SCF_GUESS_MIX, it breaks alpha-beta symmetry for UPP by mixing the alpha HOMO and LUMO orbitals according to the user-defined fraction of LUMO to add the HOMO. 100 corresponds to a 1:1 ratio of HOMO and LUMO in the mixed orbitals.\nRecommendation: : 25 often works well to break symmetry without overly impeding convergence.    '''
        if value == "":
            if "GVB_GUESS_MIX" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_GUESS_MIX"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_GUESS_MIX" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_GUESS_MIX"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_GUESS_MIX"]=value.lower()


    def gvb_local(self,value="show"):
        '''\nName: GVB_LOCAL\nType: STRING\nDefault: 1\nOptions: Boys localized:Pipek-Mezey\nDescription: Sets the localization scheme used in the initial guess wave function.\nRecommendation: : Different initial guesses can sometimes lead to different solutions. It can be helpful to try both to ensure the global minimum has been found.    '''
        if value == "":
            if "GVB_LOCAL" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_LOCAL"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_LOCAL" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_LOCAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_LOCAL"]=value.lower()


    def gvb_n_pairs(self,value="show"):
        '''\nName: GVB_N_PAIRS\nType: INTEGER\nDefault: 2\nOptions: 0:100:0:1\nDescription: Alternative to CC_REST_OCC and CC_REST_VIR for setting active space size in GVB and valence coupled cluster methods.\nRecommendation: : Use default unless one wants to study a special active space. When using small active spaces, it is important to ensure that the proper orbitals are incorporated in the active space. If not, use the $reordermo feature to adjust the SCF orbitals appropriately.    '''
        if value == "":
            if "GVB_N_PAIRS" in self.dict_of_keywords:
                del self.dict_of_keywords["GVB_N_PAIRS"]
                print "Keyword removed."
        elif value == "show":
            if "GVB_N_PAIRS" in self.dict_of_keywords:
                return self.dict_of_keywords["GVB_N_PAIRS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GVB_N_PAIRS"]=value.lower()


    def incdft(self,value="show"):
        '''\nName: INCDFT\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Toggles the use of the IncDFT procedure for DFT energy calculations.\nRecommendation: : Turning this option on can lead to faster SCF calculations, particularly towards the end of the SCF. Please note that for some systems use of this option may lead to convergence problems.    '''
        if value == "":
            if "INCDFT" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT"]
                print "Keyword removed."
        elif value == "show":
            if "INCDFT" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCDFT"]=value.lower()


    def incdft_dendiff_thresh(self,value="show"):
        '''\nName: INCDFT_DENDIFF_THRESH\nType: INTEGER\nDefault: 2\nOptions: 0:12:8:1\nDescription: Sets the threshold for screening density matrix values in the IncDFT procedure. \nRecommendation: : If the default value causes convergence problems, set this value higher to tighten the threshold.    '''
        if value == "":
            if "INCDFT_DENDIFF_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_DENDIFF_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "INCDFT_DENDIFF_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_DENDIFF_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCDFT_DENDIFF_THRESH"]=value.lower()


    def incdft_dendiff_varthresh(self,value="show"):
        '''\nName: INCDFT_DENDIFF_VARTHRESH\nType: INTEGER\nDefault: 2\nOptions: 0:12:0:1\nDescription: Sets the lower bound for the variable threshold for screening density matrix values in the IncDFT procedure. The threshold will begin at this value and then vary depending on the error in the current SCF iteration until the value specified by INCDFT_DENDIFF_THRESH is reached. This means this value must be set lower than INCDFT_DENDIFF_THRESH.\nRecommendation: : If the default value causes convergence problems, set this value higher to tighten accuracy. If this fails, set to 0 and use a static threshold.    '''
        if value == "":
            if "INCDFT_DENDIFF_VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"]
                print "Keyword removed."
        elif value == "show":
            if "INCDFT_DENDIFF_VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCDFT_DENDIFF_VARTHRESH"]=value.lower()


    def incdft_griddiff_thresh(self,value="show"):
        '''\nName: INCDFT_GRIDDIFF_THRESH\nType: INTEGER\nDefault: 2\nOptions: 0:12:8:1\nDescription: Sets the threshold for screening functional values in the IncDFT procedure\nRecommendation: : If the default value causes convergence problems, set this value higher to tighten the threshold.    '''
        if value == "":
            if "INCDFT_GRIDDIFF_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "INCDFT_GRIDDIFF_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCDFT_GRIDDIFF_THRESH"]=value.lower()


    def incdft_griddiff_varthresh(self,value="show"):
        '''\nName: INCDFT_GRIDDIFF_VARTHRESH\nType: INTEGER\nDefault: 2\nOptions: 0   :12:0:1\nDescription: Sets the lower bound for the variable threshold for screening the functional values in the IncDFT procedure. The threshold will begin at this value and then vary depending on the error in the current SCF iteration until the value specified by INCDFT_GRIDDIFF_THRESH is reached. This means that this value must be set lower than INCDFT_GRIDDIFF_THRESH.\nRecommendation: : If the default value causes convergence problems, set this value higher to tighten accuracy. If this fails, set to 0 and use a static threshold.    '''
        if value == "":
            if "INCDFT_GRIDDIFF_VARTHRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"]
                print "Keyword removed."
        elif value == "show":
            if "INCDFT_GRIDDIFF_VARTHRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCDFT_GRIDDIFF_VARTHRESH"]=value.lower()


    def incfock(self,value="show"):
        '''\nName: INCFOCK\nType: INTEGER\nDefault: 2\nOptions: 0 :100:1:1\nDescription: Iteration number after which the incremental Fock matrix algorithm is initiated.  Setting this to 0 turns INCFOCK off.\nRecommendation: : May be necessary to allow several iterations before switching on INCFOCK.    '''
        if value == "":
            if "INCFOCK" in self.dict_of_keywords:
                del self.dict_of_keywords["INCFOCK"]
                print "Keyword removed."
        elif value == "show":
            if "INCFOCK" in self.dict_of_keywords:
                return self.dict_of_keywords["INCFOCK"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INCFOCK"]=value.lower()


    def integrals_buffer(self,value="show"):
        '''\nName: INTEGRALS_BUFFER\nType: INTEGER\nDefault: 2\nOptions: 1 :128:15:1\nDescription: Controls the size of in-core integral storage buffer (in megabytes).\nRecommendation: : Use the default, or consult your systems administrator for hardware limits.    '''
        if value == "":
            if "INTEGRALS_BUFFER" in self.dict_of_keywords:
                del self.dict_of_keywords["INTEGRALS_BUFFER"]
                print "Keyword removed."
        elif value == "show":
            if "INTEGRALS_BUFFER" in self.dict_of_keywords:
                return self.dict_of_keywords["INTEGRALS_BUFFER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTEGRALS_BUFFER"]=value.lower()


    def lin_k(self,value="show"):
        '''\nName: LIN_K\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls whether linear scaling evaluation of exact exchange (LinK) is used.\nRecommendation: : Use for HF and hybrid DFT calculations with large numbers of atoms.    '''
        if value == "":
            if "LIN_K" in self.dict_of_keywords:
                del self.dict_of_keywords["LIN_K"]
                print "Keyword removed."
        elif value == "show":
            if "LIN_K" in self.dict_of_keywords:
                return self.dict_of_keywords["LIN_K"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["LIN_K"]=value.lower()


    def max_sub_file_num(self,value="show"):
        '''\nName: MAX_SUB_FILE_NUM\nType: INTEGER\nDefault: 2\nOptions: 1:64:16:1\nDescription: Sets the maximum number of sub files allowed.\nRecommendation: : Leave as default, or adjust according to your system limits.    '''
        if value == "":
            if "MAX_SUB_FILE_NUM" in self.dict_of_keywords:
                del self.dict_of_keywords["MAX_SUB_FILE_NUM"]
                print "Keyword removed."
        elif value == "show":
            if "MAX_SUB_FILE_NUM" in self.dict_of_keywords:
                return self.dict_of_keywords["MAX_SUB_FILE_NUM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MAX_SUB_FILE_NUM"]=value.lower()


    def qui_frozen_core(self,value="show"):
        '''\nName: QUI_FROZEN_CORE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_FROZEN_CORE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_FROZEN_CORE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_FROZEN_CORE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_FROZEN_CORE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_FROZEN_CORE"]=value.lower()


    def xc_smart_grid(self,value="show"):
        '''\nName: XC_SMART_GRID\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Uses SG-0 (where available) for early SCF cycles, and switches to the (larger) grid specified by XC_GRID (which defaults to SG-1, if not otherwise specified) for final cycles of the SCF.\nRecommendation: : The use of the smart grid can save some time on initial SCF cycles.    '''
        if value == "":
            if "XC_SMART_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["XC_SMART_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "XC_SMART_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["XC_SMART_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XC_SMART_GRID"]=value.lower()


    def xopt_seam_only(self,value="show"):
        '''\nName: XOPT_SEAM_ONLY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Orders an intersection seam search only, no minimization is to perform.\nRecommendation: : In systems with a large number of degrees of freedom it might be useful to locate the seam first setting this option to TRUE and use that geometry as a starting point for the minimization.    '''
        if value == "":
            if "XOPT_SEAM_ONLY" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_SEAM_ONLY"]
                print "Keyword removed."
        elif value == "show":
            if "XOPT_SEAM_ONLY" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_SEAM_ONLY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XOPT_SEAM_ONLY"]=value.lower()


    def xcis(self,value="show"):
        '''\nName: XCIS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Do an XCIS calculation in addition to a CIS calculation\n    '''
        if value == "":
            if "XCIS" in self.dict_of_keywords:
                del self.dict_of_keywords["XCIS"]
                print "Keyword removed."
        elif value == "show":
            if "XCIS" in self.dict_of_keywords:
                return self.dict_of_keywords["XCIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XCIS"]=value.lower()


    def wavefunction_analysis(self,value="show"):
        '''\nName: WAVEFUNCTION_ANALYSIS\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls the running of the default wavefunction analysis tasks.\n    '''
        if value == "":
            if "WAVEFUNCTION_ANALYSIS" in self.dict_of_keywords:
                del self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"]
                print "Keyword removed."
        elif value == "show":
            if "WAVEFUNCTION_ANALYSIS" in self.dict_of_keywords:
                return self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["WAVEFUNCTION_ANALYSIS"]=value.lower()


    def vibman_print(self,value="show"):
        '''\nName: VIBMAN_PRINT\nType: INTEGER\nDefault: 2\nOptions: 1:7:1:1\nDescription: Controls level of extra print out for vibrational analysis.\nRecommendation: : Use default.    '''
        if value == "":
            if "VIBMAN_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["VIBMAN_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "VIBMAN_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["VIBMAN_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["VIBMAN_PRINT"]=value.lower()


    def vci(self,value="show"):
        '''\nName: VCI\nType: INTEGER\nDefault: 2\nOptions: 0:10:0:1\nDescription: Specifies the number of quanta involved in the VCI calculation.\nRecommendation: : The availability depends on the memory of the machine.  For example, a machine with 1.5 GB memory and for molecules with fewer than 4 atoms, VCI(10) can be carried out, for molecule containing fewer than 5 atoms, VCI(6) can be carried out, for molecule containing fewer than 6 atoms, VCI(5) can be carried out. For molecules containing fewer than 50 atoms, VCI(2) is available. VCI(1) and VCI(3) usually overestimated the true energy while VCI(4) usually gives an answer close to the converged energy.    '''
        if value == "":
            if "VCI" in self.dict_of_keywords:
                del self.dict_of_keywords["VCI"]
                print "Keyword removed."
        elif value == "show":
            if "VCI" in self.dict_of_keywords:
                return self.dict_of_keywords["VCI"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["VCI"]=value.lower()


    def stability_analysis(self,value="show"):
        '''\nName: STABILITY_ANALYSIS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Performs stability analysis for a HF or DFT solution.\nRecommendation: : Set to TRUE when a HF or DFT solution is suspected to be unstable.    '''
        if value == "":
            if "STABILITY_ANALYSIS" in self.dict_of_keywords:
                del self.dict_of_keywords["STABILITY_ANALYSIS"]
                print "Keyword removed."
        elif value == "show":
            if "STABILITY_ANALYSIS" in self.dict_of_keywords:
                return self.dict_of_keywords["STABILITY_ANALYSIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["STABILITY_ANALYSIS"]=value.lower()


    def scf_print_frgm(self,value="show"):
        '''\nName: SCF_PRINT_FRGM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the output of Q-Chem jobs on isolated fragments.\nRecommendation: : Use TRUE if details about isolated fragments are important.    '''
        if value == "":
            if "SCF_PRINT_FRGM" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_PRINT_FRGM"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_PRINT_FRGM" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_PRINT_FRGM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_PRINT_FRGM"]=value.lower()


    def omega(self,value="show"):
        '''\nName: OMEGA\nType: INTEGER\nDefault: 2\nOptions: 0.001:5.000:0.200:0.001\nDescription: Controls the degree of attenuation of the Coulomb operator.\n    '''
        if value == "":
            if "OMEGA" in self.dict_of_keywords:
                del self.dict_of_keywords["OMEGA"]
                print "Keyword removed."
        elif value == "show":
            if "OMEGA" in self.dict_of_keywords:
                return self.dict_of_keywords["OMEGA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["OMEGA"]=value.lower()


    def chemsol_print(self,value="show"):
        '''\nName: CHEMSOL_PRINT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Increases the amount of ChemSol output.\n    '''
        if value == "":
            if "CHEMSOL_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "CHEMSOL_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHEMSOL_PRINT"]=value.lower()


    def qui_charge(self,value="show"):
        '''\nName: QUI_CHARGE\nType: INTEGER\nDefault: 2\nOptions: -100:100:0:1\nDescription: Sets the total charge of the system.\n    '''
        if value == "":
            if "QUI_CHARGE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_CHARGE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_CHARGE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_CHARGE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_CHARGE"]=value.lower()


    def thresh(self,value="show"):
        '''\nName: THRESH\nType: INTEGER\nDefault: 2\nOptions: 0:14:8:1\nDescription: Cutoff for neglect of two electron integrals. 10-THRESH (THRESH ? 14).\nRecommendation: : Should be at least three greater than the SCF convergence setting. Increase for more significant figures, at greater computational cost.    '''
        if value == "":
            if "THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["THRESH"]=value.lower()


    def unrestricted(self,value="show"):
        '''\nName: UNRESTRICTED\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the use of restricted or unrestricted orbitals.\nRecommendation: : Use default unless ROHF is desired. Note that for unrestricted calculations on systems with an even number of electrons it is usually necessary to break alpha-beta symmetry in the initial guess, by using SCF_GUESS_MIX or providing $occupied information.    '''
        if value == "":
            if "UNRESTRICTED" in self.dict_of_keywords:
                del self.dict_of_keywords["UNRESTRICTED"]
                print "Keyword removed."
        elif value == "show":
            if "UNRESTRICTED" in self.dict_of_keywords:
                return self.dict_of_keywords["UNRESTRICTED"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["UNRESTRICTED"]=value.lower()


    def write_wfn(self,value="show"):
        '''\nName: WRITE_WFN\nType: STRING\nDefault: 0\nOptions: qchem\nDescription: Specifies whether or not a wfn file is created, which is suitable for use with AIMPAC. Note that the output to this file is currently limited to f orbitals, which is the highest angular momentum implemented in AIMPAC.\n    '''
        if value == "":
            if "WRITE_WFN" in self.dict_of_keywords:
                del self.dict_of_keywords["WRITE_WFN"]
                print "Keyword removed."
        elif value == "show":
            if "WRITE_WFN" in self.dict_of_keywords:
                return self.dict_of_keywords["WRITE_WFN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["WRITE_WFN"]=value.lower()


    def aimd_moments(self,value="show"):
        '''\nName: AIMD_MOMENTS\nType: INTEGER\nDefault: 2\nOptions: 0:20:0:1\nDescription: Specifies the order of multipole moments that are output at each time step.  Setting this to 0 disables the printing of moments.\n    '''
        if value == "":
            if "AIMD_MOMENTS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_MOMENTS"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_MOMENTS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_MOMENTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_MOMENTS"]=value.lower()


    def n_frozen_core(self,value="show"):
        '''\nName: N_FROZEN_CORE\nType: INTEGER\nDefault: 2\nOptions: 0:100:0:1\nDescription: Sets the number of frozen core orbitals in a post-Hartree-Fock calculation.\nRecommendation: : While the default is not to freeze orbitals, MP2 calculations are more efficient with frozen core orbitals. Use FC if possible.    '''
        if value == "":
            if "N_FROZEN_CORE" in self.dict_of_keywords:
                del self.dict_of_keywords["N_FROZEN_CORE"]
                print "Keyword removed."
        elif value == "show":
            if "N_FROZEN_CORE" in self.dict_of_keywords:
                return self.dict_of_keywords["N_FROZEN_CORE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["N_FROZEN_CORE"]=value.lower()


    def n_frozen_virtual(self,value="show"):
        '''\nName: N_FROZEN_VIRTUAL\nType: INTEGER\nDefault: 2\nOptions: 0:500:0:1\nDescription: Sets the number of frozen virtual orbitals in a post-Hartree-Fock calculation.\n    '''
        if value == "":
            if "N_FROZEN_VIRTUAL" in self.dict_of_keywords:
                del self.dict_of_keywords["N_FROZEN_VIRTUAL"]
                print "Keyword removed."
        elif value == "show":
            if "N_FROZEN_VIRTUAL" in self.dict_of_keywords:
                return self.dict_of_keywords["N_FROZEN_VIRTUAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["N_FROZEN_VIRTUAL"]=value.lower()


    def chemsol(self,value="show"):
        '''\nName: CHEMSOL\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the use of ChemSol in Q-Chem.\n    '''
        if value == "":
            if "CHEMSOL" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL"]
                print "Keyword removed."
        elif value == "show":
            if "CHEMSOL" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHEMSOL"]=value.lower()


    def ftc(self,value="show"):
        '''\nName: FTC\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the overall use of the FTC.\nRecommendation: : Use FTC when bigger and/or diffuse basis sets are used.     '''
        if value == "":
            if "FTC" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC"]
                print "Keyword removed."
        elif value == "show":
            if "FTC" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FTC"]=value.lower()


    def qui_cfmm(self,value="show"):
        '''\nName: QUI_CFMM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_CFMM" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_CFMM"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_CFMM" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_CFMM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_CFMM"]=value.lower()


    def qui_largemol_none(self,value="show"):
        '''\nName: QUI_LARGEMOL_NONE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_LARGEMOL_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_LARGEMOL_NONE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_LARGEMOL_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_LARGEMOL_NONE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_LARGEMOL_NONE"]=value.lower()


    def qui_solvent_onsager(self,value="show"):
        '''\nName: QUI_SOLVENT_ONSAGER\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_SOLVENT_ONSAGER" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_ONSAGER"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_ONSAGER" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_ONSAGER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_ONSAGER"]=value.lower()


    def qui_plots_points(self,value="show"):
        '''\nName: QUI_PLOTS_POINTS\nType: INTEGER\nDefault: 2\nOptions: 1:10000:1:1\nDescription: \n    '''
        if value == "":
            if "QUI_PLOTS_POINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_PLOTS_POINTS"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_PLOTS_POINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_PLOTS_POINTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_PLOTS_POINTS"]=value.lower()


    def ssg(self,value="show"):
        '''\nName: SSG\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the calculation of the SSG wavefunction.\nRecommendation: : See also the UNRESTRICTED and DIIS_SUBSPACE_SIZE $rem variables.    '''
        if value == "":
            if "SSG" in self.dict_of_keywords:
                del self.dict_of_keywords["SSG"]
                print "Keyword removed."
        elif value == "show":
            if "SSG" in self.dict_of_keywords:
                return self.dict_of_keywords["SSG"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SSG"]=value.lower()


    def qui_eom_method(self,value="show"):
        '''\nName: QUI_EOM_METHOD\nType: STRING\nDefault: 0\nOptions: Spin Flip:EA (Diffuse Orbital):IP (Diffuse Orbital):IP (Proper):DIP\nDescription: Specifies the type of EOM calculation to perform\n    '''
        if value == "":
            if "QUI_EOM_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_EOM_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_EOM_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_EOM_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_EOM_METHOD"]=value.lower()


    def cc_theta_stepsize(self,value="show"):
        '''\nName: CC_THETA_STEPSIZE\nType: INTEGER\nDefault: 2\nOptions: 0:::\nDescription: Scale factor for the orbital rotation step size. The optimal rotation steps should be approximately equal to the gradient vector.\nRecommendation: : Try a smaller value in cases of poor convergence and very large orbital gradients. For example, a value of 01001 translates to 0.1    '''
        if value == "":
            if "CC_THETA_STEPSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_THETA_STEPSIZE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_THETA_STEPSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_THETA_STEPSIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_THETA_STEPSIZE"]=value.lower()


    def geom_opt_tol_displacement(self,value="show"):
        '''\nName: GEOM_OPT_TOL_DISPLACEMENT\nType: INTEGER\nDefault: 2\nOptions: 0:5000:1200:1\nDescription: Convergence on maximum atomic displacement (in micro-angstroms).\nRecommendation: : Use the default. To converge the gradient and either one of the energy or displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_DISPLACEMENT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_TOL_DISPLACEMENT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_DISPLACEMENT"]=value.lower()


    def geom_opt_tol_energy(self,value="show"):
        '''\nName: GEOM_OPT_TOL_ENERGY\nType: INTEGER\nDefault: 2\nOptions: 1:500:100:1\nDescription: Convergence on energy change of successive optimization cycles (x 10-8).\nRecommendation: : Use the default. To converge the gradient and either one of the energy or displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_TOL_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_ENERGY"]=value.lower()


    def geom_opt_tol_gradient(self,value="show"):
        '''\nName: GEOM_OPT_TOL_GRADIENT\nType: INTEGER\nDefault: 2\nOptions: 1:1000:300:1\nDescription: Convergence on maximum gradient component.\nRecommendation: : Use the default. To converge the gradient and either one of the energy and displacement tolerances must be satisfied.    '''
        if value == "":
            if "GEOM_OPT_TOL_GRADIENT" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_TOL_GRADIENT" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_TOL_GRADIENT"]=value.lower()


    def skip_cis_rpa(self,value="show"):
        '''\nName: SKIP_CIS_RPA\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Skips the solution of the CIS, RPA, TDA or TDDFT equations for wavefunction analysis.
\nRecommendation: : Set to true to speed up the generation of plot data if the same calculation has ben run previously.    '''
        if value == "":
            if "SKIP_CIS_RPA" in self.dict_of_keywords:
                del self.dict_of_keywords["SKIP_CIS_RPA"]
                print "Keyword removed."
        elif value == "show":
            if "SKIP_CIS_RPA" in self.dict_of_keywords:
                return self.dict_of_keywords["SKIP_CIS_RPA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SKIP_CIS_RPA"]=value.lower()


    def qui_coordinates(self,value="show"):
        '''\nName: QUI_COORDINATES\nType: STRING\nDefault: 0\nOptions: Cartesian:Z-matrix:Z-matrix (compact)\nDescription: Controls the format of the geometry in the output file.\n    '''
        if value == "":
            if "QUI_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_COORDINATES"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_COORDINATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_COORDINATES"]=value.lower()


    def pseudo_canonical(self,value="show"):
        '''\nName: PSEUDO_CANONICAL\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: When SCF_ALGORITHM = DM, this controls the way the initial step, and steps after subspace resets are taken.\nRecommendation: : The default is usually more efficient, but choosing TRUE sometimes avoids problems with orbital reordering.    '''
        if value == "":
            if "PSEUDO_CANONICAL" in self.dict_of_keywords:
                del self.dict_of_keywords["PSEUDO_CANONICAL"]
                print "Keyword removed."
        elif value == "show":
            if "PSEUDO_CANONICAL" in self.dict_of_keywords:
                return self.dict_of_keywords["PSEUDO_CANONICAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PSEUDO_CANONICAL"]=value.lower()


    def purecart(self,value="show"):
        '''\nName: PURECART\nType: STRING\nDefault: 0\nOptions: All pure//1111:All cartesian//2222:Pure d, f and g//2111:Pure d and f//2211:Pure d only//2221\nDescription: Controls the use of pure (spherical harmonic) or Cartesian angular forms.\n    '''
        if value == "":
            if "PURECART" in self.dict_of_keywords:
                del self.dict_of_keywords["PURECART"]
                print "Keyword removed."
        elif value == "show":
            if "PURECART" in self.dict_of_keywords:
                return self.dict_of_keywords["PURECART"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PURECART"]=value.lower()


    def moprop(self,value="show"):
        '''\nName: MOPROP\nType: STRING\nDefault: 0\nOptions: None//0:NMR Shielding Tensors//1:Static Polarizability//2:Dynamic Polarizability//100:Hyperpolarizability//101:Hyperpolarizability (Read)//102:Hyperpolarizabilty (Wigner)//103:Hyperpolarizability (Wigner+Read)//104\nDescription: Specifies the job for mopropman.  Note that for hyperpolarizabilities, the Wigner option uses the (2n+1) rule, and the Read option takes first order results from disk.\n    '''
        if value == "":
            if "MOPROP" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP"]=value.lower()


    def ecp(self,value="show"):
        '''\nName: ECP\nType: STRING\nDefault: 0\nOptions: None:----------:CRENBL:CRENBS:HWMB:HWVDZ:LACVP:LANL2DZ:SBKJC:SRLC:SRSC:---------:User-defined//gen\nDescription: Defines the effective core potential and associated basis set to be used\nRecommendation: : Pseudopotentials are recommended for first row transition metals and heavier elements. Consult the reviews for more details.    '''
        if value == "":
            if "ECP" in self.dict_of_keywords:
                del self.dict_of_keywords["ECP"]
                print "Keyword removed."
        elif value == "show":
            if "ECP" in self.dict_of_keywords:
                return self.dict_of_keywords["ECP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["ECP"]=value.lower()


    def aimd_fictitious_mass(self,value="show"):
        '''\nName: AIMD_FICTITIOUS_MASS\nType: INTEGER\nDefault: 2\nOptions: 1:500:100:1\nDescription: Specifies the value of the fictious electronic mass is atomic units.  This quantity has dimensions (energy)x(time)2. 
\nRecommendation: : Values in the range 50-200 a.u. have been employed.    '''
        if value == "":
            if "AIMD_FICTITIOUS_MASS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FICTITIOUS_MASS"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_FICTITIOUS_MASS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FICTITIOUS_MASS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_FICTITIOUS_MASS"]=value.lower()


    def anhar(self,value="show"):
        '''\nName: ANHAR\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects whether or not to perform various nuclear vibrational theory (TOSH, VPT2, VCI) calculations to obtain vibrational anharmonic frequencies. 
\nRecommendation: : Since this calculation involves the third and fourth derivatives at the minimum of the potential energy surface, it is recommended that the geometry optimization tolerances (displacement, gradient and energy) be set tighter.  Note that VPT2 calculations may fail if the system involves accidental degenerate resonances. See the VCI $rem variable for more details about increasing the accuracy of anharmonic calculations.    '''
        if value == "":
            if "ANHAR" in self.dict_of_keywords:
                del self.dict_of_keywords["ANHAR"]
                print "Keyword removed."
        elif value == "show":
            if "ANHAR" in self.dict_of_keywords:
                return self.dict_of_keywords["ANHAR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["ANHAR"]=value.lower()


    def aimd_initial_velocities(self,value="show"):
        '''\nName: AIMD_INITIAL_VELOCITIES\nType: STRING\nDefault: 0\nOptions: Read//0:ZPE:Random:Thermal\nDescription: Specifies the method for selecting initial nuclear velocities. \nRecommendation: : This variable need only be specified in the event that velocities are not specified explicitly in a $velocity  section.    '''
        if value == "":
            if "AIMD_INITIAL_VELOCITIES" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_INITIAL_VELOCITIES" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_INITIAL_VELOCITIES"]=value.lower()


    def aimd_temperature(self,value="show"):
        '''\nName: AIMD_TEMPERATURE\nType: INTEGER\nDefault: 2\nOptions: 0:1000:0:1\nDescription: Specifies a temperature (in Kelvin) for Maxwell-Boltzmann velocity sampling.\nRecommendation: : This variable is only useful in conjunction with AIMD velocit initialisation set to Thermal. Note that the simulations are run at constant energy, rather than constant temperature, so the mean nuclear kinetic energy will fluctuate in the course of the simulation.    '''
        if value == "":
            if "AIMD_TEMPERATURE" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_TEMPERATURE"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_TEMPERATURE" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_TEMPERATURE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_TEMPERATURE"]=value.lower()


    def anharmonic(self,value="show"):
        '''\nName: ANHARMONIC\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects whether or not to perform various nuclear vibrational theory (TOSH, VPT2, VCI) calculations to obtain vibrational anharmonic frequencies. 
\nRecommendation: : Since this calculation involves the third and fourth derivatives at the minimum of the potential energy surface, it is recommended that the geometry optimization tolerances (displacement, gradient and energy) be set tighter.  Note that VPT2 calculations may fail if the system involves accidental degenerate resonances. See the VCI $rem variable for more details about increasing the accuracy of anharmonic calculations.    '''
        if value == "":
            if "ANHARMONIC" in self.dict_of_keywords:
                del self.dict_of_keywords["ANHARMONIC"]
                print "Keyword removed."
        elif value == "show":
            if "ANHARMONIC" in self.dict_of_keywords:
                return self.dict_of_keywords["ANHARMONIC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["ANHARMONIC"]=value.lower()


    def basis_linear_dependence_thresh(self,value="show"):
        '''\nName: BASIS_LINEAR_DEPENDENCE_THRESH\nType: INTEGER\nDefault: 2\nOptions: 1   :10:6:1\nDescription: Sets the threshold for determining linear dependence in the basis set.  The threshold is set to 10-n.\nRecommendation: : Set to 5 or smaller if you have a poorly behaved SCF and you suspect linear dependence in you basis set. Lower values (larger thresholds) may affect the accuracy of the calculation.    '''
        if value == "":
            if "BASIS_LINEAR_DEPENDENCE_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS_LINEAR_DEPENDENCE_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "BASIS_LINEAR_DEPENDENCE_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS_LINEAR_DEPENDENCE_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["BASIS_LINEAR_DEPENDENCE_THRESH"]=value.lower()


    def basis_projection_type(self,value="show"):
        '''\nName: BASIS_PROJECTION_TYPE\nType: STRING\nDefault: 0\nOptions: Project Fock matrix//FOPPROJECTION:Project MOs//OVPROJECTION\nDescription: Determines which method to use when projecting the density matrix for the basis set projection guess.\n    '''
        if value == "":
            if "BASIS_PROJECTION_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS_PROJECTION_TYPE"]
                print "Keyword removed."
        elif value == "show":
            if "BASIS_PROJECTION_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS_PROJECTION_TYPE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["BASIS_PROJECTION_TYPE"]=value.lower()


    def cc_amplitude_response(self,value="show"):
        '''\nName: CC_AMPLITUDE_RESPONSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If set to TRUE, adds amplitude response terms to one-particle and two-particle CCSD density matrices before calculation of properties. CC_PROP must be set to TRUE. \nRecommendation: : The cost is always about the cost of an analytic gradient calculation, independent of whether or not the two-particle properties are requested. Besides, adding amplitude response terms without orbital response will unlikely improve the quality of the properties. However, it can be used for debugging purposes.    '''
        if value == "":
            if "CC_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_AMPLITUDE_RESPONSE"]=value.lower()


    def cc_properties(self,value="show"):
        '''\nName: CC_PROPERTIES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments < X2>, < Y2>, and < Z2> of electron density, and the total < R2> = < X2> +< Y2> +< Z2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.\nRecommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PROPERTIES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PROPERTIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PROPERTIES"]=value.lower()


    def auxiliary_basis(self,value="show"):
        '''\nName: AUXILIARY_BASIS\nType: STRING\nDefault: 0\nOptions: None:RIMP2-VDZ:RIMP2-TZVPP:RIMP2-cc-PVDZ:RIMP2-cc-PVTZ:RIMP2-cc-PVQZ:RIMP2-aug-cc-PVDZ:RIMP2-aug-cc-PVTZ:RIMP2-aug-cc-PVQZ\nDescription: Specifies the auxiliary basis to be used in a RI-MP2 calculation\n    '''
        if value == "":
            if "AUXILIARY_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["AUXILIARY_BASIS"]
                print "Keyword removed."
        elif value == "show":
            if "AUXILIARY_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["AUXILIARY_BASIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AUXILIARY_BASIS"]=value.lower()


    def cc_two_particle_properties(self,value="show"):
        '''\nName: CC_TWO_PARTICLE_PROPERTIES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Request for calculation of non-relaxed two-particle CCSD properties. The two-particle properties currently include. The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation is inferior compared to the cost of 2>. The variable CC_PROPERTIES must be also set to TRUE.\nRecommendation: : The two-particle properties are extremely computationally expensive, since they require calculation and use of the two-particle density matrix (the cost is approximately the same as the cost of an analytic gradient calculation). Do not request the two-particle properties unless you really need them.    '''
        if value == "":
            if "CC_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_TWO_PARTICLE_PROPERTIES"]=value.lower()


    def cc_block_tensor_buffer_size(self,value="show"):
        '''\nName: CC_BLOCK_TENSOR_BUFFER_SIZE\nType: INTEGER\nDefault: 2\nOptions: 0:8000:1000:1\nDescription: Specifies the maximum size, in Mb, of the buffers for in-core storage of block-tensors.\nRecommendation: : Larger values can give better I/O performance and are recommended for systems with large memory (add to your .qchemrc file)    '''
        if value == "":
            if "CC_BLOCK_TENSOR_BUFFER_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_BLOCK_TENSOR_BUFFER_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_BLOCK_TENSOR_BUFFER_SIZE"]=value.lower()


    def cc_diis_maximum_overlap(self,value="show"):
        '''\nName: CC_DIIS_MAXIMUM_OVERLAP\nType: INTEGER\nDefault: 2\nOptions: 0.01:1.00:1.00:0.01\nDescription: \n    '''
        if value == "":
            if "CC_DIIS_MAXIMUM_OVERLAP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS_MAXIMUM_OVERLAP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS_MAXIMUM_OVERLAP"]=value.lower()


    def cc_diis12_switch(self,value="show"):
        '''\nName: CC_DIIS12_SWITCH\nType: INTEGER\nDefault: 2\nOptions: 0:12:5:1\nDescription: When to switch from DIIS 2 to DIIS 1 procedure, or when DIIS 2 procedure is required to generate DIIS guesses less frequently. Total value of DIIS error vector must be less than 10-n, where n is the value of this option.\n    '''
        if value == "":
            if "CC_DIIS12_SWITCH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS12_SWITCH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS12_SWITCH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS12_SWITCH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS12_SWITCH"]=value.lower()


    def cc_diis_extrapolation_frequency(self,value="show"):
        '''\nName: CC_DIIS_EXTRAPOLATION_FREQUENCY\nType: INTEGER\nDefault: 2\nOptions: 0:100:2:1\nDescription: DIIS extrapolation will be attempted every n iterations. However, DIIS2 will be attempted every iteration while total error vector exceeds CC_DIIS12_SWITCH. DIIS1 cannot generate guesses more frequently than every 2 iterations. \n    '''
        if value == "":
            if "CC_DIIS_EXTRAPOLATION_FREQUENCY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_EXTRAPOLATION_FREQUENCY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS_EXTRAPOLATION_FREQUENCY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_EXTRAPOLATION_FREQUENCY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS_EXTRAPOLATION_FREQUENCY"]=value.lower()


    def cc_diis_minimum_overlap(self,value="show"):
        '''\nName: CC_DIIS_MINIMUM_OVERLAP\nType: INTEGER\nDefault: 2\nOptions: 0:12:11:1\nDescription: The DIIS procedure will be halted when the square root of smallest element of the error overlap matrix is less than 10-n, where n is the value of this option. Small values of the B matrix mean it will become near-singular, making the DIIS equations difficult to solve.\n    '''
        if value == "":
            if "CC_DIIS_MINIMUM_OVERLAP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS_MINIMUM_OVERLAP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS_MINIMUM_OVERLAP"]=value.lower()


    def cc_diis_size(self,value="show"):
        '''\nName: CC_DIIS_SIZE\nType: INTEGER\nDefault: 2\nOptions: 1:50:7:1\nDescription: Specifies the maximum size of the DIIS space.\nRecommendation: : Larger values involve larger amounts of disk storage.    '''
        if value == "":
            if "CC_DIIS_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_SIZE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_SIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS_SIZE"]=value.lower()


    def cc_diis_start(self,value="show"):
        '''\nName: CC_DIIS_START\nType: INTEGER\nDefault: 2\nOptions: 1:500:3:1\nDescription: Iteration number when DIIS is turned on. Set to a large number to disable DIIS.\nRecommendation: : Occasionally DIIS can cause optimized orbital coupled-cluster calculations to diverge through large orbital changes. If this is seen, DIIS should be disabled.    '''
        if value == "":
            if "CC_DIIS_START" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS_START"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS_START" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS_START"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS_START"]=value.lower()


    def cc_dmaxiter(self,value="show"):
        '''\nName: CC_DMAXITER\nType: INTEGER\nDefault: 2\nOptions: 1:100:30:1\nDescription: Maximum number of iteration allowed for Davidson diagonalization procedure. \nRecommendation: : Default is usually sufficient    '''
        if value == "":
            if "CC_DMAXITER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DMAXITER"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DMAXITER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DMAXITER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DMAXITER"]=value.lower()


    def cc_dmaxvectors(self,value="show"):
        '''\nName: CC_DMAXVECTORS\nType: INTEGER\nDefault: 2\nOptions: 1:500:60:1\nDescription: Specifies maximum number of vectors in the subspace for the Davidson diagonalization. \nRecommendation: : Larger values increase disk storage but accelerate and stabilize convergence.    '''
        if value == "":
            if "CC_DMAXVECTORS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DMAXVECTORS"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DMAXVECTORS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DMAXVECTORS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DMAXVECTORS"]=value.lower()


    def cc_do_disconnected(self,value="show"):
        '''\nName: CC_DO_DISCONNECTED\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Determines whether disconnected terms included in the EOM-OD equations\nRecommendation: : Inclusion of disconnected terms has very small effects and is not necessary.    '''
        if value == "":
            if "CC_DO_DISCONNECTED" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DISCONNECTED"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DO_DISCONNECTED" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DISCONNECTED"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DO_DISCONNECTED"]=value.lower()


    def cc_do_cisdt(self,value="show"):
        '''\nName: CC_DO_CISDT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the calculation of full CISDT\n    '''
        if value == "":
            if "CC_DO_CISDT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_CISDT"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DO_CISDT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_CISDT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DO_CISDT"]=value.lower()


    def cc_do_dyson_ee(self,value="show"):
        '''\nName: CC_DO_DYSON_EE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether excited state Dyson orbitals will be calculated for EOM-IP/EA-CCSD calculations.\nRecommendation: : none    '''
        if value == "":
            if "CC_DO_DYSON_EE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DYSON_EE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DO_DYSON_EE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DYSON_EE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DO_DYSON_EE"]=value.lower()


    def cc_do_dyson(self,value="show"):
        '''\nName: CC_DO_DYSON\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether ground state Dyson orbitals will be calculated for EOM-IP/EA-CCSD calculations.\nRecommendation: : none    '''
        if value == "":
            if "CC_DO_DYSON" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DO_DYSON"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DO_DYSON" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DO_DYSON"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DO_DYSON"]=value.lower()


    def cc_ea(self,value="show"):
        '''\nName: CC_EA\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If TRUE, calculates EOM-EA-CCSD excitation energies and properties using the diffuse orbital trick. A very diffuse orbital is added to the basis set, excitations from which correspond to electron attachment.\n    '''
        if value == "":
            if "CC_EA" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EA"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EA" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EA"]=value.lower()


    def cc_ip(self,value="show"):
        '''\nName: CC_IP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If TRUE, calculates EOM-IP-CCSD excitation energies and properties using the diffuse orbital trick. A very diffuse orbital is added to the basis set, excitations to which correspond to ionization. CC_IP and CC_IP_PROPER keywords are mutually exclusive. \n    '''
        if value == "":
            if "CC_IP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_IP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_IP"]=value.lower()


    def cc_ip_filter(self,value="show"):
        '''\nName: CC_IP_FILTER\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If TRUE, filters the EOM-IP-CCSD amplitudes obtained using the diffuse orbital implementation, to be used in conjunction with CC_IP keyword. Helps with convergence.\n    '''
        if value == "":
            if "CC_IP_FILTER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP_FILTER"]
                print "Keyword removed."
        elif value == "show":
            if "CC_IP_FILTER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP_FILTER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_IP_FILTER"]=value.lower()


    def cc_ip_proper(self,value="show"):
        '''\nName: CC_IP_PROPER\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If TRUE, calculates proper EOM-IP-CCSD excitation energies and properties. This implementation does not use the diffuse orbital and is the recommended method of doing EOM-IP-CCSD calculations. CC_IP and CC_IP_PROPER keywords are mutually exclusive.\n    '''
        if value == "":
            if "CC_IP_PROPER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_IP_PROPER"]
                print "Keyword removed."
        elif value == "show":
            if "CC_IP_PROPER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_IP_PROPER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_IP_PROPER"]=value.lower()


    def cc_mp2no_grad(self,value="show"):
        '''\nName: CC_MP2NO_GRAD\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If CC_MP2NO_GUESS is TRUE, what kind of one-particle density matrix is used to make the guess orbitals? \nRecommendation: : The two definitions give generally similar performance.    '''
        if value == "":
            if "CC_MP2NO_GRAD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MP2NO_GRAD"]
                print "Keyword removed."
        elif value == "show":
            if "CC_MP2NO_GRAD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MP2NO_GRAD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_MP2NO_GRAD"]=value.lower()


    def cc_mp2no_guess(self,value="show"):
        '''\nName: CC_MP2NO_GUESS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Will guess orbitals be natural orbitals of the MP1 wavefunction? Alternatively, it is possible to use an effective one-particle density matrix to define the natural orbitals. \n    '''
        if value == "":
            if "CC_MP2NO_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MP2NO_GUESS"]
                print "Keyword removed."
        elif value == "show":
            if "CC_MP2NO_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MP2NO_GUESS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_MP2NO_GUESS"]=value.lower()


    def cc_preconv_doubles(self,value="show"):
        '''\nName: CC_PRECONV_DOUBLES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: When TRUE, doubly-excited vectors are converged prior to a full excited states calculation. \nRecommendation: : Occasionally necessary to ensure a doubly excited state is found.    '''
        if value == "":
            if "CC_PRECONV_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_DOUBLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONV_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_DOUBLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONV_DOUBLES"]=value.lower()


    def cc_preconv_singles(self,value="show"):
        '''\nName: CC_PRECONV_SINGLES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: When TRUE, singly-excited vectors are converged prior to a full excited states calculation. \n    '''
        if value == "":
            if "CC_PRECONV_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONV_SINGLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONV_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONV_SINGLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONV_SINGLES"]=value.lower()


    def cc_prop(self,value="show"):
        '''\nName: CC_PROP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments < X2>, < Y2>, and < Z2> of electron density, and the total < R2> = < X2> +< Y2> +< Z2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.\nRecommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PROP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PROP"]=value.lower()


    def cc_restart(self,value="show"):
        '''\nName: CC_RESTART\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Allows an optimized orbital coupled cluster calculation to begin with an initial guess for the orbital transformation matrix U other than the unit vector. The scratch file from a previous run must be available for the U matrix to be read successfully. \nRecommendation: : Useful for restarting a job that did not converge, if files were saved.    '''
        if value == "":
            if "CC_RESTART" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTART"]
                print "Keyword removed."
        elif value == "show":
            if "CC_RESTART" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTART"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_RESTART"]=value.lower()


    def cc_restart_no_scf(self,value="show"):
        '''\nName: CC_RESTART_NO_SCF\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Should an optimized orbital coupled cluster calculation begin with optimized orbitals from a previous calculation? When TRUE, molecular orbitals are initially orthogonalized, and CC_PRECONV_T2Z and CC_CANONIZE are set to TRUE while other guess options are set to FALSE.\n    '''
        if value == "":
            if "CC_RESTART_NO_SCF" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_RESTART_NO_SCF"]
                print "Keyword removed."
        elif value == "show":
            if "CC_RESTART_NO_SCF" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_RESTART_NO_SCF"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_RESTART_NO_SCF"]=value.lower()


    def cc_sd_3(self,value="show"):
        '''\nName: CC_SD_3\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: This keyword initializes calculation of non-iterative triples corrections (fT) and (dT) for EE or SF after the EOM-CCSD is done.\n    '''
        if value == "":
            if "CC_SD_3" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SD_3"]
                print "Keyword removed."
        elif value == "show":
            if "CC_SD_3" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SD_3"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_SD_3"]=value.lower()


    def cc_spin_flip(self,value="show"):
        '''\nName: CC_SPIN_FLIP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects whether do perform a standard excited state calculation, or a spin-flip calculation. Spin multiplicity should be set to 3 for systems with an even number of electrons, and 4 for systems with an odd number of electrons.\n    '''
        if value == "":
            if "CC_SPIN_FLIP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SPIN_FLIP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_SPIN_FLIP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SPIN_FLIP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_SPIN_FLIP"]=value.lower()


    def cc_symmetry(self,value="show"):
        '''\nName: CC_SYMMETRY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the use of symmetry in coupled-cluster calculations\nRecommendation: : It is automatically turned off for any finite difference calculations, e.g. second derivatives.    '''
        if value == "":
            if "CC_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_SYMMETRY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_SYMMETRY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_SYMMETRY"]=value.lower()


    def cc_eom_amplitude_response(self,value="show"):
        '''\nName: CC_EOM_AMPLITUDE_RESPONSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If set to TRUE, adds amplitude response terms to one-particle and two-particle EOM-CCSD density matrices before calculation of properties. CC_EXSTATES_PROP must be set to TRUE.\nRecommendation: : The cost is always about the cost of an analytic gradient calculation for each state, independent of whether or not the two-particle properties are requested. Besides, adding amplitude response terms without orbital response will unlikely improve the quality of the properties. However, it can be used for debugging purposes.    '''
        if value == "":
            if "CC_EOM_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_AMPLITUDE_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_AMPLITUDE_RESPONSE"]=value.lower()


    def cc_eom_properties(self,value="show"):
        '''\nName: CC_EOM_PROPERTIES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the non-relaxed (expectation value) one-particle EOM-CCSD target state properties will be calculated. The properties currently include permanent dipole moment, the second moments 2>, 2>, and 2> of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.\nRecommendation: : Additional equations (EOM-CCSD equations for the left eigenvectors) need to be solved for properties, approximately doubling the cost of calculation for each irrep. Sometimes the equations for left and right eigenvectors converge to different sets of target states. In this case, the simultaneous iterations of left and right vectors will diverge, and the properties for several or all the target states may be incorrect! The problem can be solved by varying the number of requested states, specified with CC_NLOWSPIN and CC_NHIGHSPIN, or the number of guess vectors (CC_NGUESS_SINGLES). The cost of the one-particle properties calculation itself is low. The one-particle density of an EOM-CCSD target state can be analyzed with NBO package by specifying the state with CC_REFSYM and CC_STATE_DERIV and requesting NBO=TRUE and CC_EXSTATES_PROP=TRUE.    '''
        if value == "":
            if "CC_EOM_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROPERTIES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROPERTIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_PROPERTIES"]=value.lower()


    def cc_eom_full_response(self,value="show"):
        '''\nName: CC_EOM_FULL_RESPONSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If set to TRUE, adds both amplitude and orbital response terms to one- and two-particle EOM-CCSD density matrices before calculation of the properties. CC_EXSTATES_PROP must be set to TRUE. If both CC_EOM_AMPL_RESP=TRUE and CC_EOM_FULL_RESP=TRUE, the CC_EOM_AMPL_RESP=TRUE will be ignored.\nRecommendation: : The cost for the full response properties calculation is about the same as the cost of the analytic gradient for each state. Adding full response terms improves quality of calculated properties, but usually it is a small but expensive correction. Use it only if you really need accurate properties.    '''
        if value == "":
            if "CC_EOM_FULL_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_FULL_RESPONSE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_FULL_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_FULL_RESPONSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_FULL_RESPONSE"]=value.lower()


    def cc_canonize_frequency(self,value="show"):
        '''\nName: CC_CANONIZE_FREQUENCY\nType: INTEGER\nDefault: 2\nOptions: 1:100:50:1\nDescription: The orbitals will be semi-canonicalized every n theta resets. The thetas (orbital rotation angles) are reset every CC_RESET_THETA iterations. The counting of iterations differs for active space (VOD, VQCCD) calculations, where the orbitals are always canonicalized at the first theta-reset.\nRecommendation: : Smaller values can be tried in cases that do not converge.    '''
        if value == "":
            if "CC_CANONIZE_FREQUENCY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CANONIZE_FREQUENCY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CANONIZE_FREQUENCY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CANONIZE_FREQUENCY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CANONIZE_FREQUENCY"]=value.lower()


    def cc_eom_transition_properties(self,value="show"):
        '''\nName: CC_EOM_TRANSITION_PROPERTIES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_REFSYM and CC_STATE_DERIV must be specified for this state.\nRecommendation: : Additional equations (for the left EOM-CCSD eigenvectors plus lambda CCSD equations in case if transition properties between the CCSD reference and EOM-CCSD target states are requested) need to be solved for transition properties, approximately doubling the computational cost. The cost of the transition properties calculation itself is low.    '''
        if value == "":
            if "CC_EOM_TRANSITION_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TRANSITION_PROPERTIES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_TRANSITION_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TRANSITION_PROPERTIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_TRANSITION_PROPERTIES"]=value.lower()


    def cc_convergence_energy(self,value="show"):
        '''\nName: CC_CONVERGENCE_ENERGY\nType: INTEGER\nDefault: 2\nOptions: 0:12:10:1\nDescription: Convergence desired on the change in total energy, between iterations.\n    '''
        if value == "":
            if "CC_CONVERGENCE_ENERGY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_ENERGY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CONVERGENCE_ENERGY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_ENERGY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CONVERGENCE_ENERGY"]=value.lower()


    def cc_convergence_zvector(self,value="show"):
        '''\nName: CC_CONVERGENCE_ZVECTOR\nType: INTEGER\nDefault: 2\nOptions:  0 :12:8:1\nDescription: Convergence criterion on the RMS difference between successive doubles Z-vector amplitudes [10-n].\nRecommendation: : Use Default    '''
        if value == "":
            if "CC_CONVERGENCE_ZVECTOR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CONVERGENCE_ZVECTOR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CONVERGENCE_ZVECTOR"]=value.lower()


    def cc_convergence_amplitudes(self,value="show"):
        '''\nName: CC_CONVERGENCE_AMPLITUDES\nType: INTEGER\nDefault: 2\nOptions: 0  :12:8:1\nDescription: Convergence criterion on the RMS difference between successive sets of coupled-cluster doubles amplitudes [10-n]\nRecommendation: : Use default    '''
        if value == "":
            if "CC_CONVERGENCE_AMPLITUDES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_CONVERGENCE_AMPLITUDES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_CONVERGENCE_AMPLITUDES"]=value.lower()


    def cc_full_response(self,value="show"):
        '''\nName: CC_FULL_RESPONSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: If set to TRUE, adds both amplitude and orbital response terms to one- and two-particle CCSD density matrices before calculation of the properties. CC_PROP must be set to TRUE. If both CC_AMPL_RESP=TRUE and CC_FULL_RESP=TRUE, the CC_AMPL_RESP=TRUE will be ignored.\nRecommendation: : The cost for the full response properties calculation is about the same as the cost of the analytic gradient. Adding full response terms improves quality of calculated properties, but usually it is a small but expensive correction. Use it only if you really need accurate properties.    '''
        if value == "":
            if "CC_FULL_RESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_FULL_RESPONSE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_FULL_RESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_FULL_RESPONSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_FULL_RESPONSE"]=value.lower()


    def cc_hessian_thresh(self,value="show"):
        '''\nName: CC_HESSIAN_THRESH\nType: INTEGER\nDefault: 2\nOptions: 0.001:1.000:0.010:0.001\nDescription: Minimum alloed value for the orbital Hessian.  Smaller values are replaced with this constant.\n    '''
        if value == "":
            if "CC_HESSIAN_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_HESSIAN_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_HESSIAN_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_HESSIAN_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_HESSIAN_THRESH"]=value.lower()


    def cc_high_spin(self,value="show"):
        '''\nName: CC_HIGH_SPIN\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of high-spin excited state roots to find.  Works only for singlet reference states and triplet excited states.  The default is to not look for any excited states.  Setting this to [i,j,k...] looks for i excited states in the first irrep, j states in the second irrep and so on.\n    '''
        if value == "":
            if "CC_HIGH_SPIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_HIGH_SPIN"]
                print "Keyword removed."
        elif value == "show":
            if "CC_HIGH_SPIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_HIGH_SPIN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_HIGH_SPIN"]=value.lower()


    def cc_low_spin(self,value="show"):
        '''\nName: CC_LOW_SPIN\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of low-spin excited state roots to find.  In the cas of closed-shell reference states, excited singlet states will be found.  For any other reference state all states (e.g. both singlet and triplet) will be calculated.  The default is to not look for any excited states.  Setting this to [i,j,k...] looks for i excited states in the first irrep, j states in the second irrep and so on.\n    '''
        if value == "":
            if "CC_LOW_SPIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_LOW_SPIN"]
                print "Keyword removed."
        elif value == "show":
            if "CC_LOW_SPIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_LOW_SPIN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_LOW_SPIN"]=value.lower()


    def cc_preconverge_doubles(self,value="show"):
        '''\nName: CC_PRECONVERGE_DOUBLES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: When TRUE, doubly-excited vectors are converged prior to a full excited states calculation. \nRecommendation: : Occasionally necessary to ensure a doubly excited state is found.    '''
        if value == "":
            if "CC_PRECONVERGE_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_DOUBLES"]=value.lower()


    def cc_preconverge_fz(self,value="show"):
        '''\nName: CC_PRECONVERGE_FZ\nType: INTEGER\nDefault: 0\nOptions: 0}  {0   \nDescription: In active space methods, whether to preconverge other wavefunction variables for fixed initial guess of active space.\nRecommendation: :     '''
        if value == "":
            if "CC_PRECONVERGE_FZ" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_FZ"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_FZ" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_FZ"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_FZ"]=value.lower()


    def cc_preconverge_sd(self,value="show"):
        '''\nName: CC_PRECONVERGE_SD\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Solves the EOM-CCSD equations, prints energies, then uses EOM-CCSD vectors as initial vectors in EOM-CC(2,3). Very convenient for calculations using energy additivity schemes.\nRecommendation: : Turning this option on is recommended    '''
        if value == "":
            if "CC_PRECONVERGE_SD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_SD"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_SD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_SD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_SD"]=value.lower()


    def cc_preconverge_singles(self,value="show"):
        '''\nName: CC_PRECONVERGE_SINGLES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: When TRUE, singly-excited vectors are converged prior to a full excited states calculation. \n    '''
        if value == "":
            if "CC_PRECONVERGE_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_SINGLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_SINGLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_SINGLES"]=value.lower()


    def cc_preconverge_t2z(self,value="show"):
        '''\nName: CC_PRECONVERGE_T2Z\nType: INTEGER\nDefault: 1\nOptions: 0:0  \nDescription: Whether to pre-converge the cluster amplitudes before beginning orbital optimization in optimized orbital cluster methods.\nRecommendation: : Experiment with this option in cases of convergence failure.    '''
        if value == "":
            if "CC_PRECONVERGE_T2Z" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_T2Z"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_T2Z" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_T2Z"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_T2Z"]=value.lower()


    def cc_preconverge_t2z_each(self,value="show"):
        '''\nName: CC_PRECONVERGE_T2Z_EACH\nType: INTEGER\nDefault: 1\nOptions: 0:0   \nDescription: Whether to pre-converge the cluster amplitudes before each change of the orbitals in optimized orbital coupled-cluster methods. The maximum number of iterations in this pre-convergence procedure is given by the value of this parameter.\nRecommendation: : A very slow last resort option for jobs that do not converge.    '''
        if value == "":
            if "CC_PRECONVERGE_T2Z_EACH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRECONVERGE_T2Z_EACH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRECONVERGE_T2Z_EACH"]=value.lower()


    def qui_solvent_none(self,value="show"):
        '''\nName: QUI_SOLVENT_NONE\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Checking this disables all solvent models.\n    '''
        if value == "":
            if "QUI_SOLVENT_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_NONE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_NONE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_NONE"]=value.lower()


    def chemsol_efield(self,value="show"):
        '''\nName: CHEMSOL_EFIELD\nType: STRING\nDefault: 0\nOptions: Exact//1:Mulliken//0\nDescription: Determines how the solute charge distribution is approximated in evaluating the electrostatic field of the solute.  Either the exact solute charge distribution is used, or the charge distribution is approximated by Mulliken atomic charges. 
Recommentation: Mulliken charges are faster, but less rigorous.\n    '''
        if value == "":
            if "CHEMSOL_EFIELD" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_EFIELD"]
                print "Keyword removed."
        elif value == "show":
            if "CHEMSOL_EFIELD" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_EFIELD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHEMSOL_EFIELD"]=value.lower()


    def cis_state_derivative(self,value="show"):
        '''\nName: CIS_STATE_DERIVATIVE\nType: INTEGER\nDefault: 2\nOptions: 0:200:0:1\nDescription: Sets which CIS state to use for excited state optimizations and vibrational analysis.\nRecommendation: : Check to see that the states do no change order during an optimization.    '''
        if value == "":
            if "CIS_STATE_DERIVATIVE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_STATE_DERIVATIVE"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_STATE_DERIVATIVE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_STATE_DERIVATIVE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_STATE_DERIVATIVE"]=value.lower()


    def rpa(self,value="show"):
        '''\nName: RPA\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Do an RPA calculation in addition to a CIS calculation\n    '''
        if value == "":
            if "RPA" in self.dict_of_keywords:
                del self.dict_of_keywords["RPA"]
                print "Keyword removed."
        elif value == "show":
            if "RPA" in self.dict_of_keywords:
                return self.dict_of_keywords["RPA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPA"]=value.lower()


    def cis_guess_disk_type(self,value="show"):
        '''\nName: CIS_GUESS_DISK_TYPE\nType: STRING\nDefault: 0\nOptions: None//-1:Triplets only//0:Triplets and singlets//1:Singlets only//2\nDescription: Determines the type of guesses to be read from disk\nRecommendation: : Must be specified if a CIS guess in to be read from disk.    '''
        if value == "":
            if "CIS_GUESS_DISK_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_GUESS_DISK_TYPE"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_GUESS_DISK_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_GUESS_DISK_TYPE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_GUESS_DISK_TYPE"]=value.lower()


    def cis_ras_cutoff_occupied(self,value="show"):
        '''\nName: CIS_RAS_CUTOFF_OCCUPIED\nType: INTEGER\nDefault: 2\nOptions: 0.00:2.00:0.50:0.01\nDescription: Specifies the occupied orbital cutoff\n    '''
        if value == "":
            if "CIS_RAS_CUTOFF_OCCUPIED" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS_CUTOFF_OCCUPIED" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS_CUTOFF_OCCUPIED"]=value.lower()


    def cis_ras_cutoff_virtual(self,value="show"):
        '''\nName: CIS_RAS_CUTOFF_VIRTUAL\nType: INTEGER\nDefault: 2\nOptions: 0:1.00:0:0.01\nDescription: Specifies the virtual orbital cutoff.\n    '''
        if value == "":
            if "CIS_RAS_CUTOFF_VIRTUAL" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS_CUTOFF_VIRTUAL" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS_CUTOFF_VIRTUAL"]=value.lower()


    def cis_ras(self,value="show"):
        '''\nName: CIS_RAS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls whether reduced single excitation space is used\n    '''
        if value == "":
            if "CIS_RAS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS"]=value.lower()


    def cis_ras_n_solute_atoms(self,value="show"):
        '''\nName: CIS_RAS_N_SOLUTE_ATOMS\nType: INTEGER\nDefault: 0\nOptions: No default\nDescription: Specifies number of atoms or orbitals in solute\n    '''
        if value == "":
            if "CIS_RAS_N_SOLUTE_ATOMS" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS_N_SOLUTE_ATOMS" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS_N_SOLUTE_ATOMS"]=value.lower()


    def cis_ras_print(self,value="show"):
        '''\nName: CIS_RAS_PRINT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects whether or not to print additional output.\n    '''
        if value == "":
            if "CIS_RAS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS_PRINT"]=value.lower()


    def cis_ras_type(self,value="show"):
        '''\nName: CIS_RAS_TYPE\nType: STRING\nDefault: 0\nOptions: Localized Orbitals//1:User-defined//2\nDescription: Controls how reduced subspace is specified\n    '''
        if value == "":
            if "CIS_RAS_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_RAS_TYPE"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_RAS_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_RAS_TYPE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_RAS_TYPE"]=value.lower()


    def diis_max_cycles(self,value="show"):
        '''\nName: DIIS_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 1:100:50:1\nDescription: The maximum number of DIIS iterations before switching to (geometric) direct minimization when SCF_ALGORITHM is DIIS_GDM or DIIS_DM. See also THRESH_DIIS_SWITCH. \n    '''
        if value == "":
            if "DIIS_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "DIIS_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIIS_MAX_CYCLES"]=value.lower()


    def diis_subspace_size(self,value="show"):
        '''\nName: DIIS_SUBSPACE_SIZE\nType: INTEGER\nDefault: 2\nOptions: 1:50:15:1\nDescription: Controls the size of the DIIS and/or RCA subspace during the SCF.\n    '''
        if value == "":
            if "DIIS_SUBSPACE_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_SUBSPACE_SIZE"]
                print "Keyword removed."
        elif value == "show":
            if "DIIS_SUBSPACE_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_SUBSPACE_SIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIIS_SUBSPACE_SIZE"]=value.lower()


    def diis_error_metric(self,value="show"):
        '''\nName: DIIS_ERROR_METRIC\nType: STRING\nDefault: 0\nOptions: Maximum//false:RMS//true\nDescription: Changes the DIIS convergence metric from the maximum to the RMS error.\nRecommendation: : Use default, the maximum error provides a more reliable criterion.    '''
        if value == "":
            if "DIIS_ERROR_METRIC" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_ERROR_METRIC"]
                print "Keyword removed."
        elif value == "show":
            if "DIIS_ERROR_METRIC" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_ERROR_METRIC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIIS_ERROR_METRIC"]=value.lower()


    def dma(self,value="show"):
        '''\nName: DMA\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Specifies whether to perform Distributed Multipole Analysis.\n    '''
        if value == "":
            if "DMA" in self.dict_of_keywords:
                del self.dict_of_keywords["DMA"]
                print "Keyword removed."
        elif value == "show":
            if "DMA" in self.dict_of_keywords:
                return self.dict_of_keywords["DMA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DMA"]=value.lower()


    def raman(self,value="show"):
        '''\nName: RAMAN\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls calculation of Raman intensities. Only relevant for a frequency calculation.\n    '''
        if value == "":
            if "RAMAN" in self.dict_of_keywords:
                del self.dict_of_keywords["RAMAN"]
                print "Keyword removed."
        elif value == "show":
            if "RAMAN" in self.dict_of_keywords:
                return self.dict_of_keywords["RAMAN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RAMAN"]=value.lower()


    def aimd_time_step(self,value="show"):
        '''\nName: AIMD_TIME_STEP\nType: INTEGER\nDefault: 0\nOptions: None.\nDescription: Specifies the molecular dynamics time step, in atomic units (1 a.u. = 0.0242 fs).\nRecommendation: : Smaller time steps lead to better energy conservation; too large a time step may cause the job to fail entirely. Make the time step as large as possible, consistent with tolerable energy conservation.    '''
        if value == "":
            if "AIMD_TIME_STEP" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_TIME_STEP"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_TIME_STEP" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_TIME_STEP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_TIME_STEP"]=value.lower()


    def basis2(self,value="show"):
        '''\nName: BASIS2\nType: STRING\nDefault: 0\nOptions: None:----------:STO-3G:STO-6G:3-21G:4-21G:6-31G:6-31G(d):cc-pVDZ:----------:6-311+G**:rcc-pVTZ:rcc-pVQZ\nDescription: Selects either a small basis set to use in basis set projection for the initial guess, or a subset basis for dual basis set calculations.\n    '''
        if value == "":
            if "BASIS2" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS2"]
                print "Keyword removed."
        elif value == "show":
            if "BASIS2" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["BASIS2"]=value.lower()


    def dscf_convergence_level_1(self,value="show"):
        '''\nName: DSCF_CONVERGENCE_LEVEL_1\nType: INTEGER\nDefault: 2\nOptions: 0:10:4:1\nDescription: Sets the convergence criterion for the level-1 iterations. This preconditions the density for the level-2 calculation, and does not include any two-electron integrals. \nRecommendation: : The criterion for level-1 convergence must be less than or equal to the level-2 criterion, otherwise the D-CPSCF will not converge.    '''
        if value == "":
            if "DSCF_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"]
                print "Keyword removed."
        elif value == "show":
            if "DSCF_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_1"]=value.lower()


    def dscf_convergence_level_2(self,value="show"):
        '''\nName: DSCF_CONVERGENCE_LEVEL_2\nType: INTEGER\nDefault: 2\nOptions: 0:10:4:1\nDescription: Sets the convergence criterion for the level-2 iterations.\n    '''
        if value == "":
            if "DSCF_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"]
                print "Keyword removed."
        elif value == "show":
            if "DSCF_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DSCF_CONVERGENCE_LEVEL_2"]=value.lower()


    def dscf_diis_subspace(self,value="show"):
        '''\nName: DSCF_DIIS_SUBSPACE\nType: INTEGER\nDefault: 2\nOptions: 0:25:11:1\nDescription: Specifies the number of matrices to use in the DIIS extrapolation in the D-CPSCF.\nRecommendation: : Use the default.    '''
        if value == "":
            if "DSCF_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_DIIS_SUBSPACE"]
                print "Keyword removed."
        elif value == "show":
            if "DSCF_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_DIIS_SUBSPACE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DSCF_DIIS_SUBSPACE"]=value.lower()


    def dcpscf_pertnum(self,value="show"):
        '''\nName: DCPSCF_PERTNUM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Specifies whether to do the perturbations all together.\n    '''
        if value == "":
            if "DCPSCF_PERTNUM" in self.dict_of_keywords:
                del self.dict_of_keywords["DCPSCF_PERTNUM"]
                print "Keyword removed."
        elif value == "show":
            if "DCPSCF_PERTNUM" in self.dict_of_keywords:
                return self.dict_of_keywords["DCPSCF_PERTNUM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DCPSCF_PERTNUM"]=value.lower()


    def dcpscf_perturbations(self,value="show"):
        '''\nName: DCPSCF_PERTURBATIONS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Specifies whether to do the perturbations all together.\n    '''
        if value == "":
            if "DCPSCF_PERTURBATIONS" in self.dict_of_keywords:
                del self.dict_of_keywords["DCPSCF_PERTURBATIONS"]
                print "Keyword removed."
        elif value == "show":
            if "DCPSCF_PERTURBATIONS" in self.dict_of_keywords:
                return self.dict_of_keywords["DCPSCF_PERTURBATIONS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DCPSCF_PERTURBATIONS"]=value.lower()


    def fd_derivative_type(self,value="show"):
        '''\nName: FD_DERIVATIVE_TYPE\nType: STRING\nDefault: 0\nOptions: Energies Only//0:Gradients Only//1:Hessians Only//2:Enegies, Gradients and Hessians//3 \nDescription: Controls what types of gradient information are used to compute higher derivatives. The default uses a combination of energy, gradient and Hessian information, which makes the force field calculation faster. \nRecommendation: : When the molecule is larger than benzene with small basis set, using only Hessian information may be faster. Note that this option will be set lower if analytic derivatives of the requested order are not available.     '''
        if value == "":
            if "FD_DERIVATIVE_TYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["FD_DERIVATIVE_TYPE"]
                print "Keyword removed."
        elif value == "show":
            if "FD_DERIVATIVE_TYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["FD_DERIVATIVE_TYPE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FD_DERIVATIVE_TYPE"]=value.lower()


    def fd_step_size(self,value="show"):
        '''\nName: FD_STEP_SIZE\nType: INTEGER\nDefault: 2\nOptions: 0.0001 :0.0100:0.0010:0.0001\nDescription: Displacement used for calculating derivatives by finite difference.\nRecommendation: : Use default, unless on a very flat potential, in which case a larger value may be required.    '''
        if value == "":
            if "FD_STEP_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["FD_STEP_SIZE"]
                print "Keyword removed."
        elif value == "show":
            if "FD_STEP_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["FD_STEP_SIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FD_STEP_SIZE"]=value.lower()


    def analytic_derivative_order(self,value="show"):
        '''\nName: ANALYTIC_DERIVATIVE_ORDER\nType: INTEGER\nDefault: 2\nOptions: 0:2:2:1\nDescription: Controls the order of derivatives that are evaluated analytically. The user is not normally required to specify a value, unless numerical derivatives are desired. The derivatives will be evaluated numerically if this option is set lower than the type of job requires.\nRecommendation: : Usually set to the maximum possible for efficiency. Note that this option will be set lower if analytic derivatives of the requested order are not available.    '''
        if value == "":
            if "ANALYTIC_DERIVATIVE_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "ANALYTIC_DERIVATIVE_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["ANALYTIC_DERIVATIVE_ORDER"]=value.lower()


    def pao_algorithm(self,value="show"):
        '''\nName: PAO_ALGORITHM\nType: STRING\nDefault: 0\nOptions: Efficient//0:Conservative//1\nDescription: Algorithm used to optimize polarized atomic orbitals (see PAO_METHOD)\n    '''
        if value == "":
            if "PAO_ALGORITHM" in self.dict_of_keywords:
                del self.dict_of_keywords["PAO_ALGORITHM"]
                print "Keyword removed."
        elif value == "show":
            if "PAO_ALGORITHM" in self.dict_of_keywords:
                return self.dict_of_keywords["PAO_ALGORITHM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PAO_ALGORITHM"]=value.lower()


    def pao_method(self,value="show"):
        '''\nName: PAO_METHOD\nType: STRING\nDefault: 0\nOptions: EPAO :PAO\nDescription: Controls evaluation of polarized atomic orbitals (PAOs).\n    '''
        if value == "":
            if "PAO_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["PAO_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "PAO_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["PAO_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PAO_METHOD"]=value.lower()


    def epao_weights(self,value="show"):
        '''\nName: EPAO_WEIGHTS\nType: STRING\nDefault: 0\nOptions: 1st & 2nd order//115:1st order only//15\nDescription: Controls algorithm and weights for EPAO calculations (see PAO_METHOD).\nRecommendation: : Use default, unless convergence failure is encountered.    '''
        if value == "":
            if "EPAO_WEIGHTS" in self.dict_of_keywords:
                del self.dict_of_keywords["EPAO_WEIGHTS"]
                print "Keyword removed."
        elif value == "show":
            if "EPAO_WEIGHTS" in self.dict_of_keywords:
                return self.dict_of_keywords["EPAO_WEIGHTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EPAO_WEIGHTS"]=value.lower()


    def aimd_fock_extrapolation_order(self,value="show"):
        '''\nName: AIMD_FOCK_EXTRAPOLATION_ORDER\nType: INTEGER\nDefault: 2\nOptions: 0:20:0:1\nDescription: Specifies the polynomial order N for Fock matrix extrapolation.\n    '''
        if value == "":
            if "AIMD_FOCK_EXTRAPOLATION_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_FOCK_EXTRAPOLATION_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_ORDER"]=value.lower()


    def aimd_fock_extrapolation_points(self,value="show"):
        '''\nName: AIMD_FOCK_EXTRAPOLATION_POINTS\nType: INTEGER\nDefault: 2\nOptions: 0   :20:0:1\nDescription: Specifies the number M of old Fock matrices that are retained for use in extrapolation. \nRecommendation: : Higher-order extrapolations with more saved Fock matrices are faster and conserve energy better than low-order extrapolations, up to a point. In many cases, the scheme (N = 6, M = 12), in conjunction with SCF_CONVERGENCE = 6, is found to provide about a 50% savings in computational cost while still conserving energy.    '''
        if value == "":
            if "AIMD_FOCK_EXTRAPOLATION_POINTS" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_POINTS"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_FOCK_EXTRAPOLATION_POINTS" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_POINTS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_FOCK_EXTRAPOLATION_POINTS"]=value.lower()


    def ftc_fast(self,value="show"):
        '''\nName: FTC_FAST\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls whether or not the operator is evaluated on a large grid and stored in memory to speed up the calculation.\nRecommendation: : Use the default if possible.  Turning this option off conserves some memory, but causes a slow down in speed.    '''
        if value == "":
            if "FTC_FAST" in self.dict_of_keywords:
                del self.dict_of_keywords["FTC_FAST"]
                print "Keyword removed."
        elif value == "show":
            if "FTC_FAST" in self.dict_of_keywords:
                return self.dict_of_keywords["FTC_FAST"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["FTC_FAST"]=value.lower()


    def cis_max_cycles(self,value="show"):
        '''\nName: CIS_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 0:100:30:1\nDescription: Maximum number of CIS iterative cycles allowed\nRecommendation: : Default is usually sufficient.    '''
        if value == "":
            if "CIS_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CIS_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "CIS_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CIS_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CIS_MAX_CYCLES"]=value.lower()


    def dscf_max_cycles_level_2(self,value="show"):
        '''\nName: DSCF_MAX_CYCLES_LEVEL_2\nType: INTEGER\nDefault: 2\nOptions: 1:500:30:1\nDescription: Sets the maximum number of level-2 iterations.\nRecommendation: : Use default.    '''
        if value == "":
            if "DSCF_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"]
                print "Keyword removed."
        elif value == "show":
            if "DSCF_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_2"]=value.lower()


    def dscf_max_cycles_level_1(self,value="show"):
        '''\nName: DSCF_MAX_CYCLES_LEVEL_1\nType: INTEGER\nDefault: 2\nOptions: 0:500:100:10\nDescription: Sets the maximum number of level-1 iterations.\nRecommendation: : Use default.    '''
        if value == "":
            if "DSCF_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"]
                print "Keyword removed."
        elif value == "show":
            if "DSCF_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DSCF_MAX_CYCLES_LEVEL_1"]=value.lower()


    def geom_opt_diis_subspace(self,value="show"):
        '''\nName: GEOM_OPT_DIIS_SUBSPACE\nType: INTEGER\nDefault: 2\nOptions: -1:50:0:1\nDescription: Controls maximum size of subspace for GDIIS. 0 turns off GDIIS and -1 causes the program to select the default size.\n    '''
        if value == "":
            if "GEOM_OPT_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_DIIS_SUBSPACE"]=value.lower()


    def geom_opt_symmetry(self,value="show"):
        '''\nName: GEOM_OPT_SYMMETRY\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls the use of point-group symmetry in the optimization.\n    '''
        if value == "":
            if "GEOM_OPT_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_SYMMETRY"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_SYMMETRY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_SYMMETRY"]=value.lower()


    def geom_opt_coordinates(self,value="show"):
        '''\nName: GEOM_OPT_COORDINATES\nType: STRING\nDefault: 0\nOptions: Cartesian//0:Internal//1:Z-matrix//2\nDescription: Controls the type of optimization coordinates.\nRecommendation: : Use the default; delocalized internals are more efficient.    '''
        if value == "":
            if "GEOM_OPT_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_COORDINATES"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_COORDINATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_COORDINATES"]=value.lower()


    def geom_opt_max_step_size(self,value="show"):
        '''\nName: GEOM_OPT_MAX_STEP_SIZE\nType: INTEGER\nDefault: 2\nOptions: 0.001:0.999:0.300:0.001\nDescription: Maximum allowed step size in the geometry optimization.\n    '''
        if value == "":
            if "GEOM_OPT_MAX_STEP_SIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_MAX_STEP_SIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_MAX_STEP_SIZE"]=value.lower()


    def geom_opt_hessian_update(self,value="show"):
        '''\nName: GEOM_OPT_HESSIAN_UPDATE\nType: STRING\nDefault: 0\nOptions: Default//-1:No Update//0:Murtagh-Sargent//1:Powell//2:Powell/Murtagh-Sargent//3:BFGS//4:BFGS w/ safeguards//5\nDescription: Controls the Hessian update algorithm.  The default dpends on the type of job.\n    '''
        if value == "":
            if "GEOM_OPT_HESSIAN_UPDATE" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_HESSIAN_UPDATE" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_HESSIAN_UPDATE"]=value.lower()


    def cfmm_grain(self,value="show"):
        '''\nName: CFMM_GRAIN\nType: STRING\nDefault: 0\nOptions: Automatic//-1:Off//1:8:9:10:11:12:13:14:15:16\nDescription: Controls the number of lowest-level boxes in one dimension for CFMM.\nRecommendation: : This is an expert option; either use the default, or use a value of 1 if CFMM is not desired.    '''
        if value == "":
            if "CFMM_GRAIN" in self.dict_of_keywords:
                del self.dict_of_keywords["CFMM_GRAIN"]
                print "Keyword removed."
        elif value == "show":
            if "CFMM_GRAIN" in self.dict_of_keywords:
                return self.dict_of_keywords["CFMM_GRAIN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CFMM_GRAIN"]=value.lower()


    def plots_property(self,value="show"):
        '''\nName: PLOTS_PROPERTY\nType: STRING\nDefault: 0\nOptions: None:ESP Only//0:ESP & EFIELD//1:EFIELD Only//2\nDescription: Triggers the calculation of the electrostatic potential and/or the electric field at the points given in the file ESPGrid. \nRecommendation: : Must use this option when IGDESP is specified.    '''
        if value == "":
            if "PLOTS_PROPERTY" in self.dict_of_keywords:
                del self.dict_of_keywords["PLOTS_PROPERTY"]
                print "Keyword removed."
        elif value == "show":
            if "PLOTS_PROPERTY" in self.dict_of_keywords:
                return self.dict_of_keywords["PLOTS_PROPERTY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PLOTS_PROPERTY"]=value.lower()


    def use_case(self,value="show"):
        '''\nName: USE_CASE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects the atenuated coulomb operator (CASE approximation).\n    '''
        if value == "":
            if "USE_CASE" in self.dict_of_keywords:
                del self.dict_of_keywords["USE_CASE"]
                print "Keyword removed."
        elif value == "show":
            if "USE_CASE" in self.dict_of_keywords:
                return self.dict_of_keywords["USE_CASE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["USE_CASE"]=value.lower()


    def intracule_conserve_memory(self,value="show"):
        '''\nName: INTRACULE_CONSERVE_MEMORY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Reduce memory required in the evaluation of W(u,v). \nRecommendation: : The low memory option is slower, use default unless memory is limited.    '''
        if value == "":
            if "INTRACULE_CONSERVE_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_CONSERVE_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_CONSERVE_MEMORY"]=value.lower()


    def intracule(self,value="show"):
        '''\nName: INTRACULE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls whether intracule properties are calculated.  Setting this option causes the data in $intracule to be activated.\n    '''
        if value == "":
            if "INTRACULE" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE"]=value.lower()


    def intracule_grid(self,value="show"):
        '''\nName: INTRACULE_GRID\nType: STRING\nDefault: 10\nOptions: 6:18:26:38:50:74:86:110:146:170:194:230:266:302:350:434:590:770:974:1202:1454:1730:2030:2354:2702:3074:3470:3890:4334:4802:5294\nDescription: Specify angular Lebedev grid for Wigner intracule calculations.\nRecommendation: : Larger grids if high accuracy required.    '''
        if value == "":
            if "INTRACULE_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_GRID"]=value.lower()


    def intracule_wigner_series_limit(self,value="show"):
        '''\nName: INTRACULE_WIGNER_SERIES_LIMIT\nType: INTEGER\nDefault: 2\nOptions: 1:100:10:1\nDescription: Sets summation limit for Wigner integrals.\nRecommendation: : Increase n for greater accuracy.    '''
        if value == "":
            if "INTRACULE_WIGNER_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_WIGNER_SERIES_LIMIT"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_WIGNER_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_WIGNER_SERIES_LIMIT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_WIGNER_SERIES_LIMIT"]=value.lower()


    def intracule_j_series_limit(self,value="show"):
        '''\nName: INTRACULE_J_SERIES_LIMIT\nType: INTEGER\nDefault: 2\nOptions: 1:100:40:1\nDescription: Sets summation limit for series expansion evaluation of j_n(x).\nRecommendation: : Lower values speed up the calculation, but may affect accuracy.    '''
        if value == "":
            if "INTRACULE_J_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_J_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_J_SERIES_LIMIT"]=value.lower()


    def intracule_i_series_limit(self,value="show"):
        '''\nName: INTRACULE_I_SERIES_LIMIT\nType: INTEGER\nDefault: 2\nOptions: 1:100:40:1\nDescription: Sets summation limit for series expansion evaluation of i_n(x).\nRecommendation: : Lower values speed up the calculation, but may affect accuracy.    '''
        if value == "":
            if "INTRACULE_I_SERIES_LIMIT" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_I_SERIES_LIMIT" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_I_SERIES_LIMIT"]=value.lower()


    def intracule_method(self,value="show"):
        '''\nName: INTRACULE_METHOD\nType: STRING\nDefault: 0\nOptions: Series//0:Quadrature//1\nDescription: Use Lebedev quadrature to evaluate Wigner integrals.\n    '''
        if value == "":
            if "INTRACULE_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["INTRACULE_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "INTRACULE_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["INTRACULE_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["INTRACULE_METHOD"]=value.lower()


    def rca_max_cycles(self,value="show"):
        '''\nName: RCA_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 0:100:50:1\nDescription: The maximum number of RCA iterations before switching to DIIS when the SCF algorithm is RCA_DIIS\n    '''
        if value == "":
            if "RCA_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "RCA_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RCA_MAX_CYCLES"]=value.lower()


    def diis_switch_thresh(self,value="show"):
        '''\nName: DIIS_SWITCH_THRESH\nType: INTEGER\nDefault: 0\nOptions: 2\nDescription: The threshold for switching between DIIS extrapolation and direct minimization of the SCF energy is 10-THRESH_DIIS_SWITCH when SCF_ALGORITHM is DIIS_GDM or DIIS_DM. See also MAX_DIIS_MAX_CYCLES\n    '''
        if value == "":
            if "DIIS_SWITCH_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["DIIS_SWITCH_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "DIIS_SWITCH_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["DIIS_SWITCH_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIIS_SWITCH_THRESH"]=value.lower()


    def scf_max_cycles(self,value="show"):
        '''\nName: SCF_MAX_CYCLES\nType: INTEGER\nDefault: 2\nOptions: 0:200:50:1\nDescription: Controls the maximum number of SCF iterations permitted.\nRecommendation: : Increase for slowly converging systems such as those containing transition metals.    '''
        if value == "":
            if "SCF_MAX_CYCLES" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_MAX_CYCLES"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_MAX_CYCLES" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_MAX_CYCLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_MAX_CYCLES"]=value.lower()


    def meteco(self,value="show"):
        '''\nName: METECO\nType: STRING\nDefault: 1\nOptions: Machine precision//1:Thresh//2\nDescription: Sets the threshold criteria for discarding shell-pairs.\nRecommendation: : Use default.    '''
        if value == "":
            if "METECO" in self.dict_of_keywords:
                del self.dict_of_keywords["METECO"]
                print "Keyword removed."
        elif value == "show":
            if "METECO" in self.dict_of_keywords:
                return self.dict_of_keywords["METECO"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["METECO"]=value.lower()


    def mm_charges(self,value="show"):
        '''\nName: MM_CHARGES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Requests the calculation of multipole-derived charges (MDCs).\nRecommendation: : Set to TRUE if MDCs or the traceless form of the multipole moments are desired. The calculation does not take long.    '''
        if value == "":
            if "MM_CHARGES" in self.dict_of_keywords:
                del self.dict_of_keywords["MM_CHARGES"]
                print "Keyword removed."
        elif value == "show":
            if "MM_CHARGES" in self.dict_of_keywords:
                return self.dict_of_keywords["MM_CHARGES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MM_CHARGES"]=value.lower()


    def molden_format(self,value="show"):
        '''\nName: MOLDEN_FORMAT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Requests a $molden-formatted input file containing information from a Q-Chem job.\n    '''
        if value == "":
            if "MOLDEN_FORMAT" in self.dict_of_keywords:
                del self.dict_of_keywords["MOLDEN_FORMAT"]
                print "Keyword removed."
        elif value == "show":
            if "MOLDEN_FORMAT" in self.dict_of_keywords:
                return self.dict_of_keywords["MOLDEN_FORMAT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOLDEN_FORMAT"]=value.lower()


    def mom_print(self,value="show"):
        '''\nName: MOM_PRINT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Switches printing on within the MOM procedure.\n    '''
        if value == "":
            if "MOM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["MOM_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "MOM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["MOM_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOM_PRINT"]=value.lower()


    def moprop_convergence_level_1(self,value="show"):
        '''\nName: MOPROP_CONVERGENCE_LEVEL_1\nType: INTEGER\nDefault: 2\nOptions: 0:12:6:1\nDescription: Sets the convergence criteria for CPSCF and 1st order TDSCF.\n    '''
        if value == "":
            if "MOPROP_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_CONVERGENCE_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_1"]=value.lower()


    def moprop_convergence_level_2(self,value="show"):
        '''\nName: MOPROP_CONVERGENCE_LEVEL_2\nType: INTEGER\nDefault: 2\nOptions: 0:12:6:1\nDescription: Sets the convergence criterium for second-order TDSCF.\n    '''
        if value == "":
            if "MOPROP_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_CONVERGENCE_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_CONVERGENCE_LEVEL_2"]=value.lower()


    def moprop_diis(self,value="show"):
        '''\nName: MOPROP_DIIS\nType: INTEGER\nDefault: 2\nOptions: 0:5:5:5\nDescription: Controls the use of Pulays DIIS.\n    '''
        if value == "":
            if "MOPROP_DIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_DIIS"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_DIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_DIIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_DIIS"]=value.lower()


    def moprop_diis_subspace(self,value="show"):
        '''\nName: MOPROP_DIIS_SUBSPACE\nType: INTEGER\nDefault: 2\nOptions: 0:50:20:1\nDescription: Specified the DIIS subspace dimension.\n    '''
        if value == "":
            if "MOPROP_DIIS_SUBSPACE" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_DIIS_SUBSPACE" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_DIIS_SUBSPACE"]=value.lower()


    def moprop_max_cycles_level_2(self,value="show"):
        '''\nName: MOPROP_MAX_CYCLES_LEVEL_2\nType: INTEGER\nDefault: 2\nOptions: 1:500:50:1\nDescription: The maximal number of iterations for second-order TDSCF.\nRecommendation: : Use default.    '''
        if value == "":
            if "MOPROP_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_MAX_CYCLES_LEVEL_2" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_2"]=value.lower()


    def moprop_max_cycles_level_1(self,value="show"):
        '''\nName: MOPROP_MAX_CYCLES_LEVEL_1\nType: INTEGER\nDefault: 2\nOptions: 1:500:50:1\nDescription: The maximal number of iterations for CPSCF and first-order TDSCF.\nRecommendation: : Use default.    '''
        if value == "":
            if "MOPROP_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_MAX_CYCLES_LEVEL_1" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_MAX_CYCLES_LEVEL_1"]=value.lower()


    def moprop_perturbations(self,value="show"):
        '''\nName: MOPROP_PERTURBATIONS\nType: INTEGER\nDefault: 2\nOptions: 0:20:0:1\nDescription: Set the number of perturbed densities that will to be treated together.\nRecommendation: : Use default    '''
        if value == "":
            if "MOPROP_PERTURBATIONS" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_PERTURBATIONS"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_PERTURBATIONS" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_PERTURBATIONS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_PERTURBATIONS"]=value.lower()


    def nbo(self,value="show"):
        '''\nName: NBO\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the use of the NBO package.\n    '''
        if value == "":
            if "NBO" in self.dict_of_keywords:
                del self.dict_of_keywords["NBO"]
                print "Keyword removed."
        elif value == "show":
            if "NBO" in self.dict_of_keywords:
                return self.dict_of_keywords["NBO"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["NBO"]=value.lower()


    def print_general_basis(self,value="show"):
        '''\nName: PRINT_GENERAL_BASIS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls print out of built in basis sets in input format\nRecommendation: : Useful for modification of standard basis sets.    '''
        if value == "":
            if "PRINT_GENERAL_BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_GENERAL_BASIS"]
                print "Keyword removed."
        elif value == "show":
            if "PRINT_GENERAL_BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_GENERAL_BASIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PRINT_GENERAL_BASIS"]=value.lower()


    def qmmm_charges(self,value="show"):
        '''\nName: QMMM_CHARGES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the printing of QM charges to file.\nRecommendation: : Use default unless running calculations with $charmm where charges on the QM region need to be saved.    '''
        if value == "":
            if "QMMM_CHARGES" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_CHARGES"]
                print "Keyword removed."
        elif value == "show":
            if "QMMM_CHARGES" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_CHARGES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QMMM_CHARGES"]=value.lower()


    def qmmm_print(self,value="show"):
        '''\nName: QMMM_PRINT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls the amount of output printed from a QM/MM job.\nRecommendation: : Use default unless running calculations with $charmm.    '''
        if value == "":
            if "QMMM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "QMMM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QMMM_PRINT"]=value.lower()


    def multipole_order(self,value="show"):
        '''\nName: MULTIPOLE_ORDER\nType: INTEGER\nDefault: 2\nOptions: 0:20:4:1\nDescription: Determines the order to which the multipole expansion of the solute charge density is carried out.\nRecommendation: : Use default unless higher (or lower) precision is desired.    '''
        if value == "":
            if "MULTIPOLE_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["MULTIPOLE_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "MULTIPOLE_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["MULTIPOLE_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MULTIPOLE_ORDER"]=value.lower()


    def plots_grid(self,value="show"):
        '''\nName: PLOTS_GRID\nType: STRING\nDefault: 0\nOptions: None:User-defined//-1:Nuclear positions//0:Read from file//1\nDescription: Controls evaluation of the electrostatic potential on a grid of points. If enabled, the output is in an ACSII file, plot.esp, in the format x, y, z, esp for each point.\n    '''
        if value == "":
            if "PLOTS_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["PLOTS_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "PLOTS_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["PLOTS_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PLOTS_GRID"]=value.lower()


    def mulliken(self,value="show"):
        '''\nName: MULLIKEN\nType: STRING\nDefault: 0\nOptions: None//0:Atomic//1:Shell//2\nDescription: \n    '''
        if value == "":
            if "MULLIKEN" in self.dict_of_keywords:
                del self.dict_of_keywords["MULLIKEN"]
                print "Keyword removed."
        elif value == "show":
            if "MULLIKEN" in self.dict_of_keywords:
                return self.dict_of_keywords["MULLIKEN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MULLIKEN"]=value.lower()


    def core_character_print(self,value="show"):
        '''\nName: CORE_CHARACTER_PRINT\nType: STRING\nDefault: 0\nOptions: None//0:Occupied MOs//1:MOs and AOs//2\nDescription: Determines the print level for the CORE_CHARACTER option.\nRecommendation: : Use default, unless you are uncertain about what the core character is.    '''
        if value == "":
            if "CORE_CHARACTER_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CORE_CHARACTER_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "CORE_CHARACTER_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CORE_CHARACTER_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CORE_CHARACTER_PRINT"]=value.lower()


    def print_distance_matrix(self,value="show"):
        '''\nName: PRINT_DISTANCE_MATRIX\nType: INTEGER\nDefault: 1\nOptions: 0:15\nDescription: Controls the printing of the inter-atomic distance matrix\nRecommendation: : Use default unless distances are required for large systems    '''
        if value == "":
            if "PRINT_DISTANCE_MATRIX" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_DISTANCE_MATRIX"]
                print "Keyword removed."
        elif value == "show":
            if "PRINT_DISTANCE_MATRIX" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_DISTANCE_MATRIX"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PRINT_DISTANCE_MATRIX"]=value.lower()


    def qui_print_orbitals(self,value="show"):
        '''\nName: QUI_PRINT_ORBITALS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_PRINT_ORBITALS" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_PRINT_ORBITALS"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_PRINT_ORBITALS" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_PRINT_ORBITALS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_PRINT_ORBITALS"]=value.lower()


    def print_orbitals(self,value="show"):
        '''\nName: PRINT_ORBITALS\nType: STRING\nDefault: 2\nOptions: 0:50:0:1\nDescription: \n    '''
        if value == "":
            if "PRINT_ORBITALS" in self.dict_of_keywords:
                del self.dict_of_keywords["PRINT_ORBITALS"]
                print "Keyword removed."
        elif value == "show":
            if "PRINT_ORBITALS" in self.dict_of_keywords:
                return self.dict_of_keywords["PRINT_ORBITALS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PRINT_ORBITALS"]=value.lower()


    def qmmm(self,value="show"):
        '''\nName: QMMM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Turns on the Q-Chem/CHARMM interface.\nRecommendation: : Use default unless running calculations with $charmm.    '''
        if value == "":
            if "QMMM" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM"]
                print "Keyword removed."
        elif value == "show":
            if "QMMM" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QMMM"]=value.lower()


    def rca_print(self,value="show"):
        '''\nName: RCA_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:3:0:1\nDescription: Controls the amount of output from a RCA SCF optimization.\n    '''
        if value == "":
            if "RCA_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "RCA_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RCA_PRINT"]=value.lower()


    def rpath_coordinates(self,value="show"):
        '''\nName: RPATH_COORDINATES\nType: STRING\nDefault: 0\nOptions: Mass-weighted//1:Z-matrix//2\nDescription: Determines which coordinate system to use in the IRC search.\nRecommendation: : Mass weighted coordinates are usually more effective.    '''
        if value == "":
            if "RPATH_COORDINATES" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_COORDINATES"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_COORDINATES" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_COORDINATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_COORDINATES"]=value.lower()


    def rpath_direction(self,value="show"):
        '''\nName: RPATH_DIRECTION\nType: STRING\nDefault: 0\nOptions: Positive//1:Negative//-1\nDescription: Determines the direction of the eigen mode to follow. This will not usually be known prior to the Hessian diagonalization and thereforr both directions will have to be considered.\n    '''
        if value == "":
            if "RPATH_DIRECTION" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_DIRECTION"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_DIRECTION" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_DIRECTION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_DIRECTION"]=value.lower()


    def rpath_max_stepsize(self,value="show"):
        '''\nName: RPATH_MAX_STEPSIZE\nType: INTEGER\nDefault: 2\nOptions: 0.001:0.500:0.150:0.001\nDescription: Specifies the maximum step size to be taken (in a.u.).\n    '''
        if value == "":
            if "RPATH_MAX_STEPSIZE" in self.dict_of_keywords:
                del self.dict_of_keywords["RPATH_MAX_STEPSIZE"]
                print "Keyword removed."
        elif value == "show":
            if "RPATH_MAX_STEPSIZE" in self.dict_of_keywords:
                return self.dict_of_keywords["RPATH_MAX_STEPSIZE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RPATH_MAX_STEPSIZE"]=value.lower()


    def geom_opt_scf_guess_always(self,value="show"):
        '''\nName: GEOM_OPT_SCF_GUESS_ALWAYS\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Switch to force the regeneration of a new initial guess for each series of SCF iterations (for use in geometry optimization).\nRecommendation: : Use default unless SCF convergence issues arise    '''
        if value == "":
            if "GEOM_OPT_SCF_GUESS_ALWAYS" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_SCF_GUESS_ALWAYS" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_SCF_GUESS_ALWAYS"]=value.lower()


    def scf_final_print(self,value="show"):
        '''\nName: SCF_FINAL_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:3 :0:1\nDescription: Controls level of output from SCF procedure to Q-Chem output file at the end of the SCF.\nRecommendation: : The break-down of energies is often useful (level 1).    '''
        if value == "":
            if "SCF_FINAL_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_FINAL_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_FINAL_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_FINAL_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_FINAL_PRINT"]=value.lower()


    def scf_guess_print(self,value="show"):
        '''\nName: SCF_GUESS_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:2:0:1\nDescription: Controls printing of guess MOs, Fock and density matrices.\n    '''
        if value == "":
            if "SCF_GUESS_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_GUESS_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_GUESS_PRINT"]=value.lower()


    def scf_print(self,value="show"):
        '''\nName: SCF_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:3:0:1\nDescription: Controls level of output from SCF procedure to Q-Chem output file.\nRecommendation: : Proceed with care; can result in extremely large output files at level 2 or higher. These levels are primarily for program debugging.    '''
        if value == "":
            if "SCF_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_PRINT"]=value.lower()


    def solute_radius(self,value="show"):
        '''\nName: SOLUTE_RADIUS\nType: INTEGER\nDefault: 2\nOptions: 0.0000:99.9999:0.0000:0.0001\nDescription: Sets the Onsager solvent model cavity radius.\nRecommendation: : Use equation (\ref{eq1000}).    '''
        if value == "":
            if "SOLUTE_RADIUS" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLUTE_RADIUS"]
                print "Keyword removed."
        elif value == "show":
            if "SOLUTE_RADIUS" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLUTE_RADIUS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SOLUTE_RADIUS"]=value.lower()


    def solvent_dielectric(self,value="show"):
        '''\nName: SOLVENT_DIELECTRIC\nType: INTEGER\nDefault: 2\nOptions: 0.0000:99.9999:0.0000:0.0001\nDescription: Sets the dielectric constant of the Onsager solvent continuum.\nRecommendation: : As per required solvent.    '''
        if value == "":
            if "SOLVENT_DIELECTRIC" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLVENT_DIELECTRIC"]
                print "Keyword removed."
        elif value == "show":
            if "SOLVENT_DIELECTRIC" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLVENT_DIELECTRIC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SOLVENT_DIELECTRIC"]=value.lower()


    def symmetry_ignore(self,value="show"):
        '''\nName: SYMMETRY_IGNORE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls whether or not Q-Chem determines the point group of the molecule.\nRecommendation: : Use default unless you do not want the molecule to be reoriented. Note that symmetry usage is disabled for RIMP2 jobs.    '''
        if value == "":
            if "SYMMETRY_IGNORE" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_IGNORE"]
                print "Keyword removed."
        elif value == "show":
            if "SYMMETRY_IGNORE" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_IGNORE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SYMMETRY_IGNORE"]=value.lower()


    def symmetry_integral(self,value="show"):
        '''\nName: SYMMETRY_INTEGRAL\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls the efficiency through the use of point group symmetry for calculating integrals.\nRecommendation: : Use default unless benchmarking. Note that symmetry usage is disabled for RIMP2 jobs.    '''
        if value == "":
            if "SYMMETRY_INTEGRAL" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_INTEGRAL"]
                print "Keyword removed."
        elif value == "show":
            if "SYMMETRY_INTEGRAL" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_INTEGRAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SYMMETRY_INTEGRAL"]=value.lower()


    def rca_switch_thresh(self,value="show"):
        '''\nName: RCA_SWITCH_THRESH\nType: INTEGER\nDefault: 2\nOptions: 0:12:3:1\nDescription: The threshold for switching between RCA and DIIS when the SCF algorithm is set to RCA_DIIS\n    '''
        if value == "":
            if "RCA_SWITCH_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["RCA_SWITCH_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "RCA_SWITCH_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["RCA_SWITCH_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["RCA_SWITCH_THRESH"]=value.lower()


    def qui_angular_grid(self,value="show"):
        '''\nName: QUI_ANGULAR_GRID\nType: STRING\nDefault: 1\nOptions: SG-0:SG-1:----------:6:18:26:38:50:74:86:110:146:170:194:230:266:302:350:434:590:770:974:1202:1454:1730:2030:2354:2702:3074:3470:3890:4334:4802:5294\nDescription: Specifies the quadrature grid to be used for evaluating the exchange-correlation component of the energy.  Either a standard grid should be selected, or a Lebedev grid with the corresponding number of points.
\nRecommendation: : Use the default unless convergence difficulties arise.  Larger grids are required for calculations involving derivatives and excited states.    '''
        if value == "":
            if "QUI_ANGULAR_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_ANGULAR_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_ANGULAR_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_ANGULAR_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_ANGULAR_GRID"]=value.lower()


    def qui_radial_grid(self,value="show"):
        '''\nName: QUI_RADIAL_GRID\nType: INTEGER\nDefault: 2\nOptions: 1:200:50:1\nDescription: Specifies the number of radial point for the exchange-correlation quadrature.\n    '''
        if value == "":
            if "QUI_RADIAL_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_RADIAL_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_RADIAL_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_RADIAL_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_RADIAL_GRID"]=value.lower()


    def isotopes(self,value="show"):
        '''\nName: ISOTOPES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Specifies if non-default masses are to be used in the frequency calculation.  If this option is selected the $isotopes section is read.\n    '''
        if value == "":
            if "ISOTOPES" in self.dict_of_keywords:
                del self.dict_of_keywords["ISOTOPES"]
                print "Keyword removed."
        elif value == "show":
            if "ISOTOPES" in self.dict_of_keywords:
                return self.dict_of_keywords["ISOTOPES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["ISOTOPES"]=value.lower()


    def memory_static(self,value="show"):
        '''\nName: MEMORY_STATIC\nType: INTEGER\nDefault: 2\nOptions: 1  :512:64:1\nDescription: Sets the memory (in megabytes) for individual program modules.\nRecommendation: : For direct and semi-direct MP2 calculations, this must exceed OVN + requirements for AO integral evaluation (32-160 Mb).    '''
        if value == "":
            if "MEMORY_STATIC" in self.dict_of_keywords:
                del self.dict_of_keywords["MEMORY_STATIC"]
                print "Keyword removed."
        elif value == "show":
            if "MEMORY_STATIC" in self.dict_of_keywords:
                return self.dict_of_keywords["MEMORY_STATIC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MEMORY_STATIC"]=value.lower()


    def memory_total(self,value="show"):
        '''\nName: MEMORY_TOTAL\nType: INTEGER\nDefault: 2\nOptions: 128:8000:2000:10\nDescription: Sets the total memory available to Q-Chem, in megabytes.\nRecommendation: : Use default, or set to the physical memory of your machine.    '''
        if value == "":
            if "MEMORY_TOTAL" in self.dict_of_keywords:
                del self.dict_of_keywords["MEMORY_TOTAL"]
                print "Keyword removed."
        elif value == "show":
            if "MEMORY_TOTAL" in self.dict_of_keywords:
                return self.dict_of_keywords["MEMORY_TOTAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MEMORY_TOTAL"]=value.lower()


    def aimd_method(self,value="show"):
        '''\nName: AIMD_METHOD\nType: STRING\nDefault: 0\nOptions: BOMD:Curvy\nDescription: Selects an ab initio molecular dynamics algorithm.\nRecommendation: : Born-oppenheimer MD (BOMD) yields exact classical molecular dynamics, provided that the energy is tolerably conserved. Curvy-steps extended Lagrangian MD (Curvy) is an approximation to exact classical dynamics whose validity should be tested for the properties of interest.     '''
        if value == "":
            if "AIMD_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["AIMD_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "AIMD_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["AIMD_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AIMD_METHOD"]=value.lower()


    def cc_dthreshold(self,value="show"):
        '''\nName: CC_DTHRESHOLD\nType: INTEGER\nDefault: 2\nOptions: 0.000001:1.00000:0.00001:0.00001\nDescription: Specifies threshold for including a new expansion vector in the iterative Davidson diagonalization. Their norm must be above this threshold. \nRecommendation: : Use default unless converge problems are encountered. Should normally be set to the same values as CC_DCONVERGENCE, if convergence problems arise try setting to a value less than CC_DCONVERGENCE.    '''
        if value == "":
            if "CC_DTHRESHOLD" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DTHRESHOLD"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DTHRESHOLD" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DTHRESHOLD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DTHRESHOLD"]=value.lower()


    def cc_state_derivative(self,value="show"):
        '''\nName: CC_STATE_DERIVATIVE\nType: INTEGER\nDefault: 0\nOptions: -1 \nDescription: Selects which EOM or CIS(D) state is to be considered for optimization or property calculations.\n    '''
        if value == "":
            if "CC_STATE_DERIVATIVE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_DERIVATIVE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_STATE_DERIVATIVE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_DERIVATIVE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_STATE_DERIVATIVE"]=value.lower()


    def cc_reference_symmetry(self,value="show"):
        '''\nName: CC_REFERENCE_SYMMETRY\nType: INTEGER\nDefault: 0\nOptions: -1 \nDescription: Together with CC_STATE_DERIV, selects which EOM state is to be considered for optimization or property calculations. When transition properties are requested, the transition properties will be calculated between this state and all other EOM states.\n    '''
        if value == "":
            if "CC_REFERENCE_SYMMETRY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REFERENCE_SYMMETRY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REFERENCE_SYMMETRY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REFERENCE_SYMMETRY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REFERENCE_SYMMETRY"]=value.lower()


    def gui(self,value="show"):
        '''\nName: GUI\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Controls the output of auxiliary information for third party packages.\nRecommendation: : Use default unless the additional information is required. Please note that any existing Test.FChk file will be overwritten.    '''
        if value == "":
            if "GUI" in self.dict_of_keywords:
                del self.dict_of_keywords["GUI"]
                print "Keyword removed."
        elif value == "show":
            if "GUI" in self.dict_of_keywords:
                return self.dict_of_keywords["GUI"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GUI"]=value.lower()


    def basis(self,value="show"):
        '''\nName: BASIS\nType: STRING\nDefault: 5\nOptions: STO-3G:STO-6G:----------:3-21G:4-31G:6-31G:6-31G*:6-31+G*:6-31G**:6-31++G**:6-311G:6-311G*:6-311+G*:6-311G**:6-311++G**:6-311++G(3df,3pd):----------:cc-pVDZ:cc-pVTZ:cc-pVQZ:cc-pcVDZ:cc-pcVTZ:cc-pcVQZ:aug-cc-pVDZ:aug-cc-pVTZ:aug-cc-pVQZ:aug-cc-pcVDZ:aug-cc-pcVTZ:aug-cc-pcVQZ:----------:G3Large:G3MP2Large:----------:CRENBL:CRENBS:HWMB:HWVDZ:LACVP:LANL2DZ:SBKJC:SRLC:SRSC:----------:User-defined//gen:Mixed\nDescription: Specifies the basis sets to be used.\nRecommendation: : Consult literature and reviews to aid your selection.    '''
        if value == "":
            if "BASIS" in self.dict_of_keywords:
                del self.dict_of_keywords["BASIS"]
                print "Keyword removed."
        elif value == "show":
            if "BASIS" in self.dict_of_keywords:
                return self.dict_of_keywords["BASIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["BASIS"]=value.lower()


    def cc_eom_two_particle_properties(self,value="show"):
        '''\nName: CC_EOM_TWO_PARTICLE_PROPERTIES\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Request for calculation of non-relaxed two-particle EOM-CCSD target state properties.  The two-part properties currently include<>. The one-particle properties will also be calculated since the additional cost of these is small in comparison.  The vairable CC_EXSTATES_PROP must also be set.
\nRecommendation:  Two-particle properties are extremely computationally expensive since the require calculation and use of the two-particle density matrix (the cost being about the same as the cost of an analytic gradient calculation for each state.    '''
        if value == "":
            if "CC_EOM_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TWO_PARTICLE_PROPERTIES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_TWO_PARTICLE_PROPERTIES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TWO_PARTICLE_PROPERTIES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_TWO_PARTICLE_PROPERTIES"]=value.lower()


    def qui_section_opt(self,value="show"):
        '''\nName: QUI_SECTION_OPT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Adds the $opt section for specifying constraints in the geometry optimization\n    '''
        if value == "":
            if "QUI_SECTION_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SECTION_OPT"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SECTION_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SECTION_OPT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SECTION_OPT"]=value.lower()


    def dft_d(self,value="show"):
        '''\nName: DFT_D\nType: STRING\nDefault: 0\nOptions: None//0:Grimme//EMPIRICAL_GRIMME:Chai Head-Gordon//EMPIRICAL_CHG\nDescription: Specifies what dispersion correction to use within a DFT calculation.\n    '''
        if value == "":
            if "DFT_D" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D"]
                print "Keyword removed."
        elif value == "show":
            if "DFT_D" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DFT_D"]=value.lower()


    def dft_d_a(self,value="show"):
        '''\nName: DFT_D_A\nType: INTEGER\nDefault: 2\nOptions: 0.01:100.00:6.00:0.01\nDescription: Controls the strength of the dispersion corrections in the Chai-Head-Gordon scheme.   The default value should be apprpriate for most systems.
\n    '''
        if value == "":
            if "DFT_D_A" in self.dict_of_keywords:
                del self.dict_of_keywords["DFT_D_A"]
                print "Keyword removed."
        elif value == "show":
            if "DFT_D_A" in self.dict_of_keywords:
                return self.dict_of_keywords["DFT_D_A"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DFT_D_A"]=value.lower()


    def scf_guess_mix(self,value="show"):
        '''\nName: SCF_GUESS_MIX\nType: INTEGER\nDefault: 2\nOptions: 0:100:0:10\nDescription: Controls mixing of LUMO and HOMO to break symmetry in the initial guess. For unrestricted jobs, the mixing is performed only for the alpha orbitals.\nRecommendation: : When performing unrestricted calculations on molecules with an even number of electrons, it is often necessary to break alpha-beta symmetry in the initial guess with this option, or by specifying input for $occupied.    '''
        if value == "":
            if "SCF_GUESS_MIX" in self.dict_of_keywords:
                del self.dict_of_keywords["SCF_GUESS_MIX"]
                print "Keyword removed."
        elif value == "show":
            if "SCF_GUESS_MIX" in self.dict_of_keywords:
                return self.dict_of_keywords["SCF_GUESS_MIX"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SCF_GUESS_MIX"]=value.lower()


    def cc_print(self,value="show"):
        '''\nName: CC_PRINT\nType: INTEGER\nDefault: 0\nOptions: 1\nDescription: Controls the output from post-MP2 coupled-cluster module of QChem\nRecommendation: : Increase if you need more output and don't like trees    '''
        if value == "":
            if "CC_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "CC_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_PRINT"]=value.lower()


    def moprop_save_last_gpx(self,value="show"):
        '''\nName: MOPROP_SAVE_LAST_GPX\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Save last G[P]x when calculating dynamic polarizabilities in order to call mopropman in a second run with MOPROP = 102.\n    '''
        if value == "":
            if "MOPROP_SAVE_LAST_GPX" in self.dict_of_keywords:
                del self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"]
                print "Keyword removed."
        elif value == "show":
            if "MOPROP_SAVE_LAST_GPX" in self.dict_of_keywords:
                return self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MOPROP_SAVE_LAST_GPX"]=value.lower()


    def qui_title(self,value="show"):
        '''\nName: QUI_TITLE\nType: STRING\nDefault: 0\nOptions:  \nDescription: Sets the lable for this section of the input file.\n    '''
        if value == "":
            if "QUI_TITLE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_TITLE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_TITLE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_TITLE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_TITLE"]=value.lower()


    def smx_solvent(self,value="show"):
        '''\nName: SMX_SOLVENT\nType: STRING\nDefault: 174\nOptions: 111trichloroethane:112trichloroethane:11dichloroethane:124trimethylbenzene:14dioxane:1bromo2methylpropane:1bromopentane:1bromopropane:1butanol:1chloropentane:1chloropropane:1decanol:1fluorooctane:1heptanol:1hexanol:1hexene:1hexyne:1iodobutane:1iodopentene:1iodopropane:1nitropropane:1nonanol:1octanol:1pentanol:1pentene:1pentyne:1propanol:222trifluoroethanol:224trimethylpentane:24dimethylpentane:24dimethylpyridine:26dimethylpyridine:2bromopropane:2chlorobutane:2heptanone:2hexanone:2methylpentane:2methylpyridine:2nitropropane:2octanone:2pentanone:2propanol:2propen1ol:3methylpyridine:3pentanone:4heptanone:4methyl2pentanone:4methylpyridine:5nonanone:aceticacid:acetone:acetonitrile:aniline:anisole:benzaldehyde:benzene:benzonitrile:benzylalcohol:bromobenzene:bromoethane:bromooctane:butanal:butanoicacid:butanone:butanonitrile:butylethanoate:butylamine:butylbenzene:carbondisulfide:carbontet:chlorobenzene:chlorotoluene:cis12dimethylcyclohexane:decalin:cyclohexane:cyclohexanone:cyclopentane:cyclopentanol:cyclopentanone:decane:dibromomethane:dibutylether:dichloromethane:diethylether:diethylsulfide:diethylamine:diiodomethane:dimethyldisulfide:dimethylacetamide:dimethylformamide:dimethylpyridine:dmso:dipropylamine:dodecane:E12dichloroethene:E2pentene:ethanethiol:ethanol:ethylethanoate:ethylmethanoate:ethylphenylether:ethylbenzene:ethyleneglycol:fluorobenzene:formamide:formicacid:hexadecyliodide:hexanoic:iodobenzene:iodoethane:iodomethane:isobutanol:isopropylether:isopropylbenzene:isopropyltoluene:mcresol:mesitylene:methanol:methylbenzoate:methylethanoate:methylmethanoate:methylphenylketone:methylpropanoate:methylbutanoate:methylcyclohexane:methylformamide:methylformamide:heptane:hexadecane:hexane:nitrobenzene:nitroethane:nitromethane:methylaniline:nonane:octane:pentane:ochlorotoluene:ocresol:odichlorobenzene:onitrotoluene:oxylene:pentadecane:pentanal:pentanoicacid:pentylethanoate:pentylamine:perfluorobenzene:phenyletherphenylether:propanal:propanoicacid:propanonitrile:propylethanoate:propylamine:pxylene:pyridine:pyrrolidine:secbutanol:tbutanol:tbutylbenzene:tetrachloroethene:thf:tetrahyrothiophenedioxide:tetralin:thiophene:thiophenol:toluene:transdecalin:tribromomethane:tributylphosphate:trichloroethene:trichloromethane:triethylamine:undecane:water:Z12dichloroethene\nDescription: Specifies which solvent to use in the SM8 model.\n    '''
        if value == "":
            if "SMX_SOLVENT" in self.dict_of_keywords:
                del self.dict_of_keywords["SMX_SOLVENT"]
                print "Keyword removed."
        elif value == "show":
            if "SMX_SOLVENT" in self.dict_of_keywords:
                return self.dict_of_keywords["SMX_SOLVENT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SMX_SOLVENT"]=value.lower()


    def smx_solvation(self,value="show"):
        '''\nName: SMX_SOLVATION\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Sets whether or not to use the SM8 solvation model.\n    '''
        if value == "":
            if "SMX_SOLVATION" in self.dict_of_keywords:
                del self.dict_of_keywords["SMX_SOLVATION"]
                print "Keyword removed."
        elif value == "show":
            if "SMX_SOLVATION" in self.dict_of_keywords:
                return self.dict_of_keywords["SMX_SOLVATION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SMX_SOLVATION"]=value.lower()


    def link_atom_projection(self,value="show"):
        '''\nName: LINK_ATOM_PROJECTION\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Controls whether to perform a link-atom projection, which is necessary in a full QM/MM hessian evaluation on a system with link atoms.\n    '''
        if value == "":
            if "LINK_ATOM_PROJECTION" in self.dict_of_keywords:
                del self.dict_of_keywords["LINK_ATOM_PROJECTION"]
                print "Keyword removed."
        elif value == "show":
            if "LINK_ATOM_PROJECTION" in self.dict_of_keywords:
                return self.dict_of_keywords["LINK_ATOM_PROJECTION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["LINK_ATOM_PROJECTION"]=value.lower()


    def qmmm_full_hessian(self,value="show"):
        '''\nName: QMMM_FULL_HESSIAN\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QMMM_FULL_HESSIAN" in self.dict_of_keywords:
                del self.dict_of_keywords["QMMM_FULL_HESSIAN"]
                print "Keyword removed."
        elif value == "show":
            if "QMMM_FULL_HESSIAN" in self.dict_of_keywords:
                return self.dict_of_keywords["QMMM_FULL_HESSIAN"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QMMM_FULL_HESSIAN"]=value.lower()


    def gaussian_blur(self,value="show"):
        '''\nName: GAUSSIAN_BLUR\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Enables the use of Gaussian-delocalized external charges in a QM/MM calculation.  If set to FALSE, then regular point charges are used.\n    '''
        if value == "":
            if "GAUSSIAN_BLUR" in self.dict_of_keywords:
                del self.dict_of_keywords["GAUSSIAN_BLUR"]
                print "Keyword removed."
        elif value == "show":
            if "GAUSSIAN_BLUR" in self.dict_of_keywords:
                return self.dict_of_keywords["GAUSSIAN_BLUR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GAUSSIAN_BLUR"]=value.lower()


    def hess_proj_trm(self,value="show"):
        '''\nName: HESS_PROJ_TRM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Selects whether or not to project out the rotational and translational degrees of freedom in a frequency calculation.\n    '''
        if value == "":
            if "HESS_PROJ_TRM" in self.dict_of_keywords:
                del self.dict_of_keywords["HESS_PROJ_TRM"]
                print "Keyword removed."
        elif value == "show":
            if "HESS_PROJ_TRM" in self.dict_of_keywords:
                return self.dict_of_keywords["HESS_PROJ_TRM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["HESS_PROJ_TRM"]=value.lower()


    def jobtype(self,value="show"):
        '''\nName: JOBTYPE\nType: STRING\nDefault: 0\nOptions: Energy//SP:Forces//Force:Geometry//Optimization:Frequencies//Frequency:Reaction Path//RPath:Transition State//TS:Chemical Shifts//NMR:Ab Initio MD//AIMD:Properties//SP\nDescription: Specifies the type of calculation to run. \n    '''
        if value == "":
            if "JOBTYPE" in self.dict_of_keywords:
                del self.dict_of_keywords["JOBTYPE"]
                print "Keyword removed."
        elif value == "show":
            if "JOBTYPE" in self.dict_of_keywords:
                return self.dict_of_keywords["JOBTYPE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["JOBTYPE"]=value.lower()


    def geom_opt_iproj(self,value="show"):
        '''\nName: GEOM_OPT_IPROJ\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Allows the molecule to reorient during a geometry optimization.  Turn this option off if using external charges.\n    '''
        if value == "":
            if "GEOM_OPT_IPROJ" in self.dict_of_keywords:
                del self.dict_of_keywords["GEOM_OPT_IPROJ"]
                print "Keyword removed."
        elif value == "show":
            if "GEOM_OPT_IPROJ" in self.dict_of_keywords:
                return self.dict_of_keywords["GEOM_OPT_IPROJ"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["GEOM_OPT_IPROJ"]=value.lower()


    def qui_qchem_executable(self,value="show"):
        '''\nName: QUI_QCHEM_EXECUTABLE\nType: STRING\nDefault: 1\nOptions: qcprog.exe:qchem_s.exe\nDescription: \n    '''
        if value == "":
            if "QUI_QCHEM_EXECUTABLE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_QCHEM_EXECUTABLE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_QCHEM_EXECUTABLE"]=value.lower()


    def qui_avogadro_visualize_file(self,value="show"):
        '''\nName: QUI_AVOGADRO_VISUALIZE_FILE\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Determines which file is passed as an argument to Avogadro.  If true, then the .out file is passed, if false the .Fchk file is passed.\n    '''
        if value == "":
            if "QUI_AVOGADRO_VISUALIZE_FILE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_AVOGADRO_VISUALIZE_FILE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_AVOGADRO_VISUALIZE_FILE"]=value.lower()


    def qui_windows_directory(self,value="show"):
        '''\nName: QUI_WINDOWS_DIRECTORY\nType: STRING\nDefault: 0\nOptions: /Windows/System32\nDescription: Sets the directory for tskill or taskkill\n    '''
        if value == "":
            if "QUI_WINDOWS_DIRECTORY" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_WINDOWS_DIRECTORY" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_WINDOWS_DIRECTORY"]=value.lower()


    def qui_windows_kill_command(self,value="show"):
        '''\nName: QUI_WINDOWS_KILL_COMMAND\nType: STRING\nDefault: 0\nOptions: tskill qchem_s:taskkill /IM qchem_s.exe /F\nDescription: The command required to kill qchem jobs, only used on Windows\n    '''
        if value == "":
            if "QUI_WINDOWS_KILL_COMMAND" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_WINDOWS_KILL_COMMAND" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_WINDOWS_KILL_COMMAND"]=value.lower()


    def xc_grid(self,value="show"):
        '''\nName: XC_GRID\nType: STRING\nDefault: 0\nOptions: SG-0//0:SG-1//1:----------:6:18:26:38:50:74:86:110:146:170:194:230:266:302:350:434:590:770:974:1202:1454:1730:2030:2354:2702:3074:3470:3890:4334:4802:5294\nDescription: Specifies the quadrature grid to be used for evaluating the exchange-correlation component of the energy.  Either a standard grid should be selected, or a Lebedev grid with the corresponding number of points.
\nRecommendation: : Use the default unless convergence difficulties arise.  Larger grids are required for calculations involving derivatives and excited states.    '''
        if value == "":
            if "XC_GRID" in self.dict_of_keywords:
                del self.dict_of_keywords["XC_GRID"]
                print "Keyword removed."
        elif value == "show":
            if "XC_GRID" in self.dict_of_keywords:
                return self.dict_of_keywords["XC_GRID"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XC_GRID"]=value.lower()


    def pdb_print(self,value="show"):
        '''\nName: PDB_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:2:0:1\nDescription: Prints final coordinates at the end of the output file using the PDB format.\n    '''
        if value == "":
            if "PDB_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["PDB_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "PDB_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["PDB_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PDB_PRINT"]=value.lower()


    def aaaa(self,value="show"):
        '''\nName: AAAA\nType: STRING\nDefault: 1\nOptions: 0:New option:1:1:1\nDescription: \n    '''
        if value == "":
            if "AAAA" in self.dict_of_keywords:
                del self.dict_of_keywords["AAAA"]
                print "Keyword removed."
        elif value == "show":
            if "AAAA" in self.dict_of_keywords:
                return self.dict_of_keywords["AAAA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["AAAA"]=value.lower()


    def threads(self,value="show"):
        '''\nName: THREADS\nType: INTEGER\nDefault: 2\nOptions: 1:1024:1:1\nDescription: Number of threads in shared memory parallel calculations.\n    '''
        if value == "":
            if "THREADS" in self.dict_of_keywords:
                del self.dict_of_keywords["THREADS"]
                print "Keyword removed."
        elif value == "show":
            if "THREADS" in self.dict_of_keywords:
                return self.dict_of_keywords["THREADS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["THREADS"]=value.lower()


    def cc_max_iter(self,value="show"):
        '''\nName: CC_MAX_ITER\nType: INTEGER\nDefault: 2\nOptions: 1:1000:200:1\nDescription: Maximum number of iterations to optimize the coupled-cluster energy. \n    '''
        if value == "":
            if "CC_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MAX_ITER"]
                print "Keyword removed."
        elif value == "show":
            if "CC_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MAX_ITER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_MAX_ITER"]=value.lower()


    def cc_memory(self,value="show"):
        '''\nName: CC_MEMORY\nType: INTEGER\nDefault: 2\nOptions: 192:1000000:1500:10\nDescription: Specifies the maximum size, in Mb, of the buffers for in-core storage of block-tensors in CCMAN and CCMAN2.\n    '''
        if value == "":
            if "CC_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_MEMORY"]
                print "Keyword removed."
        elif value == "show":
            if "CC_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_MEMORY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_MEMORY"]=value.lower()


    def mem_total(self,value="show"):
        '''\nName: MEM_TOTAL\nType: INTEGER\nDefault: 2\nOptions: 200:1000000:2000:10\nDescription: Sets the total memory available to Q-Chem, in megabytes.\n    '''
        if value == "":
            if "MEM_TOTAL" in self.dict_of_keywords:
                del self.dict_of_keywords["MEM_TOTAL"]
                print "Keyword removed."
        elif value == "show":
            if "MEM_TOTAL" in self.dict_of_keywords:
                return self.dict_of_keywords["MEM_TOTAL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MEM_TOTAL"]=value.lower()


    def mem_static(self,value="show"):
        '''\nName: MEM_STATIC\nType: INTEGER\nDefault: 2\nOptions: 32:1000000:240:10\nDescription: Sets the memory (in megabytes) for individual fortran program modules.\nRecommendation: : For direct and semi-direct MP2 calculations, this must exceed OVN + requirements for AO integral evaluation (32-160 Mb).    '''
        if value == "":
            if "MEM_STATIC" in self.dict_of_keywords:
                del self.dict_of_keywords["MEM_STATIC"]
                print "Keyword removed."
        elif value == "show":
            if "MEM_STATIC" in self.dict_of_keywords:
                return self.dict_of_keywords["MEM_STATIC"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["MEM_STATIC"]=value.lower()


    def cc_incl_core_corr(self,value="show"):
        '''\nName: CC_INCL_CORE_CORR\nType: LOGICAL\nDefault: 1\nOptions: FALSE:TRUE\nDescription: Whether to include the correlation contribution from frozen core orbitals in non iterative (2) corrections, such as OD(2) and CCSD(2).\nRecommendation: : Use default unless no core-valence or core correlation is desired (e.g., for comparison with other methods or because the basis used cannot describe core correlation).    '''
        if value == "":
            if "CC_INCL_CORE_CORR" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_INCL_CORE_CORR"]
                print "Keyword removed."
        elif value == "show":
            if "CC_INCL_CORE_CORR" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_INCL_CORE_CORR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_INCL_CORE_CORR"]=value.lower()


    def cc_ref_prop(self,value="show"):
        '''\nName: CC_REF_PROP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the non-relaxed (expectation value) one-particle CCSD properties will be calculated. The properties currently include permanent dipole moment, the second moments , , and  of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.\nRecommendation: : Additional equations need to be solved (lambda CCSD equations) for properties with the cost approximately the same as CCSD equations. Use default if you do not need properties. The cost of the properties calculation itself is low. The CCSD one-particle density can be analyzed with NBO package by specifying NBO=TRUE, CC_PROP=TRUE and JOBTYPE=FORCE.    '''
        if value == "":
            if "CC_REF_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REF_PROP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REF_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REF_PROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REF_PROP"]=value.lower()


    def cc_nguess_doubles(self,value="show"):
        '''\nName: CC_NGUESS_DOUBLES\nType: INTEGER\nDefault: 0\nOptions: 0\nDescription: Specifies number of excited state guess vectors which are double excitations. \nRecommendation: : This should be set to the expected number of doubly excited states (see also EOM_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "CC_NGUESS_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_NGUESS_DOUBLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_NGUESS_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_NGUESS_DOUBLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_NGUESS_DOUBLES"]=value.lower()


    def ccman2(self,value="show"):
        '''\nName: CCMAN2\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "CCMAN2" in self.dict_of_keywords:
                del self.dict_of_keywords["CCMAN2"]
                print "Keyword removed."
        elif value == "show":
            if "CCMAN2" in self.dict_of_keywords:
                return self.dict_of_keywords["CCMAN2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CCMAN2"]=value.lower()


    def cc_t_conv(self,value="show"):
        '''\nName: CC_T_CONV\nType: INTEGER\nDefault: 2\nOptions: 0:12:8:1\nDescription: \n    '''
        if value == "":
            if "CC_T_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_T_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "CC_T_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_T_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_T_CONV"]=value.lower()


    def cc_ref_prop_te(self,value="show"):
        '''\nName: CC_REF_PROP_TE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Request for calculation of non-relaxed two-particle CCSD properties. The two-particle properties currently include . The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation is
inferior compared to the cost of . The variable CC_REF_PROP must be also set to TRUE.\n    '''
        if value == "":
            if "CC_REF_PROP_TE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_REF_PROP_TE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_REF_PROP_TE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_REF_PROP_TE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_REF_PROP_TE"]=value.lower()


    def eom_corr(self,value="show"):
        '''\nName: EOM_CORR\nType: STRING\nDefault: 0\nOptions: CIS:CIS(D):SDT:DT:SD:SD(dT):SD(fT):SD(sT)\nDescription: Specifies the correlation level\n    '''
        if value == "":
            if "EOM_CORR" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_CORR"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_CORR" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_CORR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_CORR"]=value.lower()


    def eom_sf_states(self,value="show"):
        '''\nName: EOM_SF_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of spin-flip target states roots to find.
[i, j, k...] Find i SF states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_SF_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_SF_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_SF_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_SF_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_SF_STATES"]=value.lower()


    def eom_dsf_states(self,value="show"):
        '''\nName: EOM_DSF_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of doubly spin-flip target states roots to find.
[i, j, k...] Find i doubly SF states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_DSF_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DSF_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DSF_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DSF_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DSF_STATES"]=value.lower()


    def eom_ea_states(self,value="show"):
        '''\nName: EOM_EA_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of attached target states roots to find. By default, alpha electron will be attached (see EOM_EA_ALPHA).
[i; j; k...] Find i EA states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EA_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EA_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EA_STATES"]=value.lower()


    def eom_ee_states(self,value="show"):
        '''\nName: EOM_EE_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of excited state roots to find. For closed-shell reference, defaults into EOM EE SINGLETS. For open-shell references, speccies all low-lying states.
[i, j, k...] Find i excited states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EE_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EE_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EE_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EE_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EE_STATES"]=value.lower()


    def eom_dip_states(self,value="show"):
        '''\nName: EOM_DIP_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of DIP roots to find. For closed-shell reference, defaults into EOM_DIP_SINGLETS. For open-shell references, speccies all low-lying states.
[i, j, k...] Find i DIP states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_DIP_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DIP_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DIP_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DIP_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DIP_STATES"]=value.lower()


    def eom_dip_singlets(self,value="show"):
        '''\nName: EOM_DIP_SINGLETS\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of singlet DIP roots to find. Works only for closed-shell references.
[i, j, k...] Find i DIP singlet states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_DIP_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DIP_SINGLETS"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DIP_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DIP_SINGLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DIP_SINGLETS"]=value.lower()


    def eom_ee_singlets(self,value="show"):
        '''\nName: EOM_EE_SINGLETS\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of singlet excited state roots to find. Works only for closed-shell references.
[i, j, k...] Find i excited states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EE_SINGLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EE_SINGLETS"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EE_SINGLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EE_SINGLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EE_SINGLETS"]=value.lower()


    def eom_ee_triplets(self,value="show"):
        '''\nName: EOM_EE_TRIPLETS\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of triplet excited state roots to find. Works only for closed-shell references.
[i, j, k...] Find i excited states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EE_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EE_TRIPLETS"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EE_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EE_TRIPLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EE_TRIPLETS"]=value.lower()


    def eom_dip_triplets(self,value="show"):
        '''\nName: EOM_DIP_TRIPLETS\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of triplet DIP roots to find. Works only for closed-shell references.
[i, j, k...] Find i DIP triplet states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_DIP_TRIPLETS" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DIP_TRIPLETS"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DIP_TRIPLETS" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DIP_TRIPLETS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DIP_TRIPLETS"]=value.lower()


    def eom_ip_states(self,value="show"):
        '''\nName: EOM_IP_STATES\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of ionized target states roots to find. By default, B electron will be removed (see EOM_IP_BETA).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_IP_STATES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_STATES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_IP_STATES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_STATES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_IP_STATES"]=value.lower()


    def eom_ip_alpha(self,value="show"):
        '''\nName: EOM_IP_ALPHA\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of ionized target states derived by removing alpha electron (Ms = -1/2).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_IP_ALPHA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_ALPHA"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_IP_ALPHA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_ALPHA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_IP_ALPHA"]=value.lower()


    def eom_ip_beta(self,value="show"):
        '''\nName: EOM_IP_BETA\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of ionized target states derived by removing beta electron (Ms = 1/2, default for EOM-IP).
[i, j, k...] Find i inonized states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_IP_BETA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_IP_BETA"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_IP_BETA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_IP_BETA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_IP_BETA"]=value.lower()


    def eom_ea_beta(self,value="show"):
        '''\nName: EOM_EA_BETA\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of attached target states derived by attaching beta electron (Ms = -1/2).
[i, j, k...] Find i EA states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EA_BETA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_BETA"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EA_BETA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_BETA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EA_BETA"]=value.lower()


    def eom_ea_alpha(self,value="show"):
        '''\nName: EOM_EA_ALPHA\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Sets the number of attached target states derived by attaching alpha electron (Ms = 1/2).
[i, j, k...] Find i EA states in the first irrep, j states in the second irrep etc.\n    '''
        if value == "":
            if "EOM_EA_ALPHA" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_EA_ALPHA"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_EA_ALPHA" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_EA_ALPHA"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_EA_ALPHA"]=value.lower()


    def eom_davidson_convergence(self,value="show"):
        '''\nName: EOM_DAVIDSON_CONVERGENCE\nType: INTEGER\nDefault: 2\nOptions: 0:12:5:1\nDescription: Convergence criterion for the RMS residuals of excited state vectors.
n Corresponding to 10^-n convergence criterion
Use default. Should normally be set to the same value as
EOM DAVIDSON THRESHOLD.\n    '''
        if value == "":
            if "EOM_DAVIDSON_CONVERGENCE" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DAVIDSON_CONVERGENCE" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DAVIDSON_CONVERGENCE"]=value.lower()


    def eom_davidson_maxvector(self,value="show"):
        '''\nName: EOM_DAVIDSON_MAXVECTOR\nType: INTEGER\nDefault: 2\nOptions: 0:200:60:1\nDescription: Species maximum number of vectors in the subspace for the Davidson diagonalization.
n Up to n vectors per root before the subspace is reset
RECOMMENDATION:
Larger values increase disk storage but accelerate and stabilize convergence.\n    '''
        if value == "":
            if "EOM_DAVIDSON_MAXVECTOR" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_MAXVECTOR"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DAVIDSON_MAXVECTOR" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_MAXVECTOR"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DAVIDSON_MAXVECTOR"]=value.lower()


    def eom_davidson_max_iter(self,value="show"):
        '''\nName: EOM_DAVIDSON_MAX_ITER\nType: INTEGER\nDefault: 2\nOptions: 0:100:30:1\nDescription: Maximum number of iteration allowed for Davidson diagonalization procedure.
n User-defined number of iterations\n    '''
        if value == "":
            if "EOM_DAVIDSON_MAX_ITER" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DAVIDSON_MAX_ITER" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DAVIDSON_MAX_ITER"]=value.lower()


    def eom_ngues_doubles(self,value="show"):
        '''\nName: EOM_NGUES_DOUBLES\nType: INTEGER\nDefault: 2\nOptions: 0:1000:0:1\nDescription: Specifies number of excited state guess vectors which are double excitations. 
Options: n Include n guess vectors that are double excitations\nRecommendation: : This should be set to the expected number of doubly excited states (see also CC_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "EOM_NGUES_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUES_DOUBLES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_NGUES_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUES_DOUBLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_NGUES_DOUBLES"]=value.lower()


    def eom_nguess_doubles(self,value="show"):
        '''\nName: EOM_NGUESS_DOUBLES\nType: INTEGER\nDefault: 2\nOptions: 0:1000:0:1\nDescription: Specifies number of excited state guess vectors which are double excitations. 
Options: n Include n guess vectors that are double excitations\nRecommendation: : This should be set to the expected number of doubly excited states (see also CC_PRECONV_DOUBLES), otherwise they may not be found.    '''
        if value == "":
            if "EOM_NGUESS_DOUBLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUESS_DOUBLES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_NGUESS_DOUBLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUESS_DOUBLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_NGUESS_DOUBLES"]=value.lower()


    def eom_nguess_singles(self,value="show"):
        '''\nName: EOM_NGUESS_SINGLES\nType: INTEGER\nDefault: 2\nOptions: 0:1000:0:1\nDescription: Specifies number of excited state guess vectors which are single excitations. 
Options: n Include n guess vectors that are single excitations\nRecommendation: : Should be greater or equal than the number of excited states requested.    '''
        if value == "":
            if "EOM_NGUESS_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_NGUESS_SINGLES"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_NGUESS_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_NGUESS_SINGLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_NGUESS_SINGLES"]=value.lower()


    def cc_eom_state_to_opt(self,value="show"):
        '''\nName: CC_EOM_STATE_TO_OPT\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Specifies which state to optimize.
[i,j] optimize the jth state of the ith irrep.\n    '''
        if value == "":
            if "CC_EOM_STATE_TO_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_STATE_TO_OPT"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_STATE_TO_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_STATE_TO_OPT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_STATE_TO_OPT"]=value.lower()


    def cc_eom_trans_prop(self,value="show"):
        '''\nName: CC_EOM_TRANS_PROP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_REFSYM and CC_STATE_DERIV must be specified for this state.\nRecommendation: : Additional equations (for the left EOM-CCSD eigenvectors plus lambda CCSD equations in case if transition properties between the CCSD reference and EOM-CCSD target states are requested) need to be solved for transition properties, approximately doubling the computational cost. The cost of the transition properties calculation itself is low.    '''
        if value == "":
            if "CC_EOM_TRANS_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_TRANS_PROP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_TRANS_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_TRANS_PROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_TRANS_PROP"]=value.lower()


    def cc_eom_prop(self,value="show"):
        '''\nName: CC_EOM_PROP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the non-relaxed (expectation value) one-particle EOM-CCSD target state properties will be calculated. The properties currently include permanent dipole moment, the second moments 2>, 2>, and 2> of electron density, and the total 2> = 2> +2> +2> (in atomic units). Incompatible with JOBTYPE=FORCE, OPT, FREQ.\nRecommendation: : Additional equations (EOM-CCSD equations for the left eigenvectors) need to be solved for properties, approximately doubling the cost of calculation for each irrep. Sometimes the equations for left and right eigenvectors converge to different sets of target states. In this case, the simultaneous iterations of left and right vectors will diverge, and the properties for several or all the target states may be incorrect! The problem can be solved by varying the number of requested states, specified with CC_NLOWSPIN and CC_NHIGHSPIN, or the number of guess vectors (CC_NGUESS_SINGLES). The cost of the one-particle properties calculation itself is low. The one-particle density of an EOM-CCSD target state can be analyzed with NBO package by specifying the state with CC_REFSYM and CC_STATE_DERIV and requesting NBO=TRUE and CC_EXSTATES_PROP=TRUE.    '''
        if value == "":
            if "CC_EOM_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_PROP"]=value.lower()


    def cc_eom_prop_te(self,value="show"):
        '''\nName: CC_EOM_PROP_TE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Request for calculation of non-relaxed two-particle EOM-CC properties. The two-particle properties currently include <S^2>. The one-particle properties also will be calculated, since the additional cost of the one-particle properties calculation
is inferior compared to the cost of <S^2>. The variable CC_EOM_PROP must be also set to TRUE. Alternatively, CC_CALC_SSQ can be used to request <S^2> calculation.\n    '''
        if value == "":
            if "CC_EOM_PROP_TE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_EOM_PROP_TE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_EOM_PROP_TE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_EOM_PROP_TE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_EOM_PROP_TE"]=value.lower()


    def cc_fullresponse(self,value="show"):
        '''\nName: CC_FULLRESPONSE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Fully relaxed properties (including orbital relaxation terms) will be computed. The variable CC EOM PROP must be also set to TRUE.\nRecommendation: : Not available for non-UHF/RHF references. Only available for EOM/CI methods for which analytic gradients are available.    '''
        if value == "":
            if "CC_FULLRESPONSE" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_FULLRESPONSE"]
                print "Keyword removed."
        elif value == "show":
            if "CC_FULLRESPONSE" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_FULLRESPONSE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_FULLRESPONSE"]=value.lower()


    def xopt_state_1(self,value="show"):
        '''\nName: XOPT_STATE_1\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Specify two electronic states the intersection of which will be searched.
[spin, irrep, state]
spin = 0 Addresses states with low spin, see also EOM EE SINGLETS.
spin = 1 Addresses states with high spin, see also EOM EE TRIPLETS.
irrep Species the irreducible representation to which the state belongs, for C2v point group symmetry 
irrep = 1 for A1, irrep = 2 for A2,
irrep = 3 for B1, irrep = 4 for B2.
state Species the state number within the irreducible
representation, state = 1 means the lowest excited
state, state = 2 is the second excited state, etc.
0, 0, -1 Ground state.\n    '''
        if value == "":
            if "XOPT_STATE_1" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_STATE_1"]
                print "Keyword removed."
        elif value == "show":
            if "XOPT_STATE_1" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_STATE_1"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XOPT_STATE_1"]=value.lower()


    def xopt_state_2(self,value="show"):
        '''\nName: XOPT_STATE_2\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Specify two electronic states the intersection of which will be searched.
[spin, irrep, state]
spin = 0 Addresses states with low spin, see also EOM EE SINGLETS.
spin = 1 Addresses states with high spin, see also EOM EE TRIPLETS.
irrep Species the irreducible representation to which the state belongs, for C2v point group symmetry 
irrep = 1 for A1, irrep = 2 for A2,
irrep = 3 for B1, irrep = 4 for B2.
state Species the state number within the irreducible
representation, state = 1 means the lowest excited
state, state = 2 is the second excited state, etc.
0, 0, -1 Ground state.\n    '''
        if value == "":
            if "XOPT_STATE_2" in self.dict_of_keywords:
                del self.dict_of_keywords["XOPT_STATE_2"]
                print "Keyword removed."
        elif value == "show":
            if "XOPT_STATE_2" in self.dict_of_keywords:
                return self.dict_of_keywords["XOPT_STATE_2"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["XOPT_STATE_2"]=value.lower()


    def eom_davidson_threshold(self,value="show"):
        '''\nName: EOM_DAVIDSON_THRESHOLD\nType: STRING\nDefault: 0\nOptions: 00105\nDescription: Species threshold for including a new expansion vector in the iterative Davidson diagonalization. Their norm must be above this threshold.
abcde Integer code is mapped to abc x 10^-de\n    '''
        if value == "":
            if "EOM_DAVIDSON_THRESHOLD" in self.dict_of_keywords:
                del self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"]
                print "Keyword removed."
        elif value == "show":
            if "EOM_DAVIDSON_THRESHOLD" in self.dict_of_keywords:
                return self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EOM_DAVIDSON_THRESHOLD"]=value.lower()


    def cc_diis(self,value="show"):
        '''\nName: CC_DIIS\nType: INTEGER\nDefault: 2\nOptions: 0:2:0:1\nDescription: Specify the version of Pulay's Direct Inversion of the Iterative Subspace (DIIS) convergence accelerator to be used in the coupled{cluster code.
0 Activates procedure 2 initially, and procedure 1 when gradients are smaller
than DIIS12 SWITCH.
1 Uses error vectors dened as dierences between parameter vectors from
successive iterations. Most ecient near convergence.
2 Error vectors are dened as gradients scaled by square root of the
approximate diagonal Hessian. Most ecient far from convergence.\nRecommendation: : DIIS1 can be more stable. If DIIS problems are encountered in the early stages of a calculation (when gradients are large) try DIIS 1.    '''
        if value == "":
            if "CC_DIIS" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DIIS"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DIIS" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DIIS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DIIS"]=value.lower()


    def cc_dov_thresh(self,value="show"):
        '''\nName: CC_DOV_THRESH\nType: STRING\nDefault: 0\nOptions: 2502\nDescription: Specifies the minimum allowed values for the coupled-cluster energy denominators.  Smaller values are replaced by this constant during the early iterations only, so the final results are unaffected, but initial convergence is improved when the guess is poor.OPTIONS: abcde Integer code is mapped to abc x 10^-de
RECOMMENDATION:
\nRecommendation: : Increase to 0.5 or 0.75 for non-convergent coupled-cluster calculations.    '''
        if value == "":
            if "CC_DOV_THRESH" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_DOV_THRESH"]
                print "Keyword removed."
        elif value == "show":
            if "CC_DOV_THRESH" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_DOV_THRESH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_DOV_THRESH"]=value.lower()


    def cc_nguess_singles(self,value="show"):
        '''\nName: CC_NGUESS_SINGLES\nType: INTEGER\nDefault: 2\nOptions: 0:200:0:1\nDescription: Specifies number of excited state guess vectors which are single excitations. \nRecommendation: : Should be greater or equal than the number of excited states requested.    '''
        if value == "":
            if "CC_NGUESS_SINGLES" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_NGUESS_SINGLES"]
                print "Keyword removed."
        elif value == "show":
            if "CC_NGUESS_SINGLES" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_NGUESS_SINGLES"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_NGUESS_SINGLES"]=value.lower()


    def cc_state_to_opt(self,value="show"):
        '''\nName: CC_STATE_TO_OPT\nType: undefined\nDefault: 0\nOptions: 0\nDescription: Species which state to optimize.
[i,j] optimize the jth state of the ith irrep.\n    '''
        if value == "":
            if "CC_STATE_TO_OPT" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_STATE_TO_OPT"]
                print "Keyword removed."
        elif value == "show":
            if "CC_STATE_TO_OPT" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_STATE_TO_OPT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_STATE_TO_OPT"]=value.lower()


    def cc_trans_prop(self,value="show"):
        '''\nName: CC_TRANS_PROP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Whether or not the transition dipole moment (in atomic units) and oscillator strength for the EOM-CCSD target states will be calculated. By default, the transition dipole moment is calculated between the CCSD reference and the EOM-CCSD target states. In order to calculate transition dipole moment between a set of EOM-CCSD states and another EOM-CCSD state, the CC_STATE_TO_OPT must be specified for this state.\n    '''
        if value == "":
            if "CC_TRANS_PROP" in self.dict_of_keywords:
                del self.dict_of_keywords["CC_TRANS_PROP"]
                print "Keyword removed."
        elif value == "show":
            if "CC_TRANS_PROP" in self.dict_of_keywords:
                return self.dict_of_keywords["CC_TRANS_PROP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CC_TRANS_PROP"]=value.lower()


    def symmetry_tolerance(self,value="show"):
        '''\nName: SYMMETRY_TOLERANCE\nType: INTEGER\nDefault: 2\nOptions: 0 :10:5:1\nDescription: Controls the tolerance for determining point group symmetry. Differences in atom locations less than 10-SYM_TOL are treated as zero.\nRecommendation: : Use the default unless the molecule has high symmetry which is not being correctly identified. Note that relaxing this tolerance too much may introduce errors into the calculation.    '''
        if value == "":
            if "SYMMETRY_TOLERANCE" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_TOLERANCE"]
                print "Keyword removed."
        elif value == "show":
            if "SYMMETRY_TOLERANCE" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_TOLERANCE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SYMMETRY_TOLERANCE"]=value.lower()


    def qui_multiplicity(self,value="show"):
        '''\nName: QUI_MULTIPLICITY\nType: INTEGER\nDefault: 2\nOptions: 0:20:1:2\nDescription: Sets the multiplicity of the system\n    '''
        if value == "":
            if "QUI_MULTIPLICITY" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_MULTIPLICITY"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_MULTIPLICITY" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_MULTIPLICITY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_MULTIPLICITY"]=value.lower()


    def efp(self,value="show"):
        '''\nName: EFP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: The keyword should be present if excited state calculation is requested.\n    '''
        if value == "":
            if "EFP" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP"]
                print "Keyword removed."
        elif value == "show":
            if "EFP" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EFP"]=value.lower()


    def efp_fragments_only(self,value="show"):
        '''\nName: EFP_FRAGMENTS_ONLY\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Set to true if there is no QM part to the calculation.\n    '''
        if value == "":
            if "EFP_FRAGMENTS_ONLY" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP_FRAGMENTS_ONLY"]
                print "Keyword removed."
        elif value == "show":
            if "EFP_FRAGMENTS_ONLY" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP_FRAGMENTS_ONLY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EFP_FRAGMENTS_ONLY"]=value.lower()


    def efp_input(self,value="show"):
        '''\nName: EFP_INPUT\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: True indicates the new format without a dummy atom in the $molecule section.  False indicates the old format which requires a dummy atom (e.g. He) in the $molecule section for an EFP-only calculation.\n    '''
        if value == "":
            if "EFP_INPUT" in self.dict_of_keywords:
                del self.dict_of_keywords["EFP_INPUT"]
                print "Keyword removed."
        elif value == "show":
            if "EFP_INPUT" in self.dict_of_keywords:
                return self.dict_of_keywords["EFP_INPUT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EFP_INPUT"]=value.lower()


    def exchange(self,value="show"):
        '''\nName: EXCHANGE\nType: STRING\nDefault: 0\nOptions: HF:----------:B88//B:muB88:GG99:Gill96//Gill:Slater:PW86:rPW86:PW91:mPW1PW91:mPW1LYP:mPW1PBE:mPW1B95:mPWB1K:TPSS:---------:B1B95:B3LYP:B3LYP5:B3PW91:B3P86:B5050LYP:B97:B97-1:B97-2:B97-D:BHHLYP:BMK:BOP:CAM-B3LYP:EDF1:EDF2:HC-O3LYP:HCTH:HCTH-120:HCTH-147:HCTH-407:SOGGA:SOGGA11:SOGGA11X:TPSSH:VSXC:X3LYP:----------:BR89:BR89B94h:omegaB97:omegaB97X:omegaB97X-D:omegaB97X-2(LP):omegaB97X-2(TQZ):----------:M05:M052X:M06:M06L:M06HF:M062X:M08HX:M08SO:M11:M11L:---------:PBE:PBE0//PBE1PBE:PBE50:revPBE:PBEOP:wPBE:muPBE:LRC-wPBEPBE:LRC-wPBEhPBE:---------:B05:B3tLAP:BM05:MCY2:XYG3:XYGJOS:LXYGJOS:---------:User-defined//gen\nDescription: Specifies the exchange level of theory.\nRecommendation: : Consult the literature and reviews for guidence    '''
        if value == "":
            if "EXCHANGE" in self.dict_of_keywords:
                del self.dict_of_keywords["EXCHANGE"]
                print "Keyword removed."
        elif value == "show":
            if "EXCHANGE" in self.dict_of_keywords:
                return self.dict_of_keywords["EXCHANGE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["EXCHANGE"]=value.lower()


    def chemsol_read_vdw(self,value="show"):
        '''\nName: CHEMSOL_READ_VDW\nType: STRING\nDefault: 0\nOptions: Default//false:User-defined//true\nDescription: Controls the input of user-defined atomic radii for a ChemSol calculation.\n    '''
        if value == "":
            if "CHEMSOL_READ_VDW" in self.dict_of_keywords:
                del self.dict_of_keywords["CHEMSOL_READ_VDW"]
                print "Keyword removed."
        elif value == "show":
            if "CHEMSOL_READ_VDW" in self.dict_of_keywords:
                return self.dict_of_keywords["CHEMSOL_READ_VDW"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHEMSOL_READ_VDW"]=value.lower()


    def pcm_print(self,value="show"):
        '''\nName: PCM_PRINT\nType: INTEGER\nDefault: 2\nOptions: 0:5:0:1\nDescription: Controls the print level during PCM calculations.\n    '''
        if value == "":
            if "PCM_PRINT" in self.dict_of_keywords:
                del self.dict_of_keywords["PCM_PRINT"]
                print "Keyword removed."
        elif value == "show":
            if "PCM_PRINT" in self.dict_of_keywords:
                return self.dict_of_keywords["PCM_PRINT"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["PCM_PRINT"]=value.lower()


    def qui_solvent_cosmo(self,value="show"):
        '''\nName: QUI_SOLVENT_COSMO\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: \n    '''
        if value == "":
            if "QUI_SOLVENT_COSMO" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_COSMO"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_COSMO" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_COSMO"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_COSMO"]=value.lower()


    def qui_solvent_pcm(self,value="show"):
        '''\nName: QUI_SOLVENT_PCM\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE\nDescription: Use an apparent surface charge polarizable continuum solvent model.\n    '''
        if value == "":
            if "QUI_SOLVENT_PCM" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_PCM"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_PCM" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_PCM"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_PCM"]=value.lower()


    def solvent_method(self,value="show"):
        '''\nName: SOLVENT_METHOD\nType: STRING\nDefault: 0\nOptions: SCRF:PCM:COSMO\nDescription: Sets the preferred solvent model.\n    '''
        if value == "":
            if "SOLVENT_METHOD" in self.dict_of_keywords:
                del self.dict_of_keywords["SOLVENT_METHOD"]
                print "Keyword removed."
        elif value == "show":
            if "SOLVENT_METHOD" in self.dict_of_keywords:
                return self.dict_of_keywords["SOLVENT_METHOD"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SOLVENT_METHOD"]=value.lower()


    def sol_order(self,value="show"):
        '''\nName: SOL_ORDER\nType: INTEGER\nDefault: 2\nOptions: 1:25:15:1\nDescription: Determines the order to which the multipole expansion of the solute charge density is carried out.\n    '''
        if value == "":
            if "SOL_ORDER" in self.dict_of_keywords:
                del self.dict_of_keywords["SOL_ORDER"]
                print "Keyword removed."
        elif value == "show":
            if "SOL_ORDER" in self.dict_of_keywords:
                return self.dict_of_keywords["SOL_ORDER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SOL_ORDER"]=value.lower()


    def svp(self,value="show"):
        '''\nName: SVP\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE:FALSE:TRUE\nDescription: Sets whether to perform the isodensity solvation procedure.\n    '''
        if value == "":
            if "SVP" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP"]
                print "Keyword removed."
        elif value == "show":
            if "SVP" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP"]=value.lower()


    def svp_charge_conv(self,value="show"):
        '''\nName: SVP_CHARGE_CONV\nType: INTEGER\nDefault: 0\nOptions: 7\nDescription: Determines the convergence value for the charges on the cavity. When the change in charges fall below this value, if the electron density is converged, then the calculation is considered converged.\nRecommendation: : The default value unless convergence problems arise.    '''
        if value == "":
            if "SVP_CHARGE_CONV" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_CHARGE_CONV"]
                print "Keyword removed."
        elif value == "show":
            if "SVP_CHARGE_CONV" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_CHARGE_CONV"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP_CHARGE_CONV"]=value.lower()


    def svp_guess(self,value="show"):
        '''\nName: SVP_GUESS\nType: STRING\nDefault: 0\nOptions: No Guess//0:Read//1:Specify//2\nDescription: Specifies how and if the solvation module will use a given guess for the charges and cavity points.\nRecommendation: : It is helpful to also set SCF_GUESS to READ when using a guess from a previous Q-Chem run.     '''
        if value == "":
            if "SVP_GUESS" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_GUESS"]
                print "Keyword removed."
        elif value == "show":
            if "SVP_GUESS" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_GUESS"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP_GUESS"]=value.lower()


    def svp_memory(self,value="show"):
        '''\nName: SVP_MEMORY\nType: INTEGER\nDefault: 2\nOptions: 32:2048:125:1\nDescription: Specifies the amount of memory for use by the solvation module.\nRecommendation: :     '''
        if value == "":
            if "SVP_MEMORY" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_MEMORY"]
                print "Keyword removed."
        elif value == "show":
            if "SVP_MEMORY" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_MEMORY"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP_MEMORY"]=value.lower()


    def svp_path(self,value="show"):
        '''\nName: SVP_PATH\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE:FALSE:TRUE\nDescription: Specifies whether to run a gas phase computation prior to performing the solvation procedure.\nRecommendation: : Running the gas-phase calculation provides a good guess to start the solvation stage and provides a more complete set of solvated properties.    '''
        if value == "":
            if "SVP_PATH" in self.dict_of_keywords:
                del self.dict_of_keywords["SVP_PATH"]
                print "Keyword removed."
        elif value == "show":
            if "SVP_PATH" in self.dict_of_keywords:
                return self.dict_of_keywords["SVP_PATH"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SVP_PATH"]=value.lower()


    def symmetry_decomposition(self,value="show"):
        '''\nName: SYMMETRY_DECOMPOSITION\nType: INTEGER\nDefault: 1\nOptions: None//0:Molecular Orbitals//1:One Electron Matrices//2\nDescription: Determines symmetry decompositions to calculate.\n    '''
        if value == "":
            if "SYMMETRY_DECOMPOSITION" in self.dict_of_keywords:
                del self.dict_of_keywords["SYMMETRY_DECOMPOSITION"]
                print "Keyword removed."
        elif value == "show":
            if "SYMMETRY_DECOMPOSITION" in self.dict_of_keywords:
                return self.dict_of_keywords["SYMMETRY_DECOMPOSITION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["SYMMETRY_DECOMPOSITION"]=value.lower()


    def qui_solvent_dielectric_cosmo(self,value="show"):
        '''\nName: QUI_SOLVENT_DIELECTRIC_COSMO\nType: INTEGER\nDefault: 2\nOptions: 0.0000:99.9999:0.0000:0.0001\nDescription: Sets the dielectric constant of the solvent\n    '''
        if value == "":
            if "QUI_SOLVENT_DIELECTRIC_COSMO" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_COSMO"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_DIELECTRIC_COSMO" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_COSMO"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_COSMO"]=value.lower()


    def qui_solvent_dielectric_onsager(self,value="show"):
        '''\nName: QUI_SOLVENT_DIELECTRIC_ONSAGER\nType: INTEGER\nDefault: 2\nOptions: 0.0000:99.9999:0.0000:0.0001\nDescription: Sets the dielectric constant of the solvent\n    '''
        if value == "":
            if "QUI_SOLVENT_DIELECTRIC_ONSAGER" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_ONSAGER"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_SOLVENT_DIELECTRIC_ONSAGER" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_ONSAGER"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_SOLVENT_DIELECTRIC_ONSAGER"]=value.lower()


    def cholesky_tol(self,value="show"):
        '''\nName: CHOLESKY_TOL\nType: INTEGER\nDefault: 2\nOptions: 0:16:3:1\nDescription: Tolerance for the Cholesky decomposition of two-electron integrals.
\nRecommendation: :
2 - qualitative calculations
3 - appropriate for most cases
4 - quantatative (error < 10-6 Eh)    '''
        if value == "":
            if "CHOLESKY_TOL" in self.dict_of_keywords:
                del self.dict_of_keywords["CHOLESKY_TOL"]
                print "Keyword removed."
        elif value == "show":
            if "CHOLESKY_TOL" in self.dict_of_keywords:
                return self.dict_of_keywords["CHOLESKY_TOL"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CHOLESKY_TOL"]=value.lower()


    def correlation(self,value="show"):
        '''\nName: CORRELATION\nType: STRING\nDefault: 0\nOptions: None :----------:B94:B94hyb:LYP:LYP(EDF1):PBE:PK09:PW92:PW91:PZ81:P86:TPSS:VWN:Wigner:(PBE)OP:(B88)OP:----------:MP2:MP3:MP4:MP4SDQ:ZAPT2:----------:Local_MP2:RIMP2:SOSMP2:MOSMP2:RILMP2:----------:CCD:CCD(2):CCSD:CCSD(T):CCSD(2):CCSD(fT):CCSD(dT):QCCD:QCISD:QCISD(T):OD:OD(T):OD(2):VOD:VOD(2):VQCCD:----------:CIS(D):RICIS(D):SOSCIS(D):SOSCIS(D0)\nDescription: Specifies the correlation level of theory, either DFT or wavefunction-based.\nRecommendation: : Consult the literature and reviews for guidence    '''
        if value == "":
            if "CORRELATION" in self.dict_of_keywords:
                del self.dict_of_keywords["CORRELATION"]
                print "Keyword removed."
        elif value == "show":
            if "CORRELATION" in self.dict_of_keywords:
                return self.dict_of_keywords["CORRELATION"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["CORRELATION"]=value.lower()


    def qui_integral_decomposition_none(self,value="show"):
        '''\nName: QUI_INTEGRAL_DECOMPOSITION_NONE\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE:FALSE:TRUE\nDescription: No integral decomposition\n    '''
        if value == "":
            if "QUI_INTEGRAL_DECOMPOSITION_NONE" in self.dict_of_keywords:
                del self.dict_of_keywords["QUI_INTEGRAL_DECOMPOSITION_NONE"]
                print "Keyword removed."
        elif value == "show":
            if "QUI_INTEGRAL_DECOMPOSITION_NONE" in self.dict_of_keywords:
                return self.dict_of_keywords["QUI_INTEGRAL_DECOMPOSITION_NONE"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["QUI_INTEGRAL_DECOMPOSITION_NONE"]=value.lower()


    def direct_ri(self,value="show"):
        '''\nName: DIRECT_RI\nType: LOGICAL\nDefault: 0\nOptions: FALSE:TRUE:FALSE:TRUE:FALSE:TRUE\nDescription: \nRecommendation: :
By default, all integrals are used in decomposed format allowing significant reduction of memory use.  If all integrals are transformed back (TRUE option) no memory reduction is achieved and decompostion error is introduced.  However, the integral transformation is performed significantly faster and conventional CC/EOM algorithms are used.    '''
        if value == "":
            if "DIRECT_RI" in self.dict_of_keywords:
                del self.dict_of_keywords["DIRECT_RI"]
                print "Keyword removed."
        elif value == "show":
            if "DIRECT_RI" in self.dict_of_keywords:
                return self.dict_of_keywords["DIRECT_RI"]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords["DIRECT_RI"]=value.lower()

    # ------------------------ End of keyword list ----------------------------
    
    # ------------------------ UNDOCUMENTED KEYWORDS CODE ---------------------
    
    def undoc(self,undoc_str,value="show"):
        '''\nFor rem values without documenation herein, please use the undoc keyword and add them manually'''
        if value == "":
            if undoc_str in self.dict_of_keywords:
                del self.dict_of_keywords[undoc_str]
                print "Keyword removed."
        elif value == "show":
            if undoc_str in self.dict_of_keywords:
                return self.dict_of_keywords[undoc_str]
            else:
                print "Value not set."
        else:
            self.dict_of_keywords[undoc_str]=value.lower()
    
    # ---------------------- END OF UNDOCUMENTED KEYWORDS CODE ----------------
    
        
    def add(self,keyword,value):
        self.dict_of_keywords[keyword.upper()]=value.lower()
        
    def remove(self,keyword):
        del self.dict_of_keywords[keyword.upper()]
        
    def clear(self):
        '''Removes all keywords from array.'''
        self.dict_of_keywords.clear()
        
    def __str__(self):
        str_ret =  "$rem\n"
        for key,value in self.dict_of_keywords.iteritems():
            str_ret += key.upper() + (rem_array.__tabstop-len(key))*" " + value + "\n"
        str_ret += "$end\n"
        return str_ret
    
    def info(self):
        print "Type: rem array"
        print "Keywords: " + str(len(self.dict_of_keywords))
