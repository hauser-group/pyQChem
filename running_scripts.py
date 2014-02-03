# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (C) 2014  Matthew Goldey

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
#                 pyQchem - Scripts for running Q-Chem              #
#                                                                   #
#####################################################################

from input_classes import *

def qcrun(inp_file,name='',loc53='',qchem='',nt=1,np=1,timestamp=False):
    """This is a script for running Q-Chem inside of iPython given an input file object, name, 53.0 file, location, and number of threads or processors"""
    #tag it all with current time for safety's sake
    import time,os
    curtime=time.strftime("%Y%m%d%H%M%S")
    if timestamp!=True:
        name=name.replace(".in","")
    else: 
        if name.endswith(".in"):
            name.replace(".in",curtime)
        else:
        	name=name+curtime

    #make script file
    scr=name+".sh"
    scr_out=open(scr,'w')

    #source appropriate Q-Chem
    if qchem!='':
        scr_out.write("source ~/"+qchem+"\n")

    #copy 53.0 if restarting
    if loc53!='':
        inp_file.rem.scf_guess('read')
        scr_out.write("mkdir $QCSCRATCH/"+name+".dir\n")
        if loc53.endswith(".in.53.0"):
            scr_out.write('cp '+loc53+' $QCSCRATCH/'+name+'.dir/53.0\n')
        else:
            scr_out.write('cp '+loc53+'/53.0 $QCSCRATCH/'+name+'.dir/\n')

    inp_file.write(name+".in")

    #write qchem command to script file
    if nt>1:
    	scr_out.write("qchem -nt "+str(nt)+" "+name+".in "+name+".out "+name+".dir "+"\n")
    elif np>1:
    	scr_out.write("qchem -np "+str(nt)+" "+name+".in "+name+".out "+name+".dir "+"\n")
    else:
    	scr_out.write("qchem "+name+".in "+name+".out "+name+".dir "+"\n")

    #close and run
    scr_out.close()
    os.popen("bash "+scr).read()
    return

def rmsd(a,b):
    from math import sqrt
    if type(a)!=type(b):
        return
    a=a.list_of_atoms
    b=b.list_of_atoms
    if len(a)!=len(b):
        print "Geometries incompatible"
        return
    natoms=len(a)
    tot=[0.0 for i in xrange(3)]
    for i in xrange(natoms):
        tot=[tot[j]+(float(a[i][j+1])-float(b[i][j+1]))**2.0 for j in xrange(3)]
    tot=[sqrt(tot[i]/natoms) for i in xrange(3)]
    return tot[0]+tot[1]+tot[2]