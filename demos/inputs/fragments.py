#
# This script generates a small water cluster and runs it using HF with the fragment code in Q-Chem
#
# MBG (3/2015)
#
#


# write a small file
xyz="""12

O          -0.106357    0.087598    0.127176
H           0.851108    0.072355    0.136719
H          -0.337031    1.005310    0.106947
O           2.701100   -0.077292   -0.273980
H           3.278147   -0.563291    0.297560
H           2.693451   -0.568936   -1.095771
O           2.271787   -1.668771   -2.587410
H           1.328156   -1.800266   -2.490761
H           2.384794   -1.339543   -3.467573
O          -0.518887   -1.685783   -2.053795
H          -0.969013   -2.442055   -1.705471
H          -0.524180   -1.044938   -1.342263
"""
xf=open("4water.xyz",'w')
xf.write(xyz)
xf.close()

import pyQChem as pq


#make the rem array 
rem1=pq.rem_array()
rem1.basis("6-31++G**")
rem1.exchange("hf")
rem1.thresh("14")
rem1.scf_convergence("10")
from copy import deepcopy
rem2=deepcopy(rem1)
rem2.scf_guess("fragmo")

#make a rem_frgm array
rem_frgm=pq.rem_frgm_array()
rem_frgm.thresh("7")
rem_frgm.scf_convergence("3")

#make objects for holding the molecular geometries
xyz=pq.read("4water.xyz")
frag=pq.fragment(atom_list=xyz.list_of_atoms)

#make molecule array from cartesian object
mol1=pq.mol_array(xyz)
mol2=pq.mol_array(frag)

#make input object and write to disk
job1=pq.inputfile()
job1.add(rem1)
job1.add(mol1)

job2=pq.inputfile()
job2.add(rem2)
job2.rem.scf_guess("fragmo")
job2.add(rem_frgm)
job2.add(mol2)

job=pq.multifile()
job.add(job1)
job.add(job2)

job.write("4water.in")
