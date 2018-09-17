#
# This script generates a simple H2 input and writes it to disk
#
# MBG (2/2014)
#
#

import pyqchem as qc

#make the rem array and define the basis and exchange (setting them to defaults, but showing how)
rem=qc.rem_array()
rem.basis("sto-3g")
rem.exchange("hf")

#make the cartesian object
xyz=qc.cartesian()
xyz.add_atom() #defaults to H at origin
xyz.add_atom("H","0","0",".74")

#make molecule array from cartesian object
molec=qc.mol_array(xyz)

#make input object and write to disk
job=qc.inputfile()
job.add(rem)
job.add(molec)
job.write("h2.in")
