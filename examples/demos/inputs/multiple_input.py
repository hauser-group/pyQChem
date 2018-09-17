#
# This script generates a multifile input for H2 geometry optimization and frequency calculation
#
# MBG (2/2014)
#
import pyqchem as qc


#make the rem array and fill it
rem=qc.rem_array()
rem.basis("sto-3g")
rem.jobtype("opt")

#make the cartesian object and fill it
xyz=qc.cartesian()
xyz.add_atom() #defaults to H at origin
xyz.add_atom("H","0","0",".74")

#make molecule array from cartesian object
molec=qc.mol_array(xyz)

#assemble the job from the arrays
job1=qc.inputfile()
job1.add(rem)
job1.add(molec)

#write it to file
job1.write("h2.in")

#make another rem to do the second job
rem2=qc.deepcopy(rem)
rem2.jobtype("freq")

#define the  molecule array as "READ" so we grab the right geometry
molec2=qc.mol_array()
molec2.geometry("read")

#assemble the second job
job2=qc.inputfile()
job2.add(rem2)
job2.add(molec2)

#make the multifile input and save it
job=qc.multifile()
job.add(job1)
job.add(job2)
job.write("h2_opt.in")
