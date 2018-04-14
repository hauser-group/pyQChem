import pyQChem as qc

#generate rem array, cartesian object, mol array, and job via standard procedure (see input demos)
rem=qc.rem_array()
rem.basis("sto-3g")
rem.jobtype("opt")
xyz=qc.cartesian()
xyz.add_atom()
xyz.add_atom("H","0","0",".74")
molec=qc.mol_array(xyz)
job=qc.inputfile()
job.add(rem)
job.add(molec)

#run job with name "h2" making h2.in h2.sh, and h2.out
job.run(name="h2")

#read in output
out=qc.read("h2.out")
out.opt.info()

#let's approximate how much this has changed

#grab first geometry
start=out.opt.geometries[0]

#grab last geometry
end=out.opt.geometries[-1]

#Print statistics for geometric distortions
print("\n\nApproximate change between starting and ending geometries by two metrics:\n")
print(qc.utilities.rmsd(start.xyzs,end.xyzs)," or ", qc.utilities.kabsch(start.xyzs,end.xyzs),"\n")

#Check for reordering (not necessary in this case) and print out sequential RMSDs
qc.utilities.orderset(out.opt.geometries)
