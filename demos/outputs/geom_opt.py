import pyQChem as qc

rem=qc.rem_fragment()
rem.basis("sto-3g")
rem.jobtype("opt")
xyz=qc.cartesian()
xyz.add_atom()
xyz.add_atom("H","0","0",".74")
molec=qc.mol_fragment(xyz)
job=qc.inputfile()
job.add(rem)
job.add(molec)

qc.run(job,"h2")

out=qc.read("h2.out")
out.opt.info()

#let's approximate how much this has changed

start=out.opt.geometries[0]
end=out.opt.geometries[-1]
print "\n\nApproximate change between starting and ending geometries by two metrics:\n"
print qc.utilities.rmsd(start.xyzs,end.xyzs)," or ", qc.utilities.kabsch(start.xyzs,end.xyzs)

#Check for reordering (not necessary in this case) and print out sequential RMSDs
qc.utilities.orderset(out.opt.geometries)
