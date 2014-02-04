import pyQchem as qc
import os
rem=qc.rem_fragment()
rem.basis("sto-3g")

job_list=[]
for i in os.popen("ls ../../databases/a24/*.xyz").read().splitlines():
    job=qc.inputfile()
    job.runinfo.name=i.split('/')[-1].replace(".xyz","")
    xyz=qc.read(i)
    job.add(qc.mol_fragment(xyz))
    job.add(rem)
    job_list.append(job)

qc.queue(job_list,num_workers=12)
