#
# This makes a simple list of jobs from xyz files and runs Q-Chem on all the jobs
#
# MBG (2/2014)
#
import pyQChem as qc
import os

#make a generic rem array
rem=qc.rem_array()
rem.basis("sto-3g")
rem.exchange("hf")

#make a list of jobs
job_list=[]

#for all xyzs in a database, create and append the job to the list
for i in os.popen("ls ../../databases/a24/*.xyz").read().splitlines():

    #make the jobs
    job=qc.inputfile()

    #give job a name
    job.runinfo.name=i.split('/')[-1].replace(".xyz","")

    #read the xyz
    xyz=qc.read(i)
   
    #append the molecule and rem array
    job.add(qc.mol_array(xyz))
    job.add(rem)
    job_list.append(job)

# run all jobs in list using 12 workers
qc.queue(job_list,num_workers=12)
