import pyQChem as qc

#generate rem array, cartesian object, mol array, and job via standard procedure (see input demos)
rem=qc.rem_array()
rem.basis("6-31g")
rem.jobtype("sp")
rem.method("adc(2)")
rem.ee_singlets("[2,0,0,0,0,2,0,0]")
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
out.adc.info()

#let's compute the excitation energies between the excited states
print("\nExcitation energies between excited states:")
for i in range(0,len(out.adc.list_of_excited_states)-1):
    es1 = out.adc.list_of_excited_states[i]
    for j in range(i + 1,len(out.adc.list_of_excited_states)):
        es2 = out.adc.list_of_excited_states[j]
        print(("{0:10s} -> {1:10s}: {2:12.6f}".format(es1.term_symbol, es2.term_symbol,es2.energy - es1.energy)))        

