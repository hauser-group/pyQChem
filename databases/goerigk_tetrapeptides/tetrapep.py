import os,sys

li=os.popen("ls *.in |grep -v mono | sed 's/\.in//' ").read().splitlines()

scf=[]
eos=[]
ess=[]

for i in li:
        d=float(os.popen("grep 'Total energy in the final basis set' "+i+".out | awk '{print $NF;}'").read())
        scf.append(627.5095*(d))
        d=float(os.popen("grep 'TOTAL SS RI-MP2 ENERGY' "+i+".out").read().split()[-2])
        ess.append(627.5095*(d))
        d=float(os.popen("grep 'TOTAL OS RI-MP2 ENERGY' "+i+".out").read().split()[-2])
        eos.append(627.5095*(d))

for i,j in enumerate(li):
        print(j,scf[i],eos[i],ess[i])
