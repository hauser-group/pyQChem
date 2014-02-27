import pyQChem as qc
import os
from pylab import *


a=qc.read("../../databases/s22/Water-dimer_mono1.xyz")
a.name="water-dimer_mono1"
b=qc.read("../../databases/s22/Water-dimer_mono2.xyz")
b.name="water-dimer_mono2"

#Save direction and magnitude
d=b.com-a.com

#Translate both to origin
a.move(-a.com)
b.move(-b.com)

for i in arange(.9,2.01,.1):
	a.move(-i*d/2)
	b.move(i*d/2)
	c=a+b
	c.write("water-dimer"+"_"+str(i)+".xyz")
	c=qc.cartesian(atom_list=(a.list_of_atoms+b.ghost()))
	c.write("water-dimer_mono1"+"_"+str(i)+".xyz")
	c=qc.cartesian(atom_list=(a.ghost()+b.list_of_atoms))
	c.write("water-dimer_mono2"+"_"+str(i)+".xyz")
	a.move(-a.com)
	b.move(-b.com)
