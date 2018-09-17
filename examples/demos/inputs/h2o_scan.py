#
# This sample file generates center-of-mass displaced water molecules from the standard water dimer
# and loops through different multiplicative displacements
#
# MBG (02/2014)
#
import pyqchem as qc
from pylab import *

# load some sample inputs - here the first and second monomers of the water dimer
a = qc.read("../../databases/s22/Water-dimer_mono1.xyz")
a.name = "water-dimer_mono1"
b = qc.read("../../databases/s22/Water-dimer_mono2.xyz")
b.name = "water-dimer_mono2"

# Save direction and magnitude
d = b.com - a.com

# Translate both to origin
a.move(-a.com)
b.move(-b.com)

# loop over a multiplicative factor for center of mass distance and generate jobs
for i in arange(.9, 2.01, .1):
    # translate to new coordinates
    a.move(-i * d / 2)
    b.move(i * d / 2)

    # form new dimer
    c = a + b

    # write dimer
    c.write("water-dimer" + "_" + str(i) + ".xyz")

    # write monomers with ghost functions
    c = qc.cartesian(atom_list=(a.list_of_atoms + b.ghost()))
    c.write("water-dimer_mono1" + "_" + str(i) + ".xyz")
    c = qc.cartesian(atom_list=(a.ghost() + b.list_of_atoms))
    c.write("water-dimer_mono2" + "_" + str(i) + ".xyz")

    # translate back to origin
    a.move(-a.com)
    b.move(-b.com)
