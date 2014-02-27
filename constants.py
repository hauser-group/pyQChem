# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (c) 2014, Andreas W. Hauser
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

#####################################################################
#                                                                   #
#               pyQchem - Input/Output-Tools for Q-Chem             #
#                                                                   #
#      Constants taken from:  http://physics.nist.gov/constants     #
#                           (1/14/2014)                             #
#                                                                   #
#   Atomic masses taken from http://www.nist.gov/pml/data/comp.cfm  #
#                            (2/8/2014)                             #
#                                                                   #
#####################################################################


def info():
	print "Selection of constants and conversion factors"
	print "---------------------------------------------"
	print ""
	print "Fundamental constants are given in SI units. Factors for conversion from x to y are abbreviated as x_to_y."
	print "Constants are taken from http://physics.nist.gov/constants (1/14/2014)"

# Fundamentals in SI units                                     Value                 Uncertainty           Unit
#---------------------------------------------------------------------------------------------------------------------------
Angstrom                                        =           1.00001495e-10      #  0.000 000 90 e-10        m
atomic_mass_constant                            =           1.660538921e-27     #  0.000 000 073 e-27       kg
atomic_mass_constant_energy_equivalent          =           1.492417954e-10     #  0.000 000 066 e-10       J
atomic_unit_of_1st_hyperpolarizability          =           3.206361449e-53     #  0.000 000 071 e-53       C^3 m^3 J^-2
atomic_unit_of_2nd_hyperpolarizability          =           6.23538054e-65      #  0.000 000 28 e-65        C^4 m^4 J^-3
atomic_unit_of_action                           =           1.054571726e-34     #  0.000 000 047 e-34       J s
atomic_unit_of_charge                           =           1.602176565e-19     #  0.000 000 035 e-19       C
atomic_unit_of_charge_density                   =           1.081202338e12      #  0.000 000 024 e12        C m^-3
atomic_unit_of_electric_dipole_moment           =           8.47835326e-30      #  0.000 000 19 e-30        C m
atomic_unit_of_electric_field                   =           5.14220652e11       #  0.000 000 11 e11         V m^-1
atomic_unit_of_electric_field_gradient          =           9.71736200e21       #  0.000 000 21 e21         V m^-2
atomic_unit_of_electric_polarizability          =           1.6487772754e-41    #  0.000 000 0016 e-41      C^2 m^2 J^-1
atomic_unit_of_electric_potential               =           27.21138505         #  0.000 000 60             V
atomic_unit_of_electricvquadrupole_moment       =           4.486551331e-40     #  0.000 000 099 e-40       C m^2
atomic_unit_of_energy                           =           4.35974434e-18      #  0.000 000 19 e-18        J
atomic_unit_of_force                            =           8.23872278e-8       #  0.000 000 36 e-8         N
atomic_unit_of_length                           =           0.52917721092e-10   #  0.000 000 000 17 e-10    m
atomic_unit_of_mag_dipole_moment                =           1.854801936e-23     #  0.000 000 041 e-23       J T^-1
atomic_unit_of_mag_flux_density                 =           2.350517464e5       #  0.000 000 052 e5         T
atomic_unit_of_magnetizability                  =           7.891036607e-29     #  0.000 000 013 e-29       J T^-2
atomic_unit_of_mass                             =           9.10938291e-31      #  0.000 000 40 e-31        kg
atomic_unit_of_momentum                         =           1.992851740e-24     #  0.000 000 088 e-24       kg m s^-1
atomic_unit_of_permittivity                     =           1.112650056e-10     #  (exact)                  F m^-1
atomic_unit_of_time                             =           2.418884326502e-17  #  0.000 000 000 012 e-17   s
atomic_unit_of_velocity                         =           2.18769126379e6     #  0.000 000 000 71 e6      m s^-1
Avogadro_constant                               =           6.02214129e23       #  0.000 000 27 e23         mol^-1
Bohr_magneton                                   =           927.400968e-26      #  0.000 020 e-26           J T^-1
Bohr_radius                                     =           0.52917721092e-10   #  0.000 000 000 17 e-10    m
Boltzmann_constant                              =           1.3806488e-23       #  0.000 0013 e-23          J K^-1
classical_electron_radius                       =           2.8179403267e-15    #  0.000 000 0027 e-15      m
Compton_wavelength                              =           2.4263102389e-12    #  0.000 000 0016 e-12      m
electric_constant                               =           8.854187817e-12     #  (exact)                  F m^-1
electron_charge_to_mass_quotient                =           -1.758820088e11     #  0.000 000 039 e11        C kg^-1
electron_g_factor                               =           -2.00231930436153   #  0.000 000 000 000 53     
electron_mass                                   =           9.10938291e-31      #  0.000 000 40 e-31        kg
electron_volt                                   =           1.602176565e-19     #  0.000 000 035 e-19       J
elementary_charge                               =           1.602176565e-19     #  0.000 000 035 e-19       C
Fermi_coupling_constant                         =           1.166364e-5         #  0.000 005 e-5            GeV^-2
fine_structure_constant                         =           7.2973525698e-3     #  0.000 000 0024 e-3       
lattice_parameter_of_silicon                    =           543.1020504e-12     #  0.000 0089 e-12          m
Loschmidt_constant_1_bar                        =           2.6516462e25        #  0.000 0024 e25           m^-3
Loschmidt_constant_1_atm                        =           2.6867805e25        #  0.000 0024 e25           m^-3
magnetic_constant                               =           12.566370614e-7     #  (exact)                  N A^-2
molar_gas_constant                              =           8.3144621           #  0.000 0075               J mol^-1 K^-1
molar_Planck_constant                           =           3.9903127176e-10    #  0.000 000 0028 e-10      J s mol^-1
Newtonian_constant_of_gravitation               =           6.67384e-11         #  0.000 80 e-11            m^3 kg^-1 s^-2
nuclear_magneton                                =           5.05078353e-27      #  0.000 000 11 e-27        J T^-1
Planck_constant                                 =           6.62606957e-34      #  0.000 000 29 e-34        J s
Planck_constant_over_2_pi                       =           1.054571726e-34     #  0.000 000 047 e-34       J s
Rydberg_constant                                =           10973731.568539     #  0.000 055                m^-1
speed_of_light_in_vacuum                        =           299792458           #  (exact)                  m s^-1
standard_acceleration_of_gravity                =           9.80665             #  (exact)                  m s^-2
standard_atmosphere                             =           101325              #  (exact)                  Pa
standard_state_pressure                         =           100000              #  (exact)                  Pa
Stefan_Boltzmann_constant                       =           5.670373e-8         #  0.000 021 e-8            W m^-2 K^-4
Thomson_cross_section                           =           0.6652458734e-28    #  0.000 000 0013 e-28      m^2
Wien_frequency_displacement_law_constant        =           5.8789254e10        #  0.000 0053 e10           Hz K^-1
Wien_wavelength_displacement_law_constant       =           2.8977721e-3        #  0.000 0026 e-3           m K



#  Conversion factors                  Value                 Uncertainty           Final Unit
#---------------------------------------------------------------------------------------------------------------------------
eV_to_amu                      =         1.073544150e-9      #  0.000 000 024 e-9        u
eV_to_hartree                  =         3.674932379e-2      #  0.000 000 081 e-2        E_h
eV_to_hertz                    =         2.417989348e14      #  0.000 000 053 e14        Hz
eV_to_inverse_m                =         8.06554429e5        #  0.000 000 18 e5          m^-1
eV_to_inverse_cm               =         8.06554429e3        #  0.000 000 18 e3          cm^-1
eV_to_joule                    =         1.602176565e-19     #  0.000 000 035 e-19       J
eV_to_kelvin                   =         1.1604519e4         #  0.000 0011 e4            K
eV_to_kilogram                 =         1.782661845e-36     #  0.000 000 039 e-36       kg
hartree_to_amu                 =         2.9212623246e-8     #  0.000 000 0021 e-8       u
hartree_to_eV                  =        27.21138505          #  0.000 000 60             eV
hartree_to_hertz               =         6.579683920729e15   #  0.000 000 000 033 e15    Hz
hartree_to_inverse_m           =         2.194746313708e7    #  0.000 000 000 011 e7     m^-1
hartree_to_inverse_cm          =         2.194746313708e5    #  0.000 000 000 011 e5     cm^-1
hartree_to_joule               =         4.35974434e-18      #  0.000 000 19 e-18        J
hartree_to_kelvin              =         3.1577504e5         #  0.000 0029 e5            K
hartree_to_kilogram            =         4.85086979e-35      #  0.000 000 21 e-35        kg
hertz_to_amu                   =         4.4398216689e-24    #  0.000 000 0031 e-24      u
hertz_to_eV                    =         4.135667516e-15     #  0.000 000 091 e-15       eV
hertz_to_hartree               =         1.5198298460045e-16 #  0.000 000 000 0076 e-16  E_h
hertz_to_inverse_m             =         3.335640951e-9      #  (exact)                  m^-1
hertz_to_inverse_cm            =         3.335640951e-11     #  (exact)                  cm^-1
hertz_to_joule                 =         6.62606957e-34      #  0.000 000 29 e-34        J
hertz_to_kelvin                =         4.7992434e-11       #  0.000 0044 e-11          K
hertz_to_kilogram              =         7.37249668e-51      #  0.000 000 33 e-51        kg
joule_to_amu                   =         6.70053585e9        #  0.000 000 30 e9          u
joule_to_eV                    =         6.24150934e18       #  0.000 000 14 e18         eV
joule_to_hartree               =         2.29371248e17       #  0.000 000 10 e17         E_h
joule_to_hertz                 =         1.509190311e33      #  0.000 000 067 e33        Hz
joule_to_inverse_m             =         5.03411701e24       #  0.000 000 22 e24         m^-1
joule_to_inverse_cm            =         5.03411701e22       #  0.000 000 22 e22         cm^-1
joule_to_kelvin                =         7.2429716e22        #  0.000 0066 e22           K
joule_to_kilogram              =         1.112650056e-17     #  (exact)                  kg
kelvin_to_amu                  =         9.2510868e-14       #  0.000 0084 e-14          u
kelvin_to_eV                   =         8.6173324e-5        #  0.000 0078 e-5           eV
kelvin_to_hartree              =         3.1668114e-6        #  0.000 0029 e-6           E_h
kelvin_to_hertz                =         2.0836618e10        #  0.000 0019 e10           Hz
kelvin_to_inverse_m            =        69.503476            #  0.00000063               m^-1
kelvin_to_inverse_cm           =         0.69503476          #  0.0000000063             cm^-1
kelvin_to_joule                =         1.3806488e-23       #  0.000 0013 e-23          J
kelvin_to_kilogram             =         1.5361790e-40       #  0.000 0014 e-40          kg
kilogram_to_amu                =         6.02214129e26       #  0.000 000 27 e26         u
kilogram_to_eV                 =         5.60958885e35       #  0.000 000 12 e35         eV
kilogram_to_hartree            =         2.061485968e34      #  0.000 000 091 e34        E_h
kilogram_to_hertz              =         1.356392608e50      #  0.000 000 060 e50        Hz
kilogram_to_inverse_m          =         4.52443873e41       #  0.000 000 20 e41         m^-1
kilogram_to_inverse_cm         =         4.52443873e39       #  0.000 000 20 e39         m^-1
kilogram_to_joule              =         8.987551787e16      #  (exact)                  J
kilogram_to_kelvin             =         6.5096582e39        #  0.000 0059 e39           K



# Conversions added to NIST selection
#---------------------------------------------------------------------------
bohr_to_angstrom                =         0.5291692998
angstrom_to_bohr                =         1.8897543761
 
joule_to_cal                    =         0.2390057361 
cal_to_joule                    =         4.184
joule_to_kcal                   =         0.2390057361e-3 
kcal_to_joule                   =         4184 
 
inverse_cm_to_hartree           =         1/2.194746313708e5 
inverse_cm_to_hertz             =         1/3.335640951e-11  
inverse_cm_to_joule             =         1/5.03411701e22    
inverse_cm_to_kelvin            =         1/0.69503476       
inverse_cm_to_kilogram          =         1/4.52443873e39    
inverse_cm_to_eV                =         1/8.06554429e3      
 
 
hartree_to_kcal_pro_mole        =         hartree_to_joule*joule_to_kcal*Avogadro_constant
hartree_to_kJ_pro_mole          =         hartree_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_hartree        =         1/hartree_to_kcal_pro_mole
kJ_pro_mole_to_hartree          =         1/hartree_to_kJ_pro_mole
 
 
kelvin_to_kcal_pro_mole         =         kelvin_to_joule*joule_to_kcal*Avogadro_constant
kelvin_to_kJ_pro_mole           =         kelvin_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_kelvin         =         1/kelvin_to_kcal_pro_mole
kJ_pro_mole_to_kelvin           =         1/kelvin_to_kJ_pro_mole
 
hertz_to_kcal_pro_mole          =         hertz_to_joule*joule_to_kcal*Avogadro_constant
hertz_to_kJ_pro_mole            =         hertz_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_hertz          =         1/hertz_to_kcal_pro_mole
kJ_pro_mole_to_hertz            =         1/hertz_to_kJ_pro_mole
 
inverse_cm_to_kcal_pro_mole     =         inverse_cm_to_joule*joule_to_kcal*Avogadro_constant
inverse_cm_to_kJ_pro_mole       =         inverse_cm_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_inverse_cm     =         1/inverse_cm_to_kcal_pro_mole
kJ_pro_mole_to_inverse_cm       =         1/inverse_cm_to_kJ_pro_mole
 
atomic_unit_of_time_to_picosec  =        1e-12/atomic_unit_of_time
atomic_unit_of_time_to_femtosec =        1e-15/atomic_unit_of_time
atomic_unit_of_time_to_attosec  =        1e-18/atomic_unit_of_time

picosec_to_atomic_unit_of_time  =        1/atomic_unit_of_time_to_picosec
femtosec_to_atomic_unit_of_time =        1/atomic_unit_of_time_to_femtosec
attosec_to_atomic_unit_of_time  =        1/atomic_unit_of_time_to_attosec

picosec_to_femtosec             =        1e3
picosec_to_attosec              =        1e6
femtosec_to_picosec             =        1e-3
femtosec_to_attosec             =        1e3
attosec_to_femtosec             =        1e-3
attosec_to_picosec              =        1e-6



# Relative atomic mass per most common isotope
#---------------------------------------------
# Using http://www.nist.gov/pml/data/comp.cfm 

dict_of_atomic_masses = {}

dict_of_atomic_masses['1'] = 1.00782503207
dict_of_atomic_masses['H'] = 1.00782503207
dict_of_atomic_masses['2'] = 3.0160293191
dict_of_atomic_masses['3'] = 6.015122795
dict_of_atomic_masses['Li'] = 6.015122795
dict_of_atomic_masses['4'] = 9.0121822
dict_of_atomic_masses['Be'] = 9.0121822
dict_of_atomic_masses['5'] = 10.0129370
dict_of_atomic_masses['B'] = 10.0129370
dict_of_atomic_masses['6'] = 12.0000000
dict_of_atomic_masses['C'] = 12.0000000
dict_of_atomic_masses['7'] = 14.0030740048
dict_of_atomic_masses['N'] = 14.0030740048
dict_of_atomic_masses['8'] = 15.99491461956
dict_of_atomic_masses['O'] = 15.99491461956
dict_of_atomic_masses['9'] = 18.99840322
dict_of_atomic_masses['F'] = 18.99840322
dict_of_atomic_masses['10'] = 19.9924401754
dict_of_atomic_masses['Ne'] = 19.9924401754
dict_of_atomic_masses['11'] = 22.9897692809
dict_of_atomic_masses['Na'] = 22.9897692809
dict_of_atomic_masses['12'] = 23.985041700
dict_of_atomic_masses['Mg'] = 23.985041700
dict_of_atomic_masses['13'] = 26.98153863
dict_of_atomic_masses['Al'] = 26.98153863
dict_of_atomic_masses['14'] = 27.9769265325
dict_of_atomic_masses['Si'] = 27.9769265325
dict_of_atomic_masses['15'] = 30.97376163
dict_of_atomic_masses['P'] = 30.97376163
dict_of_atomic_masses['16'] = 31.97207100
dict_of_atomic_masses['S'] = 31.97207100
dict_of_atomic_masses['17'] = 34.96885268
dict_of_atomic_masses['Cl'] = 34.96885268
dict_of_atomic_masses['18'] = 35.967545106
dict_of_atomic_masses['Ar'] = 35.967545106
dict_of_atomic_masses['19'] = 38.96370668
dict_of_atomic_masses['K'] = 38.96370668
dict_of_atomic_masses['20'] = 39.96259098
dict_of_atomic_masses['Ca'] = 39.96259098
dict_of_atomic_masses['21'] = 44.9559119
dict_of_atomic_masses['Sc'] = 44.9559119
dict_of_atomic_masses['22'] = 45.9526316
dict_of_atomic_masses['Ti'] = 45.9526316
dict_of_atomic_masses['23'] = 49.9471585
dict_of_atomic_masses['V'] = 49.9471585
dict_of_atomic_masses['24'] = 49.9460442
dict_of_atomic_masses['Cr'] = 49.9460442
dict_of_atomic_masses['25'] = 54.9380451
dict_of_atomic_masses['Mn'] = 54.9380451
dict_of_atomic_masses['26'] = 53.9396105
dict_of_atomic_masses['Fe'] = 53.9396105
dict_of_atomic_masses['27'] = 58.9331950
dict_of_atomic_masses['Co'] = 58.9331950
dict_of_atomic_masses['28'] = 57.9353429
dict_of_atomic_masses['Ni'] = 57.9353429
dict_of_atomic_masses['29'] = 62.9295975
dict_of_atomic_masses['Cu'] = 62.9295975
dict_of_atomic_masses['30'] = 63.9291422
dict_of_atomic_masses['Zn'] = 63.9291422
dict_of_atomic_masses['31'] = 68.9255736
dict_of_atomic_masses['Ga'] = 68.9255736
dict_of_atomic_masses['32'] = 69.9242474
dict_of_atomic_masses['Ge'] = 69.9242474
dict_of_atomic_masses['33'] = 74.9215965
dict_of_atomic_masses['As'] = 74.9215965
dict_of_atomic_masses['34'] = 73.9224764
dict_of_atomic_masses['Se'] = 73.9224764
dict_of_atomic_masses['35'] = 78.9183371
dict_of_atomic_masses['Br'] = 78.9183371
dict_of_atomic_masses['36'] = 77.9203648
dict_of_atomic_masses['Kr'] = 77.9203648
dict_of_atomic_masses['37'] = 84.911789738
dict_of_atomic_masses['Rb'] = 84.911789738
dict_of_atomic_masses['38'] = 83.913425
dict_of_atomic_masses['Sr'] = 83.913425
dict_of_atomic_masses['39'] = 88.9058483
dict_of_atomic_masses['Y'] = 88.9058483
dict_of_atomic_masses['40'] = 89.9047044
dict_of_atomic_masses['Zr'] = 89.9047044
dict_of_atomic_masses['41'] = 92.9063781
dict_of_atomic_masses['Nb'] = 92.9063781
dict_of_atomic_masses['42'] = 91.906811
dict_of_atomic_masses['Mo'] = 91.906811
dict_of_atomic_masses['43'] = 96.906365
dict_of_atomic_masses['Tc'] = 96.906365
dict_of_atomic_masses['44'] = 95.907598
dict_of_atomic_masses['Ru'] = 95.907598
dict_of_atomic_masses['45'] = 102.905504
dict_of_atomic_masses['Rh'] = 102.905504
dict_of_atomic_masses['46'] = 101.905609
dict_of_atomic_masses['Pd'] = 101.905609
dict_of_atomic_masses['47'] = 106.905097
dict_of_atomic_masses['Ag'] = 106.905097
dict_of_atomic_masses['48'] = 105.906459
dict_of_atomic_masses['Cd'] = 105.906459
dict_of_atomic_masses['49'] = 112.904058
dict_of_atomic_masses['In'] = 112.904058
dict_of_atomic_masses['50'] = 111.904818
dict_of_atomic_masses['Sn'] = 111.904818
dict_of_atomic_masses['51'] = 120.9038157
dict_of_atomic_masses['Sb'] = 120.9038157
dict_of_atomic_masses['52'] = 119.904020
dict_of_atomic_masses['Te'] = 119.904020
dict_of_atomic_masses['53'] = 126.904473
dict_of_atomic_masses['I'] = 126.904473
dict_of_atomic_masses['54'] = 123.9058930
dict_of_atomic_masses['Xe'] = 123.9058930
dict_of_atomic_masses['55'] = 132.905451933
dict_of_atomic_masses['Cs'] = 132.905451933
dict_of_atomic_masses['56'] = 129.9063208
dict_of_atomic_masses['Ba'] = 129.9063208
dict_of_atomic_masses['57'] = 137.907112
dict_of_atomic_masses['La'] = 137.907112
dict_of_atomic_masses['58'] = 135.907172
dict_of_atomic_masses['Ce'] = 135.907172
dict_of_atomic_masses['59'] = 140.9076528
dict_of_atomic_masses['Pr'] = 140.9076528
dict_of_atomic_masses['60'] = 141.9077233
dict_of_atomic_masses['Nd'] = 141.9077233
dict_of_atomic_masses['61'] = 144.912749
dict_of_atomic_masses['Pm'] = 144.912749
dict_of_atomic_masses['62'] = 143.911999
dict_of_atomic_masses['Sm'] = 143.911999
dict_of_atomic_masses['63'] = 150.9198502
dict_of_atomic_masses['Eu'] = 150.9198502
dict_of_atomic_masses['64'] = 151.9197910
dict_of_atomic_masses['Gd'] = 151.9197910
dict_of_atomic_masses['65'] = 158.9253468
dict_of_atomic_masses['Tb'] = 158.9253468
dict_of_atomic_masses['66'] = 155.924283
dict_of_atomic_masses['Dy'] = 155.924283
dict_of_atomic_masses['67'] = 164.9303221
dict_of_atomic_masses['Ho'] = 164.9303221
dict_of_atomic_masses['68'] = 161.928778
dict_of_atomic_masses['69'] = 168.9342133
dict_of_atomic_masses['Tm'] = 168.9342133
dict_of_atomic_masses['70'] = 167.933897
dict_of_atomic_masses['Yb'] = 167.933897
dict_of_atomic_masses['71'] = 174.9407718
dict_of_atomic_masses['Lu'] = 174.9407718
dict_of_atomic_masses['72'] = 173.940046
dict_of_atomic_masses['Hf'] = 173.940046
dict_of_atomic_masses['72'] = 175.9414086
dict_of_atomic_masses['73'] = 179.9474648
dict_of_atomic_masses['Ta'] = 179.9474648
dict_of_atomic_masses['74'] = 179.946704
dict_of_atomic_masses['W'] = 179.946704
dict_of_atomic_masses['75'] = 184.9529550
dict_of_atomic_masses['Re'] = 184.9529550
dict_of_atomic_masses['76'] = 183.9524891
dict_of_atomic_masses['Os'] = 183.9524891
dict_of_atomic_masses['77'] = 190.9605940
dict_of_atomic_masses['Ir'] = 190.9605940
dict_of_atomic_masses['78'] = 189.959932
dict_of_atomic_masses['Pt'] = 189.959932
dict_of_atomic_masses['79'] = 196.9665687
dict_of_atomic_masses['Au'] = 196.9665687
dict_of_atomic_masses['80'] = 195.965833
dict_of_atomic_masses['Hg'] = 195.965833
dict_of_atomic_masses['81'] = 202.9723442
dict_of_atomic_masses['Tl'] = 202.9723442
dict_of_atomic_masses['82'] = 203.9730436
dict_of_atomic_masses['Pb'] = 203.9730436
dict_of_atomic_masses['83'] = 208.9803987
dict_of_atomic_masses['Bi'] = 208.9803987
dict_of_atomic_masses['84'] = 208.9824304
dict_of_atomic_masses['Po'] = 208.9824304
dict_of_atomic_masses['85'] = 209.987148
dict_of_atomic_masses['At'] = 209.987148
dict_of_atomic_masses['86'] = 210.990601
dict_of_atomic_masses['Rn'] = 210.990601
dict_of_atomic_masses['87'] = 223.0197359
dict_of_atomic_masses['Fr'] = 223.0197359
dict_of_atomic_masses['88'] = 223.0185022
dict_of_atomic_masses['Ra'] = 223.0185022
dict_of_atomic_masses['89'] = 227.0277521
dict_of_atomic_masses['Ac'] = 227.0277521
dict_of_atomic_masses['90'] = 230.0331338
dict_of_atomic_masses['Th'] = 230.0331338
dict_of_atomic_masses['91'] = 231.0358840
dict_of_atomic_masses['Pa'] = 231.0358840
dict_of_atomic_masses['92'] = 233.0396352
dict_of_atomic_masses['U'] = 233.0396352
dict_of_atomic_masses['93'] = 236.046570
dict_of_atomic_masses['Np'] = 236.046570
dict_of_atomic_masses['94'] = 238.0495599
dict_of_atomic_masses['Pu'] = 238.0495599
dict_of_atomic_masses['95'] = 241.0568291
dict_of_atomic_masses['Am'] = 241.0568291
dict_of_atomic_masses['96'] = 243.0613891
dict_of_atomic_masses['Cm'] = 243.0613891
dict_of_atomic_masses['97'] = 247.070307
dict_of_atomic_masses['Bk'] = 247.070307
dict_of_atomic_masses['98'] = 249.0748535
dict_of_atomic_masses['Cf'] = 249.0748535
dict_of_atomic_masses['99'] = 252.082980
dict_of_atomic_masses['Es'] = 252.082980
dict_of_atomic_masses['100'] = 257.095105
dict_of_atomic_masses['Fm'] = 257.095105
dict_of_atomic_masses['101'] = 258.098431
dict_of_atomic_masses['Md'] = 258.098431
dict_of_atomic_masses['102'] = 259.10103
dict_of_atomic_masses['No'] = 259.10103
dict_of_atomic_masses['103'] = 262.10963
dict_of_atomic_masses['Lr'] = 262.10963
dict_of_atomic_masses['104'] = 265.11670
dict_of_atomic_masses['Rf'] = 265.11670
dict_of_atomic_masses['105'] = 268.12545
dict_of_atomic_masses['Db'] = 268.12545
dict_of_atomic_masses['106'] = 271.13347
dict_of_atomic_masses['Sg'] = 271.13347
dict_of_atomic_masses['107'] = 272.13803
dict_of_atomic_masses['Bh'] = 272.13803
dict_of_atomic_masses['108'] = 270.13465
dict_of_atomic_masses['Hs'] = 270.13465
dict_of_atomic_masses['109'] = 276.15116
dict_of_atomic_masses['Mt'] = 276.15116


