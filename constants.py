# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (C) 2014  Andreas W. Hauser

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#####################################################################
#                                                                   #
#               pyQchem - Input/Output-Tools for Q-Chem             #
#                                                                   #
#      Constants taken from:  http://physics.nist.gov/constants     #
#                                                                   #
#                           AWH  1/14/2014                          #
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

