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
    print("Selection of constants and conversion factors")
    print("---------------------------------------------")
    print("")
    print(
        "Fundamental constants are given in SI units. Factors for conversion from x to y are abbreviated as x_to_y.")
    print(
        "Constants are taken from http://physics.nist.gov/constants (1/14/2014)")


# Fundamentals in SI units                                     Value                 Uncertainty           Unit
# ---------------------------------------------------------------------------------------------------------------------------
Angstrom = 1.00001495e-10  # 0.000 000 90 e-10        m
atomic_mass_constant = 1.660538921e-27  # 0.000 000 073 e-27       kg
atomic_mass_constant_energy_equivalent = 1.492417954e-10  # 0.000 000 066 e-10       J
atomic_unit_of_1st_hyperpolarizability = 3.206361449e-53  # 0.000 000 071 e-53       C^3 m^3 J^-2
atomic_unit_of_2nd_hyperpolarizability = 6.23538054e-65  # 0.000 000 28 e-65        C^4 m^4 J^-3
atomic_unit_of_action = 1.054571726e-34  # 0.000 000 047 e-34       J s
atomic_unit_of_charge = 1.602176565e-19  # 0.000 000 035 e-19       C
atomic_unit_of_charge_density = 1.081202338e12  # 0.000 000 024 e12        C m^-3
atomic_unit_of_electric_dipole_moment = 8.47835326e-30  # 0.000 000 19 e-30        C m
atomic_unit_of_electric_field = 5.14220652e11  # 0.000 000 11 e11         V m^-1
atomic_unit_of_electric_field_gradient = 9.71736200e21  # 0.000 000 21 e21         V m^-2
atomic_unit_of_electric_polarizability = 1.6487772754e-41  # 0.000 000 0016 e-41      C^2 m^2 J^-1
atomic_unit_of_electric_potential = 27.21138505  # 0.000 000 60             V
atomic_unit_of_electricvquadrupole_moment = 4.486551331e-40  # 0.000 000 099 e-40       C m^2
atomic_unit_of_energy = 4.35974434e-18  # 0.000 000 19 e-18        J
atomic_unit_of_force = 8.23872278e-8  # 0.000 000 36 e-8         N
atomic_unit_of_length = 0.52917721092e-10  # 0.000 000 000 17 e-10    m
atomic_unit_of_mag_dipole_moment = 1.854801936e-23  # 0.000 000 041 e-23       J T^-1
atomic_unit_of_mag_flux_density = 2.350517464e5  # 0.000 000 052 e5         T
atomic_unit_of_magnetizability = 7.891036607e-29  # 0.000 000 013 e-29       J T^-2
atomic_unit_of_mass = 9.10938291e-31  # 0.000 000 40 e-31        kg
atomic_unit_of_momentum = 1.992851740e-24  # 0.000 000 088 e-24       kg m s^-1
atomic_unit_of_permittivity = 1.112650056e-10  # (exact)                  F m^-1
atomic_unit_of_time = 2.418884326502e-17  # 0.000 000 000 012 e-17   s
atomic_unit_of_velocity = 2.18769126379e6  # 0.000 000 000 71 e6      m s^-1
Avogadro_constant = 6.02214129e23  # 0.000 000 27 e23         mol^-1
Bohr_magneton = 927.400968e-26  # 0.000 020 e-26           J T^-1
Bohr_radius = 0.52917721092e-10  # 0.000 000 000 17 e-10    m
Boltzmann_constant = 1.3806488e-23  # 0.000 0013 e-23          J K^-1
classical_electron_radius = 2.8179403267e-15  # 0.000 000 0027 e-15      m
Compton_wavelength = 2.4263102389e-12  # 0.000 000 0016 e-12      m
electric_constant = 8.854187817e-12  # (exact)                  F m^-1
electron_charge_to_mass_quotient = -1.758820088e11  # 0.000 000 039 e11        C kg^-1
electron_g_factor = -2.00231930436153  # 0.000 000 000 000 53
electron_mass = 9.10938291e-31  # 0.000 000 40 e-31        kg
electron_volt = 1.602176565e-19  # 0.000 000 035 e-19       J
elementary_charge = 1.602176565e-19  # 0.000 000 035 e-19       C
Fermi_coupling_constant = 1.166364e-5  # 0.000 005 e-5            GeV^-2
fine_structure_constant = 7.2973525698e-3  # 0.000 000 0024 e-3
lattice_parameter_of_silicon = 543.1020504e-12  # 0.000 0089 e-12          m
Loschmidt_constant_1_bar = 2.6516462e25  # 0.000 0024 e25           m^-3
Loschmidt_constant_1_atm = 2.6867805e25  # 0.000 0024 e25           m^-3
magnetic_constant = 12.566370614e-7  # (exact)                  N A^-2
molar_gas_constant = 8.3144621  # 0.000 0075               J mol^-1 K^-1
molar_Planck_constant = 3.9903127176e-10  # 0.000 000 0028 e-10      J s mol^-1
Newtonian_constant_of_gravitation = 6.67384e-11  # 0.000 80 e-11            m^3 kg^-1 s^-2
nuclear_magneton = 5.05078353e-27  # 0.000 000 11 e-27        J T^-1
Planck_constant = 6.62606957e-34  # 0.000 000 29 e-34        J s
Planck_constant_over_2_pi = 1.054571726e-34  # 0.000 000 047 e-34       J s
Rydberg_constant = 10973731.568539  # 0.000 055                m^-1
speed_of_light_in_vacuum = 299792458  # (exact)                  m s^-1
standard_acceleration_of_gravity = 9.80665  # (exact)                  m s^-2
standard_atmosphere = 101325  # (exact)                  Pa
standard_state_pressure = 100000  # (exact)                  Pa
Stefan_Boltzmann_constant = 5.670373e-8  # 0.000 021 e-8            W m^-2 K^-4
Thomson_cross_section = 0.6652458734e-28  # 0.000 000 0013 e-28      m^2
Wien_frequency_displacement_law_constant = 5.8789254e10  # 0.000 0053 e10           Hz K^-1
Wien_wavelength_displacement_law_constant = 2.8977721e-3  # 0.000 0026 e-3           m K

#  Conversion factors                  Value                 Uncertainty           Final Unit
# ---------------------------------------------------------------------------------------------------------------------------
eV_to_amu = 1.073544150e-9  # 0.000 000 024 e-9        u
eV_to_hartree = 3.674932379e-2  # 0.000 000 081 e-2        E_h
eV_to_hertz = 2.417989348e14  # 0.000 000 053 e14        Hz
eV_to_inverse_m = 8.06554429e5  # 0.000 000 18 e5          m^-1
eV_to_inverse_cm = 8.06554429e3  # 0.000 000 18 e3          cm^-1
eV_to_joule = 1.602176565e-19  # 0.000 000 035 e-19       J
eV_to_kelvin = 1.1604519e4  # 0.000 0011 e4            K
eV_to_kilogram = 1.782661845e-36  # 0.000 000 039 e-36       kg
hartree_to_amu = 2.9212623246e-8  # 0.000 000 0021 e-8       u
hartree_to_eV = 27.21138505  # 0.000 000 60             eV
hartree_to_hertz = 6.579683920729e15  # 0.000 000 000 033 e15    Hz
hartree_to_inverse_m = 2.194746313708e7  # 0.000 000 000 011 e7     m^-1
hartree_to_inverse_cm = 2.194746313708e5  # 0.000 000 000 011 e5     cm^-1
hartree_to_joule = 4.35974434e-18  # 0.000 000 19 e-18        J
hartree_to_kelvin = 3.1577504e5  # 0.000 0029 e5            K
hartree_to_kilogram = 4.85086979e-35  # 0.000 000 21 e-35        kg
hertz_to_amu = 4.4398216689e-24  # 0.000 000 0031 e-24      u
hertz_to_eV = 4.135667516e-15  # 0.000 000 091 e-15       eV
hertz_to_hartree = 1.5198298460045e-16  # 0.000 000 000 0076 e-16  E_h
hertz_to_inverse_m = 3.335640951e-9  # (exact)                  m^-1
hertz_to_inverse_cm = 3.335640951e-11  # (exact)                  cm^-1
hertz_to_joule = 6.62606957e-34  # 0.000 000 29 e-34        J
hertz_to_kelvin = 4.7992434e-11  # 0.000 0044 e-11          K
hertz_to_kilogram = 7.37249668e-51  # 0.000 000 33 e-51        kg
joule_to_amu = 6.70053585e9  # 0.000 000 30 e9          u
joule_to_eV = 6.24150934e18  # 0.000 000 14 e18         eV
joule_to_hartree = 2.29371248e17  # 0.000 000 10 e17         E_h
joule_to_hertz = 1.509190311e33  # 0.000 000 067 e33        Hz
joule_to_inverse_m = 5.03411701e24  # 0.000 000 22 e24         m^-1
joule_to_inverse_cm = 5.03411701e22  # 0.000 000 22 e22         cm^-1
joule_to_kelvin = 7.2429716e22  # 0.000 0066 e22           K
joule_to_kilogram = 1.112650056e-17  # (exact)                  kg
kelvin_to_amu = 9.2510868e-14  # 0.000 0084 e-14          u
kelvin_to_eV = 8.6173324e-5  # 0.000 0078 e-5           eV
kelvin_to_hartree = 3.1668114e-6  # 0.000 0029 e-6           E_h
kelvin_to_hertz = 2.0836618e10  # 0.000 0019 e10           Hz
kelvin_to_inverse_m = 69.503476  # 0.00000063               m^-1
kelvin_to_inverse_cm = 0.69503476  # 0.0000000063             cm^-1
kelvin_to_joule = 1.3806488e-23  # 0.000 0013 e-23          J
kelvin_to_kilogram = 1.5361790e-40  # 0.000 0014 e-40          kg
kilogram_to_amu = 6.02214129e26  # 0.000 000 27 e26         u
kilogram_to_eV = 5.60958885e35  # 0.000 000 12 e35         eV
kilogram_to_hartree = 2.061485968e34  # 0.000 000 091 e34        E_h
kilogram_to_hertz = 1.356392608e50  # 0.000 000 060 e50        Hz
kilogram_to_inverse_m = 4.52443873e41  # 0.000 000 20 e41         m^-1
kilogram_to_inverse_cm = 4.52443873e39  # 0.000 000 20 e39         m^-1
kilogram_to_joule = 8.987551787e16  # (exact)                  J
kilogram_to_kelvin = 6.5096582e39  # 0.000 0059 e39           K

# Conversions added to NIST selection
# ---------------------------------------------------------------------------
bohr_to_angstrom = 0.5291692998
angstrom_to_bohr = 1.8897543761

joule_to_cal = 0.2390057361
cal_to_joule = 4.184
joule_to_kcal = 0.2390057361e-3
kcal_to_joule = 4184

inverse_cm_to_hartree = 1 / 2.194746313708e5
inverse_cm_to_hertz = 1 / 3.335640951e-11
inverse_cm_to_joule = 1 / 5.03411701e22
inverse_cm_to_kelvin = 1 / 0.69503476
inverse_cm_to_kilogram = 1 / 4.52443873e39
inverse_cm_to_eV = 1 / 8.06554429e3

hartree_to_kcal_pro_mole = hartree_to_joule * joule_to_kcal * Avogadro_constant
hartree_to_kJ_pro_mole = hartree_to_joule * Avogadro_constant / 1000
kcal_pro_mole_to_hartree = 1 / hartree_to_kcal_pro_mole
kJ_pro_mole_to_hartree = 1 / hartree_to_kJ_pro_mole

kelvin_to_kcal_pro_mole = kelvin_to_joule * joule_to_kcal * Avogadro_constant
kelvin_to_kJ_pro_mole = kelvin_to_joule * Avogadro_constant / 1000
kcal_pro_mole_to_kelvin = 1 / kelvin_to_kcal_pro_mole
kJ_pro_mole_to_kelvin = 1 / kelvin_to_kJ_pro_mole

hertz_to_kcal_pro_mole = hertz_to_joule * joule_to_kcal * Avogadro_constant
hertz_to_kJ_pro_mole = hertz_to_joule * Avogadro_constant / 1000
kcal_pro_mole_to_hertz = 1 / hertz_to_kcal_pro_mole
kJ_pro_mole_to_hertz = 1 / hertz_to_kJ_pro_mole

inverse_cm_to_kcal_pro_mole = inverse_cm_to_joule * joule_to_kcal * Avogadro_constant
inverse_cm_to_kJ_pro_mole = inverse_cm_to_joule * Avogadro_constant / 1000
kcal_pro_mole_to_inverse_cm = 1 / inverse_cm_to_kcal_pro_mole
kJ_pro_mole_to_inverse_cm = 1 / inverse_cm_to_kJ_pro_mole

atomic_unit_of_time_to_picosec = 1e-12 / atomic_unit_of_time
atomic_unit_of_time_to_femtosec = 1e-15 / atomic_unit_of_time
atomic_unit_of_time_to_attosec = 1e-18 / atomic_unit_of_time

picosec_to_atomic_unit_of_time = 1 / atomic_unit_of_time_to_picosec
femtosec_to_atomic_unit_of_time = 1 / atomic_unit_of_time_to_femtosec
attosec_to_atomic_unit_of_time = 1 / atomic_unit_of_time_to_attosec

picosec_to_femtosec = 1e3
picosec_to_attosec = 1e6
femtosec_to_picosec = 1e-3
femtosec_to_attosec = 1e3
attosec_to_femtosec = 1e-3
attosec_to_picosec = 1e-6

# Relative atomic mass per most common isotope
# ---------------------------------------------
# Using http://www.nist.gov/pml/data/comp.cfm 

dict_of_atomic_masses = {}

dict_of_atomic_masses['1'] = 1.00782503207
dict_of_atomic_masses['H'] = 1.00782503207
dict_of_atomic_masses['2'] = 4.00260325415
dict_of_atomic_masses['He'] = 4.00260325415
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

dict_of_atomic_abbr = {}
dict_of_atomic_abbr[1] = "H"
dict_of_atomic_abbr[2] = "He"
dict_of_atomic_abbr[3] = "Li"
dict_of_atomic_abbr[4] = "Be"
dict_of_atomic_abbr[5] = "B"
dict_of_atomic_abbr[6] = "C"
dict_of_atomic_abbr[7] = "N"
dict_of_atomic_abbr[8] = "O"
dict_of_atomic_abbr[9] = "F"
dict_of_atomic_abbr[10] = "Ne"
dict_of_atomic_abbr[11] = "Na"
dict_of_atomic_abbr[12] = "Mg"
dict_of_atomic_abbr[13] = "Al"
dict_of_atomic_abbr[14] = "Si"
dict_of_atomic_abbr[15] = "P"
dict_of_atomic_abbr[16] = "S"
dict_of_atomic_abbr[17] = "Cl"
dict_of_atomic_abbr[18] = "Ar"
dict_of_atomic_abbr[19] = "K"
dict_of_atomic_abbr[20] = "Ca"
dict_of_atomic_abbr[21] = "Sc"
dict_of_atomic_abbr[22] = "Ti"
dict_of_atomic_abbr[23] = "V"
dict_of_atomic_abbr[24] = "Cr"
dict_of_atomic_abbr[25] = "Mn"
dict_of_atomic_abbr[26] = "Fe"
dict_of_atomic_abbr[27] = "Co"
dict_of_atomic_abbr[28] = "Ni"
dict_of_atomic_abbr[29] = "Cu"
dict_of_atomic_abbr[30] = "Zn"
dict_of_atomic_abbr[31] = "Ga"
dict_of_atomic_abbr[32] = "Ge"
dict_of_atomic_abbr[33] = "As"
dict_of_atomic_abbr[34] = "Se"
dict_of_atomic_abbr[35] = "Br"
dict_of_atomic_abbr[36] = "Kr"
dict_of_atomic_abbr[37] = "Rb"
dict_of_atomic_abbr[38] = "Sr"
dict_of_atomic_abbr[39] = "Y"
dict_of_atomic_abbr[40] = "Zr"
dict_of_atomic_abbr[41] = "Nb"
dict_of_atomic_abbr[42] = "Mo"
dict_of_atomic_abbr[43] = "Tc"
dict_of_atomic_abbr[44] = "Ru"
dict_of_atomic_abbr[45] = "Rh"
dict_of_atomic_abbr[46] = "Pd"
dict_of_atomic_abbr[47] = "Ag"
dict_of_atomic_abbr[48] = "Cd"
dict_of_atomic_abbr[49] = "In"
dict_of_atomic_abbr[50] = "Sn"
dict_of_atomic_abbr[51] = "Sb"
dict_of_atomic_abbr[52] = "Te"
dict_of_atomic_abbr[53] = "I"
dict_of_atomic_abbr[54] = "Xe"
dict_of_atomic_abbr[55] = "Cs"
dict_of_atomic_abbr[56] = "Ba"
dict_of_atomic_abbr[57] = "La"
dict_of_atomic_abbr[58] = "Ce"
dict_of_atomic_abbr[59] = "Pr"
dict_of_atomic_abbr[60] = "Nd"
dict_of_atomic_abbr[61] = "Pm"
dict_of_atomic_abbr[62] = "Sm"
dict_of_atomic_abbr[63] = "Eu"
dict_of_atomic_abbr[64] = "Gd"
dict_of_atomic_abbr[65] = "Tb"
dict_of_atomic_abbr[66] = "Dy"
dict_of_atomic_abbr[67] = "Ho"
dict_of_atomic_abbr[68] = "Er"
dict_of_atomic_abbr[69] = "Tm"
dict_of_atomic_abbr[70] = "Yb"
dict_of_atomic_abbr[71] = "Lu"
dict_of_atomic_abbr[72] = "Hf"
dict_of_atomic_abbr[73] = "Ta"
dict_of_atomic_abbr[74] = "W"
dict_of_atomic_abbr[75] = "Re"
dict_of_atomic_abbr[76] = "Os"
dict_of_atomic_abbr[77] = "Ir"
dict_of_atomic_abbr[78] = "Pt"
dict_of_atomic_abbr[79] = "Au"
dict_of_atomic_abbr[80] = "Hg"
dict_of_atomic_abbr[81] = "Tl"
dict_of_atomic_abbr[82] = "Pb"
dict_of_atomic_abbr[83] = "Bi"
dict_of_atomic_abbr[84] = "Po"
dict_of_atomic_abbr[85] = "At"
dict_of_atomic_abbr[86] = "Rn"
dict_of_atomic_abbr[87] = "Fr"
dict_of_atomic_abbr[88] = "Ra"
dict_of_atomic_abbr[89] = "Ac"
dict_of_atomic_abbr[90] = "Th"
dict_of_atomic_abbr[91] = "Pa"
dict_of_atomic_abbr[92] = "U"
dict_of_atomic_abbr[93] = "Np"
dict_of_atomic_abbr[94] = "Pu"
dict_of_atomic_abbr[95] = "Am"
dict_of_atomic_abbr[96] = "Cm"
dict_of_atomic_abbr[97] = "Bk"
dict_of_atomic_abbr[98] = "Cf"
dict_of_atomic_abbr[99] = "Es"
dict_of_atomic_abbr[100] = "Fm"
dict_of_atomic_abbr[101] = "Md"
dict_of_atomic_abbr[102] = "No"
dict_of_atomic_abbr[103] = "Lr"
dict_of_atomic_abbr[104] = "Rf"
dict_of_atomic_abbr[105] = "Db"
dict_of_atomic_abbr[106] = "Sg"
dict_of_atomic_abbr[107] = "Bh"
dict_of_atomic_abbr[108] = "Hs"
dict_of_atomic_abbr[109] = "Mt"
dict_of_atomic_abbr[110] = "Ds"
dict_of_atomic_abbr[111] = "Rg"
dict_of_atomic_abbr[112] = "Uub"
dict_of_atomic_abbr[113] = "Uut"
dict_of_atomic_abbr[114] = "Uuq"
dict_of_atomic_abbr[115] = "Uup"
dict_of_atomic_abbr[116] = "Uuh"
dict_of_atomic_abbr[117] = "Uus"
dict_of_atomic_abbr[118] = "Uuo"

dict_of_atomic_numbers = {}
dict_of_atomic_numbers["H"] = 1
dict_of_atomic_numbers["He"] = 2
dict_of_atomic_numbers["Li"] = 3
dict_of_atomic_numbers["Be"] = 4
dict_of_atomic_numbers["B"] = 5
dict_of_atomic_numbers["C"] = 6
dict_of_atomic_numbers["N"] = 7
dict_of_atomic_numbers["O"] = 8
dict_of_atomic_numbers["F"] = 9
dict_of_atomic_numbers["Ne"] = 10
dict_of_atomic_numbers["Na"] = 11
dict_of_atomic_numbers["Mg"] = 12
dict_of_atomic_numbers["Al"] = 13
dict_of_atomic_numbers["Si"] = 14
dict_of_atomic_numbers["P"] = 15
dict_of_atomic_numbers["S"] = 16
dict_of_atomic_numbers["Cl"] = 17
dict_of_atomic_numbers["Ar"] = 18
dict_of_atomic_numbers["K"] = 19
dict_of_atomic_numbers["Ca"] = 20
dict_of_atomic_numbers["Sc"] = 21
dict_of_atomic_numbers["Ti"] = 22
dict_of_atomic_numbers["V"] = 23
dict_of_atomic_numbers["Cr"] = 24
dict_of_atomic_numbers["Mn"] = 25
dict_of_atomic_numbers["Fe"] = 26
dict_of_atomic_numbers["Co"] = 27
dict_of_atomic_numbers["Ni"] = 28
dict_of_atomic_numbers["Cu"] = 29
dict_of_atomic_numbers["Zn"] = 30
dict_of_atomic_numbers["Ga"] = 31
dict_of_atomic_numbers["Ge"] = 32
dict_of_atomic_numbers["As"] = 33
dict_of_atomic_numbers["Se"] = 34
dict_of_atomic_numbers["Br"] = 35
dict_of_atomic_numbers["Kr"] = 36
dict_of_atomic_numbers["Rb"] = 37
dict_of_atomic_numbers["Sr"] = 38
dict_of_atomic_numbers["Y"] = 39
dict_of_atomic_numbers["Zr"] = 40
dict_of_atomic_numbers["Nb"] = 41
dict_of_atomic_numbers["Mo"] = 42
dict_of_atomic_numbers["Tc"] = 43
dict_of_atomic_numbers["Ru"] = 44
dict_of_atomic_numbers["Rh"] = 45
dict_of_atomic_numbers["Pd"] = 46
dict_of_atomic_numbers["Ag"] = 47
dict_of_atomic_numbers["Cd"] = 48
dict_of_atomic_numbers["In"] = 49
dict_of_atomic_numbers["Sn"] = 50
dict_of_atomic_numbers["Sb"] = 51
dict_of_atomic_numbers["Te"] = 52
dict_of_atomic_numbers["I"] = 53
dict_of_atomic_numbers["Xe"] = 54
dict_of_atomic_numbers["Cs"] = 55
dict_of_atomic_numbers["Ba"] = 56
dict_of_atomic_numbers["La"] = 57
dict_of_atomic_numbers["Ce"] = 58
dict_of_atomic_numbers["Pr"] = 59
dict_of_atomic_numbers["Nd"] = 60
dict_of_atomic_numbers["Pm"] = 61
dict_of_atomic_numbers["Sm"] = 62
dict_of_atomic_numbers["Eu"] = 63
dict_of_atomic_numbers["Gd"] = 64
dict_of_atomic_numbers["Tb"] = 65
dict_of_atomic_numbers["Dy"] = 66
dict_of_atomic_numbers["Ho"] = 67
dict_of_atomic_numbers["Er"] = 68
dict_of_atomic_numbers["Tm"] = 69
dict_of_atomic_numbers["Yb"] = 70
dict_of_atomic_numbers["Lu"] = 71
dict_of_atomic_numbers["Hf"] = 72
dict_of_atomic_numbers["Ta"] = 73
dict_of_atomic_numbers["W"] = 74
dict_of_atomic_numbers["Re"] = 75
dict_of_atomic_numbers["Os"] = 76
dict_of_atomic_numbers["Ir"] = 77
dict_of_atomic_numbers["Pt"] = 78
dict_of_atomic_numbers["Au"] = 79
dict_of_atomic_numbers["Hg"] = 80
dict_of_atomic_numbers["Tl"] = 81
dict_of_atomic_numbers["Pb"] = 82
dict_of_atomic_numbers["Bi"] = 83
dict_of_atomic_numbers["Po"] = 84
dict_of_atomic_numbers["At"] = 85
dict_of_atomic_numbers["Rn"] = 86
dict_of_atomic_numbers["Fr"] = 87
dict_of_atomic_numbers["Ra"] = 88
dict_of_atomic_numbers["Ac"] = 89
dict_of_atomic_numbers["Th"] = 90
dict_of_atomic_numbers["Pa"] = 91
dict_of_atomic_numbers["U"] = 92
dict_of_atomic_numbers["Np"] = 93
dict_of_atomic_numbers["Pu"] = 94
dict_of_atomic_numbers["Am"] = 95
dict_of_atomic_numbers["Cm"] = 96
dict_of_atomic_numbers["Bk"] = 97
dict_of_atomic_numbers["Cf"] = 98
dict_of_atomic_numbers["Es"] = 99
dict_of_atomic_numbers["Fm"] = 100
dict_of_atomic_numbers["Md"] = 101
dict_of_atomic_numbers["No"] = 102
dict_of_atomic_numbers["Lr"] = 103
dict_of_atomic_numbers["Rf"] = 104
dict_of_atomic_numbers["Db"] = 105
dict_of_atomic_numbers["Sg"] = 106
dict_of_atomic_numbers["Bh"] = 107
dict_of_atomic_numbers["Hs"] = 108
dict_of_atomic_numbers["Mt"] = 109
dict_of_atomic_numbers["Ds"] = 110
dict_of_atomic_numbers["Rg"] = 111
dict_of_atomic_numbers["Uub"] = 112
dict_of_atomic_numbers["Uut"] = 113
dict_of_atomic_numbers["Uuq"] = 114
dict_of_atomic_numbers["Uup"] = 115
dict_of_atomic_numbers["Uuh"] = 116
dict_of_atomic_numbers["Uus"] = 117
dict_of_atomic_numbers["Uuo"] = 118

dict_of_atomic_names = {}
dict_of_atomic_names[1] = "Hydrogen"
dict_of_atomic_names[2] = "Helium"
dict_of_atomic_names[3] = "Lithium"
dict_of_atomic_names[4] = "Beryllium"
dict_of_atomic_names[5] = "Boron"
dict_of_atomic_names[6] = "Carbon"
dict_of_atomic_names[7] = "Nitrogen"
dict_of_atomic_names[8] = "Oxygen"
dict_of_atomic_names[9] = "Fluorine"
dict_of_atomic_names[10] = "Neon"
dict_of_atomic_names[11] = "Sodium"
dict_of_atomic_names[12] = "Magnesium"
dict_of_atomic_names[13] = "Aluminum"
dict_of_atomic_names[14] = "Silicon"
dict_of_atomic_names[15] = "Phosphorus"
dict_of_atomic_names[16] = "Sulfur"
dict_of_atomic_names[17] = "Chlorine"
dict_of_atomic_names[18] = "Argon"
dict_of_atomic_names[19] = "Potassium"
dict_of_atomic_names[20] = "Calcium"
dict_of_atomic_names[21] = "Scandium"
dict_of_atomic_names[22] = "Titanium"
dict_of_atomic_names[23] = "Vanadium"
dict_of_atomic_names[24] = "Chromium"
dict_of_atomic_names[25] = "Manganese"
dict_of_atomic_names[26] = "Iron"
dict_of_atomic_names[27] = "Cobalt"
dict_of_atomic_names[28] = "Nickel"
dict_of_atomic_names[29] = "Copper"
dict_of_atomic_names[30] = "Zinc"
dict_of_atomic_names[31] = "Gallium"
dict_of_atomic_names[32] = "Germanium"
dict_of_atomic_names[33] = "Arsenic"
dict_of_atomic_names[34] = "Selenium"
dict_of_atomic_names[35] = "Bromine"
dict_of_atomic_names[36] = "Krypton"
dict_of_atomic_names[37] = "Rubidium"
dict_of_atomic_names[38] = "Strontium"
dict_of_atomic_names[39] = "Yttrium"
dict_of_atomic_names[40] = "Zirconium"
dict_of_atomic_names[41] = "Niobium"
dict_of_atomic_names[42] = "Molybdenum"
dict_of_atomic_names[43] = "Technetium"
dict_of_atomic_names[44] = "Ruthenium"
dict_of_atomic_names[45] = "Rhodium"
dict_of_atomic_names[46] = "Palladium"
dict_of_atomic_names[47] = "Silver"
dict_of_atomic_names[48] = "Cadmium"
dict_of_atomic_names[49] = "Indium"
dict_of_atomic_names[50] = "Tin"
dict_of_atomic_names[51] = "Antimony"
dict_of_atomic_names[52] = "Tellurium"
dict_of_atomic_names[53] = "Iodine"
dict_of_atomic_names[54] = "Xenon"
dict_of_atomic_names[55] = "Cesium"
dict_of_atomic_names[56] = "Barium"
dict_of_atomic_names[57] = "Lanthanum"
dict_of_atomic_names[58] = "Cerium"
dict_of_atomic_names[59] = "Praseodymium"
dict_of_atomic_names[60] = "Neodymium"
dict_of_atomic_names[61] = "Promethium"
dict_of_atomic_names[62] = "Samarium"
dict_of_atomic_names[63] = "Europium"
dict_of_atomic_names[64] = "Gadolinium"
dict_of_atomic_names[65] = "Terbium"
dict_of_atomic_names[66] = "Dysprosium"
dict_of_atomic_names[67] = "Holmium"
dict_of_atomic_names[68] = "Erbium"
dict_of_atomic_names[69] = "Thulium"
dict_of_atomic_names[70] = "Ytterbium"
dict_of_atomic_names[71] = "Lutetium"
dict_of_atomic_names[72] = "Hafnium"
dict_of_atomic_names[73] = "Tantalum"
dict_of_atomic_names[74] = "Tungsten"
dict_of_atomic_names[75] = "Rhenium"
dict_of_atomic_names[76] = "Osmium"
dict_of_atomic_names[77] = "Iridium"
dict_of_atomic_names[78] = "Platinum"
dict_of_atomic_names[79] = "Gold"
dict_of_atomic_names[80] = "Mercury"
dict_of_atomic_names[81] = "Thallium"
dict_of_atomic_names[82] = "Lead"
dict_of_atomic_names[83] = "Bismuth"
dict_of_atomic_names[84] = "Polonium"
dict_of_atomic_names[85] = "Astatine"
dict_of_atomic_names[86] = "Radon"
dict_of_atomic_names[87] = "Francium"
dict_of_atomic_names[88] = "Radium"
dict_of_atomic_names[89] = "Actinium"
dict_of_atomic_names[90] = "Thorium"
dict_of_atomic_names[91] = "Protactinium"
dict_of_atomic_names[92] = "Uranium"
dict_of_atomic_names[93] = "Neptunium"
dict_of_atomic_names[94] = "Plutonium"
dict_of_atomic_names[95] = "Americium"
dict_of_atomic_names[96] = "Curium"
dict_of_atomic_names[97] = "Berkelium"
dict_of_atomic_names[98] = "Californium"
dict_of_atomic_names[99] = "Einsteinium"
dict_of_atomic_names[100] = "Fermium"
dict_of_atomic_names[101] = "Mendelevium"
dict_of_atomic_names[102] = "Nobelium"
dict_of_atomic_names[103] = "Lawrencium"
dict_of_atomic_names[104] = "Rutherfordium"
dict_of_atomic_names[105] = "Dubnium"
dict_of_atomic_names[106] = "Seaborgium"
dict_of_atomic_names[107] = "Bohrium"
dict_of_atomic_names[108] = "Hassium"
dict_of_atomic_names[109] = "Meitnerium"
dict_of_atomic_names[110] = "Darmstadtium"
dict_of_atomic_names[111] = "Roentgenium"
dict_of_atomic_names[112] = "Ununbium"
dict_of_atomic_names[113] = "Ununtrium"
dict_of_atomic_names[114] = "Ununquadium"
dict_of_atomic_names[115] = "Ununpentium"
dict_of_atomic_names[116] = "Ununhexium"
dict_of_atomic_names[117] = "Ununseptium"
dict_of_atomic_names[118] = "Ununoctium"

dict_abbr_to_name = {'Ru': 'Ruthenium', 'Re': 'Rhenium', 'Rf': 'Rutherfordium',
                     'Rg': 'Roentgenium', 'Ra': 'Radium', 'Rb': 'Rubidium',
                     'Rn': 'Radon', 'Rh': 'Rhodium', 'Be': 'Beryllium',
                     'Ba': 'Barium', 'Bh': 'Bohrium', 'Bi': 'Bismuth',
                     'Bk': 'Berkelium', 'Br': 'Bromine', 'Uuh': 'Ununhexium',
                     'H': 'Hydrogen', 'P': 'Phosphorus', 'Os': 'Osmium',
                     'Es': 'Einsteinium', 'Hg': 'Mercury', 'Ge': 'Germanium',
                     'Gd': 'Gadolinium', 'Ga': 'Gallium', 'Uub': 'Ununbium',
                     'Pr': 'Praseodymium', 'Pt': 'Platinum', 'Pu': 'Plutonium',
                     'C': 'Carbon', 'Pb': 'Lead', 'Pa': 'Protactinium',
                     'Pd': 'Palladium', 'Cd': 'Cadmium', 'Po': 'Polonium',
                     'Pm': 'Promethium', 'Hs': 'Hassium', 'Uuq': 'Ununquadium',
                     'Uup': 'Ununpentium', 'Uus': 'Ununseptium',
                     'Ho': 'Holmium', 'Hf': 'Hafnium', 'K': 'Potassium',
                     'He': 'Helium', 'Md': 'Mendelevium', 'Mg': 'Magnesium',
                     'Mo': 'Molybdenum', 'Mn': 'Manganese', 'O': 'Oxygen',
                     'Mt': 'Meitnerium', 'S': 'Sulfur', 'W': 'Tungsten',
                     'Zn': 'Zinc', 'Eu': 'Europium', 'Zr': 'Zirconium',
                     'Er': 'Erbium', 'Ni': 'Nickel', 'No': 'Nobelium',
                     'Na': 'Sodium', 'Nb': 'Niobium', 'Nd': 'Neodymium',
                     'Ne': 'Neon', 'Np': 'Neptunium', 'Fr': 'Francium',
                     'Fe': 'Iron', 'Fm': 'Fermium', 'B': 'Boron',
                     'F': 'Fluorine', 'Sr': 'Strontium', 'N': 'Nitrogen',
                     'Kr': 'Krypton', 'Si': 'Silicon', 'Sn': 'Tin',
                     'Sm': 'Samarium', 'V': 'Vanadium', 'Sc': 'Scandium',
                     'Sb': 'Antimony', 'Sg': 'Seaborgium', 'Se': 'Selenium',
                     'Co': 'Cobalt', 'Cm': 'Curium', 'Cl': 'Chlorine',
                     'Ca': 'Calcium', 'Cf': 'Californium', 'Ce': 'Cerium',
                     'Xe': 'Xenon', 'Lu': 'Lutetium', 'Cs': 'Cesium',
                     'Cr': 'Chromium', 'Cu': 'Copper', 'La': 'Lanthanum',
                     'Li': 'Lithium', 'Tl': 'Thallium', 'Tm': 'Thulium',
                     'Lr': 'Lawrencium', 'Th': 'Thorium', 'Ti': 'Titanium',
                     'Te': 'Tellurium', 'Tb': 'Terbium', 'Tc': 'Technetium',
                     'Ta': 'Tantalum', 'Yb': 'Ytterbium', 'Db': 'Dubnium',
                     'Dy': 'Dysprosium', 'Ds': 'Darmstadtium', 'I': 'Iodine',
                     'U': 'Uranium', 'Y': 'Yttrium', 'Ac': 'Actinium',
                     'Ag': 'Silver', 'Uut': 'Ununtrium', 'Ir': 'Iridium',
                     'Am': 'Americium', 'Al': 'Aluminum', 'As': 'Arsenic',
                     'Ar': 'Argon', 'Au': 'Gold', 'At': 'Astatine',
                     'In': 'Indium'}

# from http://en.wikipedia.org/wiki/Covalent_radius
covalent_radii = [0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
                  0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
                  2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.61, 1.52, 1.5,
                  1.24, 1.32, 1.22, 1.22, 1.2, 1.19, 1.2, 1.2, 1.16,
                  2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42,
                  1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.4,
                  2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98,
                  1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75,
                  1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45,
                  1.46, 1.48, 1.4, 1.5, 1.5, 2.6, 2.21, 2.15, 2.06,
                  2., 1.96, 1.9, 1.87, 1.8, 1.69]

# from http://toc.uni-muenster.de/DFTD3/index.html, times 1.1, then converted to Angstroms
vdw_radii = [1.892, 1.912, 1.559, 2.661, 2.806, 2.744, 2.640, 2.536, 2.432,
             2.349,
             2.162, 2.578, 3.097, 3.243, 3.222, 3.180, 3.097, 3.014, 2.806,
             2.785,
             2.952, 2.952, 2.952, 2.952, 2.952, 2.952, 2.952, 2.952, 2.952,
             2.952,
             3.118, 3.264, 3.326, 3.347, 3.305, 3.264, 3.076, 3.035, 3.097,
             3.097,
             3.097, 3.097, 3.097, 3.097, 3.097, 3.097, 3.097, 3.097, 3.160,
             3.409,
             3.555, 3.575, 3.575, 3.555, 3.405, 3.330, 3.251, 3.313, 3.313,
             3.313,
             3.313, 3.313, 3.313, 3.313, 3.313, 3.313, 3.313, 3.313, 3.313,
             3.313,
             3.313, 3.378, 3.349, 3.349, 3.349, 3.349, 3.349, 3.349, 3.349,
             3.322,
             3.752, 3.673, 3.586, 3.789, 3.762, 3.636]
