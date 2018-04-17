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
#                     pyQchem - Output Classes                      #
#                                                                   #
#####################################################################

# All classes here are all hidden. Their only purpose is to structure
# information for the user.


import numpy as _np
import re
from copy import deepcopy

from .input_classes import cartesian
from .utilities import _readinput
from .adc_classes import _parse_adc
from . import constants


########################## INFO CLASSES #############################

class _mm(object):
    '''
    This structure contains molecular mechanics data. Energies are in kcal/mol.
    '''

    def __init__(self, mm_list):
        if (len(mm_list[0]) > 1):
            self.etot = mm_list[0]
            self.ecoulomb = mm_list[1]
            self.evdw = mm_list[2]
            self.etorsion = mm_list[3]
            self.eimptors = mm_list[4]
            self.eureybrad = mm_list[5]
            self.eangle = mm_list[6]
            self.ebond = mm_list[7]
            self.nbonds = mm_list[8]
        else:
            self.etot = mm_list[0][0]
            self.ecoulomb = mm_list[1][0]
            self.evdw = mm_list[2][0]
            self.etorsion = mm_list[3][0]
            self.eimptors = mm_list[4][0]
            self.eureybrad = mm_list[5][0]
            self.eangle = mm_list[6][0]
            self.ebond = mm_list[7][0]
            self.nbonds = mm_list[8]

    def info(self):
        '''
        Prints an overview of the calculated MM energies.
        '''
        print("MM calculation summary")
        print("----------------------")
        print("")
        if type(self.etot) == list:
            print(str(len(self.etot)) + " energies found, printing last:")
            print("")
            print("Number of bonds:\t\t" + str(self.nbonds))
            print("Bond energy:\t\t\t" + str(self.ebond[-1]) + " kcal/mol")
            print("Angle energy:\t\t\t" + str(self.eangle[-1]) + " kcal/mol")
            print("Urey-Bradly energy:\t\t" + str(
                self.eureybrad[-1]) + " kcal/mol")
            print("Improper rotation energy:\t" + str(
                self.eimptors[-1]) + " kcal/mol")
            print(
                "Torsion energy:\t\t\t" + str(self.etorsion[-1]) + " kcal/mol")
            print(
                "van der Waals energy:\t\t" + str(self.evdw[-1]) + " kcal/mol")
            print(
                "Coulomb energy:\t\t\t" + str(self.ecoulomb[-1]) + " kcal/mol")
            print("-------------------------")
            print("Total energy:\t\t\t" + str(self.etot[-1]) + " kcal/mol (" + \
                  str(self.etot[
                          -1] * constants.kcal_pro_mole_to_hartree) + " Hartree)")
        else:
            print("Number of bonds:\t\t" + str(self.nbonds))
            print("Bond energy:\t\t\t" + str(self.ebond) + " kcal/mol")
            print("Angle energy:\t\t\t" + str(self.eangle) + " kcal/mol")
            print("Urey-Bradly energy:\t\t" + str(self.eureybrad) + " kcal/mol")
            print("Improper rotation energy:\t" + str(
                self.eimptors) + " kcal/mol")
            print("Torsion energy:\t\t\t" + str(self.etorsion) + " kcal/mol")
            print("van der Waals energy:\t\t" + str(self.evdw) + " kcal/mol")
            print("Coulomb energy:\t\t\t" + str(self.ecoulomb) + " kcal/mol")
            print("-------------------------")
            print("Total energy:\t\t\t" + str(self.etot) + " kcal/mol (" + \
                  str(
                      self.etot * constants.kcal_pro_mole_to_hartree) + " Hartree)")


class _general(object):
    '''
    This structure contains basic information about the Q-Chem jobfile.
    '''

    def __init__(self, jobtype, version, spin, basis_size, energy, status,
                 inputfile, mm_type, initial_geometry, final_geometry,
                 wall_time, cpu_time):
        self.jobtype = jobtype
        self.version = version
        self.spin = _np.float(spin)
        if mm_type != "mm":
            self.basis_size = int(basis_size)
        try:
            self.energy = _np.float(energy)
        except:
            self.energy = "undefined"
        self.status = status
        self.inputfile = inputfile
        self.mm_type = mm_type
        self.initial_geometry = initial_geometry
        self.final_geometry = final_geometry
        self.wall_time = wall_time
        self.cpu_time = cpu_time

    def info(self):
        '''
        Prints a summary of basic information. Energy is given in Hartree.
        '''
        print("About this job:")
        print("--------------")
        print("")
        print("Q-Chem version:\t\t" + self.version)

        jobstring = "Jobtype:\t\t" + self.jobtype.lower()
        if self.mm_type == "mm":
            jobstring += " (MM calculation)"
        if self.mm_type == "janus":
            jobstring += " (QM/MM calculation of type Janus)"
        if self.mm_type == "oniom":
            jobstring += " (QM/MM calculation of type ONIOM)"
        print(jobstring)

        if self.mm_type == "mm":
            print("MM energy:\t\t" + str(self.energy))
        elif self.mm_type == "":
            print("Basis functions:\t" + str(self.basis_size))
            print("Spin:\t\t\t" + str(self.spin))
            print("SCF energy:\t\t" + str(self.energy))
        else:
            print("Basis functions:\t" + str(self.basis_size))
            print("Spin:\t\t\t" + str(self.spin))
            print("QM/MM total energy:\t" + str(self.energy))

        print("Status:\t\t\t" + self.status)


class _thermo(object):
    '''
    This structure contains thermodynamics information obtained from a frequency calculation.
    '''

    def __init__(self, E, ZPE, ITE, T, p, S, H, F, G, frequencies, intensities,
                 mass, mom_inertia, rot_sym, linear_switch):
        self.E = E
        self.ZPE = ZPE
        self.ITE = ITE
        self.T = T
        self.p = p
        self.S = S
        self.H = H
        self.F = F
        self.G = G
        self.frequencies = _np.asarray(frequencies)
        self.intensities = _np.asarray(intensities)
        self.mass = _np.asarray(mass)
        self.mom_inertia = _np.asarray(mom_inertia)
        self.rot_sym = rot_sym
        self._lin_switch = linear_switch

    def info(self):
        '''
        Prints a summary of basic thermodynamics information (in Hartree).
        The internal thermal energy (ITE) is corrected for zero point energy (ZPE).
        Enthalpy (H), Helmholtz free energy (F) and Gibbs free energy (G) are
        already corrected for internal thermal energy (which includes zero point energy).
        '''
        print("Electronic energy (E):\t\t" + str(self.E))
        print("---------------------")
        print("")
        print("Corrections to E:")
        print("-----------------")
        print("Zero point energy (ZPE):\t" + str(self.ZPE))
        print("Internal thermal energy (ITE):\t" + str(self.ITE))
        print("")
        print("Temperature (T in Kelvin):\t" + str(self.T))
        print("-------------------------")
        print("")
        print("Pressure (p in Pa):\t\t" + str(self.p))
        print("------------------")
        print("")
        print("Thermodynamic potentials:")
        print("------------------------")
        print("Entropy (S):\t\t\t" + str(self.S))
        print("Enthalpy (H):\t\t\t" + str(self.H))
        print("Helmholtz free energy (F):\t" + str(self.F))
        print("Gibbs free energy (G):\t\t" + str(self.G))

    def calculate(self, temp, pressure="", loop_iso=0, loop_freq="", grimme=0,
                  grimme_thresh=100, silent=0):
        '''
        Thermodynamics calculation
        --------------------------

        This function calculates the thermodynamic potentials S, H, F and G for (T,p) pairs.
        T is mandatory, all other parameters are optional.

        Input:
        The variables 'temp' and 'pressure' may contain a single value, or a list / an array (of same length).
        The specification of pressure values is optional, standard pressure will be used if not defined.
        The 'loop_iso' variable specifies the set of isotopes to be used (default is 1).
        The 'loop_freq' variable specifies the set of frequencies to be used (default is value of 'loop_iso').
        Negative frequencies are ignored.

        If 'grimme' is set to 1, hindered rotations are approximated via Grimme's interpolation technique,
        see Grimme, S. Chemistry - A European Journal 2012, 18, 9955-9964.
        The threshold is adjusted via the 'grimme_thresh' variable. Default is 100 wavenumbers.

        To make pyQchem shut up set silent=1.

        Output:
        A numpy array of dimension N*7, with N as number of temperature/pressure
        pairs. The columns are: temperature, pressure, ITE (internal thermal energy), S, H, F and G.

        Example:
        <somefile>.thermo.calculate(373.15) calculates the potentials at 100 deg C and 1 atm.
        '''
        if loop_freq == "":
            loop_freq = loop_iso

        # Check temperature and pressure
        if type(temp) is str or type(temp) is float or type(temp) is int:
            dummy = []  # convert temperature input to list (will be iterated over)
            dummy.append(temp)
            T = _np.asarray(dummy, dtype=float)
        else:
            T = _np.asarray(temp, dtype=float)

        if pressure == "":
            p = _np.ones(_np.size(T)) * 1.01325e5  # from ATM to Pa
        else:
            if type(pressure) is str or type(pressure) is float or type(
                    pressure) is int:
                dummy = []  # convert pressure input to list
                for k in range(
                        len(T)):  # assume same pressure for set of temperature
                    dummy.append(pressure)
                p = _np.asarray(dummy, dtype=float)
            else:
                p = _np.asarray(pressure, dtype=float)

        if (_np.size(T) != _np.size(p)):
            print("Temperature and pressure input differs in length.")
        else:
            # Get frequencies
            dummy = _np.asarray(self.frequencies[loop_freq], dtype=float)
            # Remove negative entries
            freq = dummy[dummy > 0]

        # Get mass (in amu) and moment of inertia (in amu * bohr**2)
        mass = self.mass[loop_iso] * constants.atomic_mass_constant
        mom_inertia = _np.asarray(
            self.mom_inertia[loop_iso]) * constants.atomic_mass_constant * (
                              1e-10 * constants.bohr_to_angstrom) ** 2
        rot_sym = float(self.rot_sym[loop_iso])

        # Translational contributions to S and the thermal energy E
        pi = 3.141592654
        Q_trans = ((2 * pi * mass * constants.Boltzmann_constant * T / (
            constants.Planck_constant) ** 2) ** (
                               3. / 2)) * constants.Boltzmann_constant * T / p
        S_trans = constants.molar_gas_constant * (
                    _np.log(Q_trans) + 1 + 3. / 2)  # in J/mol
        E_trans = 3. / 2 * constants.molar_gas_constant * T  # in J/mol

        # Rotational contributions (for linear molecule)
        if self._lin_switch == 1:
            if silent == 0:
                print("Linear molecule detected...")
            theta_lin = (constants.Planck_constant) ** 2 / (
                    8 * pi ** 2 * _np.max(
                    mom_inertia) * constants.Boltzmann_constant)  # rotational temperature
            Q_rot_lin = T / theta_lin / rot_sym
            S_rot = constants.molar_gas_constant * (_np.log(Q_rot_lin) + 1)
            E_rot = constants.molar_gas_constant * T
        else:
            # Rotational contributions (for nonlinear molecule)
            theta = (constants.Planck_constant) ** 2 / (
                    8 * pi ** 2 * mom_inertia * constants.Boltzmann_constant)  # rotational temperature
            Q_rot = _np.sqrt(pi) / rot_sym * ((T ** (3. / 2)) / _np.sqrt(
                theta[0] * theta[1] * theta[2]))
            S_rot = constants.molar_gas_constant * (_np.log(Q_rot) + 3. / 2)
            E_rot = 3. / 2 * constants.molar_gas_constant * T

        # Vibrational contributions
        dummy = freq * constants.inverse_cm_to_hertz  # conversion from wavenumbers to Hertz
        freq = dummy * constants.Planck_constant / constants.Boltzmann_constant  # conversion to Kelvin
        S_vib = []
        E_vib = []
        S_hind = []

        for k in T:
            dum_term = (freq / k) / (_np.exp(freq / k) - 1) - _np.log(
                1 - _np.exp(-freq / k))
            S_vib.append(constants.molar_gas_constant * dum_term.sum())
            dum_term2 = (freq * (1.0 / 2 + 1.0 / (_np.exp(freq / k) - 1)))
            E_vib.append(constants.molar_gas_constant * dum_term2.sum())

            if grimme >= 1:
                B_average = 10E-44;  # kg.m^2
                if grimme == 2:  # calculate B_av from current frequencies
                    B_average = sum(constants.Planck_constant / (
                            8 * pi ** 2 * constants.speed_of_light_in_vacuum * freq)) / len(
                        freq);
                alpha = 4;
                mu = constants.Planck_constant / (
                        8 * pi ** 2 * constants.speed_of_light_in_vacuum * freq)
                mu_prime = mu * B_average / (mu + B_average)
                dum_grim = 0.5 + _np.log((
                                                 8 * pi ** 3 * constants.Boltzmann_constant * k * mu_prime / constants.Planck_constant ** 2) ** 0.5);
                weight = 1. / (1 + (grimme_thresh / freq) ** alpha);
                S_hind.append(constants.molar_gas_constant * (
                            weight * dum_term + (1 - weight) * dum_grim).sum())

        if grimme >= 1:
            S_vib = S_hind  # replace vibrational contribution with Grimme correction

        # Summmation over all contributions
        S = S_trans + S_rot + S_vib  # in J/(mol.K)
        S = S / constants.Avogadro_constant * constants.joule_to_hartree  # in Hartree
        ITE = E_trans + E_rot + E_vib
        ITE = ITE / constants.Avogadro_constant * constants.joule_to_hartree  # in Hartree

        # Calculate the thermodynamic potentials
        H = self.E + ITE + constants.Boltzmann_constant * T * constants.joule_to_hartree  # in Hartree
        F = self.E + ITE - T * S  # array issue here
        G = H - T * S

        # Print and return array of thermodynamical potentials
        if silent == 0:
            headerline = "Results based on frequencies of loop " + str(
                loop_freq) + " and isotopes of loop " + str(loop_iso)
            print(headerline)
            print("-" * len(headerline))
            print("\nT\tp\t\tITE\t\t S\t\t H\t\t F\t\t G")
            for k, l in enumerate(T):
                print("%.2f\t%.2f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t" % (
                T[k], p[k], ITE[k], S[k], H[k], F[k], G[k]))
        dummy = _np.asarray([T, p, ITE, S, H, F, G])
        data = dummy.transpose()
        return data


class _opt(object):
    '''
    This structure contains information about the geometry optimization. Energies are given in Hartree.
    '''

    def __init__(self, geometries, energies, gradient, gradient_vector,
                 displacement, change, optstat):
        self.geometries = geometries
        self.energies = energies
        self.gradient_vector = gradient_vector
        self.gradient = gradient
        self.displacement = displacement
        self.change = change
        self.status = optstat

    def info(self):
        print("Summary of geometry optimization:")
        print("--------------------------------")
        print("")
        print("Energy\t\tGradient\tDisplacement\tDeltaE")
        for i, k in enumerate(self.energies):
            # print str(k) + "\t\t" + str(self.gradient[i]) + "\t" + str(self.displacement[i]) + "\t" + str(self.change[i])
            print("%.9f\t%.9f\t%.9f\t%.9f\t" % (
            k, self.gradient[i], self.displacement[i], self.change[i]))

    def write_trajectory_xyz(self, filename='trajectory.xyz'):
        """
        This method concatenates the xyz geometries, using the current energy as title (supported by Molden).
        """
        f = open(filename, 'w')
        ret_str = ""
        for i, k in enumerate(self.geometries):
            k.title(str(self.energies[i]))
            ret_str += k.__str__()
        f.write(ret_str)
        f.close()


class _aimd(object):
    """
    This structure contains information about the AIMD steps. Time steps are given in fs, energies in Hartee.
    """

    def __init__(self, temp, N_steps, time_step, total_time, time, energies,
                 drift, kinetic_energies, geometries, aimdstat):
        self.temp = temp
        self.N_steps = N_steps
        self.time_step = time_step
        self.total_time = total_time
        self.time = time
        self.energies = energies
        self.drift = drift
        self.kinetic_energies = kinetic_energies
        self.geometries = geometries
        self.status = aimdstat

    def info(self):
        print("Summary of AIMD calculation:")
        print("--------------------------------")
        print("")
        print("Step\tenergies\ttime (fs)\tdrift")
        for i, k in enumerate(self.energies):
            print("%i\t%.9f\t%8.2f\t%.5f\t" % (
            i + 1, k, self.time[i], self.drift[i]))

    def write_trajectory_xyz(self, filename='trajectory.xyz'):
        """
        This method concatenates the xyz geometries, using the current energy as title (supported by Molden).
        """
        f = open(filename, 'w')
        ret_str = ""
        for i, k in enumerate(self.geometries):
            k.title(str(self.energies[i]))
            ret_str += k.__str__()
        f.write(ret_str)
        f.close()


class _orbitals(object):
    '''
    This structure contains information (occupation, energies) about the orbitals.
    '''

    def __init__(self, N_elec):
        self.N_elec = N_elec
        # to be continued...


####################### MULTI OUTPUTFILE  ###########################

class _multioutput(object):

    def __init__(self, jobs=[]):
        self.list_of_jobs = []
        self.list_of_content = []
        for k in jobs:
            self.add(k)

    def add(self, new_job):
        self.list_of_jobs.append(new_job)
        self.list_of_content.append(new_job.general.jobtype)

    # def remove(self,position=0): #if not specified delete last
    #    del self.list_of_content[position]
    #    del self.list_of_jobs[position]


########################## OUTPUTFILE  ##############################

class _outputfile(object):

    def __init__(self, file_input, silent=False):

        # Check input type
        if type(file_input) == list:
            content = file_input
        else:
            content = open(file_input, "r").readlines()

        spin = '0'
        energy = 'undetermined'
        jobtype = 'undetermined'
        version = 'undetermined'
        basis_size = 'undetermined'
        status = 'unfinished'
        wall_time = -99
        cpu_time = -99

        # flag for detection of basis2 job
        basis2_flag = False

        mm_type = ""
        self.aifdem = 0
        self.N_Fragments = 1
        self.N_SET = 0

        switch = 0
        for line in content:
            if "jobtype" in line.lower() or "JOB_TYPE" in line:
                jobtype = ((line.split())[-1]).lower()
            if "basis2" in line.lower():
                basis2_flag = True
            if ("QM_MM_INTERFACE" in line) or ("qm_mm_interface" in line):
                mm_type = ((line.split())[-1]).lower()
            if "aifdem" in line.lower():
                self.aifdem = ((line.split())[-1]).lower()
            if ("CIS_N_ROOTS" in line):
                self.N_SET = ((line.split())[-1]).lower()
            if "Q-Chem, Version" in line:
                version = (((line.split(","))[1]).split())[1]
            if "<S^2> =" in line:
                spin = (line.split())[2]
            if (
            "Total energy in the final basis set") in line and mm_type != "mm":
                energy = (line.split())[8]
            if ("Convergence criterion met") in line and basis2_flag:
                energy = (line.split())[1]
            if ("Etot:" in line) and (mm_type == "mm"):
                energy = (line.split())[4]
            if ("There are" in line) and ("shells" in line):
                basis_size = (line.split())[5]
            if ("Total job time:" in line):
                wall_time = float(line.split()[3].split("s")[0])
                cpu_time = float(line.split()[4].split("s")[0])
            if "MISSION" in line:
                status = 'finished'
            if "--fragment" in line:
                ifrgm = int((line.split())[-1])
                if ifrgm + 1 > self.N_Fragments:
                    self.N_Fragments = ifrgm + 1
            if "TIME STEPS COMPLETED" in line and jobtype == "aimd":
                status = 'time steps completed'
            # Create corresponding inputfile:
            if switch == 0 and "User input:" in line:
                switch = 1
                infile_content = []
            if switch == 1:
                infile_content.append(line)
            if switch == 1 and "Standard Nuclear Orientation" in line:
                switch = 2
                initial_cartesian = cartesian(
                    "sp job - initial geometry in standard orientation")
                continue
            if switch == 2 and ("Repulsion" in line or "Molecular" in line):
                switch = 3
            elif switch == 2 and "Atom" not in line and "---" not in line:
                dummy = (line.split())[1:]
                initial_cartesian.add_atom(dummy[0], dummy[1], dummy[2],
                                           dummy[3])

        inputfile = _readinput(infile_content, silent)

        # Creating geometry objects under 'general' for convenience, final geometry will be overwritten later if different

        # Final geometry is NOT read directly from inputfile, but from Q-Chem standard orientation output in order to
        # avoid issues with molecule 'read' in batch jobs

        initial_geometry = inputfile.molecule.geometry()
        final_geometry = initial_cartesian

        # The info object 'general' will be created after ALL OTHER objects are finished with parsing

        # Make another round if we have an MM or a QM/MM Janus job (just one MM per step)
        if mm_type == "mm" or mm_type == "janus":
            self._process_mm(content)

        # Make another round if we have an QM/MM ONIOM job (has two MM calulations per step)
        elif mm_type == "oniom":
            self._process_oniom(content)

        if jobtype == "sp":
            # Retrieve rem object from input file
            index = inputfile.list_of_content.index("rem")
            rem = inputfile.list_of_arrays[index]
            adc_variant = ''
            if "METHOD" in rem.dict_of_keywords:
                if "adc" in rem.dict_of_keywords["METHOD"]:
                    adc_variant = rem.dict_of_keywords["METHOD"]
            elif "ADC_ORDER" in rem.dict_of_keywords:
                adc_order = int(rem.dict_of_keywords["ADC_ORDER"])
                adc_ext = False
                if adc_order == 2 and "ADC_EXTENDED" in rem.dict_of_keywords:
                    if int(rem.dict_of_keywords["ADC_EXTENDED"] == 1):
                        adc_ext = True
                adc_variant = "adc(" + str(adc_order) + ")"
                if adc_ext:
                    adc_variant += "-x"

            if adc_variant:
                if "ADC_SOS" in rem.dict_of_keywords:
                    if int(rem.dict_of_keywords("ADC_SOS")) != 0:
                        if "adc(2)" in adc_variant and not "sos" in adc_variant:
                            adc_variant = "sos-" + adc_variant
                if "ADC_CVS" in rem.dict_of_keywords:
                    if int(rem.dict_of_keywords("ADC_CVS")) != 0:
                        if not "cvs" in adc_variant:
                            adc_variant = "cvs-" + adc_variant

                self.adc = _parse_adc(adc_variant, content, silent)

        if jobtype == "freq":
            self._process_freq(content, energy, silent)

        if jobtype == "opt" or jobtype == "optimization" or jobtype == "ts":
            final_geometry = self._process_opt(content, jobtype)

        if jobtype == "aimd":
            final_geometry = self._process_aimd(content)

        if self.aifdem != 0:
            self.EvalStrng = ""
            self.aifdem_E_Excite = 0.0
            self.aifdem_Time = 0.0
            EvalSwitch = 0
            for line in content:
                if ("AIFDEM Time:" in line):
                    self.aifdem_Time = float(line.split()[6])
                if (" EigenVectors " in line) and (EvalSwitch == 1):
                    EvalSwitch = 0
                if EvalSwitch == 1:
                    self.EvalStrng += line
                if " EigenValues \n" in line:
                    EvalSwitch = 1
                _ierr = 0
            for i in range(len(self.EvalStrng.split())):
                _ierr += 1
                if _ierr > 7:
                    print(
                        "Possible Error finding AIFDEM Excitation Energy in output_classes.py")
                if len(self.EvalStrng.split()[i]) > 1:
                    _iFirstline = i
                    break

            _EvalStrngFirstline = self.EvalStrng.split()[_iFirstline]

            self.aifdem_E_Excite = float(
                (_EvalStrngFirstline.split("-"))[1]) - float(
                (_EvalStrngFirstline.split("-"))[2])
            self.aifdem_E_Excite = round(self.aifdem_E_Excite, 7)
        if self.N_SET > 0 and self.aifdem == 0:
            self.excited_states = []
            self.cis_time = 0
            _N_SET = 0
            for line in content:
                if "Excited state" in line:
                    _E_Exc_eV = float(line.split()[-1])
                    _N_SET += 1
                if "Total energy for state" in line:
                    _E_Exc_total = float(line.split()[-1])
                if "Multiplicity:" in line:
                    _Mult = line.split()[-1]
                if "Trans. Mom.:" in line:
                    _momX = line.split()[2]
                    _momY = line.split()[4]
                    _momZ = line.split()[6]
                if "Strength" in line:
                    _osc = line.split()[-1]
                    self.excited_states.append(
                        {"Exc_eV": _E_Exc_eV, "Tot": _E_Exc_total,
                         "Mult": _Mult, "X": _momX, "Y": _momY, "Z": _momZ,
                         "Strength": _osc})
                if "CPU time" in line:
                    self.cis_time = line.split()[-1]
            self.N_SET = _N_SET

        if self.N_SET > 0 and self.aifdem == 0:
            self.excited_states = []
            self.cis_time = 0
            _N_SET = 0
            for line in content:
                if "Excited state" in line:
                    _E_Exc_eV = float(line.split()[-1])
                    _N_SET += 1
                if "Total energy for state" in line:
                    _E_Exc_total = float(line.split()[-1])
                if "Multiplicity:" in line:
                    _Mult = line.split()[-1]
                if "Trans. Mom.:" in line:
                    _momX = line.split()[2]
                    _momY = line.split()[4]
                    _momZ = line.split()[6]
                if "Strength" in line:
                    _osc = line.split()[-1]
                    self.excited_states.append(
                        {"Exc_eV": _E_Exc_eV, "Tot": _E_Exc_total,
                         "Mult": _Mult, "X": _momX, "Y": _momY, "Z": _momZ,
                         "Strength": _osc})
                if "CPU time" in line:
                    self.cis_time = line.split()[-1]
            self.N_SET = _N_SET

        # Finally, we create the global info object 'general'
        self.general = _general(jobtype, version, spin, basis_size, energy,
                                status, inputfile, mm_type, initial_geometry,
                                final_geometry, wall_time, cpu_time)

    def _process_aimd(self, content):
        drift = []
        kinetic_energies = []
        time = []
        energies = []
        geometries = []
        aimdstat = "steps not completed"
        temp = 0
        aimd_step = 0
        drift_switch = 0
        geom_switch = 0
        for line in content:
            if "Simulation temperature" in line:
                temp = float((line.split())[3])
            if "AIMD will take" in line:
                N_steps = int((line.split())[3])
            if "Time step =" in line:
                time_step = float((line.split())[7])  # use fs units
            if "Total simulation time requested" in line:
                total_time = float((line.split())[5])  # use fs units
            if "TIME STEP #" in line:  # AIMD starts
                aimd_step += 1
                time.append(float((line.split())[5]))
            if ("Drift factor =" in line) and (aimd_step > 0):
                drift_switch = 1
            if ("Total" in line) and (drift_switch == 1) and (
                    aimd_step > 0):
                drift.append(float((line.split())[2]))
                kinetic_energies.append(float((line.split())[1]))
                drift_switch = 0
            if ("Total energy in the final" in line) and (aimd_step > 0):
                dummy = float((line.split())[8])
                energies.append(dummy)
            # if "Atom           X                Y                Z" in line and (aimd_step>0):
            if "Standard Nuclear Orientation (Angstroms)" in line and (
                    aimd_step > 0):
                geom_switch = 1
                cycle_name = "time step " + str(aimd_step)
                cart_dummy = cartesian(cycle_name)
            if (geom_switch == 1) and ("------" not in line) and (
                    "Atom" not in line) and ("Nuclear" not in line):
                con = line.split()
                cart_dummy.add_atom(con[1], con[2], con[3], con[4])
            if geom_switch == 1 and "Nuclear Repulsion Energy" in line:
                geometries.append(deepcopy(cart_dummy))
                geom_switch = 0
            if "TIME STEPS COMPLETED" in line:
                aimdstat = "steps completed"
        # The geometry has changed, so let's update a variable in the 'general' info object
        final_geometry = deepcopy(geometries[-1])
        self.aimd = _aimd(temp, N_steps, time_step, total_time, time,
                          energies, drift, kinetic_energies, geometries,
                          aimdstat)
        return final_geometry

    def _process_opt(self, content, jobtype):
        energies = []
        gradient_vector = []
        gradient = [0.0]
        displacement = []
        change = []
        geometries = []
        optstat = "no convergence"
        N_step = 1
        switch = 0
        for line in content:
            if "Energy is" in line:
                dummy = float((line.split())[-1])
                energies.append(dummy)
            if "Gradient   " in line:
                try:
                    dummy = float((line.split())[1])
                except (ValueError, IndexError):
                    dummy = 0.0
                gradient.append(dummy)
            if "Displacement   " in line:
                dummy = float((line.split())[1])
                displacement.append(dummy)
            if "Energy change   " in line:
                try:
                    dummy = float((line.split())[2])
                except (ValueError, IndexError):
                    dummy = 0.0
                change.append(dummy)
            if "**  OPTIMIZATION CONVERGED  **" in line:
                optstat = "converged"
            if re.search("ATOM\s{10,15}X\s{10,17}Y\s{10,17}Z", line, flags=re.IGNORECASE):
                switch = 1
                cycle_name = "Optimization step " + str(N_step)
                cart_dummy = cartesian(cycle_name)
            if switch == 1 and "atom" not in line.lower() \
                    and "Point Group" not in line \
                    and '--------' not in line:
                con = line.split()
                cart_dummy.add_atom(con[1], con[2], con[3], con[4])
            if "Point Group" in line and switch == 1:
                geometries.append(deepcopy(cart_dummy))
                N_step += 1
                switch = 0
            if "Gradient of SCF Energy" in line:
                switch = 2
                grad_dummy = []
            elif "Max gradient component" in line and switch == 2:
                # Assuming that the array will always have a 3xN structure:
                matrix = [[], [], []]
                for i, sp in enumerate(grad_dummy):
                    if not i % 4 == 0:
                        matrix[i % 4 - 1].extend(
                            [float(si) for si in sp[1:]])
                switch = 0
                gradient_vector.append(_np.array(matrix))
            elif switch == 2:
                grad_dummy.append(line.split())
        # The geometry has changed, so let's update a variable in the 'general' info object
        try:
            final_geometry = deepcopy(geometries[-1])
        except IndexError:
            raise ValueError('Could not find the geometries from optimization output.')
        if jobtype == "opt" or jobtype == "optimization":
            self.opt = _opt(geometries, energies, gradient, gradient_vector,
                            displacement, change, optstat)
        else:
            self.ts = _opt(geometries, energies, gradient, gradient_vector,
                           displacement, change, optstat)
        return final_geometry

    def _process_freq(self, content, energy, silent):
        H2kcal = constants.hartree_to_kcal_pro_mole
        E = float(energy)
        temp = []
        press = []
        R_T = []
        enth_corr = []
        entr_corr = []
        zero_point = []
        mass = []
        mom_inertia = []
        rot_sym = []
        switch = 0
        frequencies = []
        intensities = []
        loop = 0
        linear_switch = 0
        for line in content:
            if "STANDARD THERMODYNAMIC QUANTITIES AT" in line:
                temp.append(float((line.split())[4]))
                press.append(1.01325e5 * float((line.split())[7]))
            if "Zero point vibrational energy" in line:
                zero_point.append(float((line.split())[4]) / H2kcal)
            if "gas constant (RT):" in line:
                R_T.append(float((line.split())[3]) / H2kcal)
            if "Total Enthalpy:" in line:
                enth_corr.append(float((line.split())[2]) / H2kcal)
            if "Total Entropy:" in line:
                entr_corr.append(float((line.split())[2]) / 1000 / H2kcal)
            if "We detect a D*h symmetry" in line:
                linear_switch = 1
            if "We detect a C*v symmetry" in line:
                linear_switch = 1
            if "We detect a C*h symmetry" in line:  # necessary because of wrong nomenclature in Q-chem (thermodyn.F)
                linear_switch = 1
            if "Molecular Mass:" in line:
                loop += 1
                if "*" in line:
                    if not silent:
                        print("Warning: Molecular mass in loop " + str(
                            loop) + " is unphysically large. Will use mass of first loop instead.")
                    mass.append(mass[0])
                else:
                    dummy = float((line.split())[2])
                    mass.append(dummy)
            if "Rotational Symmetry Number is" in line:
                rot_sym.append(float((line.split())[4]))
            if "Eigenvalues --" in line:
                if "*" in line:
                    if not silent:
                        print("Warning: Moment of inertia in loop " + str(
                            loop) + " is unphysically large. Will use values of first loop instead.")
                    mom_inertia.append(mom_inertia[0])
                else:
                    dummy = [float((line.split())[2]),
                             float((line.split())[3]),
                             float((line.split())[4])]
                    mom_inertia.append(dummy)

            if "VIBRATIONAL ANALYSIS" in line:
                loop1 = []
                loop2 = []
                switch = 1
                infile_content = []
            if switch == 1:
                if "Frequency:" in line:
                    dummy = line.split()
                    del dummy[0]
                    for k in dummy:
                        loop1.append(float(k))
                if "IR Intens:" in line:
                    dummy = line.split()
                    del dummy[0]
                    del dummy[0]
                    for k in dummy:
                        loop2.append(float(k))
            if switch == 1 and "STANDARD THERMODYNAMIC" in line:
                switch = 0
                frequencies.append(loop1)
                intensities.append(loop2)
        T = _np.asarray(temp)
        p = _np.asarray(press)
        ZPE = _np.asarray(zero_point)
        H = E + _np.asarray(enth_corr)
        S = _np.asarray(entr_corr)
        ITE = _np.asarray(
            enth_corr) - R_T  # RT = pV is subtracted from H to obtain the ZPE corrected ITE
        F = E + ITE - (T * S)
        G = H - (T * S)
        self.thermo = _thermo(E, ZPE, ITE, T, p, S, H, F, G, frequencies,
                              intensities, mass, mom_inertia, rot_sym,
                              linear_switch)

    def _process_oniom(self, content):
        etot = []
        ecoulomb = []
        etorsion = []
        eimptors = []
        eureybrad = []
        eangle = []
        ebond = []
        evdw = []
        etot2 = []
        ecoulomb2 = []
        etorsion2 = []
        eimptors2 = []
        eureybrad2 = []
        eangle2 = []
        ebond2 = []
        evdw2 = []
        mm_step = 0
        first_bond = 1
        first_bond2 = 1
        for line in content:
            if "Step 1: MM calculation on the entire system" in line:
                mm_step = 1
            if "Step 2: MM calculation on the model system" in line:
                mm_step = 2

            if ("MM bonds to file" in line) and (first_bond == 1) and (
                    mm_step == 1):
                nbonds = int((line.split())[1])
                first_bond = 0
            if "Ebond:" in line and (mm_step == 1):
                ebond.append(float((line.split())[1]))
            if "Eangle:" in line and (mm_step == 1):
                eangle.append(float((line.split())[1]))
            if "EUreyBrad:" in line and (mm_step == 1):
                eureybrad.append(float((line.split())[1]))
            if "Eimptors:" in line and (mm_step == 1):
                eimptors.append(float((line.split())[1]))
            if "Etorsion:" in line and (mm_step == 1):
                etorsion.append(float((line.split())[1]))
            if "Evdw:" in line and (mm_step == 1):
                evdw.append(float((line.split())[1]))
            if "Ecoulomb:" in line and (mm_step == 1):
                ecoulomb.append(float((line.split())[1]))
            if "Etot:" in line and (mm_step == 1):
                etot.append(float((line.split())[1]))

            if ("n_qm_bonds =" in line) and (first_bond2 == 1) and (
                    mm_step == 2):
                nbonds2 = int((line.split())[2])
                first_bond2 = 0
            if "Ebond:" in line and (mm_step == 2):
                ebond2.append(float((line.split())[1]))
            if "Eangle:" in line and (mm_step == 2):
                eangle2.append(float((line.split())[1]))
            if "EUreyBrad:" in line and (mm_step == 2):
                eureybrad2.append(float((line.split())[1]))
            if "Eimptors:" in line and (mm_step == 2):
                eimptors2.append(float((line.split())[1]))
            if "Etorsion:" in line and (mm_step == 2):
                etorsion2.append(float((line.split())[1]))
            if "Evdw:" in line and (mm_step == 2):
                evdw2.append(float((line.split())[1]))
            if "Ecoulomb:" in line and (mm_step == 2):
                ecoulomb2.append(float((line.split())[1]))
            if "Etot:" in line and (mm_step == 2):
                etot2.append(float((line.split())[1]))
        # Create MM info object for entire system
        self.mm_total = _mm(
            [etot, ecoulomb, evdw, etorsion, eimptors, eureybrad, eangle,
             ebond, nbonds])
        # Create MM info object for model system
        self.mm_model = _mm(
            [etot2, ecoulomb2, evdw2, etorsion2, eimptors2, eureybrad2,
             eangle2, ebond2, nbonds2])

    def _process_mm(self, content):
        etot = []
        ecoulomb = []
        etorsion = []
        eimptors = []
        eureybrad = []
        eangle = []
        ebond = []
        evdw = []
        first_bond = 1
        for line in content:
            if ("MM bonds to file" in line) and (first_bond == 1):
                nbonds = int((line.split())[1])
                first_bond = 0
            if "Ebond:" in line:
                ebond.append(float((line.split())[1]))
            if "Eangle:" in line:
                eangle.append(float((line.split())[1]))
            if "EUreyBrad:" in line:
                eureybrad.append(float((line.split())[1]))
            if "Eimptors:" in line:
                eimptors.append(float((line.split())[1]))
            if "Etorsion:" in line:
                etorsion.append(float((line.split())[1]))
            if "Evdw:" in line:
                evdw.append(float((line.split())[1]))
            if "Ecoulomb:" in line:
                ecoulomb.append(float((line.split())[1]))
            if "Etot:" in line:
                etot.append(float((line.split())[1]))
        self.mm = _mm(
            [etot, ecoulomb, evdw, etorsion, eimptors, eureybrad, eangle,
             ebond, nbonds])  # Create MM info object
