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
#                     pyQchem - Output Classes                      #
#                                                                   #
#####################################################################

# All classes here are all hidden. Their only purpose is to structure 
# information for the user.


import numpy as np 
from copy import deepcopy

from input_classes import cartesian
from utilities import _readinput
import constants


########################## INFO CLASSES #############################

class _general(object):
    '''
    This structure contains basic information about the Q-Chem jobfile.
    '''
    def __init__(self,jobtype,version,spin,basis_size,energy,status,inputfile):
        self.jobtype = jobtype
        self.version = version
        self.spin = spin
        self.basis_size = basis_size
        self.energy = energy
        self.status = status
        self.inputfile = inputfile

    def info(self):
        '''
        Prints a summary of basic information.
        '''
        print "About this job:"
        print "--------------"
        print ""
        print "Q-Chem version:\t\t" + self.version
        print "Jobtype:\t\t" + self.jobtype
        print "Basis functions:\t" + str(self.basis_size)
        print "Spin:\t\t\t" + str(self.spin)
        print "SCF energy:\t\t" + str(self.energy)
        print "Status:\t\t\t" + self.status


class _thermo(object):
    '''
    This structure contains thermodynamics information obtained from a frequency calculation.
    '''

    def __init__(self,E,ZPE,ITE,T,p,S,H,F,G,frequencies,intensities,mass,mom_inertia,rot_sym,linear_switch):
        self.E = E
        self.ZPE = ZPE
        self.ITE = ITE
        self.T = T
        self.p = p
        self.S = S
        self.H = H
        self.F = F
        self.G = G
        self.frequencies = frequencies
        self.intensities = intensities
        self.mass = mass
        self.mom_inertia = mom_inertia
        self.rot_sym = rot_sym
        self._lin_switch = linear_switch

    def info(self):
        '''
        Prints a summary of basic thermodynamics information (in Hartree).
        The internal thermal energy (ITE) is corrected for zero point energy (ZPE). 
        Enthalpy (H), Helmholtz free energy (F) and Gibbs free energy (G) are 
        already corrected for internal thermal energy (which includes zero point energy). 
        '''
        print "Electronic energy (E):\t\t" + str(self.E)
        print "---------------------"
        print ""
        print "Corrections to E:"
        print "-----------------"
        print "Zero point energy (ZPE):\t" + str(self.ZPE)
        print "Internal thermal energy (ITE):\t" + str(self.ITE)
        print ""
        print "Temperature (T in Kelvin):\t" + str(self.T)
        print "-------------------------"
        print ""
        print "Pressure (p in Pa):\t\t" + str(self.p)
        print "------------------"
        print ""
        print "Thermodynamic potentials:"
        print "------------------------"
        print "Entropy (S):\t\t\t" + str(self.S)
        print "Enthalpy (H):\t\t\t" + str(self.H)
        print "Helmholtz free energy (F):\t" + str(self.F)
        print "Gibbs free energy (G):\t\t" + str(self.G)

    def calculate(self,temp,pressure="",loop_iso=0,loop_freq="",grimme=0,grimme_thresh=100):
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

        Output:
        A numpy array of dimension N*7, with N as number of temperature/pressure 
        pairs. The columns are: temperature, pressure, ITE (internal thermal energy), S, H, F and G. 

        Example:
        <somefile>.thermo.calculate(373.15) calculates the potentials at 100 deg C and 1 atm.
        '''
        if loop_freq=="":
            loop_freq = loop_iso

        # Check temperature and pressure
        if type(temp) is str or type(temp) is float or type(temp) is int:  
            dummy = [] # convert temperature input to list (will be iterated over)
            dummy.append(temp)
            T = np.asarray(dummy,dtype=float)
        else:
            T = np.asarray(temp,dtype=float)

        if pressure == "":
            p = np.ones(np.size(T))*1.01325e5  # from ATM to Pa
        else:
            if type(pressure) is str or type(pressure) is float or type(pressure) is int:  
                dummy = [] # convert pressure input to list
                for k in range(len(T)):  # assume same pressure for set of temperature
                    dummy.append(pressure)
                p = np.asarray(dummy,dtype=float) 
            else:
                p = np.asarray(pressure,dtype=float)

        if (np.size(T) != np.size(p)):
            print "Temperature and pressure input differs in length."
        else:
            # Get frequencies
            dummy = np.asarray(self.frequencies[loop_freq],dtype=float)
            # Remove negative entries
            freq = dummy[dummy>0]

        # Get mass (in amu) and moment of inertia (in amu * bohr**2)
        mass = self.mass[loop_iso]*constants.atomic_mass_constant
        mom_inertia = np.asarray(self.mom_inertia[loop_iso])*constants.atomic_mass_constant*(1e-10*constants.bohr_to_angstrom)**2  
        rot_sym = float(self.rot_sym[loop_iso])


        # Translational contributions to S and the thermal energy E
        pi = 3.141592654
        Q_trans = ((2*pi*mass*constants.Boltzmann_constant*T/(constants.Planck_constant)**2)**(3./2))*constants.Boltzmann_constant*T/p
        S_trans = constants.molar_gas_constant*(np.log(Q_trans) + 1 + 3./2) # in J/mol
        E_trans = 3./2*constants.molar_gas_constant*T # in J/mol

        # Rotational contributions (for linear molecule)
        if self._lin_switch==1:
            print "Linear molecule detected..."
            theta_lin = (constants.Planck_constant)**2/(8*pi**2*np.max(mom_inertia)*constants.Boltzmann_constant) # rotational temperature
            Q_rot_lin = T/theta_lin/rot_sym
            S_rot = constants.molar_gas_constant*(np.log(Q_rot_lin)+1)
            E_rot = constants.molar_gas_constant*T
        else:
            # Rotational contributions (for nonlinear molecule)
            theta = (constants.Planck_constant)**2/(8*pi**2*mom_inertia*constants.Boltzmann_constant) # rotational temperature
            Q_rot = np.sqrt(pi)/rot_sym*((T**(3./2))/np.sqrt(theta[0]*theta[1]*theta[2]))
            S_rot = constants.molar_gas_constant*(np.log(Q_rot)+3./2)
            E_rot = 3./2*constants.molar_gas_constant*T

        # Vibrational contributions 
        dummy = freq*constants.inverse_cm_to_hertz # conversion from wavenumbers to Hertz
        freq = dummy*constants.Planck_constant/constants.Boltzmann_constant  # conversion to Kelvin
        S_vib = []
        E_vib = []
        S_hind = []
        
        for k in T:
            dum_term = (freq/k)/(np.exp(freq/k)-1)-np.log(1-np.exp(-freq/k))
            S_vib.append(constants.molar_gas_constant*dum_term.sum())
            dum_term2 = (freq*(1.0/2 + 1.0/(np.exp(freq/k)-1)))
            E_vib.append(constants.molar_gas_constant*dum_term2.sum())

            if grimme==1:
                B_average = 10E-44; # kg.m^2
                alpha = 4;
                mu=constants.Planck_constant/(8*pi**2*constants.speed_of_light_in_vacuum*freq)
                mu_prime=mu*B_average/(mu+B_average)
                dum_grim = 0.5+np.log((8*pi**3*constants.Boltzmann_constant*k*mu_prime/constants.Planck_constant**2)**0.5);
                weight=1./(1+(grimme_thresh/freq)**alpha);
                S_hind.append(constants.molar_gas_constant*(weight*dum_term + (1-weight)*dum_grim).sum())

        if grimme==1:
            S_vib = S_hind # replace vibrational contribution with Grimme correction

        # Summmation over all contributions
        S = S_trans + S_rot + S_vib  # in J/(mol.K)
        S = S/constants.Avogadro_constant*constants.joule_to_hartree # in Hartree
        ITE = E_trans + E_rot + E_vib
        ITE = ITE/constants.Avogadro_constant*constants.joule_to_hartree # in Hartree

        # Calculate the thermodynamic potentials 
        H = self.E + ITE + constants.Boltzmann_constant*T*constants.joule_to_hartree # in Hartree
        F = self.E + ITE - T*S # array issue here
        G = H - T*S

        # Print and return array of thermodynamical potentials
        headerline = "Results based on frequencies of loop " + str(loop_freq) + " and isotopes of loop " + str(loop_iso)
        print headerline
        print "-"*len(headerline)
        print "\nT\tp\t\tITE\t\t S\t\t H\t\t F\t\t G"
        for k,l in enumerate(T):
            print "%.2f\t%.2f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t" % (T[k],p[k],ITE[k],S[k],H[k],F[k],G[k])        
        dummy = np.asarray([T,p,ITE,S,H,F,G])
        data = dummy.transpose()
        return data



class _opt(object):
    '''
    This structure contains information about the geometry optimization.
    '''
    def __init__(self,geometries,energies,gradient,displacement,change,optstat):
        self.geometries = geometries
        self.energies = energies
        self.gradient = gradient
        self.displacement = displacement
        self.change = change
        self.optstat = optstat

    def info(self):
        print "Summary of geometry optimization:"
        print "--------------------------------"
        print ""
        print "Energy\t\t\tGradient\tDisplacement\tDeltaE"
        for i,k in enumerate(self.energies):
            print str(k) + "\t\t" + str(self.gradient[i]) + "\t" + str(self.displacement[i]) + "\t" + str(self.change[i])


class _orbitals(object):
    '''
    This structure contains information (occupation, energies) about the orbitals.
    '''
    def __init__(self,N_elec):
        self.N_elec = N_elec

####################### MULTI OUTPUTFILE  ###########################

class _multioutput(object):

    def __init__(self, jobs=[]):
        self.list_of_jobs=[]
        self.list_of_content=[]
        for k in jobs:
            self.add(k)

    def add(self,new_job):
        self.list_of_jobs.append(new_job)
        self.list_of_content.append(new_job.general.jobtype)

    #def remove(self,position=0): #if not specified delete last
    #    del self.list_of_content[position-1] 
    #    del self.list_of_jobs[position-1] 

########################## OUTPUTFILE  ##############################

class _outputfile(object):

    def __init__(self,file_input):

        #Check input type
        if type(file_input)==list:
            content = file_input
        else:
            infile = open(file_input,"r")
            content = infile.readlines()
        
        spin = 'undetermined'
        energy = 'undetermined'
        jobtype = 'undetermined'
        version = 'undetermined'
        basis_size = 'undetermined'
        status = 'unfinished'
        
        switch = 0
        for line in content:
            if "JOBTYPE" in line:
                jobtype = ((line.split())[1]).lower()
            if "Q-Chem, Version" in line:
                version = (((line.split(","))[1]).split())[1]
            if "<S^2> =" in line:
                spin = (line.split())[2]
            if "Total energy in the final basis set" in line:
                energy = (line.split())[8]
            if ("There are" in line) and ("shells" in line):
                basis_size = (line.split())[5]
            if "MISSION" in line:
                status = 'finished'
            # Create corresponding inputfile:
            if "User input:" in line:
                switch = 1
                infile_content = []
            if switch == 1:
                infile_content.append(line)
            if switch == 1 and "Standard Nuclear Orientation" in line:
                switch = 0
        inputfile = _readinput(infile_content)


        self.general = _general(jobtype,version,spin,basis_size,energy,status,inputfile)


        if jobtype=="freq":
            H2kcal=constants.hartree_to_kcal_pro_mole
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
                    press.append(1.01325e5*float((line.split())[7]))
                if "Zero point vibrational energy" in line:
                    zero_point.append(float((line.split())[4])/H2kcal)
                if "gas constant (RT):" in line:
                    R_T.append(float((line.split())[3])/H2kcal)
                if "Total Enthalpy:" in line:
                    enth_corr.append(float((line.split())[2])/H2kcal)
                if "Total Entropy:" in line:
                    entr_corr.append(float((line.split())[2])/1000/H2kcal)
                if "We detect a D*h symmetry" in line:
                    linear_switch = 1
                if "We detect a C*v symmetry" in line: 
                    linear_switch = 1
                if "We detect a C*h symmetry" in line:  # necessary because of wrong nomenclature in Q-chem (thermodyn.F)
                    linear_switch = 1
                if "Molecular Mass:" in line:
                    loop += 1
                    if "*" in line:
                        print "Warning: Molecular mass in loop " + str(loop) + " is unphysically large. Will use mass of first loop instead."
                        mass.append(mass[0])
                    else:
                        dummy = float((line.split())[2])
                        mass.append(dummy)
                if "Rotational Symmetry Number is" in line:
                    rot_sym.append(float((line.split())[4]))
                if "Eigenvalues --" in line:
                    if "*" in line:
                        print "Warning: Moment of inertia in loop " + str(loop) + " is unphysically large. Will use values of first loop instead."
                        mom_inertia.append(mom_inertia[0])
                    else:
                        dummy = [float((line.split())[2]),float((line.split())[3]),float((line.split())[4])]
                        mom_inertia.append(dummy)

                if "Vibman" in line:
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

            T = np.asarray(temp)
            p = np.asarray(press)
            ZPE = np.asarray(zero_point)
            H = E + np.asarray(enth_corr)
            S = np.asarray(entr_corr)
            ITE = np.asarray(enth_corr) - R_T # RT = pV is subtracted from H to obtain the ZPE corrected ITE
            F = E + ITE - (T*S)
            G = H - (T*S)
            self.thermo = _thermo(E,ZPE,ITE,T,p,S,H,F,G,frequencies,intensities,mass,mom_inertia,rot_sym,linear_switch)

        if jobtype=="opt":
            energies = []
            gradient = []
            displacement = []
            change = []
            geometries = []
            optstat = "no convergence"
            N_step = 1
            switch = 0
            for line in content:
                if "Energy is" in line:
                    dummy = (line.split())[2]
                    energies.append(dummy)
                if "Gradient   " in line:
                    dummy = (line.split())[1]
                    gradient.append(dummy)
                if "Displacement   " in line:
                    dummy = (line.split())[1]
                    displacement.append(dummy)
                if "Energy change   " in line:
                    dummy = (line.split())[2]
                    change.append(dummy)
                if "**  OPTIMIZATION CONVERGED  **" in line:
                    optstat = "converged"
                if "ATOM              X           Y           Z" in line:
                    switch = 1
                    cycle_name = "Optimization step " + str(N_step)
                    cart_dummy = cartesian(cycle_name)
                if switch == 1 and "ATOM" not in line and "Point Group" not in line:
                    con = line.split()
                    cart_dummy.add_atom(con[1],con[2],con[3],con[4])
                if "Point Group" in line and switch == 1:
                    geometries.append(deepcopy(cart_dummy))
                    N_step += 1
                    switch = 0


            self.opt = _opt(geometries,energies,gradient,displacement,change,optstat)


