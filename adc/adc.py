#####################################################################
#                                                                   #
#                   pyQchem - ADC Output Classes                    #
#                                                                   #
#####################################################################

# All classes here are all hidden. Their only purpose is to structure 
# information for the user.


class _state_data(object):
    '''
    State data computed for some state using ADC
    Attributes:
     - term_symbol - Term symbol of state 
     - total_energy - Absolute energy of the state (in a.u.)
     - properties - Dictionary of properties calculated
    '''
    
    def __init__(self, ts = '', te = 0.0, prop = {}):
        self.term_symbol = ts
        self.total_energy = te
        self.properties = prop
        
    def info(self)
        '''
        Print an overview of the state information 
        '''
        s = ts + ': E = ' + repr(self.total_energy);
        for name,value in self.properties.items():
            s += ', ' + name + " = " + repr(value)
        print(s)
        
class _transition_data(object):
    '''
    Transition data computed for some transition between states using ADC
    Attributes:
     - energy - Excitation energy (in eV)
     - osc_strength - Oscillator strength
     - properties - Dictionary of transition properties
    '''
    
    def __init__(self, omega = 0.0, osc = 0.0, prop = {}):
        '''
        Initialise the transition information
        '''
        self.energy = omega
        self.osc_strength = osc
        self.properties = prop
        
    def info(self)
        '''
        Print an overview of the state information 
        '''
        s = 'omega = ' + repr(self.energy)
        s += " (" + repr(self.osc_strength) + ")";
        for name,value in self.properties.items():
            s += ', ' + name + " = " + repr(value)
        print(s)

class _ground_state(_state_data):
    '''
    Ground state data obtained by ADC. 
    Additional attributes:
     - energy_contrib - Dictionary of energy contributions (in a.u.) 
    '''
    
    def __init__(self, ts = '', te = 0.0, prop = {}, contrib = {}) :
        _state_data.__init__(ts, te, prop)
        if (len(contrib) == 0) 
            self.energy_contrib = { 'hf': te }
        else  
            self.energy_contrib = contrib

    def info(self)
        '''
        Prints an overview of the ground state information
        '''
        _state_data.info()
        print("Energy composition:")
        for name,value self.energy_contrib.items():
            print(" E(" + name + ") = " + repr(value))


class _excited_state(_state_data, _transition_data):
    '''
    Excited state data obtained by ADC
    '''
    
    def __init(self, ts = '', te = 0.0, omega = 0.0, prop = {}, trprop = {}, 
            conv = true, rnorm = 0.0, vnorm = []):
        _state_data.__init__(ts, te, prop)
        _transition_data.__init__(omega, trprop)
        self.converged = conv
        self.rnorm = rnorm 
        self.vnorm = vnorm
        self.amplitudes = { }
        
    def info(self)
        '''
        Prints an overview of the ground state information
        '''
        _state_data.info()
        _transition_data.info()

class _adc(object):
    '''
    Complete information on ADC ground and excited states. 
    Attributes:
     - adc_variant - Variant of ADC used in the calculation
     - ground_state - Ground state data 
     - excited_states - List of excited state data
     - transitions - List of state-to-state transitions (only for excited states)
    '''

    def __init__(self, variant, gs = _ground_state(), eslist = [],trlist = {}):
        '''
        Initializes 
        '''
        self.adc_variant = variant
        self.ground_state = gs
        self.excited_states = eslist
        self.transitions = trlist

    def info(self):
        '''
        Prints an overview of the ADC calculation.
        '''
        print("ADC calculation summary")
        print("-----------------------")
        print("")
        print(self.ground_state)
        print("")
        for i,state in enumerate(self.excited_states):
            print("Excited state " + str(i + 1) + ":")
            print(state)
        print("")
        for i,trlist in self.transitions.items():
            if ! isinstance(trlist, dict): 
                continue
            for j,tr in trlist.items():
                print("Transition from state " + str(i + 1) + " to state " 
                    + str(j + 1) + ":")
                print(tr)
                