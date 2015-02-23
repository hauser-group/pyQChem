#####################################################################
#                                                                   #
#            pyQchem - ADC Classes and Parsing Functions            #
#                                                                   #
#####################################################################

# All classes and functions here are all hidden. Their only purpose is to 
# structure information for the user.

######################### ADC INFO CLASSES ########################## 

class _state(object):
    '''
    Information for any state computed using ADC
    Attributes:
     - term_symbol - Term symbol of state 
     - total_energy - Absolute energy of the state (in a.u.)
     - dict_of_properties - Dictionary of properties calculated
    '''
    
    def __init__(self, ts = '', te = 0.0, prop = {}):
        self.term_symbol = ts
        self.total_energy = te
        self.dict_of_properties = prop
        
    def info(self):
        '''
        Print an overview of the state information 
        '''
        if self.term_symbol != "":
            print("  Term symbol: " + self.term_symbol)
        print("  Energy = " + repr(self.total_energy))
        for name,value in self.dict_of_properties.items():
            print("  Property: " + name + " = " + repr(value))
        
        
class _transition(object):
    '''
    Information on any transition between two states computed using ADC
    Attributes:
     - energy - Excitation energy (in eV)
     - osc_strength - Oscillator strength
     - dict_of_properties - Dictionary of transition properties
    '''
    
    def __init__(self, omega = 0.0, osc = 0.0, prop = {}):
        '''
        Initialise the transition information
        '''
        self.energy = omega
        self.osc_strength = osc
        self.dict_of_properties = prop
        
    def info(self):
        '''
        Print an overview of the transition information 
        '''
        print("  Excitation energy = " + repr(self.energy) +
            " (" + repr(self.osc_strength) + ")");
        for name,value in self.dict_of_properties.items():
            print("  Property " + name + " = " + repr(value))


class _amplitude(object):
    def __init__(self, name, value):
        self.excitation = name
        self.value = value
                

class _ground_state(_state):
    '''
    Ground state computed using ADC (derived from _state) 
    Additional attributes:
     - dict_of_econtrib - Dictionary of energy contributions (in a.u.) 
    '''
    
    def __init__(self, ts = '', te = 0.0, prop = {}, contrib = {}):
        _state.__init__(self, ts, te, prop)
        if (len(contrib) == 0):
            self.dict_of_econtrib = { 'hf': te }
        else:
            self.dict_of_econtrib = contrib

    def info(self):
        '''
        Prints an overview of the ground state information
        '''
        _state.info(self)
        print("  Energy composition:")
        for name,value in self.dict_of_econtrib.items():
            print("    E(" + name + ") = " + repr(value))


class _excited_state(_state, _transition):
    '''
    Information on any excited state computed using ADC (derived form _state 
    and _transition). The transition information refers to the transition 
    from the ground into the excited state.
    Additional attributes:
     - converged - Boolean variable indicating if the calculation of the state converged
     - rnorm - Square of the residual norm
     - vnorm - Vector of the square norms of the excited state components
     - list_of_amplitudes - Vector of important amplitudes in the excited state 
    '''
    
    def __init__(self, ts = '', te = 0.0, omega = 0.0, osc = 0.0, prop = {}, 
            trprop = {}, conv = True, rnorm = 0.0, vnorm = []):
        _state.__init__(self, ts, te, prop)
        _transition.__init__(self, omega, osc, trprop)
        self.converged = conv
        self.rnorm = rnorm 
        self.vnorm = vnorm
        self.list_of_amplitudes = []
        
    def add_amplitude(self, name, value):
        self.list_of_amplitudes.append(_amplitude(name, value))
    
    def info(self):
        '''
        Prints an overview of the ground state information
        '''
        _state.info(self)
        _transition.info(self)
        s = "(R^2 = " + str(self.rnorm) + ")"
        if self.converged:
            s = "  Converged " + s
        else:
            s = "  Not converged " + s
        print(s)
        print("  Components V^2 = " + str(self.vnorm))
        if len(self.list_of_amplitudes) > 0:
            print("  Imporant amplitudes")
            for a in self.list_of_amplitudes:
                print("    " + a.excitation + ": " + str(a.value))
        

class _adc(object):
    '''
    Complete information on ADC ground and excited states at one geometry
    Attributes:
     - adc_variant - Variant of ADC used in the calculation
     - ground_state - Ground state data 
     - list_of_excited_states - List of excited state data
     - dict_of_transitions - List of state-to-state transitions (only for excited states)
    '''

    def __init__(self, variant, gs = _ground_state(), es = [], tr = {}):
        '''
        Initialize the data
        '''
        self.adc_variant = variant
        self.ground_state = gs
        self.list_of_excited_states = es
        self.dict_of_transitions = tr

    def info(self):
        '''
        Prints an overview of the ADC calculation.
        '''
        print("Ground state:")
        self.ground_state.info()
        print("")
        for i,state in enumerate(self.list_of_excited_states):
            print("Excited state " + str(i + 1) + ":")
            state.info()
            print("")
        print("")
        if isinstance(self.dict_of_transitions, dict):
            for i,trlist in self.dict_of_transitions.items():
                if not isinstance(trlist, dict): 
                    continue
                for j,tr in trlist.items():
                    print("Transition from state " + str(i + 1) + " to state " 
                        + str(j + 1) + ":")
                    tr.info()

####################### ADC PARSER FUNCTIONS ######################## 

def _parse_adc(variant, content, silent=False):
    '''
    Parses a Q-Chem job output given as array of content lines and returns 
    an _adc object. If no ADC section is found, zero is returned
    '''

    # First find the ADC section in the output
    for i in range(0,len(content)):
        if "A D C  M A N" in content[i]:
            break
    # if no ADC section is found            
    if i == len(content):
        return 0;    

    # Search the end of the ADC header
    ref_data = []
    es_data = []
    tr_data = {}
    for i in range(i + 1, len(content)):
        if content[i].find("=====") == 0: break
                 
    i += 1
    # Now search the end of the ADC section and delegate extraction to other 
    # functions
    while i < len(content):
        if content[i].find("=====") == 0: break

        if "Excited State Summary" in content[i]:
            i,es_data = _parse_es_summary(content, i + 2, silent)
        elif "Transition Summary" in content[i]:
            i,tr_data = _parse_tr_summary(content, i + 2, silent)
        elif "Summary" in content[i]:
            if "HF" in content[i] or "MP" in content[i]:
                name = content[i].strip().rstrip(" Summary")
                i,summary = _parse_ref_summary(content, i + 2, silent)
                summary['name'] = name
                ref_data.append(summary)
        i += 1
    
    data = _adc(variant, _build_gs('', ref_data), es_data, tr_data)
    
    return data
            
        
def _parse_ref_summary(content, start, silent=False):
    '''
    Parse an ADC reference state summary and extract information 
    '''

    data = { 'total_energy': 0.0, 'energy_contrib': 0.0, 'properties': {} }
    for i in range(start,len(content)):
        if content[i].find("-----") == 0: 
            break
        if "Energy:" in content[i]:
            data['total_energy'] = float(content[i].split()[-2])
            data['energy_contrib'] = float(content[i].split()[-2])
        elif "Total energy:" in content[i]:
            data['total_energy'] = float(content[i].split()[-2])
        elif "Energy contribution:" in content[i]:
            data['energy_contrib'] = float(content[i].split()[-2])
        elif "Total dipole" in content[i]:
            data['properties']['dipole'] = float(content[i].split()[-1])
        elif "Total <r^2>" in content[i]:
            data['properties']['rsq'] = float(content[i].split()[-1])
        
    return i,data

    
def _build_gs(ts, data):
    '''
    Construct the ground state from the data provided
    '''
    
    te = 0.0
    prop = {}
    contrib = {}
    for state in data: 
        te = state['total_energy']
        contrib[state['name']] = state['energy_contrib']
        for key,val in state['properties'].items():
            prop[key] = val

    return _ground_state(ts, te, prop, contrib)


def _parse_es(content, start, silent=False):
    '''
    Parse the information of a single ADC excited state
    '''

    ts=''
    te = 0.0
    ee = 0.0
    osc = 0.0
    prop = {} 
    tprop = {}
    conv=False
    rnorm=0.0
    vnorm=[1.0]
    sep=1
    
    for i in range(start,len(content)):
        if "Excited state" in content[i]:
            if sep == 0:
                break
            if content[i].split()[-1] == "[converged]":
                conv = True
            sep -= 1
        elif content[i].find("-----") == 0:
            break
        elif "Term symbol" in content[i]:
            list = content[i].split()
            ts = ' '.join(list[2:5])
            rnorm = float(list[-1])                        
        elif "Total energy:" in content[i]:
            te = float(content[i].split()[-2])
        elif "Excitation energy" in content[i]:
            ee = float(content[i].split()[-2])
        elif "Osc. strength" in content[i]:
            osc = float(content[i].split()[-1])
        elif "EFP" in content[i]:
            prop['efp'] = content[i].split()[-2]
        elif "PCM" in content[i]:
            list = content[i].split()
            if list[1] == "SS":
                prop['pcm'] = float(list[-2])
            elif list[1] == "LR":
                tprop['pcm'] = float(list[-2])
        elif "Trans. dip." in content[i]:
            tdip = content[i].strip().split(':')[1].strip('[] ').split(',')
            tprop['dipole'] = []
            for comp in tdip:
                tprop['dipole'].append(float(comp))
        elif "<i|r^2|0>" in content[i]:
            rsq = content[i].strip().split(':')[1].strip('[] ').split(',')
            prop['rsq'] = []
            for comp in rsq:
                prop['rsq'].append(float(comp))
        elif "Two-photon absorption cross-section" in content[i]:
            tprop['tpa'] = float(content[i].split()[-1])
        elif "Total dipole" in content[i]:
            prop['dipole'] = float(content[i].split()[-1])
        elif "Total <r^2>" in content[i]:
            prop['rsq'] = float(content[i].split()[-1])
        elif "V1^2" in content[i]:
            list = content[i].split(',')
            vnorm = []
            for vcomp in list:
                vnorm.append(float(vcomp.split('=')[1]))
        elif "Important amplitudes:" in content[i]:
            arange = []
            for j in range(i + 1, len(content)):
                if "----" in content[j]:
                    arange.append(j)
                if len(arange) == 2:
                    break
            i = j + 1

    state = _excited_state(ts, te, ee, osc, prop, tprop, conv, rnorm, vnorm)
    if len(arange) == 2:
        for j in range(arange[0] + 1, arange[1] - 1):
            list = content[j].replace("(","").replace(")","").split()
            value = list[-1]
            name = []
            k = 3
            while k < len(list):
                name.append(' '.join(list[k-3:k]))
                k += 3
            list=[', '.join(name[0:len(name)/2]), ', '.join(name[len(name)/2:len(name)])]
            state.add_amplitude(' -> '.join(list), value)

    return i - 1,state


def _parse_es_summary(content, start, silent=False):
    '''
    Parse the ADC excited summary section
    '''
    
    states = []
    i = start
    while (i < len(content)):
        if content[i].find("-----") == 0: break
        if "Excited state" in content[i]:
            i,state = _parse_es(content, i, silent)
            states.append(state)
        i += 1

    return i,states 
           

def _parse_tr(content, start, silent=False):
    '''
    Parses the transition data of a single ADC state-to-state transition
    '''
    
    ee = 0.0
    osc = 0.0
    prop = {}
    sep = 1
    
    for i in range(start, len(content)):
        if "Transition from" in content[i]:
            if sep == 0:
                break
            sep -= 1
        elif "Excitation energy" in content[i]:
            ee = float(content[i].split()[-2])
        elif "Osc. strength" in content[i]:
            osc = float(content[i].split()[-1])
        elif "Trans. dip." in content[i]:
            tdip = content[i].strip().split(':')[1].strip('[] ').split(',')
            prop['tdip'] = []
            for comp in tdip:
                prop['tdip'].append(float(comp))
        elif "<i|r^2|0>" in content[i]:
            rsq = content[i].strip().split(':')[1].strip('[] ').split(',')
            prop['rsq'] = []
            for comp in rsq:
                prop['rsq'].append(float(comp))
    
    return i - 1,_transition(ee, osc, prop)


def _parse_tr_summary(content, start, silent=False):
    '''
    Parses the ADC transition summary section
    '''
    
    trlist = {}
    i = start
    while (i < len(content)):
        if content[i].find("-----") == 0: break
        if "Transition from" in content[i]:
            list = content[i].strip().rstrip(':').split()
            index = [int(list[4]), int(list[8])]
            if not index[0] in trlist:
                trlist[index[0]] = {}
            i,tr = _parse_tr(content, i, silent)
            trlist[index[0]][index[1]] = tr
        i += 1

    return i,trlist 

           
