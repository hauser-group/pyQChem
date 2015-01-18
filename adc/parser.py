#####################################################################
#                                                                   #
#                    pyQchem - ADC output parser                    #
#                                                                   #
#####################################################################

# All classes here are all hidden. Their only purpose is to structure 
# information for the user.


def _parse_output(variant, content, silent=False):
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
    ref_data = {}
    es_data = []
    tr_data = []
    for i in range(i + 1, len(content)):
        if content[i].find("=====") == 0: break
                 
    i++
    # Now search the end of the ADC section and delegate extraction to other 
    # functions
    while i < len(content):
        if content[i].find("=====") == 0: break

        if "Excited State Summary" in content[i]:
            i,es_data = _parse_es_summary(content, i + 2, silent)
        else if "Transition Summary" in content[i]:
            i,tr_data = _parse_tr_summary(content, i + 2, silent)
        else if "Summary" in content[i]:
            name = content[i].lstrip().rstrip(" Summary")
            i,ref_data[name] = _parse_ref_summary(content, i + 2, silent)
        else
            i++
    
    data = _adc(variant, _build_gs('', ref_data), es_data, tr_data)
    
    return i,data
            
        
def _parse_ref_summary(content, start, silent=False):
    '''
    Parses an ADC reference state summary and extracts all important information 
    '''

    data = { 'total_energy': 0.0, 'energy_contrib': {}, 'properties': {} }
    for i in range(start,len(content)):
        if "-----------------------------------" in content[i]:
            break
        if "Energy:" in content[i]:
            data['total_energy'] = content[i].split()[2]
        else if "Total energy:" in content[i]:
            data['total_energy'] = content[i].split()[3]
        else if "contribution" in content[i]:
            data['energy_contrib']['default'] = content[i].split()[3]
        else if "EFP" in content[i]:
            data['energy_contrib']['efp'] = content[i].split()[4]
        else if "PCM" in content[i]:
            data['energy_contrib']['pcm'] = content[i].split()[4]
        else if "Total dipole" in content[i]:
            data['properties']['dipole'] = content[i].split()[4]
        else if "Total <r^2>" in content[i]
            data['properties']['rsq'] = content[i].split()[4]
        
    return i,data
    
def _build_gs(ts, data):
    gs = _ground_state(ts, )

def _parse_ref_summary(content, start, silent=False):
    '''
    Parses an ADC reference state summary and extracts all important information 
    '''

    data = { 'total_energy': 0.0, 'energy_contrib': {}, 'properties': {} }
    for i in range(start,len(content)):
        if "-----------------------------------" in content[i]:
            break
        if "Energy:" in content[i]:
            data['total_energy'] = content[i].split()[2]
        else if "Total energy:" in content[i]:
            data['total_energy'] = content[i].split()[3]
        else if "contribution" in content[i]:
            data['energy_contrib']['default'] = content[i].split()[3]
        else if "EFP" in content[i]:
            data['energy_contrib']['efp'] = content[i].split()[4]
        else if "PCM" in content[i]:
            data['energy_contrib']['pcm'] = content[i].split()[4]
        else if "Total dipole" in content[i]:
            data['properties']['dipole'] = content[i].split()[4]
        else if "Total <r^2>" in content[i]
            data['properties']['rsq'] = content[i].split()[4]
        
    return i,data

def _parse_es_summary(content, start, silent=False):
    '''
    Parses the ADC excited summary section
    '''
    
    state_list = []
    i = start
    while (i < len(content)):
        if content[i].find("----") == 0: break
        if "Excited state" in content[i]:
            i,state = _parse_es(content, i, silent)
            state_list.append(state)
        else 
            i++

    return i,state_list 
           
def _parse_es(content, start, silent=False):
    '''
    Parses the state data of a single ADC excited state
    '''
    
    return start + 2,_excited_state()

def _parse_tr_summary(content, start, silent=False):
    '''
    Parses the ADC transition summary section
    '''
    
    tr_list = []
    i = start
    while (i < len(content)):
        if content[i].find("----") == 0: break
        if "Transition from" in content[i]:
            i,tr = _parse_es(content, i, silent)
            tr_list.append(tr)
        else 
            i++

    return i,state_list 
           
def _parse_tr(content, start, silent=False):
    '''
    Parses the transition data of a single ADC state-to-state transition
    '''
    
    return start + 2,_transition_data()

