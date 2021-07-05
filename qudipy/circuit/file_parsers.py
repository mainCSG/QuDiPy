"""
Functions for handling control pulse files .ctrlp

@author: simba
"""
import os
import pandas as pd
import numpy as np
from .quantum_circuit import QuantumCircuit
from .control_pulse import ControlPulse

def _load_one_pulse(f_name):
    '''
    This function takes a single .ctrlp file and extracts the pulse information.

    Parameters
    ----------
    f_name : string
        Full file path to the control pulse file to load.

    Returns
    -------
    ctrl_pulse : ControlPulse object
        Control pulse object containing all information loaded from the file.

    '''
        
    # Check file extension
    if f_name[-6::] != ".ctrlp":
        raise ValueError("Unrecognized control pulse file type..." +
                         " Expected .ctrlp file.")
    
    # Open file and extract pulse name
    f = open(f_name, "r")
    pulse_name = os.path.basename(f.name).split('.')[0]
    
    # Count the number of lines which don't have control pulse vs time information
    line_cnt = 0
    # Initialize ideal gate to None object for if it isn't specified in the
    # pulse file
    ideal_gate = None
    for x in f:
        # Keep track of how many lines to skip for reading in the pulse table
        line_cnt += 1
        # Parse line by line (don't care about ordering except that Control 
        # pulses is the last line)
        if "control pulses:" in x.lower():
            break
        elif "pulse type:" in x.lower():
            # Remove any leading/ending white space and convert to lowercase
            pulse_type = x.split(":")[1].strip().lower()
        elif "pulse length:" in x.lower():
            pulse_length = int(x.split(":")[1].strip().split(" ")[0])
        elif "ideal gate:" in x.lower():
            ideal_gate = x.split(":")[1].strip().upper()
            
    # Initialize the control pulse object
    ctrl_pulse = ControlPulse(pulse_name, pulse_type, pulse_length, ideal_gate)
    
    # Read the rest of the pulse file to get the pulse information
    df = pd.read_csv(f_name, skiprows = line_cnt)
    
    # Loop over each column in the pulse table and add control variables to 
    # control pulse object
    ctrl_var_names = df.columns
    for ctrl_var in ctrl_var_names:
        ctrl_pulse.add_control_variable(ctrl_var, df[ctrl_var])
    
    return ctrl_pulse

def load_pulses(f_names):
    '''
    This function takes multiple .ctrlp files as inputs and constructs a
    dictionary containing every control pulse's information

    Parameters
    ----------
    f_names : list of strings
        Full file path to the control pulse file to load.

    Returns
    -------
    pulse_dict : dictionary of ControlPulse objects
        Dictionary containing many control pulse objects
        
    '''
    
    pulse_dict = {}
    
    # Check if list of file names was passed.
    # If it was not, then only a single file name was passed and we need to
    # wrap it in a list container.
    if not isinstance(f_names, list):
        f_names = [f_names]
    
    # Loop through each file, load the pulse then add to the pulse dictionary.
    for f in f_names:
        curr_pulse = _load_one_pulse(f)
        
        pulse_dict[curr_pulse.name] = curr_pulse
        
    return pulse_dict
      
def load_circuit(f_name, gate_dict=None):
    '''
    This function takes in a single quantum circuit .qcirc file as input and 
    constructs a quantum circuit object

    Parameters
    ----------
    f_name : string
        Full file path to the circuit file.
        
    Keyword Arguments
    -----------------
    gate_dict : dictionary of controlPulse objects, optional
        Contains all loaded control pulse objects to be used in the 
        quantumCircuit's circuit sequence. Default is {}.

    Returns
    -------
    q_circ : QuantumCircuit object
        Class containing information loaded from the .qcirc file 
        required for simulation.

    '''
    # If a dictionary was not passed to the function, assign gate_dict as an empty dictionary
    if not gate_dict:
        gate_dict = {}

    # Check file extension
    if f_name[-6::] != ".qcirc":
        raise ValueError("Unrecognized quantum circuit file type..." +
                         " Expected .qcirc file.")
    
    # Open the file and extract the circuit name
    f = open(f_name, "r")
    circuit_name = os.path.basename(f.name).split('.')[0]

    # Loop over every gate in the circuit file
    for x in f:
        # Parse line by line
        if "Number of qubits:" in x:
            n_qubits = int(x.split(":")[1])
            
            # Initialize the quantum circuit object
            q_circ = QuantumCircuit(circuit_name, n_qubits, gate_dict)       
        else:
            # Parse the line information
            temp = x.strip().split(" ")
            
            gate_name = temp[0]
            
            # Track if the current line in the file is an ideal gate or a
            # pulse file that was loaded
            is_ideal_gate = check_ideal_gate(gate_name)
            
            # Check that the gate name is one that was actually loaded or a
            # valid ideal gate and print an error if not
            if gate_name not in q_circ.gates.keys() and not is_ideal_gate:
                raise ValueError("Problem loading circuit file: " +
                                 f"{circuit_name}.\n" +
                                 f"Gate {gate_name} could not be loaded as " +
                                 "the corresponding pulse was not loaded or " + 
                                 "is not the gate name a ideal gate keyword.\n" +
                                 "Check .qcirc file for typos or load the " +
                                 "corresponding pulse file.")
            
            # If it's not an ideal gate, get the ideal_gate equivalent from
            # the pulse object dictionary in quantumCircuit object 
            if not is_ideal_gate:
                # Update tracker that checks if .qcirc is all ideal gates or not
                q_circ.ideal_circuit = False
                
                # Get the corresponding ideal gate
                ideal_gate = q_circ.gates[gate_name].ideal_gate
                # Check that it's a valid ideal gate keyword and that the 
                # ideal_gate isn't None (i.e. not specified in the pulse file)
                if not check_ideal_gate(ideal_gate) and ideal_gate:
                    raise ValueError("Problem loading circuit file: " +
                                     f"{circuit_name}.\n" +
                                     f"Gate {gate_name}.ctrlp could not be loaded " +
                                     "as the ideal gate keyword was not recognized.")
            # If it is an ideal gate
            elif is_ideal_gate:
                ideal_gate = gate_name
                gate_name = "IDEAL_GATE"
            
            gate_aff_qubits = temp[1:]
            # int() always rounds down so as a precaution we convert to float,
            # then round, then convert to int.
            gate_aff_qubits = [int(round(float(qubit_idx))) 
                               for qubit_idx in gate_aff_qubits]
            q_circ.add_gate(gate_name, ideal_gate, gate_aff_qubits)

    # Check now if every gate in the circuit file loaded had a ideal gate 
    # correctly specified or not.  If not, then warn user that we cannot print
    # the ideal circuit nor do simulations of the ideal circuit as well.
    if not q_circ.specified_all_ideal:
        print(f'WARNING: Problem loading circuit file: {circuit_name}.')
        print('Not every gate loaded had a corresponding ideal specified')
        print('correctly or specified at all in the file. This will prevent')
        print(f'you from simulating the ideal quantum circuit for {circuit_name}.\n')
    
    return q_circ    
    
def check_ideal_gate(gate_name, qubit_idx=None):
    '''
    This function checks if the supplied gate_name is a valid ideal gate
    keyword used to simulate an ideal quantum circuit. 
    Current supported keywords are
    I, RX###, RY###, RZ###, H, CTRLX, CTRLY, CTRLZ, SWAP, RSWAP
    where ### in the R gates indicates the gate's rotation angle in degrees.

    Parameters
    ----------
    gate_name : string
        Gate name to be tested.
        
    Keyword Arguments
    -----------------
    qubit_idx : list of ints, optional
        Indices of qubits used by gate. Default is [].

    Returns
    -------
    boolean

    '''

    # If qubit_idx was not passed to function, assign an empty list
    if not qubit_idx:
        qubit_idx = []
    
    # If gate name is None type then that means there was no ideal gate line
    # specified in the corresponding pulse file
    # Also, quick check by looking at gate name length
    if gate_name == None or len(gate_name) not in [1,4,5]:
        return False
    
    # Check for an R gate first
    if gate_name[0] == "R" and gate_name[1] in ["X","Y","Z"]:
        # Now check that the next three characters are ints
        for idx in range(2,5):
            try:
                int(gate_name[idx])
            except (ValueError, IndexError):
                return False
        return True
    
    # Check I, H, and CTRL gates, and (if CTRL, SWAP, or RSWAP gate)
    # check qubit number (must be an even number of used qubits)
    if gate_name in ("I", "H") \
            or (len(gate_name) == 5 \
            and gate_name[:4] == "CTRL" \
            and gate_name[4] in ["Z","X","Y"] \
            or gate_name in ("SWAP", "RSWAP")) \
            and np.mod(len(qubit_idx),2) == 0:
        return True
    
    # Otherwise
    return False

    
    

