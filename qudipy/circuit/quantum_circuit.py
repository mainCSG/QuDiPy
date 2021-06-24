"""
Class for a quantum circuit

@author: simba
"""

class QuantumCircuit:
    
    def __init__(self, circuit_name, n_qubits, pulse_dict):
        '''
        Initialize QuantumCircuit object

        Parameters
        ----------
        circuit_name : string
            Name of quantum circuit.
        n_qubits : int
            Number of qubits in the circuit.
        pulse_dict : dictionary of ControlPulse objects
            Dictionary containing all control pulses that will be used in the
            quantum circuit.

        Returns
        -------
        None.

        '''
        
        # Name of circuit file
        self.name = circuit_name
        # Number of qubits in the circuit
        self.n_qubits = n_qubits
        # Loaded gates 
        self.gates = pulse_dict
        
        # The circuit sequence for this quantum circuit
        self.circuit_sequence = []
        
        # Index to track which gate in sequence we are on
        self.curr_gate_idx = 0
        
        # Flag to determine is every gate in the circuit has a correctly
        # specified ideal gate
        self.specified_all_ideal = True 
        # Flag to determine if the .qirc file that was loaded is comprised 
        # ONLY of ideal gates. Default assumes it is
        self.ideal_circuit = True
        
    def reset_circuit_index(self):
       '''
       Reset the current gate index for the circuit sequence back to the 
       begining of the circuit.

       Returns
       -------
       None.

       '''
       
       self.curr_gate_idx = 0
        
    def add_gate(self, gate_name, ideal_gate, used_qubits):
        '''
        This function adds a gate into the quantum circuit sequence.

        Parameters
        ----------
        gate_name : string
            Name of gate being added to sequence.
        used_qubits : int list
            Indices of qubits acted on by the gate.

        Returns
        -------
        None.

        '''
        
        # Check if a single qubit index was loaded or if it was a list.
        # If not a list, make it one.
        if not isinstance(used_qubits, list):
            used_qubits = [used_qubits]
        
        # Make sure used_qubits contains only ints
        used_qubits = [int(qubit_idx) for qubit_idx in used_qubits]
            
        # Check that the gate we are adding does not have an invalid qubit
        # index (i.e. outside of the allowable values)
        if (any(qubit_idx > self.n_qubits for qubit_idx in used_qubits) or
            any(qubit_idx < 1 for qubit_idx in used_qubits)):
            raise ValueError("Problem loading circuit file: " +
                            f"{self.name}.\n" +
                            f"Gate {gate_name} could not be loaded as the " +
                            f"affected qubit indices {used_qubits}\n are " +
                            " greater than the number of qubits in the circuit " +
                            f"({self.n_qubits}) or is <= 0.\n")
            
        # If the ideal_gate is None type, then there was an issue reading the
        # gate when the .ctrlp file was loaded
        if ideal_gate is None:
            self.specified_all_ideal = False
        
        # Add the gate to the circuit sequence
        self.circuit_sequence.append([gate_name, ideal_gate, used_qubits])
        
    def get_next_gate(self):
        '''
        Get the next gate in the quantum circuit sequence. If no more exist, 
        return None.

        Returns
        -------
        next_gate : [string, int list]
            Returns a list containing the next gate name and the affected 
            qubits by the gate.

        '''
        
        try:
            next_gate = self.circuit_sequence[self.curr_gate_idx]
            self.curr_gate_idx += 1
        except IndexError:
            next_gate = None
        
        return next_gate
    
    def load_more_gates(self, pulse_dict):
        '''
        Adds more controlPulse objects to the gate dictionary.

        Parameters
        ----------
        pulse_dict : dict of controlPulse objects
            Dictionary containing all the pulse objects to be added.

        Returns
        -------
        None.

        '''
        
        for pulse_key, pulse_value in pulse_dict.items():
            self.gates[pulse_key] = pulse_value
        
    
    def print_ideal_circuit(self):
        '''
        Prints out an ascii display of the loaded circuit sequence for the user.
    
        Parameters
        ----------
        None.
    
        Returns
        -------
        None.
    
        '''
        
        # Check that every gate has an ideal gate specified
        if not self.specified_all_ideal:
            print(f'Cannot print ideal circuit for {self.name}.')
            print('Some or all of the gates in the circuit do not have an')
            print('ideal gate specified. Please check the .ctrlp or .qcirc')
            print('files for errors.')
            return
        
        
        # Initialize the circuit to print
        circ_str = []
        for idx in range(self.n_qubits):
            circ_str.append('Q' + str(idx+1) + ' --')
            if idx != self.n_qubits:
                if idx < 10:
                    circ_str.append('     ')
                else:
                    circ_str.append('      ')
        
        # Each odd idx in circ_str corresponds to a qubit in the circuit
        # Each even idx correspond to gaps between qubit lines.
        
        # Store the current gate index to change back to later
        initial_gate_idx = self.curr_gate_idx
        # Now reset index
        self.curr_gate_idx = 0
        
        # Now loop through each gate in circuit and add the respective strings
        curr_gate = 0
        gate_flag = -1
        curr_gate = self.get_next_gate()
        while curr_gate is not None:
            
            # Extract ideal gate and affected qubits
            ideal_gate = curr_gate[1]
            aff_qubits = curr_gate[2]
            
            # Build the strings for a qubit affected by gate, nont affected by 
            # gate, and empty space between qubit lines
            if ideal_gate in ['H','I']:
                gate_flag = 1
                used_str = ideal_gate + '-'
                non_used_str = '--'
                empty_space = '  '
                
            if ideal_gate[:2] in ['RX','RY','RZ']:
                gate_flag = 1
                used_str = ideal_gate + '-'
                non_used_str = '------'
                empty_space = '      '
                                                
            # Now append respective strings as appropriate for each qubit
            # Single qubit gate
            if gate_flag == 1:
                for idx in range(1,self.n_qubits+1):
                    # Is qubit affected by gate
                    if idx in aff_qubits:
                        circ_str[2*(idx-1)] += used_str
                    else:
                        circ_str[2*(idx-1)] += non_used_str
                        
                    # Update empty space
                    circ_str[2*idx-1] += empty_space
            
            # Double qubit gates are trickier
            # SWAP gates
            if ideal_gate in ['RSWAP','SWAP']:
                used_str = ideal_gate + '-'
                non_used_str = ''.join(['-']*(len(ideal_gate)+1))
                
                # Edit the qubit lines first
                for idx in range(1,self.n_qubits+1):
                    # Is qubit affected by gate
                    if idx in aff_qubits:
                        circ_str[2*(idx-1)] += used_str
                    else:
                        circ_str[2*(idx-1)] += non_used_str
                        
                # Now fill in the empty spaces
                for idx in range(1,self.n_qubits):
                    if idx in range(min(aff_qubits),max(aff_qubits)):
                        if ideal_gate == 'SWAP':
                            circ_str[2*idx-1] += '  |  '
                        elif ideal_gate == 'RSWAP':
                            circ_str[2*idx-1] += '  |   '
                    else:
                        circ_str[2*idx-1] += ''.join([' ']*(len(ideal_gate)+1))
                        
            # CTRL gates        
            if ideal_gate[:4] == 'CTRL':
                targ_str = ideal_gate + '-'
                ctrl_str = '--o---'
                non_used_str = '------'
                
                # Edit the qubit lines first
                for idx in range(1,self.n_qubits+1):
                    # First qubit indices are always the ctrl qubits
                    if idx in aff_qubits[:-1]:
                        circ_str[2*(idx-1)] += ctrl_str
                    # Last qubit index is always the target qubit
                    elif idx == aff_qubits[-1]:
                        circ_str[2*(idx-1)] += targ_str
                    else:
                        circ_str[2*(idx-1)] += non_used_str
    
                # Now fill in the empty spaces
                for idx in range(1,self.n_qubits):
                    if idx in range(min(aff_qubits),max(aff_qubits)):
                        circ_str[2*idx-1] += '  |   '
                    else:
                        circ_str[2*idx-1] += ''.join([' ']*(len(ideal_gate)+1))
            
            # Reset things for next loop        
            curr_gate = self.get_next_gate()
            gate_flag = -1
    
        # Tidy up output by adding an extra padding at the end of each qubit
        # line
        for idx in range(len(circ_str)):
            if idx%2 == 0:
                circ_str[idx] += '-'
            else:
                circ_str[idx] += ' '
     
        # Print the circuit
        print(f'Ideal circuit: {self.name}\n')
        for idx in range(len(circ_str)):
            print(circ_str[idx])
            
        self.curr_gate_idx = initial_gate_idx
