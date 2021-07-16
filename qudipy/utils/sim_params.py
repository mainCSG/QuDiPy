"""
SimulationParameters class

@author: simba
"""

class SimulationParameters:
    
    def __init__(self):
        
        # TODO: Eventually I want all of this stuff parsed from an inputted
        # file.
        
        # [PATH/TO/WHERE/POTENTIAL/FILES/ARE/STORED]
        self.path_to_pot = ''
        # [PATH/TO/WHERE/PULSE/FILES/ARE/STORED]
        self.path_to_pulses = ''
        
        self.gate_voltages_to_load = {
              "V1": [0.1],
              "V2": [0.2, 0.22, 0.24, 0.26, 0.27, 0.28, 0.29, 0.3],
              "V3": [0.2, 0.22, 0.24, 0.26, 0.27, 0.28, 0.29, 0.3],
              "V4": [0.2, 0.22, 0.24, 0.26, 0.27, 0.28, 0.29, 0.3],
              "V5": [0.1]
            }
        
        self.number_of_gates = len(self.gate_voltages_to_load.keys())
        self.z_2DEG = -1E-9 # z-axis location of 2DEG
