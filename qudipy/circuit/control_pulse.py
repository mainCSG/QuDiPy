"""
Class for a control pulse

@author: simba
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from copy import deepcopy

class ControlPulse:
    
    def __init__(self, pulse_name, pulse_type, pulse_length=-1, ideal_gate=None):
        '''
        Initialize a ControlPulse object.

        Parameters
        ----------
        pulse_name : string
            Name of pulse.
        pulse_type : string
            Specify whether pulse is described with "experimental" or 
            "effective" control variables.
            
        Keyword Arguments
        -----------------
        pulse_length : int, optional
            Total length of pulse in [s]. The default is -1.
            This is an optional keyword because 
            sometimes you may wish to vary the pulse length across
            different simulations for a given pulse. Pulse length 
            can be set later using the set_pulse_length() method.
        ideal_gate : string, optional
            The ideal gate keyword for this control pulse. The default
            is None.

        Returns
        -------
        None.

        '''
        
        self.pulse_type = pulse_type
        
        self.name = pulse_name
        self.length = pulse_length # Units are s
        self.n_pts = 0 # Number of points in pulse (default 0), must be same among control variables
        
        self.ctrl_pulses = {}
        self.ctrl_names = []
        self.n_ctrls = 0
        
        self.ideal_gate = ideal_gate
       
        self.ctrl_time = None

        # Attribute to be initialized later in _generate_ctrl_interpolators
        self.ctrl_interps = None
        
        
    def __call__(self, time_pts):
                
        '''
        Call method for class.

        Parameters
        ----------
        time_pts : 1D float array
            Time points where we want to obtain the interpolated pulse
            values.

        Returns
        -------
        interp_pulse : 2D float array
            The interpolated pulse values at the inputted time points.

        '''
        
        # Check that we actually have a valid pulse length set
        if self.length == -1:
            print('Cannot call control pulse object.\nPulse length has not been specified.'\
                            + '\nPlease set using the .set_pulse_length() method.')
            return None
        
        # Check if the interpolators have been constructed. If not, then 
        # make them.
        if not hasattr(self,'ctrl_interps'):
            self._generate_ctrl_interpolators()
        
        # Loop through each control variable and get the interpolated pulse
        # for each time point.
        interp_pulse = np.zeros((len(time_pts),len(self.ctrl_names)))
        for ctrl_idx, ctrl in enumerate(self.ctrl_names):
            interp_pulse[:,ctrl_idx] = self.ctrl_interps[ctrl](time_pts)

        return interp_pulse    
    
    def copy(self):
        '''
        This method will do a deep copy of the current class object.

        Returns
        -------
        ControlPulse
            Deep copy of current ControlPulse object.

        '''
        return deepcopy(self)
    
    def plot(self, plot_ctrls='all', time_int='full', n=250):
        '''
        Plot the control pulse. Can plot a subset of the control
        variables and within some time interval.

        Keyword Arguments
        -----------------
        plot_ctrls : list of strings, optional
            Specify the name of each control variable pulse you wish
            to plot or plot 'all'. The default is 'all'.
        time_int : list of floats
            Specify the time interval over which to plot the pulse.
            Cannot be less than 0 or greater than the current pulse
            length. The default is the full time interval.
        n : int, optional
            Number of points to use in the pulse when plotting. The 
            default is 250.

        Returns
        -------
        None.

        '''
        
        # Get the time interval and time points to plot
        if time_int == 'full':
            min_time = 0
            max_time = self.length
        else:
            min_time = time_int[0]
            max_time = time_int[1]
            
        t_pts = np.linspace(min_time,max_time,n)
        
        # Get the actual pulse
        pulse = self(t_pts)
        
        # Check the plot_ctrls input
        if not isinstance(plot_ctrls,(list,tuple,set)):
            if plot_ctrls.lower() == 'all':
                plot_ctrls = self.ctrl_names
            # A single control variable was inputted but wasn't wrapped in a
            # list, tuple, or set, so we will wrap it in a list
            elif plot_ctrls in self.ctrl_names:
                plot_ctrls = [plot_ctrls]
            else:
                raise ValueError('Unrecognized input for control variables '+
                                 f'to plot {plot_ctrls}.\nAllowed inputs are '+
                                 'either ''all'' or a list of allowed names:\n'+
                                 f' {self.ctrl_names}')
        # If a list of names was specified, check that all are valid ctrl 
        # variable names
        elif not set(plot_ctrls).issubset(self.ctrl_names):
            raise ValueError('Unrecognized input for control variables '+
                             f'to plot {plot_ctrls}.\nAllowed inputs are '+
                             'either ''all'' or a list of allowed names:\n'+
                             f'{self.ctrl_names}')
           
        # For each pulse specified by ctrl_vars, plot those pulses
        ctrl_idxs = []
        for ctrl in plot_ctrls:
            ctrl_idxs.append(self.ctrl_names.index(ctrl))
            
        # Figure out scale to apply (if any) on time axis to make more
        # readable as default length units are ps
        if 1E-15 < max(t_pts) <= 1E-12:
            scale = 1E15
            units = '[fs]'
        elif 1E-12 < max(t_pts) <= 1E-9:
            scale = 1E12
            units = '[ps]'
        elif 1E-9 < max(t_pts) <= 1E-6:
            scale = 1E9
            units = '[ns]'
        elif 1E-6 < max(t_pts) <= 1E-3:
            scale = 1E6
            units = '[us]'
        elif 1E-3 < max(t_pts):
            scale = 1E3
            units = '[ms]'
            
        # Generate figure
        plt.figure()
        plt.plot(t_pts*scale, pulse[:,ctrl_idxs])   
        lgd_names = [self.ctrl_names[idx] for idx in ctrl_idxs]
        plt.legend(lgd_names,loc='best')
        plt.xlabel('Time ' + units)
        plt.show()
        
        
    def set_pulse_length(self, pulse_length):
        '''
        Change the pulse length by scaling the time axis of the pulse with 
        respect to the previously specified pulse.

        Parameters
        ----------
        pulse_length : int
            Pulse legnth for control pulse in [s].

        Returns
        -------
        None.

        '''
                
        # Keep track of old length
        old_length = self.length
        self.length = pulse_length
        
        # If previous length is -1, then pulse length has not been specified
        # yet so we can just set the attribute as we did above and then return
        if old_length == -1:
            return
            
        # If self.ctrl_time is defined, then we need to just scale all the
        # current time points
        if self.ctrl_time is not None:
            self.ctrl_time *= pulse_length/old_length
            
        # Check if the interpolators have been constructed. If they have, 
        # then we need to update them with the new ctrl_time attribute
        if hasattr(self,'ctrl_interps'):
            self._generate_ctrl_interpolators()
        
    def add_control_variable(self, var_name, var_pulse):
        '''
        Adds a control variable to the pulse. If the variable name is
        time, then it will store the time points in the ctrl_time
        variable.
        
        Parameters
        ----------
        var_name : string
            Name of new control variable being added.
        var_pulse : 1D array
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        # Check that the new control variable pulse has the same number of 
        # points as all the other control variables (if this is first one, 
        # then save the number of points)
        if self.n_pts != 0:
            if len(var_pulse) != self.n_pts:
                raise ValueError('Number of pulse points is not equal to the '+
                                 'current\nnumber of pulse points in previously '+
                                 'loaded control variable pulses.\n'+
                                 f'Expected {self.n_pts}, got {len(var_pulse)}.')
        else:
            self.n_pts = len(var_pulse)
        
        if var_name.lower() == 'time':
            self.ctrl_time = np.array(var_pulse)
            # Double check that the time points start at 0.
            if not np.isclose(self.ctrl_time[0], 0):
                raise ValueError(f'Cannot load time variable for {self.name}'+
                                 ' because time points do not start at 0.')

            self.set_pulse_length(var_pulse[-1])
        else:
            self.ctrl_pulses[var_name] = np.array(var_pulse)
            # self.ctrl_pulses.keys() should be in the order items are added 
            # to the dict if using python >3.7. But to keep things backwards
            # compatible, we will manually append our own list.
            self.ctrl_names.append(var_name)
            self.n_ctrls = len(self.ctrl_names)
            
    def _generate_ctrl_interpolators(self):
        '''
        Loop through every control variable and make a 1D time 
        interpolation object.

        Returns
        -------
        None.

        '''
        
        # Initialize an empty dictionary to store all the interpolators
        self.ctrl_interps = {}
        
        # Check if the time array is set, if not, then assume each point in
        # the pulse is linearly spaced in time
        if self.ctrl_time is None:
            # Build a linear time array
            self.ctrl_time = np.linspace(0, self.length, self.n_pts)
        
        # For each pulse, find the time axis points and then make the 1D
        # interpolator. All ctrl pulses will be linearly interpolated.
        for ctrl in self.ctrl_names:
            curr_pulse = self.ctrl_pulses[ctrl]
                
            self.ctrl_interps[ctrl] = interp1d(self.ctrl_time, curr_pulse)        
            
        
        
        
        
        
        
        
        
        
        
