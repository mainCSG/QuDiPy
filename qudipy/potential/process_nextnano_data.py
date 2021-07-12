"""
Functions for covnerting 3D nextnano data to 2D data along user defined z-coordinate.

@author: Kewei, Zach
"""

import numpy as np
import os
import re


def load_file(filename):
    '''
    This function loads data files from a user defined nextnano simulation directory.

    Parameters
    ----------
    filename: relative path for the files of interest
    
    Returns
    -------
    x,y,z: a single array ordered by the coordinates for potential.dat files or
        a tuple of 3 element, x, y, z for potential.coord files
    '''

    # import .dat data
    if filename[-4:] == '.dat':

        data = np.genfromtxt(filename,dtype=float)

        return data

    # import .coord data
    if filename[-6:] == '.coord':

         # read in xyz dimensions from .fld file for extracting .coord data
        with open(filename.replace('.coord','.fld'), 'r') as f:
            d = f.readlines()

            xdim = [int(i) for i in d[3].split() if i.isdigit()][0]
            ydim = [int(i) for i in d[4].split() if i.isdigit()][0]
            zdim = [int(i) for i in d[5].split() if i.isdigit()][0]

        # extract xyz coordinate data
        with open(filename, 'r') as f:
            d = f.readlines()

            # convert list of strings to list of floats
            data = []
            for i in list(filter(lambda x: x != '\n', d)):
                data.append(float(i.strip()))

            x = data[:xdim]
            y = data[xdim:xdim+ydim]
            z = data[xdim+ydim:xdim+ydim+zdim]

            return x, y, z

def parse_ctrl_items(filename,ctrl_type = 'value'):
    '''
    This function collects all of the control items in the filename string and
    potential files are assumed to follow the syntax: 
    'TYPE_C1NAME_C1VAL_C2NAME_C2VAL_..._CNNAME_CNVAL.txt'
    where TYPE = 'Uxy' or 'Ez'. 
    Refer to tutorial for a more explicit example.

    Parameters
    ----------
    filename: String which is the relative path to the filename and containts control items (either name/value) information
    
    Returns
    -------
    ctrl_type: List of control items of the specific control type
    '''

    # parse string via _,\, /
    parsed_filename = re.split(r'[_/\\]',filename)

    if ctrl_type.lower() in ['value','values']:

        ctrl_vals = []

        # search through list of strings for floats from the parsed file name
        for i in parsed_filename:
            try:
                if float(i) < 100:
                    ctrl_vals.append(float(i))
            except ValueError:
                pass

        return ctrl_vals

    elif ctrl_type.lower() in ['name','names']:

        ctrl_names = []

        # search through list of strings for control name attached to float and seperated by an _
        for idx, strg in enumerate(parsed_filename):
            try:
                if float(strg) < 100:
                    ctrl_names.append(parsed_filename[idx-1])
            except ValueError:
                pass

        return ctrl_names

def import_folder(folder):
    '''
    Parameters
    ----------
    folder: String, name of the folder where nextnano files are stored
    
    Returns
    -------
    data: List, where each element is a list of voltages, potentials, and coordinates

    nextnano file structure:

    /simulation_runs
        /simulation_run_#_with_gate_voltages
            /output directory
                /data files

    '''

    # list which holds all of the simulation run data grouped by run
    data = []

    # return lists of all subdirectories, base directories (ignored), files in folder
    for subdir, _, files in os.walk(folder):


        # list containing run data i.e. voltages, potential, 3-tuple of coordinates
        data_per_run = []

        # parse voltage information for the directory one level higher than /output
        if subdir != folder and subdir[-6:] == 'output':

            # print('Importing .coord and .dat data files from {}:'.format(subdir.replace(str(folder),'')), end = '\r')
            print('Importing .coord and .dat data files from {}:'.format(subdir.replace(str(folder),'')))
            voltage = parse_ctrl_items(subdir,'value')
            data_per_run.append(voltage)

            # first append potential data
            for file in files:
                filename = os.path.join(subdir, file)

                if filename[-13:] == 'potential.dat':
                    # print(file)
                    data_per_run.append(load_file(filename))
                    

            # second append coordinate data
            for file in files:
                filename = os.path.join(subdir, file)

                if filename[-15:] == 'potential.coord':
                    # print(file)
                    data_per_run.append(load_file(filename))

            data.append(data_per_run)

    return data

def retrieve_ctrl_vals(potential):
    '''
    Parameters
    ----------
    potential: List containing lists of voltages, potentials, and coordinates
    
    Returns
    -------
    voltages: List of all unique control values per control name for all simulation runs in data directory

    ex:

        simulation runs:
        [0.1, 0.2, 0.2, 0.20, 0.1]
        [0.1, 0.2, 0.2, 0.22, 0.1]
        [0.1, 0.2, 0.2, 0.24, 0.1]
        [0.1, 0.2, 0.2, 0.25, 0.1]
        [0.1, 0.2, 0.2, 0.26, 0.1]

        unique voltages per gate:
        V1 = [0.1]
        V2 = [0.2]
        V3 = [0.2]
        V4 = [0.2, 0.22, 0.24, 0.25, 0.26]
        V5 = [0.1]

        output:
        [V1, V2, V3, V4, V5]

    '''

    # store unique voltages per gate
    voltages = []

    # loop over all gate voltages
    for i in range(len(potential[0][0])):

        # initialize with voltage of the ith gate's first simulation run
        gate = [potential[0][0][i]]

        # loop over all simulation runs
        for j in range(len(potential)-1):

            # store voltage if unique
            if potential[j][0][i] != potential[j+1][0][i]:
                gate.append(potential[j+1][0][i])

        # convert between set/list to remove any duplicate voltages per gate and sort floats in ascending order for later interpolation
        gate = set(gate)
        gate = sorted(list(gate))

        voltages.append(gate)
    return voltages

def reshape_potential(potential, x, y, z, slice, f_type):
    '''
    input:  1d potential array,
            lists of x, y ,z coordinates
            the z coordinate indicating the slice of x-y plane
    output: a 2d array of the potentials in the x-y plane
    '''

    '''
    Parameters
    ----------
    potential: 1d potential array
    x,y,z: coordinate data for potential
    slice: z coordinate value to generate a 2D potential interpolation object for from all simulation runs
        
    Keyword Arguments
    ----------
    f_type: field type identifier
    
    Returns
    -------

    
    '''

    # find the index for the desired z-coordinate
    index = z.index(slice)

    # number of data points per axis provided from simulations 
    xsize = len(x)
    ysize = len(y)
    zsize = len(z)

    ''' 
        Original setup to reshape the potential data. However, the format of the
        preprocessed potential files in /QuDiPy tutorial data/Pre-precessed potentials
        have the x/y coordinate potential data running along the columns/rows where as
        the code below formats the potnetial data with y/x coordinate potential data 
        running along the columns/rows. 2D potential data files will be transposed prior 
        to writting to data file.
    '''
    pot3DArray = np.reshape(potential,(xsize,ysize,zsize), order='F')



    # test
    # for i in range(zsize):
    #     pot3DArray[:,:,i] = pot3DArray[:,:,i] * -1

    print('size of 3d array {}'.format(np.shape(pot3DArray))) # delete

    if f_type in ['field', 'electric', 'Ez']:
        gradient = np.gradient(pot3DArray,x,y,z)[-1]
        field2DArray = gradient[:, :, index]
    elif f_type in ['pot', 'potential', 'Uxy']:
        field2DArray = pot3DArray[:, :, index]

    return field2DArray

def xy_potential(potentialL, gates, slice, f_type, dir_path):
    '''
    Parameters
    ----------
    potentialL: 
    slice:
    gates:
        
    Keyword Arguments
    ----------
    f_type:
    dir_path:
    
    Returns
    -------
    Potential or electric field XY-plane data is saved for slice
    '''

    potentialL_copy = potentialL.copy()
    # loop through each combination of gate voltages
    for i in potentialL_copy:

        if f_type in ['pot', 'potential', 'Uxy']:
            f_name = 'Uxy'
            # slice an x-y plane of the potentials
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, f_name)
        elif f_type in ['field', 'electric', 'Ez']:
            f_name = 'Ez'
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, f_name)

        # create an array of zeros with the demension of the potential 2D slice and x/y coordinate axis
        coords_and_pot = np.zeros((np.shape(potential2D)[0]+1, np.shape(potential2D)[1]+1),dtype=float)

        # insert x,y, and potential 2D slice into array
        coords_and_pot[0,1:] = i[2][1]
        coords_and_pot[1:,0] = i[2][0]
        # coords_and_pot[1:,1:] = potential2D
        coords_and_pot[1:,1:] = -1*potential2D # delete

        # transpose array to align with load_data.load_potential()
        coords_and_pot = np.transpose(coords_and_pot)

        for j in range(len(i[0])):
            f_name = f_name+ '_' + gates[j] + '_' + "{:.3f}".format(i[0][j])

        f_name +='.txt'

        # create directory for preprocessed data
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

        # join file name to directory path
        f_path = os.path.join(dir_path,f_name)

        # save potential data for xy slice
        np.savetxt(f_path, coords_and_pot, delimiter=',')
        
    return 0

def write_data(input_nextnano,output_preprocessed,slice,data_forms):
    '''
    Parameters
    ----------
    input_nextnano: 
    output_preprocessed:
    slice:
        
    Keyword Arguments
    ----------
    data_forms:
    
    Returns
    -------
    A data directory containing field data based on the data_forms argument.
    '''
    potentialL = import_folder(input_nextnano)

    for subdir, _, _ in os.walk(input_nextnano):

        # parse voltage information for the directory one level higher than /output
        if subdir != input_nextnano and subdir[-6:] == 'output':
            gates = parse_ctrl_items(subdir,'name')
            break

    # write xy potential files
    for i in data_forms:
        print('Converting 3D nextnano simulation data too 2D XY-plane {} data slice for z = {}.'.format(i,slice))
        xy_potential(potentialL, gates, slice, i,output_preprocessed)


