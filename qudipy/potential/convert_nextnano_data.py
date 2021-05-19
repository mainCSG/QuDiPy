import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os
import re
from pathlib import Path
import itertools
import errno

def load_file(filename):
    '''
    returns a single array ordered by the coordinates for potential.dat
            a tuple of 3 element, x, y, z for coord files
    '''

    # import .dat data
    if filename[-4:] == '.dat':

        data = np.genfromtxt(filename,dtype=float)
        print(data)
        return data

    # import .coord data
    if filename[-6:] == '.coord':

         # read in xyz dimensions from .fld file for extracting .coord data
        with open(filename.replace('.coord','.fld'), 'r') as f:
            d = f.readlines()

            xdim = [int(i) for i in d[3].split() if i.isdigit()]
            ydim = [int(i) for i in d[4].split() if i.isdigit()]
            zdim = [int(i) for i in d[5].split() if i.isdigit()]

        # extract xyz coordinate data
        with open(filename, 'r') as f:
            d = f.readlines()

            x = []
            y = []
            z = []
            # x-coordinate data
            for i in range(xdim[0]):
                x.append(float(d[i]))

            # y-coordinate data
            for i in range(ydim[0]):
                y.append(float(d[i]))

            # z-coordinate data
            for i in range(zdim[0]):
                z.append(float(d[i]))

            return x, y, z

def parse_ctrl_names(filename):
    '''

    input: a string which is the relative path to the filename and containts control name information
    output: a list of control names

    '''

    # parse string via _,\, /
    parsed_filename = re.split(r'[_/\\]',str(filename))
    ctrl_names_list = []

    # search through list of strings for control name attached to float and seperated by an _
    for idx, strg in enumerate(parsed_filename):
        try:
            if float(strg) < 100:
                ctrl_names_list.append(parsed_filename[idx-1])
        except ValueError:
            pass

    return ctrl_names_list

def parse_ctrl_vals(filename):
    '''
    input: a string which is the relative path to the filename and containts control value information
    output: a list of control values
    '''

    # parse string via _,\, /
    parsed_filename = re.split(r'[_/\\]',filename)
    ctrl_vals_list = []

    # search through list of strings for floats from the parsed file name
    for i in parsed_filename:
        try:
            if float(i) < 100:
                ctrl_vals_list.append(float(i))
        except ValueError:
            pass

    return ctrl_vals_list

def import_folder(folder):
    '''
    input: a string, name of the folder where nextnano++ files are stored
    output: a list, where each element is a list of voltages, potentials, and coordinates

    nextnano++ file structure:

    /simulation_runs
        /simulation_run_#_with_gate_voltages
            /output directory
                /data files

    '''

    # list which holds all of the simulation run data grouped by run
    data_list = []

    # return lists of all subdirectories, base directories (ignored), files in folder
    for subdir, _, files in os.walk(folder):

        # list containing run data i.e. voltages, potential, 3-tuple of coordinates
        data_per_run = []

        # parse voltage information for the directory one level higher than /output
        if subdir != folder and subdir[-6:] == 'output':
            voltage = parse_ctrl_vals(subdir)
            data_per_run.append(voltage)

            # first append potential data
            for file in files:
                filename = os.path.join(subdir, file)
                if filename[-13:] == 'potential.dat':
                    data_per_run.append(load_file(filename))

            # second append coordinate data
            for file in files:
                filename = os.path.join(subdir, file)
                if filename[-15:] == 'potential.coord':
                    data_per_run.append(load_file(filename))

            data_list.append(data_per_run)

    return data_list

def retrieve_ctrl_vals(potential):
    '''
        input: potential is list of voltages, potentials, and coordinates
        output: list of all unique control values per control name for all simulation runs in data directory

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
    return(voltages)

def reshape_potential(potential, x, y, z, slice, option):
    '''
    input:  1d potential array,
            lists of x, y ,z coordinates
            the z coordinate indicating the slice of x-y plane
    output: a 2d array of the potentials in the x-y plane
    '''
    index = z.index(slice)

    xsize = len(x)
    ysize = len(y)
    zsize = len(z)

    pot3DArray = np.reshape(potential,(xsize,ysize,zsize))

    print('shape of array: {}'.format(pot3DArray.shape))

    if option == "field":
        gradient = np.gradient(pot3DArray,x,y,z)[-1]
        pot2DArray = gradient[:, :, index]
    else:
        pot2DArray = pot3DArray[:, :, index]

    return pot2DArray

def xy_potential(potentialL, gates, slice, option, dir_path):
    '''
    input:  
    output: 
    '''
    
    overwrite_trig = False
    potentialL_copy = potentialL.copy()
    # loop through each combination of gate voltages
    for i in potentialL_copy:
        if option == "potential":
            
            f_name = 'Uxy'

            # slice an x-y plane of the potentials
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, option)
        elif option == "field":

            f_name = 'Ez'

            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, option)


        for j in range(len(i[0])):
            f_name = f_name+ '_' + gates[j] + '_' + "{:.3f}".format(i[0][j])

        f_name +='.txt'

        # create directory for preprocessed data
        p = Path(dir_path)
        p.mkdir(exist_ok=True)

        # join file name to directory path
        f_path = Path(os.path.join(dir_path,f_name))

        # check if user already decided to overwrite existing files in the directory for a given type of data i.e. field/potential
        if overwrite_trig == True:
            # save potential data for xy slice
            np.savetxt(f_path, potential2D, delimiter=',')
        else:
            # check if there are files existing in the data directory with the same file name as the data file to be created
            if f_path.is_file():

                print('Files already exist in directory: {}'.format(dir_path))
                overwrite = input('Overwrite existing files? <y/n>: ')

                if overwrite == 'y':
                    overwrite_trig = True
                    np.savetxt(f_path, potential2D, delimiter=',')
                else:
                    print('File creation terminated.')
                    break
            else:
                # save potential data for xy slice
                np.savetxt(f_path, potential2D, delimiter=',')
        
    return 0

def write_data(input_nextnano,output_preprocessed,slice,data_forms):
    
    potentialL = import_folder(input_nextnano)

    for subdir, _, _ in os.walk(input_nextnano):

        # parse voltage information for the directory one level higher than /output
        if subdir != input_nextnano and subdir[-6:] == 'output':
            gates = parse_ctrl_names(subdir)
            break

    # write xy potential files
    for i in data_forms:
        print(i)
        xy_potential(potentialL, gates, slice, i,output_preprocessed)


########## Tests ##########

# take data source and converted data destinations paths for next nano data
input_nextnano = Path(r"qudipy/potential/nextnanoSims_Small")
output_preprocessed = Path(r"qudipy/potential/test_data")

potentialL = import_folder(input_nextnano)

# Enter the z coordinate of the x-y plane
coord = potentialL[0][2]
z = coord[2][0]

write_data(input_nextnano,output_preprocessed, z, ['potential','field'])


