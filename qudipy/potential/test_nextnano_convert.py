import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os
import re

########## Global Variables ##########

# number of gates
numberOfGates = 5

# gate voltages
V1 = [0.1]
V2 = [0.2]
V3 = [0.2]
V4 = [0.2, 0.22, 0.24, 0.25, 0.26]
V5 = [0.1]

########## Helper Functions ##########

def is_float(string):
    """ True if given string is float else False"""
    try:
        return float(string)
    except ValueError:
        return False

def is_int(string):
    """ True if given string is int else False"""
    try:
        return int(string)
    except ValueError:
        return False

def load_file(filename):
    """
    returns a single array ordered by the coordinates for potential.dat
            a tuple of 3 element, x, y, z for coord files
    """
    data = []
    x = []
    y = []
    z = []
    counter = 0
    with open(filename, 'r') as f:
        d = f.readlines()
        if filename[-4:] == '.dat':
            for i in d:
                k = i.rstrip().split(" ")
                data.append(float(k[0]))     
            data = np.array(data, dtype='O')
            return data
        else:
            for i in d:
                k = i.rstrip().split(" ")
                if is_float(i)==False:
                    # append number list if the element is an int but not float
                    try:
                        int(i)
                        if counter == 0:
                            x.append(float(k[0]))
                        elif counter == 1:
                            y.append(float(k[0]))
                        else:
                            z.append(float(k[0]))
                    # ValueError happens when it hits an empty line
                    except ValueError:
                        # print(i)
                        counter+=1
                # counter keeps track of which coord the data belong to
                elif counter == 0:
                    x.append(float(k[0]))
                elif counter == 1:
                    y.append(float(k[0]))
                else:
                    z.append(float(k[0]))
            x = np.array(x, dtype='O')
            y = np.array(y, dtype='O')
            z = np.array(z, dtype='O')
            return x, y, z

def reshape_potential(potential, x, y, z, slice, option):
    """
    input:  1d potential array, 
            lists of x, y ,z coordinates
            the z coordinate indicating the slice of x-y plane
    output: a 2d array of the potentials in the x-y plane
    """
    index = np.where(z==slice)[0]
    N = len(x)
    M = len(y)
    Q = len(z)
    pot3DArray = np.reshape(potential,(N,M,Q))
    if option == "field":
        gradient = np.gradient(pot3DArray,x,y,z)[-1]
        pot2DArray = gradient[:, :, index]
    else:
        pot2DArray = pot3DArray[:, :, index]
    return pot2DArray
    
def parse_voltage(filename):
    """
    input: a string, the filename 
           an int, number of gates
    output: a list of voltages of each gate
    """
    org = re.split("[_/]",filename)
    s = []
    delete = []
    for i in org:
        try:
            if float(i) < 100:
                s.append(float(i))
        except ValueError:
            delete.append(i)
    return s

def import_folder(folder):
    """
    input: a string, name of the folder where nextnano++ files are stored 
    output: a list, where each element is a list of voltages, potentials, and coordinates
    """
    L = []                  # each element in L would be a list of voltages, potentials, and coordinates
    counter = 0             # track which subdirectory 
    for subdir, dirs, files in os.walk(folder):
        if subdir != folder and subdir[-7:] != '/output':
            counter += 1
            voltage = parse_voltage(subdir)
            L.append([voltage])
        for file in files:
            filename = os.path.join(subdir, file)
            # always first .dat then .coord
            if filename[-4:] == '.dat' or filename[-6:] == '.coord':
                L[counter-1].append(load_file(filename))
    return L

def group_2D_potential(potentialL, voltages, coord, slice, option):
    """
    input:  a list, where each element is a list of voltages, potentials, and coordinates
            a list of gate voltages
            a float indicating the x-y plane
    output: an n-dimensial potential file, where n = number of gates + 2
    """
    potentialL_copy = potentialL.copy()
    # loop through each combination of gate voltages
    for i in potentialL_copy:
        if option == "potential":
            # slice an x-y plane of the potentials
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, option)
        elif option == "field":
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, option)
        i[1] = potential2D
        # reverse the list of voltages for sorting purpose
        i[0].reverse()
    potentialL_copy.sort()

    # stack up the potential arrays in the correct order
    potential_elmt = ()
    for i in range(len(potentialL_copy)):
        potential_elmt = potential_elmt + (potentialL_copy[i][1],) 
    potential_overall = np.stack(potential_elmt, axis = 0)

    # get the shape of the potential based on the number of gates and the voltages of each gate
    shape = ()
    for v in voltages:
        if len(v) > 1:
            shape = shape + (len(v),)
    shape = shape+ (len(coord[0]), len(coord[1]))
    
    potential_reshaped = np.reshape(potential_overall,shape)
    return potential_reshaped

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
        coords_and_pot[1:,1:] = potential2D

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
    
    potentialL = import_folder(input_nextnano)

    for subdir, _, _ in os.walk(input_nextnano):

        # parse voltage information for the directory one level higher than /output
        if subdir != input_nextnano and subdir[-6:] == 'output':
            gates = parse_ctrl_names(subdir)
            break

    # write xy potential files
    for i in data_forms:
        print('Converting 3D nextnano++ simulation data too 2D XY-plane {} data slice for z = {}.'.format(i,slice))
        xy_potential(potentialL, gates, slice, i,output_preprocessed)
        group_2D_potential(potentialL, voltages, coord, slice, option)
        


