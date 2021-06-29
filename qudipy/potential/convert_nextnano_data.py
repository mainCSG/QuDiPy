import numpy as np
import os
import re

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
            # x = np.array(x, dtype='O')
            # y = np.array(y, dtype='O')
            # z = np.array(z, dtype='O')
            # x = np.array(x, dtype='float')
            # y = np.array(y, dtype='float')
            # z = np.array(z, dtype='float')

            x = [float(i) for i in x]
            y = [float(i) for i in y]
            z = [float(i) for i in z]

            return x, y, z

def load_file_test(filename):
    '''
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

def parse_ctrl_names(filename):
    '''
    Parameters
    ----------
    filename: String which is the relative path to the filename and containts control name information
    
    Returns
    -------
    ctrl_names: List of control names
    '''

    # parse string via _,\, /
    parsed_filename = re.split(r'[_/\\]',str(filename))
    ctrl_names = []

    # search through list of strings for control name attached to float and seperated by an _
    for idx, strg in enumerate(parsed_filename):
        try:
            if float(strg) < 100:
                ctrl_names.append(parsed_filename[idx-1])
        except ValueError:
            pass

    return ctrl_names

def parse_ctrl_vals(filename):
    '''
    Parameters
    ----------
    filename: String which is the relative path to the filename and containts control value information
    
    Returns
    -------
    ctrl_names: List of control values
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
    Parameters
    ----------
    folder: String, name of the folder where nextnano++ files are stored
    
    Returns
    -------
    data: List, where each element is a list of voltages, potentials, and coordinates

    nextnano++ file structure:

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
            voltage = parse_ctrl_vals(subdir)
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
    pot3DArray = np.reshape(potential,(xsize,ysize,zsize))

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

def xy_potential_test(potentialL, gates, slice, f_type, dir_path):
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


    potentialL_copy = potentialL.copy()
    # loop through each combination of gate voltages
    for i in potentialL_copy:
        if f_type == "potential":
            # slice an x-y plane of the potentials
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, f_name)
        elif f_type == "field":
            potential2D = reshape_potential(i[1], i[2][0], i[2][1], i[2][2], slice, f_name)
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
    for v in gates:
        if len(v) > 1:
            shape = shape + (len(v),)
      
    shape = shape+ (len(i[2][0]), len(i[2][1]))

    potential_reshaped = np.reshape(potential_overall,shape)

    # # transpose array to align with load_data.load_potential()
    # potential_reshaped = np.transpose(potential_reshaped)

    for j in range(len(i[0])):
        f_name = f_name+ '_' + gates[j] + '_' + "{:.3f}".format(i[0][j])

    f_name +='.txt'

    # create directory for preprocessed data
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    # join file name to directory path
    f_path = os.path.join(dir_path,f_name)

    # save potential data for xy slice
    np.savetxt(f_path, potential_reshaped, delimiter=',')

    return 0


    # for j in range(len(i[0])):
    #     f_name = f_name+ '_' + gates[j] + '_' + "{:.3f}".format(i[0][j])

    # f_name +='.txt'

    # # create directory for preprocessed data
    # if not os.path.exists(dir_path):
    #     os.mkdir(dir_path)

    # # join file name to directory path
    # f_path = os.path.join(dir_path,f_name)

    # # save potential data for xy slice
    # np.savetxt(f_path, coords_and_pot, delimiter=',')


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


