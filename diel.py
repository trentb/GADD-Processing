# -*- coding: utf-8 -*-
"""
Python-GADD module for reading .diel files
Author: Jared Carter, Clive Randall Group
Version 1.0
Last updated: 2015-01-20
"""
import numpy as np
import pandas as pd

def getdata(filename):
    '''
    This function reads the .diel file and returns a dictionary with the values
    of each sweep as well as a dictionary with a list of each measurement.
    
    Parameters
    ----------
    filename : (string)
        Name of .diel file to be imported

    Returns
    -------
    s : (dict)
        A dictionary where the keys are the sweeps and their values are lists
        of the swept parameter
    m : (dict)
        A dictionary where the keys are the thing being measured and the values
        are lists of the measured values
    '''
    # initialize values
    az = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    letter = 0
    s = {}
    m = {}
    H = True
    # open the file
    with open(filename,'r+b') as f:
        # iterate over every line in the file
        for line in f:
            line = line.replace('\r\n','')
            # What to do in the header
            if H is True:
                # Define a new sweep value
                if '_set' in line or 'COUNTER' in line:
                    sline = line.split('\t')
                    xname = az[letter]+' '+sline[0]
                    line = f.next().replace('+','').replace(' ','')
                    sline = line.split('\t')
                    del sline[-1]   # Last value in list is empty
                    xdata = [float(x) for x in sline]
                    s[xname] = xdata
                    letter += 1
                # Make a key in the measurement dictionary
                elif 'REAL' in line or 'CMPLX' in line:
                    xname = line.replace('  ',' ')
                    m[xname] = []
            # What to do when outside of the header
            elif H is False:
                # This line starts a measurement
                if line in m:
                    key = line
                    line = f.next()
                    # complex
                    if '\t' in line:
                        sline = line.split('\t')
                        realline = float(sline[0])
                        imagline = float(sline[1])
                        xdata = complex(realline,imagline)
                        m[key].append(xdata)
                    # real
                    else:
                        m[key].append(float(line))
                        
                        
            # Are we in the header?
            if '*' in line:
                H = False
    return s,m

def createdf(data, rows, cols):
    '''
    Takes 1-D lists of data, rows, and cols to create a 2-D dataframe with the
    rows and columns labeled. The the rows, cols, and data should be in this
    order in the .diel file:
    
        SWEEP (rows)
            SWEEP (cols)
                MEAS (data)
            ENDSWEEP
        ENDSWEEP
    
    Parameters
    ----------
    data : (list)
        List of measurement values (output from diel.getdata)
    rows : (list)
        List of values that were swept. These should have been swept before the
        sweep going in the row field.
    cols : (list)
        List of values that were swept. These should be from the same sweep
        that the measurement was taken.

    Returns
    -------
    df : (pandas.DataFrame)
        A 2 dimensional array of the data with rows and columns labeled. Use
        df.T if the data needs to be transposed.
    '''
    tempdata = np.array(data)
    tempdata = tempdata.reshape(len(rows),len(cols))
    df = pd.DataFrame(data=tempdata, index=rows, columns=cols)
    return df

def createabsdf(data, rows, cols):
    '''
    Takes 1-D lists of data, rows, and cols to create a 2-D dataframe with the
    rows and columns labeled. This function returns the absolute value of the
    input data.
    
    Parameters
    ----------
    data : (list)
        List of measurement values (output from diel.getdata)
    rows : (list)
        List of values that were swept. These should have been swept before the
        sweep going in the row field.
    cols : (list)
        List of values that were swept. These should be from the same sweep
        that the measurement was taken.

    Returns
    -------
    df : (pandas.DataFrame)
        A 2 dimensional array of the data with rows and columns labeled. Use
        df.T if the data needs to be transposed.
    '''
    tempdata = np.array(data)
    tempdata = tempdata.reshape(len(rows),len(cols))
    tempdata = np.abs(tempdata)
    df = pd.DataFrame(data=tempdata, index=rows, columns=cols)
    return df

def zprime(eps, freq):
    '''
    This function takes a one dimensional list of complex permittivity data and
    the corresponding list of frequencies and returns a one dimensional numpy 
    array of complex impedance.
    
    Parameters
    ----------
    eps : (list)
        A one dimensional list of complex permittivity data
    freq: (list)
        Frequencies that correspond 1:1 with the complex permittivity data

    Returns
    -------
    result : (numpy.array)
        A 1 dimensional numpy array of complex impedance (instead of complex
        permittivity)
    '''
    eps = np.array(eps)
    freq = np.array(freq)
    sumsquares = np.add(np.square(eps.imag), np.square(eps.real))
    denom = np.multiply(2*np.pi*freq, sumsquares)
    re = np.divide(eps.imag, denom)
    im = np.divide(-eps.real, denom)*1j
    return np.add(re, im)
    
def zprimeprime(eps, freq): # Depreciated
    '''
    This function takes a one dimensional list of complex permittivity data and
    the corresponding list of frequencies and returns a one dimensional numpy 
    array of complex impedance.
    
    Parameters
    ----------
    eps : (list)
        A one dimensional list of complex permittivity data
    freq: (list)
        Frequencies that correspond 1:1 with the complex permittivity data

    Returns
    -------
    result : (numpy.array)
        A 1 dimensional numpy array of complex impedance (instead of complex
        permittivity)
    '''
    result = np.zeros(len(freq),dtype=np.complex128)
    for i in range(len(freq)):
        re = eps[i].imag/(2*np.pi*freq[i]*(np.square(eps[i].imag)+np.square(eps[i].real)))
        im = -eps[i].real/(2*np.pi*freq[i]*(np.square(eps[i].imag)+np.square(eps[i].real)))
        result[i] = complex(re, im)
    return result
    
def impedance(filename):
    '''
    This function takes a .diel if impedance data and returns a dictionary of
    dataframes with columns of frequency, real impedance, and imaginary
    impedance for each temperature.
    
    Parameters
    ----------
    filename : (string)
        Name of .diel file to be imported

    Returns
    -------
    f : (dict)
        A dictionary where the keys are the filename_temperature and the values
        are pandas.DataFrame objects with columns of frequency, real impedance,
        and imaginary impedance.
    '''
    f = {}
    s,m = getdata(filename)
    data = np.array(m['CMPLX EPS1'])
    data.resize(len(s['A TEMPERATURE_set']),len(s['B FREQ_set']))
    data = data.T
    cmplxz = np.apply_along_axis(zprime, 0, data, freq=s['B FREQ_set'])
    df = pd.DataFrame(cmplxz, index=s['B FREQ_set'], columns=s['A TEMPERATURE_set'])

    for i in range(len(s['A TEMPERATURE_set'])):
        points = df.iloc[:,i]
        d = {'Real at '+str(s['A TEMPERATURE_set'][i]) : pd.Series(np.real(points), index=points.index),
             'IMAG' : pd.Series(np.imag(points), index=points.index)}
        xname = filename[:-5]+'_%d' %s['A TEMPERATURE_set'][i]
        newdf = pd.DataFrame(d, columns=['Real at '+str(s['A TEMPERATURE_set'][i]), 'IMAG'])
        f[xname] = newdf
    return f

def ivloop(files):
    '''
    Reads .diel files and converts them to spreadsheets for plotting abs(I)-V
    curves.
    
    Parameters
    ----------
    files : (list)
        List of of .diel files to be imported

    Returns
    -------
    f : (dict)
        A dictionary where the keys are the filenames and the values
        are pandas.DataFrame objects with rows of voltage and columns of
        current at time t.
    '''
    f = {}
    for filename in files:
        s,m = getdata(filename)
        agefirst = False
        # Was the sample aged before the measurement?
        if 'REAL PA0' in m.keys():
            agefirst = True
        
        if agefirst is True:
            df1 = createabsdf(m['REAL PA1'], s['B DC_set'], s['C TIME_set'])
            df2 = createabsdf(m['REAL PA3'], s['E DC_set'], s['F TIME_set'])
            df3 = createabsdf(m['REAL PA5'], s['H DC_set'], s['I TIME_set'])
        
        elif agefirst is False:
            try:
                df1 = createabsdf(m['REAL PA1'], s['A DC_set'], s['B TIME_set'])
            except KeyError:
                print filename+' is the wrong kind of file for this operation.'
                break
            df2 = createabsdf(m['REAL PA3'], s['D DC_set'], s['E TIME_set'])
            df3 = createabsdf(m['REAL PA5'], s['G DC_set'], s['H TIME_set'])
            
        final = pd.concat([df1, df2, df3])
        xname = filename[:-5]+'_IVloop'
        f[xname] = final
    return f
    
def ageprogression(files):
    '''
    Pulls current during aging from .diel files
    Reads .diel files and extracts the leakage current during aging. If 
    multiple files are listed, aging data will be appended to a 
    pandas.DataFrame where the index is the time in minutes and the columns
    contain current values. The second column starts at the index immediately
    followng the previous column and so on.
    
    Parameters
    ----------
    files : (list)
            List of .diel files to be imported (should be in order)
            
    Returns
    -------
    DataFrame : (pandas.DataFrame)
                Dataframe that contains current during the aging process
    '''
    adj = 0
    d = {}
    for filename in files:
        s,m = getdata(filename)
        if 'REAL PA0' in m.keys():
            values = np.array(m['REAL PA0'])
            index = np.arange(adj, len(values)+adj)
            d[filename] = pd.Series(values, index=index)
            adj += len(values)
    return pd.DataFrame(d)
    
def rayleigh(filename,freq,thickness,diameter,minpts,maxpts,maxi):
        '''
        Reads .diel file and extracts CD0/Voltage field for requested 
        frequency for all 3 runs. Converts to dielectric constant and electric
        field. Reshapes data and performs linear regressions on all ranges i:j
        where i ranges from point 1 to max i. And j ranges from 1+minpts to
        maxi maxpts. Averages linear regressions for counter 2 and 3 and checks
        if average R^2 is better than existing average R^2. If this is the case,
        A new best R^2 is saved along with the best i and j values, slope,
        intercept and confidence intervals. Averages counter 2 and 3 data
        and exports that data
        
        Paramters
        ---------
        filename : (string)
                    Name of file of interest
        freq : (integer)
                    Frequency of interest in Hz
        thickness : (float)
                    Thickness of dielectric in microns
        diameter : (float)
                    Diameter of contact in microns
        minpts : (integer)
                    Minimum number of points for fitting
        maxpts : (integer)
                    Maximum number of points for fitting
        maxi : (integer)
                    Maximum value of i (start of fit)
                    
        Returns
        -------
        average : (numpy.array)
                    Contains average dielectric constant and electric field data
        besti : (integer)
                    Best value for i by maximizing R^2
        bestj : (integer)
                    Best value for j by maximizing R^2
        bestslope : (float)
                    Best value of slope by maximizing R^2
        bestint : (float)
                    Best value of intercept by maximizing R^2
        bestrsq : (float)
                    Best value of R^2
        bestconf : (list)
                    95% confidence interval on slope and intercept
        '''
        
        s,m = getdata(filename)