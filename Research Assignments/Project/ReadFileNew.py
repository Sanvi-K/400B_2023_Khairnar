#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import astropy.units as u


# In[2]:


def ReadNew(filename, merger = False):
    """ This function reads text file and return the time and total number of particles. 
        It also gives a data array with information of particle type, mass(in 1e10 Msun), 'x', 'y', 'z' position(in kpc) 
        and the velocities in x,y,z direction (in km/s)
    
        Inputs:
        filename: 'string'
            name of text file from which data need to be extracted
        
        Outputs:
        data: 'numpy array'
            information with mass of each particle and their positions, and velocities in x,y,z directions 
    
    """
    if merger == False:
        "when reading file for individual galaxy"
        file = open(filename,'r') # open the file
        line1 = file.readline() # read first line in the file
        label, value = line1.split() # split the line to read the label and the numerical value
        time = float(value)*u.Myr # add units to the time

        line2 = file.readlines()[0] # read second line in the file
        label, value = line2.split() # split the line to read the label and the numerical value
        number_of_particles = int(value) # specify that the value is an integer and rename it

        file.close()

        data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3) # get remaining of the data for each particle
    elif merger == True:
        "When reading file for merger remnant"
        data = np.genfromtxt(filename, dtype=None, names=True) # get remaining of the data for each particle
   
    return(data)

