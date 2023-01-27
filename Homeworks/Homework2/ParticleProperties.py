#!/usr/bin/env python
# coding: utf-8

# In[21]:


import numpy as np
import astropy.units as u
from ReadFile import Read


# In[22]:


def ParticleInfo(filename, particle_type, particle_number):
    """ Given the type of the particle (Dark matter, disk star, bulge star) and its number within the type,
        this function returns the properties of any given particle from a given txt file.
        
        Inputs:
        filename: 'string'
            name of text file from which data need to be extracted
            
        particle_type: 'float'
            a float value specifies the type of the particle
                1.0 = dark matter
                2.0 = disk star
                3.0 = bulge star
        
        particle_number: 'int'
            the number of particle for a specific particle type
            
        Outputs:
            Distance: 'astropy quantity'
                The 3D distance of the particle from the Galactic center (in kpc)
                
            Velocity: 'astropy quantity'
                The 3D velocity of the partcile in v direction (in km/s)
                
            Mass: 'astropy quantity'
                The mass of the particle (in Msun)
    """
    # extract the data array from a predefined function
    data = (Read(filename)[2])
    
    # create new indexing based on the particle type 
    if (particle_type == 1.0): 
        index = np.where(data["type"]==1.0)
    if (particle_type == 2.0):    
        index = np.where(data["type"]==2.0)
    if (particle_type == 3.0):    
        index = np.where(data["type"]==3.0)
    
    # define new arrays with reqiured position values based on new index
    xnew = data['x'][index]
    ynew = data['y'][index]
    znew = data['z'][index]
    # get the position value in each direction (x,y,z), and round it to 3 decimals
    x = np.around(xnew[particle_number-1], 3) # index of [particle_number-1] as 1 particle in the array has index 0
    y = np.around(ynew[particle_number-1], 3) 
    z = np.around(znew[particle_number-1], 3) 
    
    dist = (x**2+y**2+z**2)**(1/2) # compute the 3 dimensional distance of the particle
    Distance = float(dist) * u.kpc # convert to astropy quantity with units in kpc
   
    # define new arrays with reqiured velocity values based on new index
    vxnew = data['vx'][index]
    vynew = data['vy'][index]
    vznew = data['vz'][index]
    # get the velocity value in each direction (x,y,z), and round it to 3 decimals
    vx = np.around(vxnew[particle_number-1], 3) 
    vy = np.around(vynew[particle_number-1], 3) 
    vz = np.around(vznew[particle_number-1], 3)
    
    vel = (vx**2+vy**2+vz**2)**(1/2) # compute the 3 dimensional velocity of the particle
    Velocity = float(vel)*u.kilometer/u.second # convert to astropy quantity with units in km/s
    
    M_sun = 1.988e30 # mass of the Sun
    mnew = data['m'][index] # define new arrays with reqiured velocity values based on new index
    converted_mass = mnew[particle_number-1]*1e-10*M_sun # Convert mass of the particle from 1e10 kg to Msun
    Mass = (converted_mass) * u.M_sun # convert to astropy quantity with units in Msun
    
    return(Distance, Velocity, Mass)


# In[25]:


solution = ParticleInfo("MW_000.txt",1.0, 100) 

Distance_in_kpc = ParticleInfo("MW_000.txt",1.0, 100)[0] # get the distance of required particle
Distance_in_lyr = Distance_in_kpc.to(u.lyr) # convert units from kpc to lyrs for the distance
distance = np.around(Distance_in_lyr, 3) # approximate the distance to 3 decimal places

print(solution[0:3]) # solution for question 5 : 
                     # prints the distance(in kpc), velocity(in km/s) and mass (in Msun) of given particle 
print(distance) # distance for the same particle in units light years

