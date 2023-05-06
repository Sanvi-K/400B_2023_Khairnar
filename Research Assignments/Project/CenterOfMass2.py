#!/usr/bin/env python
# coding: utf-8

# In[43]:


# Modified CenterOfMass.ipynb for Homework 6
# Center of Mass Position and Velocity
# Sanvi Khairnar


# In[44]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[50]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]

        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''

        # xcomponent Center of mass
        a_com = np.sum(a*m)/np.sum(m) 
        # ycomponent Center of mass
        b_com = np.sum(b*m)/np.sum(m)
        # zcomponent Center of mass
        c_com = np.sum(c*m)/np.sum(m)
        
        # return the 3 components separately (use ';')
        return a_com, b_com, c_com; 
    
    
    def COM_P(self, volDec, delta=0.1):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        volDec : `float`
            fraction by which r_max(3D distance of particles from COM position) is reduced 
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        r_COM = (x_COM**2 + y_COM**2 + z_COM**2)**0.5


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = (x_new**2 + y_new**2 + z_new**2)**0.5

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    

        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new <= r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = (x_COM2**2 + y_COM2**2 + z_COM2**2)**0.5

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)                                                                                   

            # Before loop continues, reset : r_max, particle separations and COM                                        
            # reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = (x_new**2 + y_new**2 + z_new**2)**0.5

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below        
        p_COM = np.round(p_COM*u.kpc, 2)
        
        return p_COM
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        xV = (self.x - x_COM/u.kpc) # remove units to perform operations
        yV = (self.y - y_COM/u.kpc)
        zV = (self.z - z_COM/u.kpc)
        rV = ((xV**2 + yV**2 + zV**2)**0.5) * u.kpc # add units back
        
        # determine the index for those particles within the max radius
        indexV = np.where(rV <= rv_max)
        
        # determine the velocity and mass of those particles within the max radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new = self.m[indexV] 
        
        # compute the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        v_COM = np.array([vx_COM, vy_COM, vz_COM])

        # return the COM vector
        # set the correct units using astropy
        # round all values                                                                                        
        return (np.round(v_COM*u.km/u.s, 2))
    


# # 6
# 
# ### 1] 
# Position and Velocity COM for MW, M31 and M33 in kpc and km/s respectively
# 

# In[51]:


# Create a Center of mass object for the MW, M31 and M33

# for all particle type selected is disk type (2.0), as it works best for COM determination
MW_COM = CenterOfMass("MW_000.txt", 2.0)
M31_COM = CenterOfMass("M31_000.txt", 2.0)
M33_COM = CenterOfMass("M33_000.txt", 2.0)


# In[52]:


# Store the position and velocity COM for MW, M31, M33
# print the values

MW_COM_p = MW_COM.COM_P(2.0, 0.1)
print("Position COM for MilkyWay is,", MW_COM_p)
MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
print("Velocity COM for MilkyWay is,", MW_COM_v,"\n")

M31_COM_p = M31_COM.COM_P(2.0,0.1)
print("Position COM for MilkyWay is,", M31_COM_p)
M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
print("Velocity COM for MilkyWay is,", M31_COM_v,"\n")

M33_COM_p = M33_COM.COM_P(2.0,0.1)
print("Position COM for MilkyWay is,", M33_COM_p)
M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
print("Velocity COM for MilkyWay is,", M33_COM_v,"\n")


# ### 2]

# In[56]:


# Compute the distance and relative velocity between MW and M31 in kpc and km/s respectively

# magnitude of position COM
Seperation_MW_M31 = np.sum((MW_COM_p - M31_COM_p)**2)**0.5
print("The seperation between MW and M31 is,", np.round((Seperation_MW_M31),3))


# magnitude of velocity COM 
relvel_MW_M31 = np.sum((MW_COM_v - M31_COM_v)**2)**0.5
print("The relative velocity between MW and M31 is,", np.round((relvel_MW_M31),3))


# ### 3]

# In[57]:


# Compute the distance and relative velocity between M33 and M31 in kpc and km/s respectively

# magnitude of position COM
Seperation_M33_M31 = np.sum((M33_COM_p - M31_COM_p)**2)**0.5
print("The seperation between M33 and M31 is,", np.round((Seperation_M33_M31),3))


# magnitude of velocity COM 
relvel_M33_M31 = np.sum((MW_COM_v - M31_COM_v)**2)**0.5
print("The relative velocity between M33 and M31 is,", np.round((relvel_M33_M31),3))


# ### 4]
# Since the merger of the two galaxies wil happen about their common COM, we first need to know the COM of individual galaxy. Also, it is neccessary to iterate it down to partciles within a certain radius as during the merger, bunch of stars in outer regions of both the galaxy will fling out making it difficult to determine the COM then. Thus calculating COM for paratciles that are relatively stil stataionary about the COM is more effecient.
