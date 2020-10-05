import numpy as np
import astropy.units as u

def ipos(PA_disk,ri,PAi):

    ############################################################
    #
    # This routine returns the x,y values for the position of a 
    # planet to be included in MCMax3D input file such that the
    # planet gets correclty positionated in the output image. 
    # 
    # IMPORTANT: the phi value in the Image.out file must be zero
    #
    # PA_disk: the position angle of the disk measured from 
    # north to east in deg.
    # ri: radial separation (projected) of the planet in AU.
    # PAi: position angle (projected) of the planet in deg.
    # 
    ############################################################
    
    # Poision vector of the object
    thetai=((PAi+90)*u.deg).to(u.rad).value
    posi=np.array([ri*np.cos(thetai),ri*np.sin(thetai)])

    # Rotation matrix
    theta=((PA_disk-90)*u.deg).to(u.rad).value
    M=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    
    # Clockwise rotation by PAi
    pos_rotated=np.dot(M,posi)
    
    # Converting to MCMax3D coordinates
    x_mcmax=-1.0*pos_rotated[1]
    y_mcmax=pos_rotated[0]


    return (x_mcmax,y_mcmax)

#print(ipos(158.6,21.0,146.8))
#print()
#print(ipos(158.6,38.6,277.0))
 

def CPD_dust_mass(m_p,f):
    # m_p: planet mass (M_jup)
    # f: fraction of the planet mass (%). f x m_p is the disc's mass (gas) 
    
    m_gas=(f/100.0)*m_p
    m_dust=((0.01*m_gas)*u.M_jup).to(u.M_sun)
    
    return m_dust.value

    
    

print(CPD_dust_mass(5,10))
    
