import numpy as np
import astropy.units as u
import astropy.constants as c

R=2*u.jupiterRad
Mdot=(1e-8)*u.jupiterMass/u.a


def properties_calculator(M):
    M=M*u.jupiterMass
    L=c.G*M*Mdot/R
    T=(L/(4*np.pi*c.sigma_sb*R**2))**0.25
    return (T.to(u.K),L.to(u.solLum))

print(properties_calculator(5))
