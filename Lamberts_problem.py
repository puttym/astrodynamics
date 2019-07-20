
# coding: utf-8

from math import sin, radians, cos, sqrt, acos, exp
import numpy as np
import Stumpff
import OrbitalElements

tol = exp(-8)

mu = 398600

r1 = [5000, 10000, 2100]
r2 = [-14600, 2500, 7000]
delta_t = 3600

r1_mag = np.linalg.norm(r1, ord=2)
r2_mag = np.linalg.norm(r2, ord=2)

#Compute delta_theta
def Lambert(r1, r2, delta_t):
    if np.cross(r1, r2)[2] >= 0: #z-component of np.cross(r1, r2)
        delta_theta = np.degrees(acos((np.dot(r1,r2)/(r1_mag*r2_mag))))
    else:
        delta_theta = 360 - np.degrees(cos((np.dot(r1,r2)/(r1_mag*r2_mag))))
    
    A = sin(radians(delta_theta)) * sqrt(r1_mag*r2_mag/(1-cos(radians(delta_theta))))
    
    z = 1.5
    
    def F1(z):
        y = r1_mag + r2_mag + A * (z * Stumpff.stumpS(z) - 1)/sqrt(Stumpff.stumpC(z))
        F = (y/Stumpff.stumpC(z))**(1.5) * Stumpff.stumpS(z) + (A * sqrt(y)) - (delta_t * sqrt(mu))
        x1 = (1/(2*z)) * (Stumpff.stumpC(z) - (3 * Stumpff.stumpS(z)/(2 * Stumpff.stumpC(z))))
        x2 = (3 * (Stumpff.stumpS(z))**2)/(4 * Stumpff.stumpC(z))
        x3 = 3 * (Stumpff.stumpS(z)/Stumpff.stumpC(z)) * sqrt(y)
        x4 = A * sqrt(Stumpff.stumpC(z)/y)
        F_prime = (y/Stumpff.stumpC(z))**(1.5) * (x1 + x2) + ((A/8) * (x3 + x4))
        return F/F_prime
    
    ratio = 1
    
    while abs(ratio) > tol: #Newton's Method
        ratio = F1(z)
        z = z - ratio
        
    y = r1_mag + r2_mag + A * (z * Stumpff.stumpS(z) - 1)/sqrt(Stumpff.stumpC(z))

    f = 1 - (y/r1_mag)
    g = A * sqrt(y/mu)
    g_dot = 1 - (y/r2_mag)

    v1 = (1/g) * (r2 - f*np.asarray(r1))
    v2 = (1/g) * (g_dot*np.asarray(r2) - r1)
        
    angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly = OrbitalElements.computeOrbitalElements(r1, v1)
        
    return angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly
    
def main():
    angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly = Lambert(r1, r2, delta_t)
    print("Angular Momentum (km^2/s): ", round(angMomentum, 3))
    print("Semi Major Axis (km): ", round(SemiMajorAxis, 3))
    print("Preiod (s): ", round(Period, 3))
    print("Inclination (deg): ", round(inclination, 3))
    print("RAAN (deg): ", round(RAAN, 3))
    print("Eccentricity: ", round(eccentricity, 3))
    print("Argument of Perigee (deg): ", round(omega, 3))
    print("True Anomaly (deg): ", round(trueAnomaly, 3))
    
if __name__ == '__main__':
    main()