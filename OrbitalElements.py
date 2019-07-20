import numpy as np
import math
mu = 398600
K = [0, 0, 1] #Unit vector in the z-direction

r = [-6045, -3490, 2500]
v = [-3.457, 6.618, 2.533]

def computeOrbitalElements(r, v):
    radius = math.sqrt(np.dot(r,r))
    velocity = math.sqrt(np.dot(v,v))
    radialVelocity = np.dot(r,v)/radius
    
    #Angular momentum
    h = np.cross(r, v)
    angMomentum = math.sqrt(np.dot(h,h))
    
    #Inclination
    inclination = np.degrees(math.acos(h[2]/angMomentum))
    
    #RAAN
    N = np.cross(K, h)
    nodeLine = math.sqrt(np.dot(N,N))
    
    if N[1] >= 0:
        RAAN = np.degrees(math.acos(N[0]/nodeLine))
    else:
        RAAN = 360 - np.degrees(math.acos(N[0]/nodeLine))
        
    #Eccentricity
    e = (((velocity)**2 - (mu/radius))*np.asarray(r) - (radius*radialVelocity*np.asarray(v)))/mu
    eccentricity = math.sqrt(np.dot(e,e))    
    
    #Argument of perigee, omega
    if e[2] >= 0:
        omega = np.degrees(math.acos(np.dot(N,e)/(nodeLine*eccentricity)))
    else:
        omega = 360 - np.degrees(math.acos(np.dot(N,e)/(nodeLine*eccentricity)))

    #True anomaly
    if radialVelocity >= 0:
        trueAnomaly = np.degrees(math.acos(np.dot(e,r)/(eccentricity*radius)))
    else:
        trueAnomaly = 360 - np.degrees(math.acos(np.dot(e,r)/(eccentricity*radius)))
        
    #Semi Major Axis
    SemiMajorAxis = (angMomentum**2/mu) * (1/(1 - eccentricity**2))
    
    #Period
    Period = (2 * math.pi/math.sqrt(mu)) * (SemiMajorAxis)**(3/2)   
    
    return angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly 

def main():
    angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly = computeOrbitalElements(r, v)
    print("Angular Momentum (km^2/s): ", round(angMomentum, 3))
    print("Semi Major Axis (km): ", round(SemiMajorAxis, 3))
    print("Preiod (s): ", round(Period, 3))
    print("Inclination (deg): ", round(inclination, 3))
    print("RAAN (deg): ", round(RAAN, 3))
    print("Eccentricity: ", round(eccentricity, 3))
    print("Argument of Perigee (deg): ", round(omega, 3))
    print("True Anomaly (deg): ", round(trueAnomaly, 3))
    
if __name__ == "__main__":
    main()