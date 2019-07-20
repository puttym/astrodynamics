from numpy import linalg, cross, asarray, dot
import math
import OrbitalElements

mu = 398600
tolerance = 10**-4

r1 = [-294.32, 4265.1, 5986.7] 
r2 = [-1365.5, 3637.6, 6346.8] 
r3 = [-2940.3, 2473.7, 6555.8]

def GibbsMethod(r1, r2, r3):
    r1_mag = linalg.norm(r1)
    r2_mag = linalg.norm(r2)
    r3_mag = linalg.norm(r3)
    
    r1_unitvector = r1/r1_mag
    c23_unitvector = cross(r2,r3)/linalg.norm(cross(r2,r3))

    if dot(r1_unitvector, c23_unitvector) <= abs(tolerance):
        N = r1_mag*cross(r2, r3) + r2_mag*cross(r3, r1) + r3_mag*cross(r1, r2)
        D = cross(r1, r2) + cross(r2, r3) + cross(r3, r1)
        S = asarray(r1)*(r2_mag - r3_mag) + asarray(r2)*(r3_mag - r1_mag) + asarray(r3)*(r1_mag - r2_mag)
        v2 = math.sqrt(mu/(linalg.norm(N)*linalg.norm(D))) * ((cross(D, r2)/r2_mag) + S)
        angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly =   OrbitalElements.computeOrbitalElements(r2, v2)
    
    else:
        print("The vectors are not coplar")
        
    return angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly

def main():
    angMomentum, SemiMajorAxis, Period, inclination, RAAN, eccentricity, omega, trueAnomaly = GibbsMethod(r1, r2, r3)
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
