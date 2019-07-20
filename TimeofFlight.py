import math

mu = 398600

#rp = 9600
#ra = 21000
#theta = 120
rp = input("Radius of perigee (in km): ")
ra = input("Radius of apogee (in km): ")
theta = input("Theta (in deg): ")

def ecc(rp, ra):
    e = (ra - rp)/(ra + rp)
    return e

def angmom(rp, e):
    h = math.sqrt(rp * mu * (1 + e))
    return h

def period(h, e):
    T = (2 * math.pi/(mu)**2) * (h/math.sqrt(1 - e**2))**3
    return T

def e_anom(e, theta):
    E = 2 * math.atan(math.sqrt((1 - e)/(1 + e)) * math.tan(math.radians(theta/2)))
    return E

def mean_anom(E, e):
    Me = E - (e*math.sin(E))
    return Me

def tof(Me, T):
    ToF = Me*T/(2*math.pi)
    return ToF

def main():
#   from argparse import ArgumentParser
#   x = ArgumentParser()
#   x.add_argument("r_p", type = float)
#   x.add_argument("r_a", type = float)
#   x.add_argument("theta", type = float)
    e = ecc(rp, ra)
    h = angmom(rp, e)
    E = e_anom(e, theta)
    Me = mean_anom(E, e)
    T = period(h, e)
    ToF = tof(Me, T)
    print(ToF/3600)
    
main()
