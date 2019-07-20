from math import sin, radians, cos, sqrt, acos, exp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mu = 398600

r1 = [5000, 10000, 2100]
r2 = [-14600, 2500, 7000]
delta_t = 3600

z = np.arange(-2, 5.5, 0.2)

r1_mag = np.linalg.norm(r1, ord=2)
r2_mag = np.linalg.norm(r2, ord=2)

def Compute_A(r1, r2):
    if np.cross(r1, r2)[2] >= 0: #z-component of np.cross(r1, r2)
        delta_theta = np.degrees(acos((np.dot(r1,r2)/(r1_mag*r2_mag))))
    else:
        delta_theta = 360 - np.degrees(cos((np.dot(r1,r2)/(r1_mag*r2_mag))))
    
    A = sin(radians(delta_theta)) * sqrt(r1_mag*r2_mag/(1-cos(radians(delta_theta))))
    return A

def Stumpff_S(z):
    S = (1/6) - (z/120) + (z**2/5040) - (z**3/362880) + (z**4/39916800)
    return S

def Stumpf_C(z):
    C = (1/2) - (z/24) + (z**2/720) - (z**3/40320) + (z**4/3628800)
    return C

def main():
    A = Compute_A(r1, r2)
    y1 = r1_mag + r2_mag + A * (((z * Stumpff_S(z)) - 1)/np.sqrt(Stumpf_C(z)))
    F = (y1/Stumpf_C(z))**(1.5) * Stumpff_S(z) + (A * np.sqrt(y1)) - (delta_t * np.sqrt(mu))
    
    fig, ax = plt.subplots()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True, useOffset=False)
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('F(z)')
    ax.plot(z, F)
    plt.show()
    fig.savefig('Lambert.pdf')
    
if __name__ == '__main__':
    main()