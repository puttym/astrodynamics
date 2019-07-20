import numpy as np
import UniversalAnomaly
import LagrangeCoefficients

mu = 398600

def compute_StateVector(r0, v0, delta_t):
    r0_mag = np.linalg.norm(r0, ord = 2)
    v0_mag = np.linalg.norm(v0, ord = 2)

    radial_velocity = np.dot(r0, v0)/r0_mag
    alpha = (2/r0_mag) - ((v0_mag)**2/mu)

    chi = UniversalAnomaly.compute_chi(r0=r0_mag, vr0=radial_velocity, a=1/alpha, delta_t=delta_t)
    z = alpha * chi**2
    f = LagrangeCoefficients.lagrange_f(r0_mag=r0_mag, chi=chi, z=z)
    g = LagrangeCoefficients.lagrange_g(delta_t=delta_t, chi=chi, z=z)

    r = f*r0 + g*v0
    r_mag = np.linalg.norm(r, ord = 2)

    fdot = LagrangeCoefficients.lagrange_fdot(r_mag=r_mag, r0_mag=r0_mag, alpha= alpha, chi=chi, z=z)
    gdot = LagrangeCoefficients.lagrange_gdot(r_mag=r_mag, chi=chi, z=z)

    v = fdot*r0 + gdot*v0
    return r, v

def main():
    r, v = compute_StateVector(r0=np.array([1600, 5310, 3800]),
                               v0=np.array([-7.350, 0.4600, 2.470]), delta_t=3200)
    print(r)
    print(v)

if __name__ == "__main__":
    main()


