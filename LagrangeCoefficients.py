from math import sqrt
import Stumpff

mu = 398600

def lagrange_f(r0_mag, chi, z):
    f = 1 - (chi**2 * Stumpff.stumpC(z)/r0_mag)
    return f

def lagrange_g(delta_t, chi, z):
    g = delta_t - (chi**3 * Stumpff.stumpS(z)/sqrt(mu))
    return g

def lagrange_fdot(r_mag, r0_mag, alpha, chi, z):
    fdot = (sqrt(mu)/(r_mag * r0_mag)) * (alpha * chi**3 * Stumpff.stumpS(z) - chi)
    return fdot

def lagrange_gdot(r_mag, chi, z):
    gdot = 1 - (chi**2 * Stumpff.stumpC(z)/r_mag)
    return gdot


def main():
    f = lagrange_f(r0_mag, chi, z)
    g = lagrange_g(delta_t, chi, z)
    fdot = lagrange_fdot(r_mag, r0_mag, alpha, chi, z)
    gdot = lagrange_gdot(r_mag, chi, z)

if __name__ == "__main__":
    main()