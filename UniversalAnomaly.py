from math import sqrt
import Stumpff

mu = 398600

def compute_chi(r0, vr0, a, delta_t):
    alpha = 1 / a
    chi = sqrt(mu) * abs(alpha) * delta_t
    z = alpha * (chi) ** 2

    ratio = 1
    tolerance = 10**(-6)

    def compute_fchi(r0, vr0, delta_t):
        A = (r0 * vr0 / sqrt(mu)) * (chi) ** 2 * Stumpff.stumpC(z)
        B = (1 - alpha * r0) * (chi) ** 3 * Stumpff.stumpS(z)
        C = r0 * chi - (sqrt(mu) * delta_t)
        return A + B + C

    def compute_fchi_derivative(r0, vr0):
        D = (r0 * vr0 / sqrt(mu)) * (chi) * (1 - (alpha * (chi) ** 2 * Stumpff.stumpS(z)))
        E = (1 - alpha * r0) * (chi) ** 2 * Stumpff.stumpC(z) + r0
        return D + E

    while abs(ratio) > tolerance:
        ratio = compute_fchi(r0, vr0, delta_t) / compute_fchi_derivative(r0, vr0)
        chi = chi - ratio

    return chi

def main():
    chi = compute_chi(r0=14000, vr0=-2.6679, a=14000, delta_t=3600)
    #chi = compute_chi(r0=10000, vr0=3.0752, a=-19655, delta_t=3600)
    print(chi)

if __name__ == "__main__":
    main()
