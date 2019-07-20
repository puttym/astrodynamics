from math import pi, sin, cos, exp

"""This programe implements Newton's method for solving the Kepler equation.
   The program takes eccentricity and mean anomaly as inputs and gives eccentric
   anomaly as output."""

tolerance = exp(-8)

def compute_E(e, Me):
    if Me > pi:
        E = Me - (e/2)
    if Me < pi:
        E = Me + (e/2)
        
    ratio = 1;
    
    while abs(ratio) > tolerance:
        ratio = (E - e * sin(E) - Me)/(1 - e * cos(E))
        E = E - ratio
        
    return E

def main():
    from argparse import ArgumentParser

    x = ArgumentParser()
    x.add_argument("eccentricity", type = float)
    x.add_argument("mean_anomaly", type = float)
    parameters = x.parse_args()
    E = compute_E(parameters.eccentricity, parameters.mean_anomaly)
    print('Eccentricity (given) = ', parameters.eccentricity)
    print('Mean Anomaly (given) = ', parameters.mean_anomaly)
    
    print('Eccentric Anomaly (computed) = ', E)
    
if __name__ == '__main__':
    main()
