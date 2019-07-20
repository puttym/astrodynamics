from math import sqrt, sin, cos, sinh, cosh

def stumpS(z):
    if z > 0:
        S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))**3
    if z < 0:
        S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))**3
    if z == 0:
        S = 1/6
    return S

def stumpC(z):
    if z > 0:
        C = (1 - cos(sqrt(z)))/z
    if z < 0:
        C = (cosh(sqrt(-z)) - 1)/(-z)
    if z == 0:
        C = 1/2
    return C

def main():
    print(stumpS(4.5911))
    print(stumpC(4.5911))

if __name__ == "__main__":
    main()
