import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'

def StumpNeg(z):
    S_Neg = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
    C_Neg = (np.cosh(np.sqrt(-z)) - 1)/(-z)
    return S_Neg, C_Neg

def StumpPos(z):
    S_Pos = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    C_Pos = (1 - np.cos(np.sqrt(z)))/z
    return S_Pos, C_Pos

def main():
    z_Neg = np.linspace(-50, -1, 100)
    z_Pos = np.linspace(0.001, 30, 100)
    z_Pos1 = np.linspace(0.001, 500, 1000)
    
    S_Neg, C_Neg = StumpNeg(z_Neg)
    S_Pos, C_Pos = StumpPos(z_Pos)
    S_Pos1, C_Pos1 = StumpPos(z_Pos1)
    
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    
    ax1.set_xlabel('z')
    ax2.set_xlabel('z')
    ax3.set_xlabel('z')

    ax3.set_ylim(0, 0.04)
    ax1.plot(z_Neg, S_Neg, label = "S(z)")
    ax1.plot(z_Neg, C_Neg, label = "C(z)")
    ax1.legend()
    plt.tight_layout()

    ax2.plot(z_Pos, S_Pos, label = "S(z)")
    ax2.plot(z_Pos, C_Pos, label = "C(z)")
    ax2.legend()
    plt.tight_layout()
    
    ax3.plot(z_Pos1, S_Pos1, label = "S(z)")
    ax3.plot(z_Pos1, C_Pos1, label = "C(z)")
    ax3.legend()
    plt.tight_layout()

    fig1.savefig('stumpff1.png')
    fig2.savefig('stumpff2.png') 
    fig3.savefig('stumpff3.png')

if __name__ == "__main__":
    main()
