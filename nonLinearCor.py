### This example program shows how non-linear a gas is relative to an ideal gas at a specific T,P

import CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
from plotData import addPlot
import csv
#temp in K, press in Pa
def calcFugacity(st,p,T):
    st.update(CP.PT_INPUTS, p, T)
    return st.fugacity_coefficient(0)



if __name__ == '__main__':
    gasName='Krypton'
    st = CP.AbstractState('HEOS', gasName)

    press=np.arange(100, 10e6, 100)

    st.specify_phase(CP.CoolProp.iphase_gas)
    st.set_mole_fractions([1])
    temps=np.arange(333.3,233.15,-10)
    fug={}
    chem={}
    plt.figure()
    outfile = open("Krypton.csv",'w')
    out = csv.writer(outfile)
    out.writerow(['Press (MPa)','Temp (K)','Fugacity','Chem Potential (J/mol)'])
    for T in temps:

        fug[T]=[]
        chem[T]=[]
        for p in press:
            st.update(CP.PT_INPUTS, p, T)
            fug[T].append(st.fugacity_coefficient(0))
            chem[T].append(st.chemical_potential(0))
            a=[p,T,fug[T][-1],chem[T][-1]]
            out.writerow(a)
        addPlot(press/1e6, fug[T],T,'-',minTemp=270)
        plt.draw()
    outfile.close()
    plt.ylabel("Fugacity Coefficient")
    plt.xlabel("Press (MPa)")
    plt.title(gasName)
    plt.legend()
    plt.savefig(gasName+"Fugacity")

    plt.figure()
    for T in chem:
        addPlot(press/1e6,chem[T],T,'-')
    plt.ylabel("Chemical Potential (J/mol)")
    plt.xlabel("Press (MPa)")
    plt.title(gasName)
    plt.legend()
    plt.savefig(gasName+"chemPot")





    plt.show()


