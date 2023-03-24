
import CoolProp as CP
from CoolProp.CoolProp import PropsSI

def critT(gas):
    return PropsSI("Tcrit",gas)
#Input: T in K
#output: P in MPa
def calcSatPress(t,gas):
    return PropsSI('P','T',t,'Q',1,gas)/1e6
#inputs: P in MPa, T in K
#outputs: D in mol/m^3
def calcDensity(p,t,gas):
    #change HEOS to refprop to use refprop instead
    name="HEOS::"+gas
    #convert MPa to Pa
    p=p*1e6
    #P in Pa, T in K, D in mol/m^3
    return PropsSI('Dmolar', 'T', t, 'P', p, name)
def calcFugacity(st,P,T):
    P=P*1e6
    st.update(CP.PT_INPUTS, P, T)
    return st.fugacity_coefficient(0)
#output in kJ/mol
def calcU(st,P,T):
    P=P*1e6
    st.update(CP.PT_INPUTS, P, T)
    st.chemical_potential(0)/1000
def initFugacity(gasName):
    st = CP.AbstractState('HEOS', gasName)
    st.specify_phase(CP.CoolProp.iphase_gas)
    st.set_mole_fractions([1])
    return st
