
import CoolProp as CP
from CoolProp.CoolProp import PropsSI
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
def initFugacity(gasName):
    st = CP.AbstractState('HEOS', gasName)
    st.specify_phase(CP.CoolProp.iphase_gas)
    st.set_mole_fractions([1])
    return st
