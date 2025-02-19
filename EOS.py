
import CoolProp as CP
from CoolProp.CoolProp import PropsSI
from userConfig import useIdealGas
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
    # P in Pa, T in K, D in mol/m^3
    if useIdealGas:
        r=8.314462 #Pa*m^3/(mol*K)
        #mol/m^3
        return p/(r*t)

    # P in Pa, T in K, D in mol/m^3
    return PropsSI('Dmolar', 'T', t, 'P', p, name)
def calcFugacity(st,P,T):
    P=P*1e6
    st.update(CP.PT_INPUTS, P, T)
    return st.fugacity_coefficient(0)
    
#output in Z, P in MPa, T in K
def calcZ(P,T,gas):
    P=P*1e6
    CP.CoolProp.set_config_bool(CP.CoolProp.OVERWRITE_BINARY_INTERACTION, True)
    if (T<=40):
        pP=0.887
    elif(T<=50):
        pP=0.77    
    elif(T<=77):   
        pP=0.5
    elif(T<=90):
        pP=.4288
    elif(T<=100):
        pP=.386
    elif(T<=120):
        pP=.3296
    else:
        pP=.25
        
    return pP*PropsSI('Z','P',P,'T',T,"parahydrogen")+(1-pP)*PropsSI('Z','P',P,'T',T,"orthohydrogen")
    
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
def molMass(gasName):
    return PropsSI("M",gasName)*1000
