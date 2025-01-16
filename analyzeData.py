import numpy as np
from EOS import initFugacity,calcFugacity,molMass
from plotData import plotIsoHeat,plotAbsUptake,plotAdsLayer
from EOS import calcDensity,calcSatPress,critT, calcZ
from parseData import molToWt
from userConfig import *
def sqrt(x):
    return np.sqrt(x)
def findAdsDensity(coef):
    #nmax in mol/g, vmax in m3/g
    #if vmax is zero divide by zero error! fix by manually setting adsrho to -1 if this is the case
    nmax=coef["nmax"]#coef['SSA']/(6.022E23*15.7E-20)
    try:
        adsrho = nmax*1e3/coef['vmax']/1e6 #coef["ρAds"]
    except:
        adsrho=-1
    print("Adsorbed Phase Density = {} mmol/ml".format(adsrho))
    print("Adsorbed Molar Volume = {} m\u00b3/mol".format(0.001/adsrho))
    return adsrho
def AdsLayerGrowth(Press,actualPress,coef,isoModel,name='test'):
    fv1={}
    vmax = coef["vmax"]
    for temp in Press:
        P=Press[temp]
        T=float(temp)
        P=np.array(P)
        fv1[temp]= 1e6 * vmax *isoModel.theta(P,T,coef)
    plotAdsLayer(actualPress, fv1, name)
    return fv1
#press is a dictionary of temperatures to arrays of pressures
def absAds(Press,actualPress,coef,isoModel,name="test",newFigure=False):
    nmax = coef["nmax"]#coef['SSA']/(6.022E23*15.7E-20)
    fa1={}
    for temp in Press:
        T=float(temp)
        P=np.array(Press[temp])
        fa1[temp]=1000*(nmax)*isoModel.theta(P,T,coef)#((K1*P*(1+a6)/2+intTerm)/(1+(1+a6)*K1*P+intTerm))
    plotAbsUptake(actualPress,fa1,name,newFigure)
    return fa1

#adjustPressure for fugacity
def adjustPressure(expP,normalPress, gasName):
    actualPress={}
    calcPress={}
    for temp in expP:
            
            if(float(temp)<critT(gasName) or CUTOFF_THETA):
                actualPress[temp]=[]
                if(CUTOFF_THETA):
                    maxPress=max(expP[temp])
                else:
                    gasSatPress= calcSatPress(float(temp),gasName)
                for press in normalPress:
                    if(CUTOFF_THETA):
                        if( press<=maxPress):
                            actualPress[temp].append(press)
                    elif(press<gasSatPress):            
                        actualPress[temp].append(press)
            else:
                actualPress[temp]=normalPress
            calcPress[temp]=[]
            st=initFugacity(gasName)
            if(useFugacity):
                for press in actualPress[temp]:
                    calcPress[temp].append(press*calcFugacity(st,press,float(temp)))
            elif(useCompress):
                for press in actualPress[temp]:
                    calcPress[temp].append(press*calcZ(press,float(temp),gasName))
            else:
                calcPress[temp]=actualPress[temp]
    return actualPress,calcPress


#find pressure at a certain uptake
def findPress(press,uptake,value):
    array = np.asarray(uptake)
    idx = (np.abs(array - value)).argmin()
    return press[idx]
##calculate isoexcess heat of adsorption, fitExcess in wt%, press in bar, uptake in mmol/g
def calcIsoExcess(fitExcess,press,uptake,T2=87,T1=77):
    #first find name of relevant isotherm
    T2Name=''
    T1Name=''
    isoHeat=[]
    r = 8.314462 /1000
    #find name that is within 2K of desired temp
    for names in press:

        if(abs(float(names)-T1)<2):
            T1Name=names
            t1n=float(T1Name)
        elif(abs(float(names)-T2)<2):
            T2Name=names
            t2n=float(T2Name)
    wtPer=molToWt(uptake,"hydrogen")
    #find pressure index of max excess uptake
    wt1Idx=np.argmax(wtPer[T1Name])
    wt2Idx=np.argmax(wtPer[T2Name])
    for ex in fitExcess:
        #we need to cut off the press, excess data at max excess so the function is well defined
        P1=findPress(press[T1Name][:wt1Idx], wtPer[T1Name][:wt1Idx], ex)
        P2=findPress(press[T2Name][:wt2Idx], wtPer[T2Name][:wt2Idx], ex)
        #print("T1: {}, T2: {}, P1: {}, P2: {}".format(t1n,t2n,P1,P2))
        isoHeat.append(r*np.log(P1/P2)*(t2n*t1n/(t1n-t2n)))
    return isoHeat
##calculate volumetric uptake (total and excess)
def calcUptake(ads,coef,gasName,sampBulkDens,tempPress,Xpore,den,isoModel):
    rssr=coef['rssr']
    y_fit={}
    fitDens={}
    actualPress,fitPress=adjustPressure(ads,tempPress,gasName=gasName)
    for temp in ads:
        tempFitDen=[]
        for i,press in enumerate(actualPress[temp]):
            tempFitDen.append(calcDensity(press,float(temp),gasName))
        fitDens[temp]=tempFitDen
        y_fit[temp] = isoModel.genExcess(fitPress[temp], t=float(temp),a=coef,den=fitDens[temp])
    measExcessVolUptake={}
    measTotalVolUptake={}
    calcExcessVolUptake={}
    calcTotalVolUptake={}
    press5bar=-1
    for p in range(len(tempPress)):
                if(round(tempPress[p],2)==0.50):
                    press5bar=p
    totVol5bar={}
    for temp in ads:
        measExcessVolUptake[temp],measTotalVolUptake[temp]=calcVolUptake(ads[temp],den[temp],sampBulkDens,Xpore,gasName)
        calcExcessVolUptake[temp],calcTotalVolUptake[temp]=calcVolUptake(y_fit[temp],fitDens[temp],sampBulkDens,Xpore,gasName)

        #calculate 5 bar total volume based on calculated

        calcVolUptakeTemp=calcTotalVolUptake[temp]
        totVol5bar[temp]= calcVolUptakeTemp[press5bar]

    return y_fit,actualPress,fitDens,fitPress,rssr,measExcessVolUptake,measTotalVolUptake, calcTotalVolUptake,totVol5bar

# calculate volumetric excess and volumetric total:
def calcVolUptake(exUptake,gasRho,bulkDens,Xpore,gasName="Methane"):
    exUptake=np.asarray(exUptake)
    gasRho=np.asarray(gasRho)
    #we need to convert mol/m^3-> mmol/ml
    gasRho=gasRho/1000
    gasSTP=1/calcDensity(0.101325,273.15,gasName)/1000 #22.413969545014 #L/mol

    #excessUptake in mmol/g,excessVolUptake in v/v
    #excessVolUptake= exUptake*bulkDens*gasSTP
    #total density is in v/v
    #totalVolUptake=excessVolUptake+(gasRho*Xpore)*gasSTP
    
    #excessUptake in mmol/g,excessVolUptake in g/L
    excessVolUptake= exUptake*bulkDens*molMass(gasName)#*1000
    #total density is in g/L
    totalVolUptake=excessVolUptake+(gasRho*Xpore)*molMass(gasName)#*1000
    
    return excessVolUptake, totalVolUptake

def calcQstInt(fitPress,actualPress,coef,adsRho,fa1,name,closeFig,isoModel,gasName):
    Qst={}
    theta={}

    for temp in fitPress:
        T=float(temp)
        dens=[]


        for press in actualPress[temp]:
            dens.append(calcDensity(press,T,gasName))
        dens=np.array(dens)
        P=np.array(fitPress[temp])
        theta[temp] = isoModel.theta(P, T, coef)
        nmax = coef["nmax"]
        vmax=coef['vmax']
        if(ADS_VOL== "Static"): #static volume means dynamic adsorbed density (gas-like)
            n = nmax * theta[temp]
        elif (ADS_VOL== "Dynamic"): #dynamic volume with constant adsorbed density (liquid-like)
            n = nmax
        try:
            adsrho = n * 1e3 / vmax / 1e6  # coef["ρAds"]
        except:
            adsrho = -1
        liqmolarv = 0.001 / adsRho

        #Qst=-TdP/dT(vgas-vads)
        #-dP/dT=(dTheta/dP)^-1*dtheta/dT
        Qst[temp]= -T*(1/dens-liqmolarv)*1000*isoModel.dPdT(P,T,coef)

    plotIsoHeat(actualPress, fa1, theta, Qst, name)
    return Qst,theta

