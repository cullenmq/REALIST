import numpy as np
from EOS import initFugacity,calcFugacity
from plotData import plotIsoHeat,plotAbsUptake,plotAdsLayer
from EOS import calcDensity
def sqrt(x):
    return np.sqrt(x)
def findAdsDensity(coef):
    #nmax in mol/g, vmax in m3/g
    #if vmax is zero divide by zero error! fix by manually setting adsrho to -1 if this is the case
    try:
        adsrho = coef['nmax']*1e3/coef['vmax']/1e6
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
    nmax = coef['nmax']
    fa1={}
    for temp in Press:
        T=float(temp)
        P=np.array(Press[temp])
        fa1[temp]=1000*(nmax)*isoModel.theta(P,T,coef)#((K1*P*(1+a6)/2+intTerm)/(1+(1+a6)*K1*P+intTerm))
    plotAbsUptake(actualPress,fa1,name,newFigure)
    return fa1

#adjustPressure for fugacity
def adjustPressure(expP,normalPress, useFugacity,gasName):
    actualPress={}
    calcPress={}
    for temp in expP:
            actualPress[temp]=normalPress#
            calcPress[temp]=[]
            st=initFugacity(gasName)
            if(useFugacity):
                for press in actualPress[temp]:
                    calcPress[temp].append(press*calcFugacity(st,press,float(temp)))
            else:
                calcPress[temp]=actualPress[temp]
    return actualPress,calcPress
##calculate volumetric uptake (total and excess)
def calcUptake(ads,coef,gasName,sampBulkDens,tempPress,Xpore,den,useFugacity,isoModel):
    rssr=coef['rssr']
    y_fit={}
    fitDens={}
    actualPress,fitPress=adjustPressure(ads,tempPress,useFugacity=useFugacity,gasName=gasName)
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
        measExcessVolUptake[temp],measTotalVolUptake[temp]=calcVolUptake(ads[temp],den[temp],sampBulkDens,Xpore)
        calcExcessVolUptake[temp],calcTotalVolUptake[temp]=calcVolUptake(y_fit[temp],fitDens[temp],sampBulkDens,Xpore)

        #calculate 5 bar total volume based on calculated

        calcVolUptakeTemp=calcTotalVolUptake[temp]
        totVol5bar[temp]= calcVolUptakeTemp[press5bar]

    return y_fit,actualPress,fitPress,rssr,measTotalVolUptake, calcTotalVolUptake,totVol5bar

# calculate volumetric excess and volumetric total:
def calcVolUptake(exUptake,gasRho,bulkDens,Xpore):
    exUptake=np.asarray(exUptake)
    gasRho=np.asarray(gasRho)
    #we need to convert mol/m^3-> mmol/ml
    gasRho=gasRho/1000
    gasSTP=22.413969545014 #methane at STP, mL/mmol
    #excessUptake in mmol/g,excessVolUptake in v/v
    excessVolUptake= exUptake*bulkDens*gasSTP
    #total density is in v/v
    totalVolUptake=excessVolUptake+(gasRho*Xpore)*gasSTP
    return excessVolUptake, totalVolUptake

def calcQstInt(Press,actualPress,coef,adsRho,fa1,name,closeFig,isoModel,gasName):
    Qst={}
    theta={}
    liqmolarv = 0.001 / adsRho
    for temp in Press:
        T=float(temp)
        dens=[]
        for press in Press[temp]:
            dens.append(calcDensity(press,T,gasName))
        dens=np.array(dens)
        P=np.array(Press[temp])
        #Qst=-TdP/dT(vgas-vads)
        #-dP/dT=(dTheta/dP)^-1*dtheta/dT
        Qst[temp]= -T*(1/dens-liqmolarv)*1000*isoModel.dPdT(P,T,coef)
        theta[temp]= isoModel.theta(P,T,coef)
    plotIsoHeat(actualPress, fa1, theta, Qst, name)
    return Qst,theta

