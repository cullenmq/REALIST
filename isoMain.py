import numpy as np
from parseData import loadData,grabAdsorption,saveFile,convBartoMPa
from analyzeData import  absAds,calcQstInt,findAdsDensity,AdsLayerGrowth,calcUptake,adjustPressure
from plotData import plotTotal,plotVolTotal,plotExcessUptake
from pathlib import Path
import os
import matplotlib.pyplot as plt
from isotherm import isotherm
def startRun(names,bulkDens,skelDens,Vpore,gasName,model="dL",useFugacity=True,RECALC_FITS=False,CLOSE_FIGS=True,numThreads=48):
    #excess data, T,P,n
    TAll={}
    PAll={}
    adsAll={}
    adjPressAll={}
    #analysis data
    #fitted data
    yFitAll={}
    adjFitPressAll={}
    rssrAll={}
    #volumetric data
    totalVolAdsAll={}
    totalVolAdsFitAll={}
    #save the 5 bar value
    totalVol5bar={}
    adsRhoAll={}
    totVolZeroC={}
    totVolZeroCfit={}
    pZeroC={}
    totVolRT={}
    totVolRTfit={}
    pRT={}
    absUptake={}
    deltaH={}
    thetaAll={}
    Xpore={}
    #fitting Data
    coefAll={}
    denAll={}
    homeDir=os.getcwd()
    for name in names:
        os.chdir(homeDir)
        #Step 0: Make a separate folder for each sample results
        Path(name+'Results'+"/"+gasName).mkdir(parents=True, exist_ok=True)
        os.chdir(name+'Results'+"/"+gasName)
        isBar= True
        #####Step 1: Load Excess Uptake Data from csv file########
        TAll[name], PAll[name], adsAll[name]=loadData(fileName=homeDir+'/rawData/'+name+'Excess.csv', isBar=isBar)
        #PAll[name], adsAll[name] =grabAdsorption(adsAll[name], PAll[name])

        #####Step 2: Run fit (or grab existing fit from file)########
        isoTest=isotherm(model,numThreads=numThreads)
        coefAll[name],denAll[name],adjPressAll[name]=isoTest.runFit(PAll[name], adsAll[name],gasName=gasName,useFugacity= useFugacity,isRunFit=RECALC_FITS, name=name)

        ######Step 3: Run analysis on fitted data
        #3a: find adsorption density, void fraction
        adsRhoAll[name]=findAdsDensity(coefAll[name])
        #Xpore Common
        Xpore[name]=Vpore[name]*bulkDens[name]
        #Xpore,best (skeletal)
        #Xpore[name]=1-bulkDens[name]/skelDens[name]
        #Xpore swollen
        #Xpore[name]=Vpore[name]/(Vpore[name]+1/skelDens[name])
        #3b: calculate volumetric uptake
        tempPress=np.arange(0.0001, 10.0, 0.01)
        #tempPress=np.arange(0.0001, 0.1, 0.001)
        yFitAll[name],fitPress,adjFitPressAll[name],rssrAll[name],totalVolAdsAll[name], totalVolAdsFitAll[name],totalVol5bar[name]=calcUptake(
                                                                                              ads=adsAll[name],coef=coefAll[name],gasName=gasName,
                                                                                              sampBulkDens=bulkDens[name],tempPress=tempPress,
                                                                                              Xpore=Xpore[name],den=denAll[name],
                                                                                              useFugacity=useFugacity,isoModel=isoTest)

        plotExcessUptake(adsAll[name],PAll[name],yFitAll[name],fitPress,name+model)
        #3c: save RT, 0C data for DOE Targets
        totalVolAds=totalVolAdsAll[name]
        totalVolAdsFit=totalVolAdsFitAll[name]
        newP=PAll[name]
        for tempName in totalVolAds:
            if(abs(int(float(tempName))-273)<5):
                totVolZeroC[name]=totalVolAds[tempName]
                totVolZeroCfit[name]=totalVolAdsFit[tempName]
                pZeroC[name]=newP[tempName]
            if(abs(int(float(tempName))-298)<5):
                totVolRT[name]=totalVolAds[tempName]
                totVolRTfit[name]=totalVolAdsFit[tempName]
                pRT[name]=newP[tempName]

        ######Step 4: Calculate Thermodynamic Quantities (Absolute Uptake. Qst,Theta) TODO: Add Entropy
        #Step 4a: set up calcPress and actualPress based on a log scale
        #calcPress is usually fugacity ("input" used to correct for nonideal gas behaviour)
        #we plot against "actual press". If our simplistic models took nonidealities into account we wouldnt need calcPress

        tempPress=np.logspace(-4,1,50)#np.arange(0.0001, 10.0, 0.05)
        actualPress,calcPress=adjustPressure(PAll[name],tempPress,useFugacity=useFugacity,gasName=gasName)


        #Step 4b: Calculate absolute uptake
        plotExcessUptake(adsAll[name],PAll[name],yFitAll[name],fitPress,name+model)
        fa1=absAds(calcPress,actualPress,coefAll[name],isoModel=isoTest,name=name+model)
        absUptake[name]=fa1
        #adsorbed layer growth
        #AdsLayerGrowth(calcPress,coefAll[name],isSqRoot=USE_SQRT,name=name)
        #Step 4c: calculate Isosteric Heat of Adsorption
        deltaH[name],thetaAll[name]=calcQstInt(calcPress,actualPress,coefAll[name], adsRhoAll[name],fa1,isoModel=isoTest,closeFig=CLOSE_FIGS,name=name+model,gasName=gasName)

        ########Step 5: Save data to be read by plotting code
        params=isoTest.names
        info={"Sample Name":name+"_"+model}
        for a in params:
            info.update({a:coefAll[name][a]})
        info.update({'adsRho (mmol/ml)':adsRhoAll[name]})
        info.update({'RSSR':rssrAll[name]})
        data={'Raw Press (MPa)':newP,'Raw Adjusted Press (MPa)':adjPressAll[name],'Excess (mmol/g)':adsAll[name],'Total Vol (V/V)':totalVolAds,'fitPress (MPa)':fitPress,'Adjusted Fit Press (MPa)':adjFitPressAll[name],'Excess Fit (mmol/g)':yFitAll[name],'Total Vol Fit (V/V)':totalVolAdsFit, 'adjustedModelPress':calcPress,'actualPress':actualPress,
              'Absolute (mmolg)':fa1,'Isosteric Heat (kJ/mol)':deltaH[name],'fill factor':thetaAll[name]}

        saveFile(info,data)
        if (CLOSE_FIGS):
            plt.close('all')
        else:
            plt.show(block=True)
    os.chdir(homeDir)
    return PAll, adsAll, yFitAll, absUptake, deltaH, thetaAll, totalVolAdsAll, totVolZeroC, totVolZeroCfit, pZeroC, totVolRT, totVolRTfit, pRT, totalVol5bar


if __name__ == '__main__':
    names=[]
    gasName="Hydrogen"
    #gasName="CarbonDioxide"
    useFugacity= True
    #parameter to force recalculating of fits, even if fit was already done
    RECALC_FITS=False
    #Closes figures automatically, otherwise you have to manually close them to continue 
    CLOSE_FIGS=False
    #names.append("ET095")
    #names.append("ET094")
    names.append("P2H2")
    bulkDens={"MSC30":0.170,"ET094":0.170,"ET095":0.163,"EH046":0.142,"P2H2":0.142}
    Vpore={"ET094":1.13,"ET095":1.22,"EH046":1.66,"MSC30":1.5,"P2H2":1.63}
    SSA={"ET094":1994,"ET095":2108,"EH046":3806,"MSC30":3350,"P2H2":3300}
    skelDens={"MSC30":2.1,"ET094":2.03,"ET095":2.05,"EH046":1.73,"P2H2":1.73}
    models = []
    #models.append("dL")
    #models.append("coAds")
    models.append("unilan")
    #models.append("unilanPurewal")
    #models.append("toth")
    #models.append("sips")
    for model in models:
        print("Using model: {}".format(model))
        startRun(names=names,bulkDens=bulkDens,skelDens=skelDens,Vpore=Vpore,useFugacity=useFugacity,RECALC_FITS=RECALC_FITS,CLOSE_FIGS=CLOSE_FIGS,gasName=gasName,model=model,numThreads=6)