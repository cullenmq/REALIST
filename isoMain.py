import numpy as np
from parseData import loadData,grabAdsorption,saveFile,convBartoMPa,loadSampleData, wtToMols, molToWt
from analyzeData import  absAds,calcQstInt,findAdsDensity,AdsLayerGrowth,calcUptake,adjustPressure,calcIsoExcess
from plotData import plotTotal,plotVolTotal,plotExcessUptake,plotExcessUptakeParams
from pathlib import Path
import os
import matplotlib.pyplot as plt
from userConfig import *
from isotherm import isotherm
from scipy.interpolate import CubicSpline
def startRun(names,gasName,model="dL",sampleParams=None,RECALC_FITS=False,CLOSE_FIGS=True,numThreads=48):
    #excess data, T,P,n
    TAll={}
    PAll={}
    adsAll={}
    adjPressAll={}
    #analysis data
    #fitted data
    yFitAll={}
    yFitWtAll={}
    adjFitPressAll={}
    rssrAll={}
    #volumetric data
    exVolAdsAll={}
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
    denFitAll={}
    wtAll={}
    nAbsAll={}
    isoEx={}
    homeDir=os.getcwd()
    for name in names:
        os.chdir(homeDir)
        #Step 0: Make a separate folder for each sample results
        Path(name+'Results'+"/"+gasName).mkdir(parents=True, exist_ok=True)
        os.chdir(name+'Results'+"/"+gasName)
        isBar= True
        #####Step 1a: Load Sample Data from JSON file########
        #only load from Json if sample params has not been passed in as a dict
        if sampleParams is None:
            sampleParams=loadSampleData(fileName=homeDir+'/rawData/'+name+'.json')
        bulkDens=sampleParams["bulkDens"]
        skelDens=sampleParams["skelDens"]
        SSA=sampleParams["SSA"]
        Vpore=sampleParams["Vpore"]
        #####Step 1b: Load Excess Uptake Data from csv file########
        TAll[name], PAll[name], adsAll[name],isAbsolute=loadData(fileName=homeDir+'/rawData/'+name+gasName, isBar=isBar)
        if (OVERRIDE_ABSOLUTE):
            isAbsolute=True
        PAll[name], adsAll[name] =grabAdsorption(adsAll[name], PAll[name])
        if(CONVERT_WT_MOL):
            wtAll[name]=adsAll[name]
            adsAll[name]=wtToMols(adsAll[name],gasName)
        else:
            wtAll[name]=molToWt(adsAll[name],gasName)

        #####Step 2: Run fit (or grab existing fit from file)########
        isoTest=isotherm(model,numThreads=numThreads,isAbsolute=isAbsolute)
        coefAll[name],denAll[name],adjPressAll[name],adsAll[name]=isoTest.runFit(PAll[name], adsAll[name],gasName=gasName,isRunFit=RECALC_FITS, name=name)

        ######Step 3: Run analysis on fitted data
        #3a: find adsorption density, void fraction
        adsRhoAll[name]=findAdsDensity(coefAll[name])
        #Xpore Common
        #Xpore[name]=Vpore*bulkDens
        #XporeName="XCommon"
        #Xpore,best (skeletal)
        Xpore[name]=1-bulkDens/skelDens
        XporeName="XBest"
        #Xpore swollen
        #Xpore[name]=Vpore/(Vpore+1/skelDens)
        #XporeName="XSwollen"
        #3b: calculate volumetric uptake
        tempPress=np.arange(0.0001, MAX_PLOT_PRESS, STEP_PRESS)
        #tempPress=np.arange(0.0001, 0.1, 0.001)
        yFitAll[name],fitPress,denFitAll[name],adjFitPressAll[name],rssrAll[name],exVolAdsAll[name],totalVolAdsAll[name],totalVolAdsFitAll[name],totalVol5bar[name]=calcUptake(
                                                                                              ads=adsAll[name],coef=coefAll[name],gasName=gasName,
                                                                                              sampBulkDens=bulkDens,tempPress=tempPress,
                                                                                              Xpore=Xpore[name],den=denAll[name],
                                                                                              isoModel=isoTest)
        


        yFitWtAll[name]=molToWt(yFitAll[name],gasName)
        plotExcessUptakeParams(adsAll[name],PAll[name],coefAll[name],adjFitPressAll[name],denFitAll[name],isoTest, fitPress,name)
        #plotExcessUptake(adsAll[name],PAll[name],yFitAll[name],fitPress,name+model)
        #3c: save RT, 0C data for DOE Targets
        totalVolAds=totalVolAdsAll[name]
        totalVolAdsFit=totalVolAdsFitAll[name]
        newP=PAll[name]
        
        
        
        #Code for calcuating total volumetric uptake at 0C ant 25C, generally for methane
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
        #New code to calc absolute adsorption
        actualPress, calcPress = adjustPressure(PAll[name], tempPress, gasName=gasName)
        nAbsAll[name] = absAds(calcPress, actualPress, coefAll[name], isoModel=isoTest, name=name + model)
        nAbsWtAll=molToWt(nAbsAll[name],gasName)

        #calc isoexcess
        nExCalc={}
        nAbCalc = {}
        prevTemp=-1
        isoEx={}
        isoSt={}
        for temp in fitPress:
            if prevTemp==-1:
                nExCalc[temp]=[]
                nAbCalc[temp]=[]
                isoEx[temp]=[]
                isoSt[temp]=[]
                prevTemp=temp
                continue
            else:
                #only fit up to the maximum measured excess uptake (where the isotherm curves over)
                nExCalc[temp]=np.arange(0.1, min(max(wtAll[name][temp]),max(wtAll[name][prevTemp])), .1)

                nAbCalc[temp]=np.arange(0.1, min(max(nAbsWtAll[temp]),max(nAbsWtAll[prevTemp])), .1)
                cs={}
                absAd={}
                if(USE_SPLINE):
                        cs[temp]=CubicSpline(PAll[name][temp],wtAll[name][temp])(fitPress[temp])
                        cs[prevTemp]=CubicSpline(PAll[name][prevTemp],wtAll[name][prevTemp])(fitPress[prevTemp])
                        
                else:
                    cs=yFitAll[name]
                    absAd=nAbsAll[name]
                isoEx[temp]=calcIsoExcess(nExCalc[temp],fitPress,cs,float(temp),float(prevTemp))
                isoSt[temp] = calcIsoExcess(nAbCalc[temp], fitPress, absAd, float(temp), float(prevTemp))
        

        
        
#Step 4b: Calculate absolute uptake
        tempPress=np.logspace(-4,np.log10(MAX_PLOT_PRESS),50)#np.arange(0.0001, 10.0, 0.05)
        actualPress,calcPress=adjustPressure(PAll[name],tempPress,gasName=gasName)
        plotExcessUptake(adsAll[name],PAll[name],yFitAll[name],fitPress,name+model)
        fa1=absAds(calcPress,actualPress,coefAll[name],isoModel=isoTest,name=name+model)
        absUptake[name]=fa1
        
        
        
        #adsorbed layer growth
        #AdsLayerGrowth(calcPress,coefAll[name],isSqRoot=USE_SQRT,name=name)
        #Step 4c: calculate Isosteric Heat of Adsorption
        deltaH[name],thetaAll[name]=calcQstInt(calcPress,actualPress,coefAll[name], adsRhoAll[name],fa1,isoModel=isoTest,closeFig=CLOSE_FIGS,name=name+model,gasName=gasName)
        #convert absolute uptake to wt%
        abswt=molToWt(fa1,gasName)
        ########Step 5: Save data to be read by plotting code
        params=isoTest.names
        info={"Sample Name":name+"_"+model}
        info.update({"Xpore Model":XporeName})
        for i,a in enumerate(params):
            info.update({a:[coefAll[name][a],str(isoTest.bounds[i])]})
        info.update({'adsRho (mmol/ml)':adsRhoAll[name]})
        info.update({'RSSR':rssrAll[name]})
        info.update({'T_Pre':preFact})
        info.update({'Vol': ADS_VOL})
        data={'Raw Press (MPa)':newP,'Raw Adjusted Press (MPa)':adjPressAll[name],'Excess (mmol/g)':adsAll[name],'Excess (wt%)':wtAll[name],'Excess Vol (g/L)':exVolAdsAll[name],
              'Total Vol (g/L)':totalVolAdsAll[name],'fitPress (MPa)':fitPress,'Adjusted Fit Press (MPa)':adjFitPressAll[name],'Gas Phase Density':denFitAll[name],'Excess Fit (mmol/g)':yFitAll[name],
              'Excess Fit (wt%)':yFitWtAll[name],'Total Vol Fit (g/L)':totalVolAdsFit, 'adjustedModelPress':calcPress,'actualPress':actualPress,'Absolute (mmolg)':fa1,'Absolute (wt%)':abswt,
              'fill factor':thetaAll[name],'Isosteric Heat (kJ/mol)':deltaH[name], 'excess uptake (wt%)':nExCalc, 'isoExcess (kJ/mol)':isoEx,'absolute uptake (wt%)':nAbCalc,'isosteric (kJ/mol)':isoSt}

        saveFile(info,data)
        if (CLOSE_FIGS):
            plt.close('all')
        else:
            plt.show(block=True)
    os.chdir(homeDir)
    return PAll, adsAll, yFitAll, absUptake, deltaH, thetaAll, totalVolAdsAll, totVolZeroC, totVolZeroCfit, pZeroC, totVolRT, totVolRTfit, pRT, totalVol5bar


if __name__ == '__main__':

    for model in models:
        print("Using model: {}".format(model))
        startRun(names=names,RECALC_FITS=RECALC_FITS,CLOSE_FIGS=CLOSE_FIGS,gasName=gasName,model=model,numThreads=numThreads)
