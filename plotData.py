import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from userConfig import *
#from analyzeData import calcDensity,calcVolUptake
matplotlib.use("TkAgg")
def getMinMax(press):
    min=1000
    max=-1
    for temp in press:
        T=float(temp)
        if T>max:
            max=T
        if T<min:
            min=T
    return min,max
def plotIsoHeat(actualPress,fa1,theta,Qst,name):
    plt.figure()
    TMin,TMax=getMinMax(actualPress)
    for temp in actualPress:
        T = float(temp)
        actualP=actualPress[temp]
        addPlot(actualP, Qst[temp], T, '-',minTemp=TMin,maxTemp=TMax)
    plt.ylabel("\u0394Hads (kJ/mol)", fontsize=12)
    plt.xlabel("Pressure (MPa)", fontsize=12)
    plt.ylim(DELTAH_MIN,DELTAH_MAX)
    plt.title(name)
    plt.draw()
    plt.pause(0.001)
    plt.savefig(name + '_deltaHvsPress.png')
    plt.figure()
    for temp in actualPress:
        T = float(temp)
        addPlot(fa1[temp], Qst[temp], T, '-',minTemp=TMin,maxTemp=TMax)
    plt.ylabel("\u0394Hads (kJ/mol)", fontsize=12)
    plt.xlabel("Absolute Uptake (mmol/g)", fontsize=12)
    plt.ylim(DELTAH_MIN, DELTAH_MAX)
    plt.title(name)
    plt.draw()
    plt.pause(0.001)
    plt.savefig(name + '_deltaHvsAbsUptake.png')

    plt.figure()
    for temp in actualPress:
        T = float(temp)
        addPlot(theta[temp], Qst[temp], T, '-',minTemp=TMin,maxTemp=TMax)
    plt.ylabel("\u0394Hads (kJ/mol)", fontsize=12)
    plt.xlabel("Site Occupancy Fraction", fontsize=12)
    plt.ylim(DELTAH_MIN, DELTAH_MAX)
    plt.title(name)
    plt.draw()
    plt.pause(0.001)
    plt.savefig(name + '_deltaHvsTheta.png')
def plotAbsUptake(actualPress,fa1,name,newFigure):
    Tmin,Tmax=getMinMax(actualPress)
    if (newFigure):
        plt.figure()
    for temp in actualPress:
        T=float(temp)
        actualP=actualPress[temp]
        addPlot(actualP, fa1[temp],T, '-.',minTemp=Tmin,maxTemp=Tmax)
    plt.ylabel("Ads(mmol/g)")
    plt.xlabel("Pressure (MPa)")
    plt.ylim(0, ABSOLUTE_MAX)
    plt.title(name)
    plt.draw()
    plt.pause(0.01)
    plt.savefig(name + '_AbsUptakevsPress.png')
def plotExcessUptake(ads,p,y_fit,tempPress,name):
    plt.figure()
    Tmin,Tmax=getMinMax(ads)
    for temp in ads:
        T=float(temp)
        addPlot(p[temp], ads[temp],T, 'o',minTemp=Tmin,maxTemp=Tmax)
        addPlot(tempPress[temp], y_fit[temp],T, '-',labels= str(("%.1f" % T))+'K',minTemp=Tmin,maxTemp=Tmax)
    plt.ylabel("Excess Uptake (mmol/g)",fontsize=12)
    plt.ylim(0, EXCESS_MAX)
    plt.xlim(0, MAX_PLOT_PRESS)
    plt.xlabel("Pressure (MPa)",fontsize=12)
    plt.title(name)
    plt.legend()
    plt.draw()
    plt.pause(0.001)
    plt.savefig(name+'_ExcessvsPressFit.png')
def plotAdsLayer(actualPress,fv1,name):
    plt.figure()
    for temp in actualPress:
        P = actualPress[temp]
        T = float(temp)
        addPlot(P,fv1,T,'-')
    plt.ylabel("Adsorbed Volume (mL/g)")
    plt.xlabel("Pressure (MPa)")
    plt.title(name)
    plt.draw()
    plt.pause(0.01)
    plt.savefig(name+'_AdsorbedVolvsPress.png')
def plotVolTotal(totVol,totVolFit,absUptake,lineType,xlbl,ylbl,pltname,gasName,T=273):

    fig=plt.figure()
    #returns number of dict elements in deltaH i.e temperatures or samples
    fact=len(totVol)
    i=0
    isFirst=True
    sampNum=0
    press5bar=-1
    fitPress=np.arange(0.0001, 10.0, 0.01)
    for name in totVol:
        tVol=totVol[name]
        tVolFit=totVolFit[name]
        abs=absUptake[name]

        line=lineType[name]
        # For the 3D plots of Enthalpy, Uptake
        if(isFirst):
            for p in range(len(fitPress)):
                if(round(fitPress[p],2)==0.50):
                    print(fitPress[p])
                    press5bar=p
                    if press5bar ==-1:
                        print("5 bar not found!")
            isFirst=False
            press=np.arange(0.5, MAX_PLOT_PRESS, 0.01)
            gasRho=calcDensity(press, T, gasName)
            excess,total=calcVolUptake(0,gasRho,0,1)
            #subtract away 5 bar
            total=total-total[0]
            plt.plot(press, total, '--', color='k',label='gas')
        if name=="EH046":
            color='k'
        else:
            color=(i/fact,0,1-i/fact)
        i+=1
        tVol=tVol-tVolFit[press5bar]
        tVolFit=tVolFit-tVolFit[press5bar]
        plt.plot(abs, tVol, 'o', color=color,label=name)
        #plt.plot(fitPress, tVolFit, '-', color=color,label=None)

        plt.ylabel(ylbl,fontsize=12)
        plt.xlabel(xlbl,fontsize=12)
        loc='lower right'
        ncol=fact+1
        x=1
        y=1.03
    #plt.ylim(bottom=0)
    plt.legend(bbox_to_anchor=(x, y), loc=loc,ncol=ncol)

    plt.title(pltname)
    plt.draw()

    plt.pause(0.01)
    plt.savefig(pltname)
    plt.close('all')

def plotTotal(deltaH,absUptake,lineType,xlbl,ylbl,pltname,gasName,T=273,n=1,maxTemp=400,minTemp=230):

    fig=plt.figure()
    #returns number of dict elements in deltaH i.e temperatures or samples
    fact=len(deltaH)
    i=0
    isFirst=True
    sampNum=0
    volUptake5bar=0

    #middle temp for color map
    midTemp=(maxTemp+minTemp)/2

    for name in deltaH:
        Qst=deltaH[name]
        abs=absUptake[name]
        line=lineType[name]
        # For the 3D plots of Enthalpy, Uptake
        if type(Qst) is dict:
            if(isFirst):
                ax = fig.add_subplot(111, projection='3d')
                isFirst=False
            for temp in Qst:
                T=float(temp)
                if(T<=midTemp):
                    fact=(midTemp-minTemp)
                    i=T-minTemp
                    color=(i/fact,0,1-i/fact)
                else:
                    fact=(maxTemp-midTemp)
                    i=T-midTemp
                    color=(1,i/fact,0)
                    #return plt.plot(x, y, sym, color=,label=labels)
                #grab sample number
                sampNum = int(''.join(filter(str.isdigit, name)))
                label=None
                if sampNum==46:
                    sampNum=0
                    label=str(("%.1f" % T))+'K'
                y=abs[temp]
                z=Qst[temp]
                ax.scatter(np.asarray(sampNum),y[::n], z[::n], line,alpha=0.6, color=color,label=label)
                ax.set_zlabel(ylbl,fontsize=12)
                ax.set_ylabel(xlbl,fontsize=12)
                ax.set_xlabel("sample num",fontsize=12)
                loc='upper left'
                ncol=1
                x=1.03
                y=1

                    #addPlot(abs[temp],Qst[temp], T, line,labels=name+' '+str(("%.1f" % T))+'K',zs=T)
        else:
            if(isFirst):
                for p in range(len(abs)):
                    if(round(abs[p],1)==0.5):
                        print(abs[p])
                        volUptake5bar=Qst[p]
                if volUptake5bar ==0:
                    print("5 bar not found!")
                isFirst=False
                press=np.arange(0.5, 10.0, 0.01)
                gasRho=calcDensity(press, T, gasName)
                excess,total=calcVolUptake(0,gasRho,0,1)
                #subtract away 5 bar
                total=total-total[0]
                plt.plot(press, total, '-', color='k',label='gas')
            if name=="EH046":
                color='k'
            else:
                color=(i/fact,0,1-i/fact)
            i+=1
            Qst=Qst-volUptake5bar
            plt.plot(abs, Qst, 'o', color=color,label=name)
            plt.ylabel(ylbl,fontsize=12)
            plt.xlabel(xlbl,fontsize=12)
            loc='lower right'
            ncol=fact+1
            x=1
            y=1.03

    plt.legend(bbox_to_anchor=(x, y), loc=loc,ncol=ncol)

    plt.title(pltname)
    plt.draw()

    plt.pause(0.01)
    plt.savefig(pltname)

def addPlot(x,y,T,sym,minTemp=230,maxTemp=400,labels=None):
    #add in some "cushion" so that the last isotherm isnt completely yellow
    maxTemp=maxTemp*1.05
    midTemp=(minTemp+maxTemp)/2
    if(T<=midTemp):
            fact=(midTemp-minTemp)
            i=T-minTemp
            return plt.plot(x, y, sym, color=(i/fact,0,1-i/fact),label=labels)
    else:
            fact=(maxTemp-midTemp)
            i=T-midTemp
            return plt.plot(x, y, sym, color=(1,i/fact,0),label=labels)
