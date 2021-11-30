import numpy as np
from parseData import saveFitData,loadFitData
from EOS import calcDensity,calcFugacity,initFugacity
from multiprocessing import Process,Queue
from scipy.optimize import differential_evolution

coefNames={} #list of names for coefficients for specific model
bounds={} #list of bounds for each coef (must equal coefNames)
thetas={} #theta function to call for specific model
dPdTs={} #dp/dt function to call for specific model (used for enthalpy)


### Calculating K, denominator either no T, T, or square root T dependence!
def calcK(A,E,T,x):
    r = 8.314462 / 1000
    denom=T**x
    return A/denom*np.exp(E/r/T)
def addModel(name,coefs,bound,theta,dPdT):
    coefNames[name]=coefs
    bounds[name]=bound
    thetas[name]=theta
    dPdTs[name]=dPdT
## DUAL Langmuir
dlNames=["nmax","a","vmax","A1","E1","A2","E2"]

dlBounds=[(0.0,0.1),(0.29,1),(0.0, 1e-5), (0.0,10),(0.0,30),(0.0,10),(0.0, 30)]
def theta_dL(p,t,coef):
        a = coef["a"]
        A1 = coef["A1"]
        E1 = coef["E1"]
        A2 = coef["A2"]
        E2 = coef["E2"]
        #x=coef["x"]
        K1=calcK(A1,E1,t,1)#K1=a4*np.exp(a5/r/t)/denom
        #assume both sites are equal
        #K2=K1
        K2=calcK(A2,E2,t,1)#a6*np.exp(a7/r/t)/denom
        return ((1-a)*(K1*p/(1+K1*p))+a*(K2*p/(1+K2*p)))
def dPdT_dL(P,T,coef):
    r = 8.314462 / 1000
    a = coef["a"]
    A1 = coef["A1"]
    E1 = coef["E1"]
    A2 = coef["A2"]
    E2 = coef["E2"]
    #x=coef["x"]
    K1=calcK(A1,E1,T,1)
    K2=calcK(A2,E2,T,1)
    #plt.figure()
    #-dP/dT=(dTheta/dP)^-1*dtheta/dK*dK/dT
    #X=dTheta/dP
    X= (1-a)*(K1/((1+K1*P)**2))+a*(K2/((1+(K2)*P)**2))
    #Y=dTheta/dK
    Y1=(1-a)*(P/((1+K1*P)**2))
    Y2=a*(P/((1+(K2)*P)**2))
    #Z=dk/dT
    #note:square root has an extra 0.5 in front
    #if(isSqRoot):
     #   mult=0.5
    #else:
     #   mult=1
    mult=1
    Z1=-K1*((mult*r*T+E1)/(r*T**2))
    Z2=-K2*((mult*r*T+E2)/(r*T**2))
    return (Y1*Z1+Y2*Z2)/X
addModel('dL',dlNames,dlBounds,theta_dL,dPdT_dL)

dlDisNames=["nmax","a","vmax","A1","E1","A2","E2"]

dlDisBounds=[(0.0,0.1),(0,1),(5e-7, 1e-5), (0.0,10),(0.0,30),(0.0,10),(0.0, 30),(0,3)]
def theta_dLDis(p,t,coef):
    a = coef["a"]
    A1 = coef["A1"]
    E1 = coef["E1"]
    A2 = coef["A2"]
    E2 = coef["E2"]
    x=1
    K1=calcK(A1,E1,t,x)#K1=a4*np.exp(a5/r/t)/denom
    #assume both sites are equal
    #K2=K1
    K2=calcK(A2,2*E2,t,x)#a6*np.exp(a7/r/t)/denom
    return ((1-a)*(K1*p/(1+K1*p))+a*(np.sqrt(K2*p)/(1+np.sqrt(K2*p))))
def dPdT_dLDis(P,T,coef):
    r = 8.314462 / 1000
    a = coef["a"]
    A1 = coef["A1"]
    E1 = coef["E1"]
    A2 = coef["A2"]
    E2 = coef["E2"]
    x=1
    K1=calcK(A1,E1,T,x)
    K2=calcK(A2,2*E2,T,x)
    #plt.figure()
    #-dP/dT=(dTheta/dP)^-1*dtheta/dK*dK/dT
    #X=dTheta/dP
    num= a*P(2*K1*T**2)*(E1+r*T)*(a*K2)
    return (Y1*Z1+Y2*Z2)/X
addModel('dLDis',dlDisNames,dlDisBounds,theta_dLDis,dPdT_dLDis)

#UnilanPurewal Model
unilanPurewalNames=["nmax","deltaS","vmax","Emin","Emax"]
unilanPurewalBounds=[(0.0,0.1),(0,10),(1e-7, 1e-5), (0.0,30),(0.0,30)]
def theta_unilanPurewal(p,t,a):
        r = 8.314462 / 1000
        deltaS = a["deltaS"]
        Emin = a["Emin"]
        Emax = a["Emax"]
        c=np.exp(deltaS)
        return (r*t)/(Emax-Emin)*np.log((c+p*np.exp(Emax/r/t))/(c+p*np.exp(Emin/r/t)))
def dPdT_unilanPurewal(p,t,a):
    r = 8.314462 / 1000
    deltaS = a["deltaS"]
    Emin = a["Emin"]
    Emax = a["Emax"]
    EminExp=p*np.exp(Emin/r/t)
    EmaxExp=p*np.exp(Emax/r/t)
    c = np.exp(deltaS)
    num= p/t*(Emin-c*Emin/(c+EminExp)-Emax+Emax*c/(c+EmaxExp)) +r*p*np.log((c + EmaxExp) / (c + EminExp))
    denom=c*r*t*(1/(c+EminExp)-1/(c+EmaxExp))
    return num/denom
addModel('unilanPurewal',unilanPurewalNames,unilanPurewalBounds,theta_unilanPurewal,dPdT_unilanPurewal)
#Unilan Model
unilanNames=["nmax","deltaS","vmax","Emin","Emax"]
unilanBounds=[(0.0,0.1),(1e-4,1),(1e-7, 1e-5), (0,30),(0.1,30)]
def theta_unilan(p,t,a):
        r = 8.314462 / 1000
        deltaS = a["deltaS"]
        Emin = a["Emin"]
        Emax = a["Emax"]
        c=deltaS
        try:
            x=(r*t)/(Emax-Emin)*np.log((r*t/c+p*np.exp(Emax/r/t))/(r*t/c+p*np.exp(Emin/r/t)))
        except:
            pass
        return x
def dPdT_unilan(p,T,a):
    r = 8.314462 / 1000
    deltaS = a["deltaS"]
    Emin = a["Emin"]
    Emax = a["Emax"]
    c = deltaS
    EminExp=c*p*np.exp(Emin/r/T)
    EmaxExp=c*p*np.exp(Emax/r/T)
    num=(r*(r*T+Emax))/(EmaxExp+r*T)-(r*(r*T+Emin))/(EminExp+r*T)+(r*T*np.log((EmaxExp+r*T)/(EminExp+r*T))-Emax+Emin)/T
    denom=1/(EminExp+r*T)-1/(EmaxExp+r*T)
    return p/(r*T)**2*num/denom
addModel('unilan',unilanNames,unilanBounds,theta_unilan,dPdT_unilan)
#cooperative adsorption model
coAdsNames=["nmax","Eint","vmax","A1","E1","A2","E2"]
coAdsBounds=[(0.0,0.1),(-7,7),(1e-7, 1e-5), (0.001,100),(0.0,30),(0.001,100),(0.0, 30)]
def theta_coAds(p,t,a):
    r = 8.314462 / 1000
    Eint = a["Eint"]
    A1 = a["A1"]
    E1 = a["E1"]
    A2 = a["A2"]
    E2 = a["E2"]
    K1=calcK(A1,E1,t)#a4*np.exp(a5/r/t)/denom
    K2=calcK(A2,E2,t)#a6*np.exp(a7/r/t)/denom
    #assume both sites are equal
    #K2=K1
    intTerm=np.square(p)*K1*K2*np.exp(-Eint/r/t)
    return 0.5*(((K1+K2)*p+2*intTerm)/(1+(K1+K2)*p+intTerm))
def dPdT_coAds(P,T,coef):
    r = 8.314462 / 1000
    Ei = coef['Eint']
    A1 = coef['A1']
    E1 = coef['E1']
    A2 = coef['A2']
    E2 = coef['E2']
    B=1/(r*T)
    K1=calcK(A1,E1,T)
    K2=calcK(A2,E2,T)
    Ki=np.exp(Ei*B)
    return (P*(K2*Ki*T**2*(-E2-1/B)+ K1**2*K2*np.square(P)*(-E2+Ei-1/B)+K1*(Ki*T**2*(-E1-1/B)+K2*P*(T*(-2*E1-2*E2+2*Ei-4/B)+K2*P*(-E1+Ei-1/B)))))\
    /(T/B*(K1*K2*np.square(P)*(K1+K2)+4*K1*K2*P*T+T**2*Ki*(K1+K2)))
addModel('coAds',coAdsNames,coAdsBounds,theta_coAds,dPdT_coAds)

#toth
tothNames=["nmax","t","vmax","A1","E1"]
tothBounds= [(0.0,0.1),(.1,3),(0, 1e-2), (1e-3,1e3),(0.0,20)]
def theta_toth(p,T,a):
    t = a["t"]
    A1 = a["A1"]
    E1 = a["E1"]
    K=calcK(A1,E1,T)

    return (K*p)/(1+(K*p)**t)**(1/t)
def dPdT_toth(P,T,a):
    r = 8.314462 / 1000
    E1 = a["E1"]
    return -P*(E1+r*T)/(r*T**2)
addModel('toth',tothNames,tothBounds,theta_toth,dPdT_toth)

##Sips
sipsNames=["nmax","c","vmax","A1","E1"]
sipsBounds= [(0.0,0.1),(0.01,3),(0, 1e-2), (1e-3,1e3),(0.0,30)]
def theta_sips(p,t,a):
    c = a["c"]
    A1 = a["A1"]
    E1 = a["E1"]
    K=calcK(A1,E1,t)
    return (K*p)**c/(1+(K*p)**c)
def dPdT_sips(P,T,coef):
    r = 8.314462 / 1000
    E1 = a["E1"]
    return -P*(E1+r*T)/(r*T**2)
addModel('sips',sipsNames,sipsBounds,theta_sips,dPdT_sips)




class isotherm():
    def __init__(self,name,numThreads=48,isAbsolute=False):
        self.name=name
        self.numThreads=numThreads
        self.bounds=bounds[name]
        self.names=coefNames[name]
        self.isAbsolute=isAbsolute
    def diff_evol(self,bounds,p,ads,den,i,queue):
            strat='randtobest1exp'#'randtobest1bin'
            queue.put(differential_evolution(func=self.objective,bounds=bounds,args=(p,ads,den),seed=i,disp=False,tol=.006,maxiter=10000,
                                                strategy=strat,init='random'))
    def calcExcess(self,calcPress,ads,den,t,a):
        p=np.array(calcPress)
        den=np.array(den)
        vmax = a["vmax"]
        #dynamic adsorbed phase volume assumption
        return ads-1000*(vmax*den)*self.theta(p,t,a)
    def genExcess(self,p,t,a,den,isAbsolute=False):
            #density in mol/m^3
            #nmax in mols
            #vmax in m^3/g
            #convert lists to arrays
            p=np.array(p)
            den=np.array(den)
            nmax = a["nmax"]
            vmax = a["vmax"]
            if (isAbsolute):
                #absolute adsorption
                return 1000*(nmax)*self.theta(p,t,a)
            #dynamic adsorbed phase volume assumption
            return 1000*(nmax-vmax*den)*self.theta(p,t,a)
            #stagnant adsorbed phase volume assumption
            #return 1000*nmax*self.theta(p,t,a)-vmax*den
    def objective(self,params, p, ads,den):
            """Calculate total residual for fits of General Langmuir Isotherms to several data sets."""
            resid =np.array([])
            #convert array into dict with coefficient names
            params=dict(zip(self.names,params))
            # make residual per data set
            for temp in ads:
                modelData = self.genExcess(p[temp],float(temp),params,den[temp],self.isAbsolute)
                curAds=np.array(ads[temp])
                curResid=np.square(curAds - modelData)
                if resid.size==0:
                    resid=curResid
                else:
                    resid=np.append(resid,curResid)
            # now flatten this to a 1D array, as minimize() needs
            return np.sum(resid)#resid.flatten()

    def runFit(self,p,ads,gasName,useFugacity=True, isRunFit=True,name='test'):

        print("calculating densities")
        den={}
        if(useFugacity):
            adjPress={}
            st=initFugacity(gasName)
        for temp in p:
            T=float(temp)
            tempDen=[]
            tempPress=[]
            for press in p[temp]:
                tempDen.append(calcDensity(press, T, gasName))
                if(useFugacity):
                    tempPress.append(press*calcFugacity(st,press,T))
            den[temp]=tempDen
            if useFugacity:
                adjPress[temp]=tempPress
        coef=None
        if useFugacity:
            calcPress=adjPress
            print("Using fugacity")
        else:
            print("Using pressure")
            calcPress=p
        if isRunFit==False:
            coef=loadFitData(name+self.name)
            if coef !=None:
                rssr=coef['rssr']
        if isRunFit or coef==None:
            print("running fit")


            #coef,rssr=leastsq(calcPress,ads,den,isSqRoot)
            coef,rssr=self.diffEV(calcPress,ads,den)
            print("done running fit")
            #test=optimize.minimize(self.objective, coef, args=(p,ads,den), method="Nelder-Mead", tol=1e-30)
            coef=dict(zip(self.names,coef))

            coef['rssr']=rssr
            saveFitData(name+self.name,coef)
        #we need to calculate excess uptake if given absolute uptake
        if self.isAbsolute:
            for temp in ads:
                ads[temp] = self.calcExcess(calcPress[temp],ads[temp],den[temp],float(temp),coef)
        return coef, den,calcPress,ads


    def theta(self,p,t,a):
        thetaModel= thetas.get(self.name)
        if thetaModel is not None:
            return thetaModel(p,t,a)
        raise ValueError("Need to define theta function for this model!")
    def dPdT(self,p,t,a):
        dPdTModel= dPdTs.get(self.name)
        if dPdTModel is not None:
            return dPdTModel(p,t,a)
        raise ValueError("Need to define theta function for this model!")


    def diffEV(self,p,ads,den):
        minimized_function=1e9
        numRuns=self.numThreads
        #lets do this in parallel
        jobs=[]
        queue = Queue()
        bounds=self.bounds
        for i in range(numRuns):
            proc=(Process(target=self.diff_evol,args=(bounds,p,ads,den,i,queue)))
            jobs.append(proc)
            proc.start()
        for proc in jobs:
            proc.join()
        for i in range(numRuns):
            temp=queue.get()
            print("{}= {}".format(i,temp.fun))
            if minimized_function > temp.fun:
                minimized_function = temp.fun
                fin = temp.x
                ssr=temp.fun
            #resultt = result
        print("coefs= {}".format(fin))
        rssr=np.sqrt(ssr)
        print("RSSR :{}".format(rssr))
        #for i,key in enumerate(a):
            #a[key]=fin[i]
        return fin,rssr
