
############################# Choose Gas Type ##########################################
gasName="Methane"
#gasName="CarbonDioxide"
#gasName="Hydrogen"
#gasName="Nitrogen"
#gasName="Krypton"
#flag for non-ideal gas corrections
useFugacity= True
#flag for using ideal gas density instead of real gas density
useIdealGas=False

########################### Choose Fitting Parameters ##########################################
#parameter to force recalculating of fits, even if fit was already done
RECALC_FITS=False
#adsorbed volume definition, only uncomment one of the two
#ADS_VOL="Static"
ADS_VOL="Dynamic"

#parameter to force or fit prefactor dependences. True means that code will fit prefactor 1/T^x
FIT_PREFACT=False
#Denominator Temperature Dependence prefactor for Langmuir (FIT_PREFACT must be false) should be 0-3
preFact=0.5
#Number of differential evolution instances to run (multiplies by PC core count)
numThreads=8

######################### Define format for raw data ##########################################
#Assumes that all data is absolute uptake, default is false (if excess file exists will assume excess uptake)
OVERRIDE_ABSOLUTE=False
#Convert data from weight percent to mmol/g, for fitting (set to false if using data in mmol/g)
CONVERT_WT_MOL=False

#######################PLOTTING PARAMETERS#########################################################################

#Closes figures automatically, otherwise you have to manually close them to continue
CLOSE_FIGS=False
# parameter to cut off absolute isosteric enthalpy calcuation at thetamax
CUTOFF_THETA=True
#parameter to fit isoexcess enthalpy with smoothing spline
#TODO: fix or remove spline, does not work currently
USE_SPLINE=False
#max pressure to plot on X axis
MAX_PLOT_PRESS=10
STEP_PRESS=0.0001
DELTAH_MAX=10
DELTAH_MIN=0
ABSOLUTE_MAX=35
EXCESS_MAX=35
###################################################################################################################

######################### Choose Sample Name ##########################################
#Names of samples that you want to run
names=[]
#names.append("EE2")
names.append("EH046")


######################### Choose Adsorption Model ##########################################
#append models to run more than one model for a given dataset
models = []
#models.append("Hill")
#models.append("sL")
models.append("dL")
#models.append("tL")
#models.append("coAds")
#models.append("unilan")
#models.append("unilanPurewal")
#models.append("MDA")
#models.append("toth")
#models.append("sips")


