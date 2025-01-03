#gasName="Methane"
gasName="CarbonDioxide"
#gasName="Hydrogen"
#gasName="Nitrogen"
#gasName="Krypton"
#flags for non-ideal gas corrections
useFugacity= True
useCompress=False
#Assumes that all data is absolute uptake, default is false (if excess file exists will assume excess uptake)
OVERRIDE_ABSOLUTE=False
#Convert data from weight percent to mmol/g, for fitting
CONVERT_WT_MOL=False
#######################PLOTTING PARAMETERS#########################################################################
#parameter to force recalculating of fits, even if fit was already done
RECALC_FITS=True
#Closes figures automatically, otherwise you have to manually close them to continue
CLOSE_FIGS=False
# parameter to cut off absolute isosteric enthalpy calcuation at thetamax
CUTOFF_THETA=False
#parameter to fit isoexcess enthalpy with smoothing spline
USE_SPLINE=True
#max pressure to plot on X axis
MAX_PLOT_PRESS=.1
STEP_PRESS=0.001
DELTAH_MAX=12
DELTAH_MIN=0
ABSOLUTE_MAX=6
EXCESS_MAX=6
###################################################################################################################


#Names of samples that you want to run
names=[]
names.append("EE2")
#names.append("Nuchar")
#names.append("ET094")
#names.append("EH046")
#names.append("HKUST1_BAM")
#names.append("HKUST1_NREL")
#names.append("HKUST1_LBNL")
#names.append("HKUST1_UCB")
#names.append("HKUST1_ICMPE")
#names.append("HKUST1_Hiden")
#names.append("HKUST1_NWU")
#names.append("HKUST1_HU")
#names.append("HKUST1_Max")
#names.append("HKUST1_A2ML")
#names.append("HKUST1_UNT")
#names.append("HKUST1_MSU")

#Names of Langmuir Models that you want to run
models = []
models.append("Hill")
#models.append("sL")
#models.append("dL")
#models.append("tL")
#models.append("coAds")
#models.append("unilan")
#models.append("unilanPurewal")
#models.append("MDA")
#models.append("toth")
#models.append("sips")

#Number of populationSize to run multiplies by core count
numThreads=8
