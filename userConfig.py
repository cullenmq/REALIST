gasName="Methane"
#gasName="CarbonDioxide"
#gasName="Hydrogen"
#gasName="Nitrogen"
#gasName="Krypton"
useFugacity= True

#Assumes that all data is absolute uptake, default is false (if excess file exists will assume excess uptake)
OVERRIDE_ABSOLUTE=False


#######################PLOTTING PARAMETERS#########################################################################
#parameter to force recalculating of fits, even if fit was already done
RECALC_FITS=False
#Closes figures automatically, otherwise you have to manually close them to continue
CLOSE_FIGS=True
#max pressure to plot on X axis
MAX_PLOT_PRESS=15
STEP_PRESS=0.001
DELTAH_MAX=20
DELTAH_MIN=12
ABSOLUTE_MAX=30
EXCESS_MAX=22
###################################################################################################################


#Names of samples that you want to run
names=[]
#names.append("ET095")
#names.append("ET094")
names.append("EH046")

#Names of Langmuir Models that you want to run
models = []
models.append("sL")
models.append("dL")
#models.append("tL")
#models.append("coAds")
#models.append("unilan")
#models.append("unilanPurewal")
#models.append("toth")
#models.append("sips")

#Number of threads that you want to run, this should be a multiple of the number of logic threads your computer has for max efficiency
numThreads=12
