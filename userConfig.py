gasName="Methane"
#gasName="CarbonDioxide"
#gasName="Hydrogen"
#gasName="Nitrogen"
#gasName="Krypton"
useFugacity= True
#parameter to force recalculating of fits, even if fit was already done
RECALC_FITS=False
#Closes figures automatically, otherwise you have to manually close them to continue
CLOSE_FIGS=False


#Names of samples that you want to run
names=[]
#names.append("ET095")
#names.append("ET094")
names.append("EH046")

#Names of Langmuir Models that you want to run
models = []
models.append("dL")
#models.append("coAds")
models.append("unilan")
#models.append("unilanPurewal")
#models.append("toth")
#models.append("sips")

#Number of threads that you want to run, this should be a multiple of the number of logic threads your computer has for max efficiency
numThreads=24
