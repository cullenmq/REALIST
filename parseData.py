import csv
import xlsxwriter
import json as json
import os.path

def loadFitData(name):
    file=name+'_fit.json'
    if not os.path.isfile(file):
        print("Fit does not exist! run fitting algorithm")
        return None
    print("loading fit from file")
    coef= json.load(open(file))
    print("RSSR: {}".format(coef["rssr"]))
    return coef

def saveFitData(name,coef):
    json.dump(coef, open(name+"_fit.json", 'w'))


"""Returns data from an excel spreadsheet in the form of T(K),P(MPa), and ads (mmol/g)"""
#data is in mmol/g, MPa, and K
def loadData(fileName='sampleData.csv',isBar=True):
    # 2d array of excess uptake at different temps
    ads = {}
    # 2d array of pressure at different temps
    P = {}
    # 1D array of temperatures
    T = []
    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        if(isBar):
            conv=10
        else:
            conv=1
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                T=list(filter(None, row))
                #make array elements for pressure and composition
                for temp in T:
                    ads[str(temp)]=[]
                    P[str(temp)]=[]

            elif line_count == 1:
                pass
            else:
                for i, temp in enumerate(T):
                    if row[3*i] and row[3*i+1]:

                        #convert bar to MPa
                        P[str(temp)].append(float(row[3*i])/conv)
                        ads[str(temp)].append(float(row[3*i+1]))
            line_count+=1
    return T,P,ads
def grabAdsorption(ads,press):
    newAds= {}
    newPress= {}
    for temp in press:
        #find index where pressure starts decreasing, this is the end of adsorption data
        prevPress=-1
        for i,curPress in enumerate(press[temp]):
            if(curPress<prevPress):
                newAds[temp]=ads[temp][:i]
                newPress[temp]=press[temp][:i]
                break
            else:
                prevPress=curPress
            if (curPress ==press[temp][-1]):
                print("no desorption found!")
                newAds[temp]=ads[temp]
                newPress[temp]=press[temp]
    return newPress,newAds
def convBartoMPa(press):
    newPress={}
    for temp in press:
        newPress[temp]=[x*.1 for x in press[temp]]
    return newPress
def saveFits(fits,name):
    # opening the csv file in 'w' mode
    file = open('name'+'.csv', 'w', newline ='')
def saveFile(fit,data):
    name=fit["Sample Name"]+'.xlsx'
    workbook = xlsxwriter.Workbook(name)
    worksheet={}
    worksheet['fit'] = workbook.add_worksheet('Fit Parameters')
    for i,name in enumerate(fit):
        worksheet['fit'].write(i,0, name)
        worksheet['fit'].write(i,1,fit[name])
    for i, names in enumerate(data):
        curData=data[names]
        for temp in data[names]:
            if temp not in worksheet:
                worksheet[temp]=workbook.add_worksheet(temp)
            worksheet[temp].write(0,i,names)
            worksheet[temp].write_column(1,i,curData[temp])
    workbook.close()
if __name__ == '__main__':
    T,P,ads=loadData()
    newP,newAds=grabAdsorption(ads,P)
