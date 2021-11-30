import os
import json
if __name__ == '__main__':
    sampleName = input("Sample Name: ")
    data={}
    data["bulkDens"] = float(input("Bulk Density (g/cc): "))
    data["Vpore"] = float(input("Nitrogen Pore Volume (cc/g): "))
    data["SSA"] = float(input("Specific Surface Area (m2/g): "))
    data["skelDens"] = float(input("Skeletal Density (g/cc): "))
    home=os.getcwd()
    os.chdir(home+'/rawData')
    with open(sampleName+'.json','w') as outfile:
        json.dump(data,outfile)
    os.chdir(home)
