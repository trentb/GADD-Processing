import os
import numpy as np
import pprint
import pickle

#import Tkinter, tkFileDialog
#root = Tkinter.Tk()
#root.withdraw()
#root.overrideredirect(True)
#root.geometry('0x0+0+0')
#root.deiconify()
#root.lift()
#root.focus_force()
#filename = tkFileDialog.askopenfilename(parent=root)
#root.destroy()
currentLineNumber = 0

filename = 'C:/Dropbox/Research/Data Sync/Leakage/8-27-14/X-E1/X-E1_1-81V-4Vsteps.diel'

### path variables. load_path = parent folder, load_name = file name itself
load_path, load_name = os.path.split(filename)

### HERE, declare variables and read files

fileID = open(filename, "r")
lines =  file.read(fileID).splitlines()
dat = {}
dat["name"] = load_name
dat["swp"] = {}
dat["meas"] = {}
names = {}
names["tot"] = list()
names["meas"] = {}
val = {}

### fgetl mimicks MatLab's fgetl() function. Don't need to pass in fileID since we only care about one file at the moment
def fgetl():
        global currentLineNumber
        
        returnLine = lines[currentLineNumber]
        currentLineNumber = currentLineNumber+1
        return returnLine

def str2num(lineOfText):
        stringsInLine = lineOfText.split("\t")
        numbersInLine = list()
        for x in range(len(stringsInLine)):
                if stringsInLine[x] != "":
                        numbersInLine.append(float(stringsInLine[x]))
        return numbersInLine

textt = fgetl()

for i in range(0, 3):
        textt = fgetl()

i=0
j=0

n = list()
n.append(1)
dat["swp"]["none"] = 1
name = list()
name.append("none")

while textt != "*****":
        if textt == "SWEEP":
                j += 1
                textt = fgetl()
                name.append(textt.rsplit("\t")[0])
                textt = fgetl()
                sweep_vals = str2num(textt)
                dat["swp"][name[j]] = sweep_vals
                n.append(len(dat["swp"][name[j]]))
        elif textt == "ENDSWEEP":
                j -=1
        else:
                names["tot"].append(textt)
                i+=1
                tok = textt.rsplit()[0]
                remain = textt.rsplit()[1]
                dat["meas"][remain.strip()] = {}
                dat["meas"][remain.strip()]["type"] = tok
                if j>0:
                        nr = n[1:j+1]
                        #nr.reverse()
                        namer = name[1:j+1]
                        namer.reverse()
                        dat["meas"][remain.strip()]["size"] = nr
                        dat["meas"][remain.strip()]["swp"] = namer
                else:
                        dat["meas"][remain.strip()]["size"] = n[0:j+1]
                        dat["meas"][remain.strip()]["swp"] = name[0:j+1]
                val[remain.strip()] = {}
        textt = fgetl()

### Import Data???
while textt != "**********":
        textt = fgetl()

linesRemaining = lines[currentLineNumber:len(lines)]

for i in range(0, len(linesRemaining)/2):
        textt = fgetl()
        temptextt = fgetl()
        remain = textt.rsplit(" ")[1]
        if i==0:
                val[remain.strip()] = np.array(str2num(temptextt))
        elif i==len(linesRemaining)/4:
                val[remain.strip()] = np.array(str2num(temptextt))
        else:
                val[remain.strip()] = np.vstack([val[remain.strip()], str2num(temptextt)])

names["meas"] = dat["meas"].keys()
num_var = len(names["meas"])

for i in range(0, num_var):
        j = len(val[names["meas"][i]])
        k = dat["meas"][names["meas"][i]]["size"]
        
        # ### Not sure why it's there?
        # # if numel(k) == 1:
        # #     k.append[1]
        
        if dat["meas"][names["meas"][i]]["type"]=='REAL':
            if np.prod(k) == j:                 
                dat["meas"][names["meas"][i]]["real"] = np.reshape(np.transpose(val[names["meas"][i]])[0],k)
            else:
                dat["meas"][names["meas"][i]]["real"] = np.transpose(val[names["meas"][i]])[0]
#        dat["meas"][names["meas"][i]] = {}
        
        elif np.prod(k) == j:
                dat["meas"][names["meas"][i]]["real"] = np.reshape(np.transpose(val[names["meas"][i]])[0],k)
                dat["meas"][names["meas"][i]]["imag"] = np.reshape(np.transpose(val[names["meas"][i]])[1],k)
        else:
                dat["meas"][names["meas"][i]]["real"] = np.transpose(val[names["meas"][i]])[0]
                dat["meas"][names["meas"][i]]["imag"] = np.transpose(val[names["meas"][i]])[1]
pickle.dump(dat, open(os.path.join(load_path,os.path.splitext(load_name)[0]+".pkl"), "wb"))