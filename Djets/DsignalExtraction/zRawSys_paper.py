import os
import subprocess
import shlex


boundSigma=0
fsigmafactor=1
fixedMass=0
bkgType=0
minfSys=1.71
maxfSys=2.1
fMassBinWidthFactor=2
sysNum=1


#boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor, sysNum = 0, 1, 0, 0, 1.71, 2.1, 2, 1
trial4eachvariation = 20

var_settings = [[4,1   ,0,0,1.71,2.1 ,2],
                [4,1.15,0,0,1.71,2.1 ,2],
                [4,0.85,0,0,1.71,2.1 ,2],
                [0,1   ,1,0,1.71,2.1 ,2],
                [0,1   ,0,0,1.71,2.1 ,2],
                [0,1   ,0,0,1.72,2.1 ,2],
                [0,1   ,0,0,1.71,2.0 ,2],
                [0,1   ,0,0,1.71,2.03,2],
                [0,1   ,0,0,1.71,2.1 ,4] 
               ]

#-------------
for i in range(trial4eachvariation):
    for j in range(len(var_settings)):
        [boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor] = var_settings[j]
        subprocess.call(['bash', 'zrun_main.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
        sysNum += 1
#-------------

print(sysNum)

