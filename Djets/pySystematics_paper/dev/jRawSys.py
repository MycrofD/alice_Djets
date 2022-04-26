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
trial4eachvariation = 1

#1
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 4, 1, 0, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1


exit()

################################################################
#2
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 4, 1.15, 0, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1

#3
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 4, 0.85, 0, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1

#4 #fixed sigma fixed mass
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 4,  1, 1, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1

#5 free sigma fixed mass
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 1, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1


#6
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.71, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1

#7
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.72, 2.1, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1

#9
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.71, 2.0, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1
#10
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.71, 2.03, 2
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1
#11
for i in range(trial4eachvariation):
    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.71, 2.1, 4
    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
    sysNum += 1
##8
#for i in range(trial4eachvariation):
#    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 0, 1.74, 2.1, 2
#    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
#    sysNum += 1
##12
#for i in range(trial4eachvariation):
#    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 1, 1.71, 2.1, 2
#    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
#    sysNum += 1
##13
#for i in range(trial4eachvariation):
#    boundSigma, fsigmafactor, fixedMass, bkgType, minfSys, maxfSys, fMassBinWidthFactor = 0, 1, 0, 2, 1.71, 2.1, 2
#    subprocess.call(['bash', 'zrun_mainsys.csh', str(boundSigma), str(fsigmafactor), str(fixedMass), str(bkgType), str(minfSys), str(maxfSys), str(fMassBinWidthFactor), str(sysNum)])
#    sysNum += 1

print(sysNum)

