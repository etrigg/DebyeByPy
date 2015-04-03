#THIS SCRIPT WRAPS newDebye.py

import debyeByPy as nD
from datetime import datetime


############## User-specified parameters ################################

inputFile = 'inputFile.xyz' # specify an input file (XYZ format)
outputFile = 'output.txt'   # name the output file
cromerMannFile = 'cromer_mann.tab'  # must be tab-delimited
qMin = 0.8          # X-ray scattering vector in inverse Angstroms
qMax = 15.
qStep=0.01
dwfB=.1             #Debye Waller factor "B" parameter (use dwfB=None if not desired)
pdfStepSize=.01     #Defines the resolution of real space atomic positions (Angstrom).
                    # Note: be sure that pdfStepSize < 1/qMax
maxAtomDistance=100 #Must be greater than the maximum separation between two atoms in your file (Angstrom)
                    # Note: overestimating this number by a factor of 2 or 3 won't affect computation time

# Specify oxidation states. If not specified, the zero oxidation
# state will be used. Atom types are case-sensitive, use 'Na' not 'na' 
oxidationStates = {
    'V': 2,   # for example, this specifies all vanadium atoms should be taken in oxidation state +2
    'O': -2   # and all oxygen atoms should be taken in oxidation state -2
    }

# Additional Note:
# -Computation time is mostly a function of number of atoms, and number of
#   different elements. Changing the parameters above won't affect computation
#   time by a whole lot.


#########################################################################























######### get Cromer Mann values from file ##############################
cromerMannDict = {}
usedOxidationStates = {}
cromerMannFile = 'cromer_mann.tab'
f = open(cromerMannFile, 'r')
for line in f.readlines():
    if line:
        line = line.split('\t')
        atomType = line[0]
        oxidationState = line[1]
        if (atomType not in cromerMannDict.keys()) and (int(oxidationState) == 0):
            cromerMannDict[atomType] = [float(ea) for ea in line[2:]]
            usedOxidationStates[atomType] = oxidationState
        if (atomType in oxidationStates.keys() and int(oxidationState) == oxidationStates[atomType]):
            cromerMannDict[atomType] = [float(ea) for ea in line[2:]]
            usedOxidationStates[atomType] = oxidationState
for atom in oxidationStates.keys():
    if int(usedOxidationStates[atom]) != oxidationStates[atom]:
        print 'WARNING: Atom "%s" oxidation state "%d" not found in Cromer Mann file.' % (atom, oxidationStates[atom])
        print ' -defaulting to zero oxidation state.'
#########################################################################
            

########### get start time, get diffraction #############################
startTime = datetime.now()
# get diffraction
q,I,pairTypes,pdfEachPairType = nD.calculateDiffraction(inputFile, qMin, qMax, qStep,
                                                        cromerMannDict, pdfStepSize, maxAtomDistance,
                                                        dwfB, nAtomLimit=None)
#########################################################################


#### report oxidation states used #######################################
atomTypes = []
for p in pairTypes:
    for ea in p:
        if ea not in atomTypes:
            atomTypes.append(ea)
for atomType in atomTypes:
    print 'Atom "%s" used oxidation state "%s"' % (atomType, usedOxidationStates[atomType])
#########################################################################


#### write intensity vs. q to file ######################################
f = open(outputFile,'w')
s = '\n'.join([ '\t'.join([str(qea),str(Iea)]) for qea,Iea in zip(q,I) ])
f.write(s)
f.close()
print 'Wrote file.'
print 'Time to run script (hh:mm:ss) ', str(datetime.now()-startTime).split('.')[0]
#########################################################################


# write pair distribution functions to files (uncomment if desired) #####
'''
dVector = np.arange(0,maxAtomDistance,pdfStepSize)
for i in range(len(pairTypes)):
    pdf_filename = '%s-%s pdf.txt' % (pairTypes[i][0],pairTypes[i][1])
    f = open(pdf_filename,'w')
    s = '\n'.join([ '\t'.join([str(d),str(pd)]) for d,pd in zip(dVector,pdfEachPairType[i])])
    f.write(s)
    f.close()
'''
