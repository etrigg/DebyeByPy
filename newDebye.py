#THIS SCRIPT SIMULATES X-RAY DIFFRACTION
import math
import numpy as np
from datetime import datetime
import os

# parameters defined at bottom of script #

def getCromerMannFactor(q,cromerMannDict,atom1Type,atom2Type,dwfB=None):
    if atom1Type not in cromerMannDict.keys():
        print 'Error: atom type "%s" Cromer Mann parameters not found. Aborted.' % (atom1Type)
        os.sys.exit()
    if atom2Type not in cromerMannDict.keys():
        print 'Error: atom type "%s" Cromer Mann parameters not found. Aborted.' % (atom2Type)
        os.sys.exit()
    cmVars1 = cromerMannDict[atom1Type]
    cmVars2 = cromerMannDict[atom2Type]
    factor1, factor2 = 0,0
    qOver4Pi = q/(4*math.pi)
    for i in [0,2,4,6]: 
        factor1 += cmVars1[i]*math.exp(-cmVars1[i+1]*qOver4Pi*qOver4Pi)
        factor2 += cmVars2[i]*math.exp(-cmVars2[i+1]*qOver4Pi*qOver4Pi)
    factor1 += cmVars1[8]
    factor2 += cmVars2[8]

    if dwfB:
        debyeWallerFactor = math.exp(-qOver4Pi*qOver4Pi*dwfB)
    else:
        debyeWallerFactor = 1.
    return factor1*factor2*debyeWallerFactor*debyeWallerFactor

def getAtomsAndCoordinatesFromFile(filename):
    atomTypes = []
    positions = []
    f = open(filename, 'r')
    for line in f.readlines()[2:]: #skip first two lines
        if line:
            line = line.split()
            atom, coords = line[0], [float(ea) for ea in line[1:]]
            if atom in atomTypes:
                index = atomTypes.index(atom)
                positions[index].append(coords)
            else:
                atomTypes.append(atom)
                positions.append([coords])
    f.close()
    return atomTypes, positions

def getHomoPDF(positions,pdfStepSize,maxAtomDistance):
    counts = np.zeros(math.ceil(maxAtomDistance/pdfStepSize))
    pdfStartTime = datetime.now() # time estimate
    timeEstimatePrinted = False
    for i in range(len(positions)):
        for j in range(i+1,len(positions)):
            d = 0
            for k in range(3):
                diff = positions[i][k]-positions[j][k]
                d += diff*diff
            d = math.sqrt(d)
            if d < 1:
                print 'Warning: <1A separation between atom pair %d-%d' % (i, j)
            counts[math.floor(d/pdfStepSize)] += 1

        # time estimate
        if not timeEstimatePrinted and (datetime.now()-pdfStartTime).total_seconds() >= 10: # time estimate
            timeRemaining = .5 * float(len(positions)) * (datetime.now()-pdfStartTime).total_seconds()/float(i)
            print 'Current PDF time remaining: %s seconds' % str(timeRemaining).split('.')[0]
            timeEstimatePrinted = True

    # determine actual max atom-atom distance and print
    maxIndex = 0
    for i,c in enumerate(counts):
        if c > 0:
            maxIndex = i
    print 'Actual maximum atom-atom distance is %f.' % (float(maxIndex)*pdfStepSize)
    return counts

def getHeteroPDF(positions1,positions2,pdfStepSize,maxAtomDistance):
    counts = np.zeros(math.ceil(maxAtomDistance/pdfStepSize))
    pdfStartTime = datetime.now() # time estimate
    timeEstimatePrinted = False
    for i in range(len(positions1)):
        for j in range(len(positions2)):
            d = 0
            for k in range(3):
                diff = positions1[i][k]-positions2[j][k]
                d += diff*diff
            d = math.sqrt(d)
            counts[math.floor(d/pdfStepSize)] += 1
            
        # time estimate
        if not timeEstimatePrinted and (datetime.now()-pdfStartTime).total_seconds() >= 10: 
            timeRemaining = float(len(positions1)) * (datetime.now()-pdfStartTime).total_seconds()/float(i)
            print 'Current PDF time remaining: %s seconds' % str(timeRemaining).split('.')[0]
            timeEstimatePrinted = True

    # determine actual max atom-atom distance and print
    maxIndex = 0
    for i,c in enumerate(counts):
        if c > 0:
            maxIndex = i
    print 'Actual maximum atom-atom distance is %f.' % (float(maxIndex)*pdfStepSize)
    return counts
              
def getPDFForEachPairType(atomTypes, positions, pdfStepSize, maxAtomDistance):
    r = np.arange(0.0,maxAtomDistance,pdfStepSize)
    pairTypes = []
    pdfs = []
    for i in range(len(atomTypes)):
        for j in range(i,len(atomTypes)):
            print 'Calculating %s-%s PDF...' % (atomTypes[i],atomTypes[j])
            pairTypes.append([atomTypes[i],atomTypes[j]])
            if i==j:
                pdfs.append(getHomoPDF(positions[i],pdfStepSize,maxAtomDistance))
            else:
                pdfs.append(getHeteroPDF(positions[i],positions[j],pdfStepSize,maxAtomDistance))
            print 'Finished %s-%s PDF.' % (atomTypes[i],atomTypes[j])
    return pairTypes, pdfs

def getPairIntensity(q,pdf,pdfStepSize):
    I = 0
    for index in range(len(pdf)):
        if pdf[index]:
            d = index*pdfStepSize
            I += math.sin(q*d)/(q*d) * pdf[index] * 2. # multiply by 2 because we half-counted the atom pairs.
    return I

def calculateDiffraction(inputFile, qMin, qMax, qStep, cromerMannDict, pdfStepSize,
                         maxAtomDistance, dwfB=None, nAtomLimit=None):
    
    atomTypes, positions = getAtomsAndCoordinatesFromFile(inputFile)
    # positions is a 3-d list: positions[atomtype][atomnumber][coordinatenumber]
    print 'Finished reading input file.'
    # now calculate the PDF for each pair type
    pairTypes, pdfEachPairType = getPDFForEachPairType(atomTypes, positions, pdfStepSize, maxAtomDistance)
    print 'Finished calculating PDFs. Calculating intensities...'
    # get diffraction
    qVector = np.arange(qMin, qMax, qStep)
    iVector = []
    startDiffractionTime = datetime.now()
    for q in qVector:
        I = 0
        for pairType,pdf in zip(pairTypes,pdfEachPairType):
            cromerMannFactor = getCromerMannFactor(q,cromerMannDict,pairType[0],pairType[1],dwfB)
            I += cromerMannFactor*getPairIntensity(q,pdf,pdfStepSize)
            # get single-atom scattering if relevant
            if pairType[0] == pairType[1]:
                nAtomsOfThisType = len(positions[atomTypes.index(pairType[0])])
                I += cromerMannFactor*float(nAtomsOfThisType)
        iVector.append(I)
        if len(iVector)%200. == 0:
            prog,tot = float(len(iVector)), float(len(qVector))
            seconds_remaining = (datetime.now()-startDiffractionTime).total_seconds()*((tot-prog)/prog)
            print 'Finished %d/%d intensity calculations (%d seconds to go)' % (len(iVector),len(qVector), int(seconds_remaining))
    print 'Finished %d/%d intensity calculations.' % (len(qVector),len(qVector))
    return qVector, iVector, pairTypes, pdfEachPairType


