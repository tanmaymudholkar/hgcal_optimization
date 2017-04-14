#!/usr/bin/env python

from __future__ import print_function, division

import optparse, numpy, sys

parser = optparse.OptionParser()
parser.add_option("-f", "--inputGeometryDir", dest="inputGeometryDir", type='string', default="", help="Name of input geometry directory")

(options, args) = parser.parse_args()

if options.inputGeometryDir == "":
    print ("Please enter input geometry directory!")
    sys.exit(1)

# def isCloseToExistingBoundary(existingBoundaries, tolerance, testValue):
#     for existingBoundary in existingBoundaries:
#         if (numpy.fabs(testValue - existingBoundary) < tolerance):
#             return True
#     return False

def minimumDifference(existingBoundaries, testValue):
    # minDiff = numpy.fabs(existingBoundaries[0]-testValue)
    minDiff = -1
    if (len(existingBoundaries) > 0):
        minDiff = numpy.fabs(testValue - existingBoundaries[0])
    for existingBoundary in existingBoundaries:
        if (numpy.fabs(testValue - existingBoundary) < minDiff):
            minDiff = numpy.fabs(testValue - existingBoundary)
    return minDiff

def testAndAdd(existingBoundaries, toleranceZero, tolerance, testValue):
    # print ("testAndAdd called for existingBoundaries = {sortedExistingBoundaries}, testValue = {testValue}".format(sortedExistingBoundaries=sorted(existingBoundaries), **locals()))
    # if (isCloseToExistingBoundary(existingBoundaries, tolerance, testValue)):
    #     if (minimumDifference(existingBoundaries, testValue) > 0.01):
    #         print ("minR={minR} is close but not equal to existing value!".format(**locals))
    #         print ("Currently, boundaryValues = {existingValues}".format(**locals))
    # else:
    #     boundaryValues += [testValue]
    minimumDiff = minimumDifference(existingBoundaries, testValue)
    if (minimumDiff < tolerance and minimumDiff > toleranceZero):
        print ("testValue = {testValue} is close but not equal to an existing value!".format(**locals()))
        sortedExistingBoundaries = sorted(existingBoundaries)
        print ("Currently, sortedExistingBoundaries = {sortedExistingBoundaries}".format(**locals()))
        sys.exit(0)
    elif(minimumDiff >= tolerance or minimumDiff == -1):
        existingBoundaries += [testValue]

boundaryValues = []
for layer in range(36, 52):
    print ("Fetching information for layer {layer}".format(**locals()))
    prefix = ''
    offsetLayerNumber = 0
    if (layer <= 39):
        prefix = 'FH'
        offsetLayerNumber = layer-27
    else:
        prefix = 'BH'
        offsetLayerNumber = layer-39
    minRArray, maxRArray = numpy.genfromtxt("{options.inputGeometryDir}/geometry_{prefix}{offsetLayerNumber}.txt".format(**locals()), usecols=(10,9), unpack=True)
    for radialLayerIndex in range(len(minRArray)):
        minR = minRArray[radialLayerIndex]
        maxR = maxRArray[radialLayerIndex]
        # print ("Min r is {minR}, max R is {maxR}".format(**locals()))
        
        # if (isCloseToExistingBoundary(boundaryValues, 10, maxR)):
        #     if (minimumDifference(boundaryValues, max) > 0.01):
        #         print ("maxR={maxR} is close but not equal to existing value!".format(**locals))
        #         print ("Currently, boundaryValues = {boundaryValues}".format(**locals))
        testAndAdd(boundaryValues, 10, 0.01, minR)
        testAndAdd(boundaryValues, 10, 0.01, maxR)

sortedBoundaryValues = sorted(boundaryValues)
print ("Final boundary values: {sortedBoundaryValues}".format(**locals()))
