#!/usr/bin/env python

from numpy import histogram, loadtxt, inf
from math import log
from sys import argv

############################################################################################################################

if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser()
    parser.add_option('-f', '--datafile', help='data file', dest='inputfile', action='store', type='string')
    parser.add_option('-i', '--column', help='column', dest='icolumn', action='store', type='int')
    parser.add_option('-w', '--weightfile', help='weight file', dest='weightfile', action='store', type='string')
    parser.add_option('-j', '--weightcolumn', help='weight column', dest='wcolumn', action='store', type='int')
    parser.add_option('-n', '--minbin', help='the lefe edge of the first bin', dest='minbin', action='store', type='float', default=None)
    parser.add_option('-x', '--maxbin', help='the right edge of the last bin', dest='maxbin', action='store', type='float', default=None)
    parser.add_option('-b', '--nbins', help='the total number of bins', dest='nbins', action='store', type='int')
    parser.add_option('-e', '--freeE', help='calcualte the free energy', dest='freeE_flag', action='store_true', default=False)
    parser.add_option('-o', '--unit', help='unit, 0:kBT, 1:kcal/mol, 2:kj/mol, default=0', dest='ounit', action='store', type='int', default=0)
    parser.add_option('-t', '--temperature', help='temperature, default=300K', dest='temp', action='store', type='float', default=300.0)
    parser.add_option('-a', '--alldata', help='coorect results including data outside of [min:max]', dest='all_flag', action='store_true', default=False)
    parser.add_option('-r', '--nonormalization', help='do not normalize the histogram', dest='normal_flag', action='store_false', default=True)

    (opts, args) = parser.parse_args()

    # commandline
    command = " ".join(argv)
    command = "#" + command
    print(command)

    # read data
    datapoints = loadtxt(opts.inputfile, comments=['#', '@'], usecols=[(opts.icolumn-1)], unpack=True)

    # read weight
    if ((opts.weightfile is not None) and (opts.wcolumn is not None)):
        myweights = loadtxt(opts.weightfile, comments=['#', '@'], usecols=[(opts.wcolumn-1)], unpack=True)        
    else:
        myweights = None

    # min and max
    if ((opts.minbin is not None) or (opts.maxbin is not None)):
        datarange = [opts.minbin, opts.maxbin]
    else:
        datarange = None
        
    # calculate the histogram
    hist_array, bin_edges_array = histogram(datapoints, weights=myweights, bins=opts.nbins, range=datarange, density=opts.normal_flag)
    hist = hist_array.tolist()
    bin_edges = bin_edges_array.tolist()

    # x coordinate
    bin_mid = []
    for dummyi in range (len(bin_edges)-1):
        bin_mid.append((bin_edges[dummyi]+bin_edges[dummyi+1])/2.0)

    # correct histogram
    if (opts.all_flag and (datarange is not None)):
        leftsum = 0.0
        rightsum = 0.0
        wsum = 0.0
        for dummyi in range(len(datapoints)):
            if (myweights is not None):
                tmpweight = myweights[dummyi]
            else:
                tmpweight = 1.0
            wsum += tmpweight
            if (datapoints[dummyi] < opts.minbin):
                leftsum += tmpweight
            elif (datapoints[dummyi] > opts.maxbin):
                rightsum += tmpweight
        leftsum /= wsum
        rightsum /= wsum
        moresum = leftsum + rightsum        
        if (moresum > 0.0):
            binsize = bin_mid[1] - bin_mid[0]            
            factor = 1.0 - moresum
            for dummyi in range(len(bin_mid)):
                hist[dummyi] *= factor
            if (leftsum > 0.0):
                tmpmid = bin_mid[0] - binsize
                tmphist = leftsum/binsize
                bin_mid.insert(0, tmpmid)
                hist.insert(0, tmphist)
            if (rightsum > 0.0):
                tmpmid = bin_mid[-1] + binsize
                tmphist = rightsum/binsize
                bin_mid.append(tmpmid)
                hist.append(tmphist)
        
    # free energy calcualtion
    if (opts.freeE_flag):
        # unit
        kB = 0.0019872041
        cal2j = 4.184
        if (opts.ounit == 0):
            unitfactor = 1.0
            unitname = 'kBT'
        elif (opts.ounit == 1):
            unitfactor = kB*opts.temp
            unitname = 'kcal/mol'
        elif (opts.ounit == 2):
            unitfactor = kB*opts.temp*cal2j
            unitname = 'kJ/mol'
        else:
            exit("Undefined unit!")
        # free energy
        base = 1.0e1000
        freeElsit = []
        for dummyi in range(len(bin_mid)):
            if (hist[dummyi] == 0.0):
                tmpfe = inf
            else:
                tmpfe = -log(hist[dummyi])*unitfactor # in unit of kB T
            if (tmpfe < base):
                base = tmpfe
            freeElsit.append(tmpfe)
        
    # print 
    for dummyi in range (len(bin_mid)):
        printline = "%20g\t%20g"%(bin_mid[dummyi], hist[dummyi])
        # printline += "%20g\t%20g"%(bin_edges[dummyi], bin_edges[dummyi+1])                
        if (opts.freeE_flag):
            printline += "\t%20g"%(freeElsit[dummyi] - base)
        print(printline)
