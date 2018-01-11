#!/usr/bin/env python

from numpy import histogram2d, loadtxt, inf
from math import log
from sys import argv

############################################################################################################################

if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser()
    parser.add_option('-f', '--datafile', help='data file', dest='inputfile', action='store', type='string')
    parser.add_option('-i', '--column', help='columns', dest='columns', action='store', type='string')
    parser.add_option('-w', '--weightfile', help='weight file', dest='weightfile', action='store', type='string')
    parser.add_option('-j', '--weightcolumn', help='weight column', dest='wcolumn', action='store', type='int')
    parser.add_option('-n', '--minmid', help='the lefe edge of the first bin', dest='minmid', action='store', type='string', default=None)
    parser.add_option('-x', '--maxmid', help='the right edge of the last bin', dest='maxmid', action='store', type='string', default=None)
    parser.add_option('-b', '--nbins', help='the total number of bins', dest='nbins', action='store', type='string')
    parser.add_option('-e', '--freeE', help='calcualte the free energy', dest='freeE_flag', action='store_true', default=False)
    parser.add_option('-o', '--unit', help='unit, 0:kBT, 1:kcal/mol, 2:kj/mol, default=0', dest='ounit', action='store', type='int', default=0)
    parser.add_option('-t', '--temperature', help='temperature, default=300K', dest='temp', action='store', type='float', default=300.0)    

    (opts, args) = parser.parse_args()

    # commandline
    command = " ".join(argv)
    command = "#" + command
    print(command)

    # get options
    xcolumn, ycolumn = [int(item) for item in opts.columns.split(',')]
    xbins, ybins = [int(item) for item in opts.nbins.split(',')]
    if ((opts.minmid is not None) or (opts.maxmid is not None)):    
        xminmid, yminmid = [float(item) for item in opts.minmid.split(',')]
        xmaxmid, ymaxmid = [float(item) for item in opts.maxmid.split(',')]    
        xyrange = [[xminmid, xmaxmid], [yminmid, ymaxmid]]
    else:
        xyrange = None
        
    # read data
    xvalues, yvalues = loadtxt(opts.inputfile, comments=['#', '@'], usecols=[(xcolumn-1), (ycolumn-1)], unpack=True)

    # read weight
    if ((opts.weightfile is not None) and (opts.wcolumn is not None)):
        myweights = loadtxt(opts.weightfile, comments=['#', '@'], usecols=[(opts.wcolumn-1)], unpack=True)
    else:
        myweights = None

    # calculate the histogram
    hist, xedges, yedges = histogram2d(xvalues, yvalues, weights=myweights, bins=[xbins, ybins],
                                       range=xyrange, normed=True)

    # x and y coordinates:
    xbin_mid = []
    ybin_mid = []
    for dummyi in range (len(xedges)-1):
        xbin_mid.append((xedges[dummyi]+xedges[dummyi+1])/2.0)
    for dummyj in range (len(yedges)-1):
        ybin_mid.append((yedges[dummyj]+yedges[dummyj+1])/2.0)
    
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
        for dummyi in range(len(xbin_mid)):
            tmplist = []
            for dummyj in range(len(ybin_mid)):
                if (hist[dummyi][dummyj] == 0.0):
                    tmpfe = inf
                else:
                    tmpfe = -log(hist[dummyi][dummyj])*unitfactor # in unit of kB T
                if (tmpfe < base):
                    base = tmpfe
                tmplist.append(tmpfe)
            freeElsit.append(tmplist)

    #print out
    for dummyi in range(len(xbin_mid)):
        for dummyj in range(len(ybin_mid)):
            printline = "%12g\t%12g\t%12g"%(xbin_mid[dummyi], ybin_mid[dummyj], hist[dummyi][dummyj])
            if (opts.freeE_flag):
                printline += "\t%12g"%(freeElsit[dummyi][dummyj] - base)
            print(printline)
        print()
