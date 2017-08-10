#!/usr/bin/python
'''Plot the cross-over density with different temperatures
'''
import numpy as np
import matplotlib.pyplot as plt
import re, os, glob, io, sys
from scipy import stats

# Read critical rho values
def readData(file_list):
    datacollect = {};
    summary_title = re.compile('^T\tbeta\*p')
    for file in file_list:
        basename = os.path.splitext(os.path.basename(file))[0]
        fin = open(file, 'r')
        lines = fin.readlines()
        fin.close()
        for index, line in enumerate(lines):
            if re.search(summary_title, line):
                buff = '\n'.join(lines[index+1:])
                data = np.genfromtxt(io.BytesIO(buff))
                xi = re.search('_(\d+\.\d+)$', basename).group(1)
                datacollect[xi] = data
                break
    return datacollect


# Read not only critical rho values, but also the neighbor points above
# critical rho, to calculate the B3(cluster)
def readDataPlus(file_list):
    criticalcollect = {}
    rawcollect = {}
    summary_title = re.compile('^T\tbeta\*p')
    for file in file_list:
        basename = os.path.splitext(os.path.basename(file))[0]
        fin = open(file, 'r')
        lines = fin.readlines()
        fin.close()
        for index, line in enumerate(lines):
            if re.search(summary_title, line):
                xi = re.search('_(\d+\.\d+)$', basename).group(1)
                rawcollect[xi] = lines[:index]
                criticalcollect[xi] = lines[index+1:]
                break
    return criticalcollect, rawcollect


# Compute cluster B3 by the given transition point and other points
# return: a dict, with keys of each xi, values of
#  a list, with row correspond to every temperature,
#    every row include T, B3(cluster), corelation coefficient r of fitting, p_value, std_err
def clusterB3searcher(criticalcollect, rawcollect):
    B3collect = {}
    Xis = sorted(criticalcollect.iterkeys())
    for xi in Xis:
        criticallist = criticalcollect[xi]
        rawlist = rawcollect[xi]
        buff = '\n'.join(criticallist)
        cdata = np.genfromtxt(io.BytesIO(buff))
        # Filter out the points when no transition detected
        cutlow, cuthigh = Tfilter(cdata[:, 0], cdata[:,2], returnchoice=1)
        criticallist = criticallist[cutlow:cuthigh]
        rawlist_splitted = [line.split() for line in rawlist]
        Ts = cdata[cutlow:cuthigh, 0]
        result = []
        for criticalline in criticallist:
            T_str, rho_str, _ = criticalline.split(None, 2)
            rhoc_float = float(rho_str)
            # Extract raw data from this particular temperature
            # Include values: T, p*beta, rho, h
            # We need rho and h
            # filter the [rho, h] within the range (rho_critical, 10*rho_critical) to make linear fitting
            rhos = []
            hs = []
            for aT, _, arho, ah in rawlist_splitted:
                if aT==T_str:
                    rho_float = float(arho)
                    if rho_float>rhoc_float and rho_float<10*rhoc_float:
                        rhos.append(rho_float)
                        hs.append(float(ah))
            if rhos:
                B3s, intercept, r_value, p_value, std_err = stats.linregress(rhos,hs)
            result.append([T_str, B3s, intercept, r_value, p_value, std_err])
        B3collect[xi] = result
    return B3collect


# filter abnormal datapoints
# return: Filted data, if returnchoice is not 1;
# lower and upper bound, if returnchoice is 1
def Tfilter(Ts, rhos, returnchoice=0):
    # first filt low temperatures
    statusflag = 0
    cutlow = 0
    cuthigh = len(rhos)
    highbound = Ts[cuthigh-1]
    for index, val in enumerate(rhos):
        if 0==statusflag:  # first element
            statusflag = 1
        elif 1==statusflag:
            if previusrho>val:
                cutlow = index
            else:
                statusflag = 2
        elif 2==statusflag:
            if previusrho > val:
                cuthigh = index
                highbound = Ts[cuthigh-1]
                break
        previusrho = val
    if 1==returnchoice:
        return cutlow, cuthigh
    else:
        return Ts[cutlow:cuthigh], rhos[cutlow:cuthigh], highbound


def printB3s(B3collect):
    for xi in sorted(B3collect.iterkeys()):
        print("Xi="+xi)
        print('T\tSlope\tIntercept\tr')
        for Ts, B3s, inter, r_value, _, _ in B3collect[xi]:
            print('%s\t%.3f\t%f\t%f' % (Ts, B3s, inter, r_value))


def statCritical(datacollect):
    data_out = {}
    for xi in sorted(datacollect.iterkeys()):
        data = datacollect[xi]
        Ts = data[:, 0]
        rhos = data[:, 2]
        Ts, rhos, criticalT = Tfilter(Ts, rhos)
        data_out[xi] = [Ts, rhos, criticalT]
    return data_out


def printRaw(datacollect):
    data_out = statCritical(datacollect)
    for xi in sorted(data_out):
        content = data_out[xi]
        for rho, y in zip(content[0], content[1]):
            print "%s\t%.3f\t%.3e" % (xi, rho, y)


def draw(datacollect):
    data_out = statCritical(datacollect)
    ccoeff = []  # critical T at a zeta
    cTs = []     # critical T
    plt.figure(figsize=(12,5))
    ax1 = plt.subplot(121)

    for key in sorted(data_out):
        content = data_out[key]
        plt.plot(content[0], content[1], label=key)  # Ts, rhos
        ccoeff.append(key)
        cTs.append(content[2])
    plt.yscale('log')
    plt.xlabel('$T$')
    plt.ylabel(r'$\rho_{CMC}$')
    plt.legend(loc=2,prop={'size':10})

    ax2 = plt.subplot(122)
    plt.plot(ccoeff, cTs)
    floatzetas = [float(tmp) for tmp in ccoeff]
    ax2.set_xlim([min(floatzetas) - 0.1, max(floatzetas) + 0.1])
    filtedcTs = filter(lambda x: x is not None, cTs)
    ax2.set_ylim([min(filtedcTs)-0.01, max(filtedcTs)+0.01])
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$T_{CMC}$')
    plt.show()

if __name__ == '__main__':
    path = 'data_min'
    file_list = [f for f in glob.glob(path+"/*.dat") if re.search(r'1.0_2.2_4.0_1.0_\d+\.\d+\.dat$', f)]
    datacollect = readData(file_list)
    if len(sys.argv)<2 or sys.argv[1]=='draw':
        draw(datacollect)
    elif sys.argv[1]=='raw':
        # ./plotmin.py raw > pyout/rhoc-1,2.2,4,1.dat
        printRaw(datacollect)
    elif sys.argv[1]=='critical':
        # ./plotmin.py critical > pyout/Tc-1,2.2,4,1.dat
        data_out = statCritical(datacollect)
        for key in sorted(data_out):
            content = data_out[key]
            print('%s\t%.3f' % (key, content[2]))
    elif sys.argv[1]=='B3':
        criticalcollect, rawcollect = readDataPlus(file_list)
        B3collect = clusterB3searcher(criticalcollect, rawcollect)
        printB3s(B3collect)
    else: print '''
Usage: plotmin.py [options]
Options:
* draw (default)
    Draw T - rho_critical, and T_critical - Xi
* raw
    Print raw data
* critical
    print critical temperature of different xi
* B3
    Fitting the function and obtain the third virial coefficent B3(cluster)
'''
