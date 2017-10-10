from __future__ import print_function
from bisect import bisect_left
from sklearn import metrics
import numpy as np
import sys
import os

def removeNaN(a):
    return a[~np.isnan(a).any(axis=1)]

def createDir(directory):
    if not os.path.exists(directory):
        print("Creating ", directory, " . . . ",)
        os.makedirs(directory)
        print("Done.")

def formatDir(directory):
    if not os.path.isdir(directory):
        print(directory, " not found!")
        exit(1)
    if directory[-1] != '/':
        return directory + '/'
    return directory

def plotROC(ax, trueY, scores, pos_label, type):
    print("Plotting ROC curve for ", type, " . . . ",)
    fpr, tpr, thresholds = metrics.roc_curve(trueY, scores, pos_label=pos_label)
    roc_auc = metrics.auc(fpr, tpr)

    #fpr = fpr[fpr <= 0.15]
    #thresholds = thresholds[:50]
    ax.plot(fpr, tpr, label=type + ' (area = ' + str(round(roc_auc, 2)) + ')')
    print("Done.")

def getExons(chr, exonDir):
    exons = []
    exonFile = open(exonDir + chr + "_exons.out", "r")
    for line in exonFile:
        exons.append(int(line.strip().split(",")[0]))
        exons.append(int(line.strip().split(",")[1]))
    return exons

def isInInterval(a, x):
    pos = bisect_left(a, x, 0, len(a))
    if pos < len(a) and a[pos] == x: # x is the very end or beginning of an interval
        return True
    else:
        if pos % 2 == 0: # x is after the end of an interval
            return False
        return True # x is after the start of an interval

def removeExons(bases, exons):
    basesNoExons = []
    for base in bases:
        if not isInInterval(exons, base):
            basesNoExons.append(base)
    return basesNoExons
