# Christa Caggiano, Arya Boudaie, 2015
# this program calculates ANOVA statistics and Tukey HSD post-hoc
# for multiple groups in a many line file. Takes a csv file as input and outputs
# to a csv file.


# imports necessary python libraries to complete program
from pandas import *
import csv
from scipy import stats
from pvalue import rprb

# asks user to input the number of lines they hav
lineNumber=raw_input("How many lines of data do you have for analysis? (i.e. number of genes, replicates etc")

data=read_csv("alldata.csv",header=None)
newData=data.transpose()

def sumofsquares(l):
    sum=0
    for i in l:
        sum+=i[0]**2
    return sum

with open('one_way_anova.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')

    header=["symbol","f-value","p-value","difference DvE", "difference DvM", "difference DvT", "difference EvM", "difference EvT", "difference MvT", "HSD DvE", "HSD DvM", "HSD DvT", "HSD EvM", "HSD EvT", "HSD MvT", "pvalue DvE","pvalue DvM", "pvalue DvT", "pvalue EvM", "pvalue EvT", "pvalue Mvt" ]

    writer.writerow(header)

    for i in range (1,lineNumber):
        symbols=data[0][i]

        morning=newData.ix[2:13,i:i]
        evening=newData.ix[14:25,i:i]
        dn1=newData.ix[26:37,i:i]
        th=newData.ix[38:49,i:i]

        F, p=stats.f_oneway(morning,evening,dn1,th)
        try:
            meanMorning=DataFrame.mean(morning)[1]
        except (IndexError, KeyError):
            meanMorning = float(DataFrame.mean(morning))
        try:
            meanEvening=DataFrame.mean(evening)[1]
        except (IndexError, KeyError):
            meanEvening = float(DataFrame.mean(evening))
        try:
            meanDn1=DataFrame.mean(dn1)[1]
        except (IndexError, KeyError):
            meanDn1 = float(DataFrame.mean(dn1))
        try:
            meanTh=(DataFrame.mean(th))[1]
        except (IndexError, KeyError):
            meanTh = float(DataFrame.mean(th))


        morningPH=DataFrame.as_matrix(morning)
        eveningPH=DataFrame.as_matrix(evening)
        dn1PH=DataFrame.as_matrix(dn1)
        thPH=DataFrame.as_matrix(th)

        for item in morningPH:
            item[0]-=meanMorning
        for item in eveningPH:
            item[0]-=meanEvening
        for item in dn1PH:
            item[0]-=meanDn1
        for item in thPH:
            item[0]-=meanTh

        dn1Evening=meanDn1-meanEvening
        dn1Morning=meanDn1-meanMorning
        dn1Th=meanDn1-meanTh
        eveningMorning=meanEvening-meanMorning
        eveningTh=meanEvening-meanTh
        morningTh=meanMorning-meanTh

        sumsqMorning = (sumofsquares(morningPH))
        sumsqEvening = (sumofsquares(eveningPH))
        sumsqDn1 = (sumofsquares(dn1PH))
        sumsqTh = (sumofsquares(thPH))

        MS=(sumsqMorning+sumsqEvening+sumsqDn1+sumsqTh)/44

        tukeyDE=dn1Evening/(MS/12)**0.5
        tukeyDM=dn1Morning/(MS/12)**0.5
        tukeyDT=dn1Th/(MS/12)**0.5
        tukeyEM=eveningMorning/(MS/12)**0.5
        tukeyET=eveningTh/(MS/12)**0.5
        tukeyMT=morningTh/(MS/12)**0.5

        pvalueDE = rprb(abs(tukeyDE), 4, 44)
        pvalueDM = rprb(abs(tukeyDM), 4, 44)
        pvalueDT = rprb(abs(tukeyDT), 4, 44)
        pvalueEM = rprb(abs(tukeyEM), 4, 44)
        pvalueET = rprb(abs(tukeyET), 4, 44)
        pvalueMT = rprb(abs(tukeyMT), 4, 44)

        anova=[symbols,F,p,dn1Evening,dn1Morning,dn1Th,eveningMorning,eveningTh,morningTh, tukeyDE, tukeyDM, tukeyDT, tukeyEM, tukeyET, tukeyMT, pvalueDE, pvalueDM, pvalueDT, pvalueEM, pvalueET, pvalueMT]


        writer.writerow(anova)