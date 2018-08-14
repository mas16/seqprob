'''
This script is for generating peptides of arbitrary length
for testing the accuracy of my statistical chemical shift amino acid
typing program (seqprob) used for assigning resonance frequencies

by MAS 08/2018
'''

#Import libraries
from __future__ import division
import numpy as np
import math
import pylab as plt
import time

start_time=time.time()

#For exporting figures
plt.rcParams['ps.fonttype']=42
plt.rcParams['ps.useafm']=True
plt.rcParams['axes.linewidth'] = 2.0

#################

#Input data path to CA CB shifts
DATAPATH='C:/Users/matt/Desktop/seqprobtest/seqprobtest.txt'

#Length of Peptides to test
PLENGTH=[2,3,4,5,6,7,8,9,10,11,12]

#Easiet to copy FASTA sequence from ExPASy 
SEQUENCE='QSLLIGVATSSLALHAPSQIVAAIKSRADQLGASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIALLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDDTEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLA'

##################

#Read Text File
def read_file(datafile):
    data=file(datafile,'r')
    contents=data.readlines()
    data.close()
    return contents

#Extract relevant information from input
#Convert input discrete vectors
def parse_file(data):
    f,b,a,n=[],[],[],[]
    for x in range(len(data)):
        split=data[x].split()
        if split[0]!='-':
            f.append(float(split[0]))
        else:
            f.append(split[0])
        if split[1]!='-':
            b.append(float(split[1]))
        else:
            b.append(split[1])
        a.append(split[2])
        n.append(int(split[3]))
    return f,b,a,n,

#Generate peptides by iteratively incrememnet by one amino acid
def gen_peptides(freq,beta,assign,number,length):
    pfreq,pb,passign,pnumber=[],[],[],[]
    f,b,a,n,=[],[],[],[]
    i=0
    for x in range(len(freq)-length):
        while i < length:
            pfreq.append(freq[x+i])
            pb.append(beta[x+i])
            passign.append(assign[x+i])
            pnumber.append(number[x+i])
            i+=1
        i=0
        f.append(pfreq)
        b.append(pb)
        a.append(passign)
        n.append(pnumber)
        pfreq,pb,passign,pnumber=[],[],[],[]
    return f,b,a,n

#Remove peptides that have missing data
#Missing data is denoted by '-'
#Note glycine CB shifts are assigned to 1.0
def remove_baddata(freq,beta,assign,number):
    f,b,a,n,=[],[],[],[]
    for x in range(len(freq)):
        c=0
        for y in range(len(freq[x])):
            if freq[x][y]=='-' or beta[x][y]=='-':
                c+=1
        if c==0:
            f.append(freq[x])
            b.append(beta[x])
            a.append(assign[x])
            n.append(number[x])
    return f,b,a,n

#This is the seqprob scoring function
def ScoreL(sequence, CA, CB):

#chemshift is a dictionary of average chemical shifts and standard deviations
#from the BMRB (calculated from a total of 3615400 chemical shifts current as of 10/01/2009)
    chemshift = {'A':[[53.168, 2.892], [19.066, 3.235]],         
                 'R':[[56.791, 3.507], [30.711, 2.612]],
                 'D':[[54.689, 2.831], [40.896, 2.691]],
                 'N':[[53.560, 3.732], [38.738, 3.827]],
                 'C':[[58.058, 3.460], [33.232, 6.442]],
                 'E':[[57.337, 3.432], [30.030, 3.352]],
                 'Q':[[56.566, 2.716], [29.208, 2.627]],
                 'H':[[56.530, 3.522], [30.320, 3.296]],
                 'I':[[61.624, 3.449], [38.607, 3.051]],
                 'L':[[55.647, 2.250], [42.268, 2.038]],
                 'K':[[56.945, 3.354], [32.811, 3.083]],
                 'M':[[56.138, 2.340], [32.973, 3.376]],
                 'F':[[58.113, 3.992], [39.980, 3.804]],
                 'P':[[63.345, 4.030], [31.895, 3.364]],
                 'S':[[58.707, 2.900], [63.714, 5.345]],
                 'T':[[62.211, 2.785], [69.586, 6.055]],
                 'W':[[57.747, 5.210], [30.123, 5.065]],
                 'Y':[[58.145, 3.178], [39.321, 3.200]],
                 'V':[[62.505, 3.236], [32.737, 2.430]],
                 'G':[[45.378, 2.289], [1.00, 1.00]]
                }

    incr = len(CA)
    i = 0
    k = 0
    seqstore=[]
    seqincr=[]
    list1=[]
    list2=[]
    sequence.upper()

    while i <= (len(sequence)-int(incr)):   
        while k < int(incr):
            seqstore.append(sequence[k+i])              
            list1.append(k+i)
            k = k + 1
        seqincr.append(seqstore)
        list2.append(list1)
        list1 = []          
        seqstore = []
        k = 0           
        i = i + 1
    i = 0
    k = 0
    j = 0
    prodlist=[]
    scorelist=[]
    testlist=[] 

    while (k < len(seqincr)):
        while (i < int(incr)):
            comCA = (1/(chemshift[(seqincr[k][i])][0][1]*math.sqrt(2*math.pi)))*math.exp(-1*(abs((CA[int(i)] - chemshift[(seqincr[k][i])][0][0]))**2)/(2*((float(chemshift[(seqincr[k][i])][0][1]))**2)))
            comCB = (1/(chemshift[(seqincr[k][i])][1][1]*math.sqrt(2*math.pi)))*math.exp(-1*(abs((CB[int(i)] - chemshift[(seqincr[k][i])][1][0]))**2)/(2*((float(chemshift[(seqincr[k][i])][1][1]))**2)))
            prod = comCA*comCB
            prodlist.append(prod)
            i = i + 1
        scorelist.append([sum(prodlist)])
        prodlist=[]
        i = 0
        k = k+1

    #This creates an ordered list of scores (highest to lowest) and each score's respective index.
    #Once stored in the list, scores were set to a very very small number instead of removed.
    #This is because identical scores are given the SAME index by default.
    #By doing this you can separate identical scores by index.
    scolist=[]
    k=0
    while k < len(scorelist):
        scolist.append([(max(scorelist)), ((scorelist.index(max(scorelist))))])
        scorelist[scorelist.index(max(scorelist))]=-10000000000000000000
        k = k + 1
    return scolist,list2


def scoreIndex(scolist,list2):
    i=0
    scolistseq=[]
    final=[]
    while i < len(scolist):
        scolistseq.append(list2[scolist[i][1]])
        i+=1
    for x in range(len(scolistseq[0])):
        final.append(scolistseq[0][x]+60)
    return final

#Main
def main(a,b,c,d,plength,sequence):
    clist=[]
    samples=[]
    for y in range(len(plength)):
        print 'looking for peptides with length = '+str(plength[y])
        a1,b1,c1,d1=gen_peptides(a,b,c,d,plength[y])
        a1,b1,c1,d1=remove_baddata(a1,b1,c1,d1)
        correct=0
        for x in range(len(a1)):
            test1,test2=ScoreL(SEQUENCE, a1[x], b1[x])
            test3=scoreIndex(test1,test2)
            if d1[x]==test3:
                correct+=1
            else:
                pass
        percent=correct/len(a1)*100
        samples.append(len(a1))
        clist.append(percent)
        a1,b1,c1,d1=[],[],[],[]
    return clist,samples

if __name__=='__main__':

    DATA=read_file(DATAPATH)
    A,B,C,D=parse_file(DATA)
    DATA=[]
    accuracy,n=main(A,B,C,D,PLENGTH,SEQUENCE)

    print time.time()-start_time

    #Plotting accuracy test results
    yy=np.ones(len(PLENGTH))*90
    xx=np.linspace(1,13,len(PLENGTH))
    plt.plot(PLENGTH,accuracy,'ok',markersize=12)
    plt.plot(xx,yy,'--r',linewidth=4.0)
    plt.axis([1,13,0,110])
    plt.tick_params(axis='both',which='major',labelsize=20,width=2.0,top="off",right="off",direction='out')
    plt.ylabel('Accuracy',fontsize=20)
    plt.xlabel('Length',fontsize=20)
    plt.show()
    
