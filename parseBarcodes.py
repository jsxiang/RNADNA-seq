#!/usr/bin/python
import sys
from numpy import *
import re
import matplotlib as mpl
mpl.use ('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import scipy.io
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 
from itertools import repeat
import time

class inputdata:
    def __init__(self,filedata):
        self.data=filedata
        print 'File was loaded...'
        
    def getcol(self,col):
        colword=[]
        for i in xrange(len(self.data)):
            params=re.split('\s',self.data[i])
            colword.append(params[col])
        return colword
    
    def getrow(self,row):
        rowword=[]
        params=re.split('\s',self.data[row])
        rowword.append(params)
        return rowword
    

class seqscount:
    # a struct that puts unique sequences at the center of the universe
    # tallies up the counts for each unique sequence that occurred more than once
    # seqs here are unique sequences, processed and prepared in bash
    # barcoded is a 2d list with sequences assigned to each bin
    # 
    def __init__(self,seqRead,outputfile):
        self.seqRead=seqRead
        self.outputfile=outputfile
        
    def barcode(self,bc):
        # this function outputs a sort of 2D m-by-n array of strings,
        # m being length of barcode pairs, n being number of seqs
        # so that the next step can use re.findall to count each seq in each barcode group
        print 'Binning sequences...'
        bcRead1=bc.getcol(0)
        bcRead2=bc.getcol(1)
        
        barcoded=[[] for i in range(len(bcRead1))]
        barcodearray=[]
        for j in xrange(len(bcRead1)):
            print 'Finding barcode ',j
            for s in self.seqRead:
                if bcRead1[j] in s and bcRead2[j] in s:
                    pat=re.compile('^('+bcRead1[j]+')(\w+)('+bcRead2[j]+'$)')
                    match=re.match(pat,s)
                    if match:
                        barcoded[j].append(match.group(2))

            f=open(self.outputfile+"_barcoded.txt","a")
            out=" ".join(barcoded[j])
            f.write(out+"\n")
            f.close()
            barcodearray.append(out)
        
        self.barcoded=barcodearray



def main():
    start = time.time()

    barcodes_file = sys.argv[1]
    seq_file = sys.argv[2]
    # uniqseq_file = sys.argv[3]
    outputfile = sys.argv[3]
    
    f=open(barcodes_file)
    barcodes=f.readlines()
    f.close()
    
    f=open(seq_file)
    seqs=f.readlines()
    f.close()
    
    # f=open(uniqseq_file)
    # num_uniqseq=len(list(enumerate(f)))
    # f.close()
    
    seqReads=inputdata(seqs)
    bc=inputdata(barcodes)
    
    ubin=seqscount(seqReads.getcol(0),outputfile)
    barcoded = ubin.barcode(bc)
    
    
    
    print 'All done!\n'
    end = time.time()
    print "Elapsed time is",end - start,"\n"
if __name__ == '__main__':
    main()




