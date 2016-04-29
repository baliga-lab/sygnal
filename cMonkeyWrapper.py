#################################################################
# @Program: cMonkeyWrapper.py                                   #
# @Version: 2 (python-cMonkey)                                  #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  2/18/2013                      #
#################################################################

import os
import sqlite3 as lite
import gzip
from bicluster import Bicluster
from subprocess import *
from sys import stdout

class cMonkeyWrapper:
    """
    A class designed to hold the information from a cMonkey RData object
    to facilitate downstream analyses."""

    def __init__(self, sqliteDb, maxEValue='NA',
                 meme_upstream=False, weeder_upstream=False,
                 weeder_3pUTR=False, tfbs_db=False,
                 pita_3pUTR=False, targetscan_3pUTR=False,
                 geneConv=False,promoterSeq='NA',p3utrSeq='NA'):
        # What has been run on this cMonkey run, legend = [0: not run, 1: run]
        de_novo_method_upstream = None
        de_novo_method_3pUTR = None
        if meme_upstream and weeder_upstream:
            raise RuntimeError('You trained the same run on both MEME and Weeder! Are you stupid or something?')
        elif meme_upstream:
            de_novo_method_upstream = 'meme'
        elif weeder_upstream:
            de_novo_method_upstream = 'weeder'

        self.meme_upstream = meme_upstream
        self.weeder_upstream = weeder_upstream

        if weeder_3pUTR:
            de_novo_method_3pUTR = 'weeder'

        self.weeder_3pUTR = weeder_3pUTR
        self.tfbs_db = tfbs_db
        self.pita_3pUTR = pita_3pUTR
        self.targetscan_3pUTR = targetscan_3pUTR

        # Attach to the database
        con = lite.connect(sqliteDb)
        con.row_factory = lite.Row
        cur = con.cursor()
        # Get the number of biclusters in run
        q1 = 'SELECT * FROM run_infos'
        cur.execute(q1)
        data = cur.fetchall()
        ks = data[0]['num_clusters']
        # Get the final iteration number
        q1 = 'SELECT max (iteration) from row_members'
        cur.execute(q1)
        self.maxIter = cur.fetchall()[0][0]
        con.close()
        print 'Found '+str(ks)+' clusters.'
        self.biclusters = {}
        for k in range(1,ks+1):
            self.biclusters[k] = Bicluster(k, self.maxIter, de_novo_method_upstream=de_novo_method_upstream, de_novo_method_3pUTR=de_novo_method_3pUTR, sqliteDb=sqliteDb)
            if k%10==0:
                stdout.write(str(k))
            else:
                stdout.write('.')
            stdout.flush()
        # Now read in the upstream sequences
        upstreamSeqsFile = gzip.open(promoterSeq,'rb')
        upstreamSeqsFile.readline() # Skip header
        self.seqsUpstream = {}
        for line in upstreamSeqsFile.readlines():
            splitUp = line.strip().split(',')
            tmp = splitUp[0].strip('"')
            if geneConv==False:
                self.seqsUpstream[tmp] = splitUp[1].strip('"')
            else:
                if tmp in geneConv:
                    for gene in geneConv[tmp]:
                        self.seqsUpstream[gene] = splitUp[1].strip('"')
        upstreamSeqsFile.close()
        # Now read in the 3' UTR sequences
        p3utrSeqsFile = gzip.open(p3utrSeq,'rb')
        p3utrSeqsFile.readline() # Skip header
        self.seqs3pUTR = {}

        for line in p3utrSeqsFile.readlines():
            splitUp = line.strip().split(',')
            tmp = splitUp[0].strip('"')
            if geneConv == False:
                self.seqs3pUTR[tmp] = splitUp[1].strip('"')
            else:
                if tmp in geneConv:
                    for gene in geneConv[tmp]:
                        self.seqs3pUTR[gene] = splitUp[1].strip('"')
        p3utrSeqsFile.close()
        # Now read in nucleotide frequencies
        nucFreqsFile = open('seqs/nucFreqs.csv','r')
        nucFreqsFile.readline()  # Skip header
        upFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqsUpstream = {'A': upFreq[1], 'C': upFreq[2],
                                 'G': upFreq[2], 'T': upFreq[1]}
        p3utrFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqs3pUTR = {'A': p3utrFreq[1], 'C': p3utrFreq[2],
                              'G': p3utrFreq[2], 'T': p3utrFreq[1]}
        nucFreqsFile.close()
        # Close database connection
        con.close()
        print '\nDone loading.\n'

    # Get all Upstream pssms
    def pssms_upstream(self, maxNormResid='NA', maxEValue='NA', maxSurv='NA',
                       de_novo_method='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = False
            if (maxNormResid == 'NA' or
                float(self.biclusters[bi].getNormResidual()) <= float(maxNormResid)):
                if maxSurv == 'NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue']) <= float(maxSurv):
                    biOk = True
            if biOk:
                tmpPssms = self.biclusters[bi].pssms_upstream
                for pssm in tmpPssms:
                    if de_novo_method == 'NA' or de_novo_method == pssm.de_novo_method:
                        # Only add it if it is less than an E-Value threshold
                        if (maxEValue == 'NA' or
                            float(pssm.evalue) <= float(maxEValue)):
                            pssms.append(pssm)
                            pssmsNames.append(pssm.name)
        return dict(zip(pssmsNames, pssms))

    # Get all 3' UTR pssms
    def pssms_3putr(self, maxNormResid='NA', maxEValue='NA', maxSurv='NA',
                    de_novo_method='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = False
            if (maxNormResid == 'NA' or
                float(self.biclusters[bi].getNormResidual()) <= float(maxNormResid)):
                if maxSurv == 'NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue']) <= float(maxSurv):
                    biOk = True
            if biOk:
                tmpPssms = self.biclusters[bi].pssms_3putr
                for pssm in tmpPssms:
                    # Only add it if it is less than an E-Value threshold
                    if de_novo_method == 'NA' or de_novo_method == pssm.de_novo_method:
                        if maxEValue == 'NA' or float(pssm.evalue) <= float(maxEValue):
                            pssms.append(pssm)
                            pssmsNames.append(pssm.name)
        return dict(zip(pssmsNames,pssms))

    def bicluster_seqs_upstream(self, k):
        """returns the upstream sequences for a bicluster as a dictionary of
        {<gene_name>: <seqeunce>, ...}"""
        genes = self.biclusters[k].genes
        seqs = dict(zip([gene for gene in genes if gene in self.seqsUpstream],
                        [self.seqsUpstream[gene] for gene in genes
                         if gene in self.seqsUpstream]))
        return seqs

    def bicluster_seqs_3putr(self, k):
        """returns the 3' UTR sequences for a bicluster as a dictionary of
        {<gene_name>: <seqeunce>, ...}"""
        genes = self.biclusters[k].genes
        seqs = dict(zip([gene for gene in genes if gene in self.seqs3pUTR],
                        [self.seqs3pUTR[gene] for gene in genes
                         if gene in self.seqs3pUTR]))
        return seqs
