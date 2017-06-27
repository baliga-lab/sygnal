#################################################################
# @Program: miRvestigator.py                                    #
# @Version: 1                                                   #
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
# Copyrighted by Chris Plaisier  8/12/2009                      #
#################################################################

from copy import deepcopy
from sys import stdout
import os, cPickle
import pssm as pssm_mod


class miRvestigator:
    """A class designed to compute and hold the information from analyzing
    miRNA seeds against motifs from 3' UTRs."""
    def __init__(self, pssms, seqs3pUTR, seedModel=[6, 7, 8], minor=True,
                 p5=True, p3=True, textOut=True, wobble=True, wobbleCut=0.25,
                 baseDir='', outName='', species='hsa', miRNA_dir='miRNA'):
        print '\nmiRvestigator analysis started...'
        self.pssms = pssms
        self.species = species
        self.miRNA_dir = miRNA_dir  # must come before setMiRNAs
        self.miRNAs = self.setMiRNAs(0,8,minor,p5,p3)

        # Trim sequences down
        self.miRNAs_6mer_1 = self.trimSeqs(deepcopy(self.miRNAs),0,6)
        self.miRNAs_6mer_2 = self.trimSeqs(deepcopy(self.miRNAs),1,7)
        self.miRNAs_6mer_3 = self.trimSeqs(deepcopy(self.miRNAs),2,8)
        self.miRNAs_7mer_m8 = self.trimSeqs(deepcopy(self.miRNAs),1,8)
        self.miRNAs_7mer_a1 = self.trimSeqs(deepcopy(self.miRNAs),0,7)
        self.miRNAs_8mer = self.trimSeqs(deepcopy(self.miRNAs),0,8)
        p3utrSeqs = 'X'.join(seqs3pUTR)

	if not baseDir=='':
            dirName = baseDir+'/miRNA'
	else:
            dirName = 'miRNA'

        if not os.path.exists(dirName):
            os.mkdir(dirName)
        if 6 in seedModel:
            print 'Screening out 6mers not present in 3\' UTRs...'
            if not os.path.exists(dirName+'/permKMers_6mers.pkl'):
                permKMers_6mer = self.allKmers(6)
                tmpKMers = []
                for i in permKMers_6mer:
                    if not p3utrSeqs.find(i)==-1:
                        tmpKMers.append(i)
                self.permKMers_6mer = tmpKMers
                pklFile = open(dirName+'/permKMers_6mers.pkl','wb')
                cPickle.dump(self.permKMers_6mer,pklFile)
            else:
                pklFile = open(dirName+'/permKMers_6mers.pkl','rb')
                self.permKMers_6mer = cPickle.load(pklFile)
            pklFile.close()

        if 7 in seedModel:
            print 'Screening out 7mers not present in 3\' UTRs...'
            if not os.path.exists(dirName+'/permKMers_7mers.pkl'):
                permKMers_7mer = self.allKmers(7)
                tmpKMers = []
                for i in permKMers_7mer:
                    if not p3utrSeqs.find(i)==-1:
                        tmpKMers.append(i)
                self.permKMers_7mer = tmpKMers
                pklFile = open(dirName+'/permKMers_7mers.pkl','wb')
                cPickle.dump(self.permKMers_7mer,pklFile)
            else:
                pklFile = open(dirName+'/permKMers_7mers.pkl','rb')
                self.permKMers_7mer = cPickle.load(pklFile)
            pklFile.close()

        if 8 in seedModel:
            print 'Screening out 8mers not present in 3\' UTRs...'
            if not os.path.exists(dirName+'/permKMers_8mers.pkl'):
                permKMers_8mer = self.allKmers(8)

                tmpKMers = []
                for i in permKMers_8mer:
                    if not p3utrSeqs.find(i)==-1:
                        tmpKMers.append(i)
                self.permKMers_8mer = tmpKMers
                pklFile = open(dirName+'/permKMers_8mers.pkl','wb')
                cPickle.dump(self.permKMers_8mer,pklFile)
            else:
                pklFile = open(dirName+'/permKMers_8mers.pkl','rb')
                self.permKMers_8mer = cPickle.load(pklFile)
            pklFile.close()
        print 'Done.\n'
        miRNAScores = {}
        cur = 1
        # Building HMM Model
        outMe = []
        for pssm in pssms:
            print '\n%s' % pssm.name
            print 'Building HMM model for ' + pssm_mod.consensus_motif(pssm) + '...'
            miRNAScores[pssm.name] = ['NA','NA']
            # Then setup the HMM
            ## States ##
            ## and Starting Probabilities ##
            # NM1 = no match 1
            # NM2 = no match 2
            # PSSMi = PSSM at spot i
            maxPSSMi = len(pssm.matrix)
            states = ['NM1', 'NM2']
            sp = {'NM1': float(1)/float(maxPSSMi+1), 'NM2': 0}
            # Add the PSSM states
            for i in range(maxPSSMi):
                states += ['PSSM'+str(i)]
                sp['PSSM'+str(i)] = float(1)/float(maxPSSMi+1)
                if wobble==True:
                    states += ['WOBBLE'+str(i)]
                    sp['WOBBLE'+str(i)] = 0
            ## Transition probabilities
            tp = {}
            # NM1
            nm1_2_nm1 = 0.01
            tp['NM1'] = { 'NM1': nm1_2_nm1, 'NM2': 0 }
            leftOver1 = float(1-nm1_2_nm1)/float(maxPSSMi)
            for i in range(maxPSSMi):
                tp['NM1']['PSSM'+str(i)] = leftOver1
                if wobble==True:
                    tp['NM1']['WOBBLE'+str(i)] = 0 # Don't start a seed with a wobble
            # NM2
            tp['NM2'] = { 'NM1': 0, 'NM2': 1 }
            for i in range(maxPSSMi):
                tp['NM2']['PSSM'+str(i)] = 0
                if wobble==True:
                    tp['NM2']['WOBBLE'+str(i)] = 0
            # PSSMis
            for i in range(maxPSSMi):
                tp['PSSM'+str(i)] = { 'NM1': 0, 'NM2': 0.01 }
                if wobble==True:
                    tp['WOBBLE'+str(i)] = { 'NM1': 0, 'NM2': 0.01 }
                if i==(maxPSSMi-1):
                    tp['PSSM'+str(i)]['NM2'] = 1
                    if wobble==True:
                        tp['WOBBLE'+str(i)]['NM2'] = 1
                for j in range(maxPSSMi):
                    if j == i+1:
                        if wobble==True:
                            # Allow wobbly matches if T is >= wobbleCut
                            if float(pssm.matrix[i+1][2])>=float(wobbleCut) or float(pssm.matrix[i+1][3])>=float(wobbleCut):
                                tp['PSSM'+str(i)]['PSSM'+str(j)] = 0.80
                                tp['PSSM'+str(i)]['WOBBLE'+str(j)] = 0.19
                            # Otherwise don't allow wobbly matches
                            else:
                                tp['PSSM'+str(i)]['PSSM'+str(j)] = 0.99
                                tp['PSSM'+str(i)]['WOBBLE'+str(j)] = 0
                            tp['WOBBLE'+str(i)]['PSSM'+str(j)] = 1
                            tp['WOBBLE'+str(i)]['WOBBLE'+str(j)] = 0
                        else:
                            tp['PSSM'+str(i)]['PSSM'+str(j)] = 0.99
                    else:
                        tp['PSSM'+str(i)]['PSSM'+str(j)] = 0
                        if wobble==True:
                            tp['PSSM'+str(i)]['WOBBLE'+str(j)] = 0
                            tp['WOBBLE'+str(i)]['PSSM'+str(j)] = 0
                            tp['WOBBLE'+str(i)]['WOBBLE'+str(j)] = 0
            ## Emission probabilities
            ep = {}
            # NM1
            ep['NM1'] = { 'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25 }
            # S1 - None
            # S2 - None
            # NM2
            ep['NM2'] = { 'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25 }
            # PSSMis
            for i in range(maxPSSMi):
                ep['PSSM'+str(i)] = { 'A': pssm.matrix[i][0], 'C': pssm.matrix[i][1], 'G': pssm.matrix[i][2], 'T': pssm.matrix[i][3] }
                if wobble==True:
                    # If motif has both G and U probability greater than wobblecut or random (0.25)
                    if float(pssm.matrix[i][2])>=float(wobbleCut) and float(pssm.matrix[i][3])>=float(wobbleCut):
                        ep['WOBBLE'+str(i)] = { 'A': 0.5, 'C': 0.5, 'G': 0, 'T': 0 }
                    # If motif has G greater than wobblecut or random (0.25)
                    elif float(pssm.matrix[i][2])>=float(wobbleCut):
                        ep['WOBBLE'+str(i)] = { 'A': 1, 'C': 0, 'G': 0, 'T': 0 }
                    # If motif has U greater than wobblecut or random (0.25)
                    elif float(pssm.matrix[i][3])>=float(wobbleCut):
                        ep['WOBBLE'+str(i)] = { 'A': 0, 'C': 1, 'G': 0, 'T': 0 }
                    # Otherwise be random (0.25 x 4)
                    else:
                        ep['WOBBLE'+str(i)] = { 'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25 }
            #print sp, ep, tp
            print 'Done.\n'

            print 'Starting miRNA detection for ' + pssm_mod.consensus_motif(pssm) + ':'
            # First try for perfect 8mer
            #  a. Do the Viterbi for all miRNAs
            vitMirnas = []
            maxVitP = 0
            testVitWiki = []
            for miRNA in self.miRNAs_6mer_1:
                #forVit_8mer = self.forwardViterbi(list(self.miRNAs_8mer[miRNA]), states, sp, tp, ep)
                forVit_8mer = self.viterbi(list(self.miRNAs_8mer[miRNA]), states, sp, tp, ep)
                forVit_8mer = [0,0,forVit_8mer[0]]
                if float(forVit_8mer[2]) > float(maxVitP):
                    maxVitP = float(forVit_8mer[2])
                    vitMirnas = [miRNA]
                elif float(forVit_8mer[2])==float(maxVitP):
                    vitMirnas.append(miRNA)
            print maxVitP, vitMirnas
            totPs_8mer = []
            vitPs_8mer = []
            hits_8mer = 0
            for kMer in self.permKMers_8mer:
                #permVit = self.forwardViterbi(list(kMer), states, sp, tp, ep)
                permVit = self.viterbi(list(kMer), states, sp, tp, ep)
                permVit = [0,0,permVit[0]]
                if float(permVit[2]) > float(maxVitP):
                    hits_8mer = 2
                    break
                elif float(permVit[2]) == float(maxVitP):
                    hits_8mer += 1
                    if hits_8mer > 1:
                        break
                totPs_8mer.append(permVit[0])
                vitPs_8mer.append(permVit[2])
            if hits_8mer<=1:
                print '8mer match!'
                outMe.append(pssm.name + ',' + '_'.join(vitMirnas) + ',8mer')

            if hits_8mer>1:
                # Then try for perfect 7mer-m8
                vitMirnas = []
                maxVitP = 0
                for miRNA in self.miRNAs_6mer_1:
                    #forVit_7mer_m8 = self.forwardViterbi(list(self.miRNAs_7mer_m8[miRNA]), states, sp, tp, ep)
                    forVit_7mer_m8 = self.viterbi(list(self.miRNAs_7mer_m8[miRNA]), states, sp, tp, ep)
                    forVit_7mer_m8 = [0,0,forVit_7mer_m8[0]]
                    if float(forVit_7mer_m8[2]) > float(maxVitP):
                        maxVitP = float(forVit_7mer_m8[2])
                        vitMirnas = [miRNA]
                    elif float(forVit_7mer_m8[2])==float(maxVitP):
                        vitMirnas.append(miRNA)
                print maxVitP, vitMirnas
                totPs_7mer = []
                vitPs_7mer = []
                hits_7mer_m8 = 0
                for kMer in self.permKMers_7mer:
                    #permVit = self.forwardViterbi(list(kMer), states, sp, tp, ep)
                    permVit = self.viterbi(list(kMer), states, sp, tp, ep)
                    permVit = [0,0,permVit[0]]
                    if float(permVit[2]) > float(maxVitP):
                        hits_7mer_m8 = 2
                        break
                    elif float(permVit[2]) >= float(maxVitP):
                        hits_7mer_m8 += 1
                        if hits_7mer_m8 > 1:
                            break
                    totPs_7mer.append(permVit[0])
                    vitPs_7mer.append(permVit[2])
                if hits_7mer_m8<=1:
                    print '7mer-m8 match!'
                    outMe.append(pssm.name + ',' + '_'.join(vitMirnas) + ',7mer_m8')

                # Finally try for perfect 7mer-a1
                vitMirnas = []
                maxVitP = 0
                for miRNA in self.miRNAs_6mer_1:
                    #forVit_7mer_a1 = self.forwardViterbi(list(self.miRNAs_7mer_a1[miRNA]), states, sp, tp, ep)
                    forVit_7mer_a1 = self.viterbi(list(self.miRNAs_7mer_a1[miRNA]), states, sp, tp, ep)
                    forVit_7mer_a1 = [0,0,forVit_7mer_a1[0]]
                    if float(forVit_7mer_a1[2]) > float(maxVitP):
                        maxVitP = float(forVit_7mer_a1[2])
                        vitMirnas = [miRNA]
                    elif float(forVit_7mer_a1[2])==float(maxVitP):
                        vitMirnas.append(miRNA)
                print maxVitP, vitMirnas
                totPs_7mer = []
                vitPs_7mer = []
                hits_7mer_a1 = 0
                for kMer in self.permKMers_7mer:
                    #permVit = self.forwardViterbi(list(kMer), states, sp, tp, ep)
                    permVit = self.viterbi(list(kMer), states, sp, tp, ep)
                    permVit = [0,0,permVit[0]]
                    if float(permVit[2]) > float(maxVitP):
                        hits_7mer_a1 = 2
                        break
                    elif float(permVit[2]) >= float(maxVitP):
                        hits_7mer_a1 += 1
                        if hits_7mer_a1 > 1:
                            break
                    totPs_7mer.append(permVit[0])
                    vitPs_7mer.append(permVit[2])
                if hits_7mer_a1<=1:
                    print '7mer-a1 match!'
                    outMe.append(pssm.name + ',' + '_'.join(vitMirnas) + ',7mer_a1')
            if hits_8mer > 1 and hits_7mer_m8 > 1 and hits_7mer_a1 > 1:
                print 'No match!'
                outMe.append(pssm.name + ',NA,NA')

        print 'miRvestigator analysis completed.\n'
        outFile = open(dirName+'/scores'+str(outName)+'.csv','w')
        outFile.write('pssm,miRNAs,match_type')
        outFile.write('\n'+'\n'.join(outMe))
        outFile.close()

    def getScore(self, forVit, vitPs):
        return (float(forVit[0])/float(forVit[2]))*float(self.getPValue(forVit[2],vitPs))

    # Generates all possible sequences with letters to the length of depth.
    def allKmers(self,depth,letters=['A','C','G','T'],seqs=[''],curdepth=0):
        newseqs = []
        for seq in seqs:
            for letter in letters:
                newseqs.append(seq + letter)
        if depth > curdepth:
            return(self.allKmers(depth,letters,newseqs,curdepth + 1))
        else:
            return(seqs)

    # Get the miRNAs to compare against
    def setMiRNAs(self,seedStart,seedEnd, minor=True, p5=True, p3=True):
        mature_fa_gz_path = os.path.join(self.miRNA_dir, 'mature.fa.gz')
        if not os.path.exists(mature_fa_gz_path):
            print '\nDownloading miRNA seeds from miRBase.org...'
            # Grab down the latest miRNA data from mirbase.org:
            #  ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
            from ftplib import FTP
            ftp1 = FTP('mirbase.org')
            ftp1.login()
            ftp1.cwd('/pub/mirbase/CURRENT/')
            # Get the miRBase.org version number for reference.
            self.miRNAver = (ftp1.pwd().split('/'))[-1]
            outFile = open(mature_fa_gz_path,'wb')
            ftp1.retrbinary('RETR mature.fa.gz', outFile.write)
            outFile.close()
            ftp1.quit()
            print 'Done.\n'
        else:
            print '\nUsing already downloaded miRNA seeds.\n'

        # Read in miRNAs: miRNAs are labeled by the hsa-* names and grabbing 2-8bp
        ### Could merge these as they come in so that don't do redundant, and also so that the labels are together
        import gzip
        miRNAFile = gzip.open(mature_fa_gz_path, 'r')
        miRNAs = {}
        while 1:
            miRNALine = miRNAFile.readline()
            seqLine = miRNAFile.readline()
            if not miRNALine:
                break
            # Get the miRNA name
            if miRNALine[0:4]=='>hsa':
                # Now grab out the 2-8bp and do reverse complement on it
                curMiRNA = miRNALine.split(' ')[1]
                miRNAs[curMiRNA] = self.reverseComplement((seqLine.strip())[seedStart:seedEnd])
        miRNAFile.close()

        # How many distinct kMers in miRNAs
        miRNAuniq = {}
        for miRNA in miRNAs:
            if not miRNAs[miRNA] in miRNAuniq:
                miRNAuniq[miRNAs[miRNA]] = [miRNA]
            else:
                miRNAuniq[miRNAs[miRNA]].append(miRNA)
        # Then merge them
        miRNAsMerged = {}
        for seed in miRNAuniq:
            tmpName = '_'.join(miRNAuniq[seed])
            miRNAsMerged[tmpName] = seed
        return miRNAsMerged

    def trimSeqs(self, miRNAs, start, stop):
        tmp = {}
        for miRNA in miRNAs:
            tmp[miRNA] = miRNAs[miRNA][start:stop]
            #print tmp[miRNA]
        return tmp

    # Complement
    def complement(self,seq):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'U':'A'}
        complseq = [complement[base] for base in seq]
        return complseq

    # Reverse complement
    def reverseComplement(self,seq):
        seq = list(seq)
        seq.reverse()
        return ''.join(self.complement(seq))

    # Reverse complement
    def reverseMe(self,seq):
        seq = list(seq)
        seq.reverse()
        return ''.join(seq)

    # Modified from From wikipedia to do both forward calculation and viterbi
    def forwardViterbi(self, obs, states, start_p, trans_p, emit_p):
       T = {}
       for state in states:
           ##          prob.           V. path  V. prob.
           T[state] = (start_p[state], [state], start_p[state])
       for output in obs:
           U = {}
           for next_state in states:
               total = float(0)
               argmax = None
               valmax = float(0)
               for source_state in states:
                   (prob, v_path, v_prob) = T[source_state]
                   p = emit_p[source_state][output] * trans_p[source_state][next_state]
                   prob *= p
                   v_prob *= p
                   total += prob
                   if v_prob > valmax:
                       argmax = v_path + [next_state]
                       valmax = v_prob
               U[next_state] = (total, argmax, valmax)
           T = U
       ## apply sum/max to the final states:
       total = float(0)
       argmax = None
       valmax = float(0)
       for state in states:
           (prob, v_path, v_prob) = T[state]
           total += prob
           if v_prob > valmax:
               argmax = v_path
               valmax = v_prob
       return (total, argmax, valmax)

    def viterbi(self, obs, states, start_p, trans_p, emit_p):
        V = [{}]
        path = {}

        # Initialize base cases (t == 0)
        for y in states:
            V[0][y] = start_p[y] * emit_p[y][obs[0]]
            path[y] = [y]

        # Run Viterbi for t > 0
        for t in range(1,len(obs)):
            V.append({})
            newpath = {}

            for y in states:
                (prob, state) = max([(V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states])
                V[t][y] = prob
                newpath[y] = path[state] + [y]

            # Don't need to remember the old paths
            path = newpath

        #print_dptable(V)
        (prob, state) = max([(V[len(obs) - 1][y], y) for y in states])
        return (prob, path[state])


    # Decide whether to add and add at correct position
    # Strucuture of elements:
    # [pssm, getConsensus(pssms[pssm]), miRNA, miRNAs[miRNA], forVit]
    def addSorted(self,all,new):
        inserted = 0
        if len(all)>0:
            for i in range(len(all)):
                if new['vitPValue']<all[i]['vitPValue']:
                    all.insert(i,new)
                    inserted = 1
                    break
            if inserted==0:
                all.append(new)
        else:
            all.append(new)
        return all

    # Get p-value from a probability and a distribution
    def getPValue(self,prob,dist):
        a = float(len([i for i in dist if i>=prob]))
        b = float(len(dist))
        if not a==0:
            return a/b
        else:
            return 0

    # Get score list for a PSSM
    def getScoreList(self,pssmName):
        return self.miRNAScores[pssmName]

    # Get top hit(s) for a PSSM.
    # Returns either one or more miRNAs based on whether a clear winner or a tie.
    def getTopHit(self,pssmName):
        scoreList = self.getScoreList(pssmName)
        if scoreList[0]['vitPValue']<scoreList[1]['vitPValue']:
            return [scoreList[0]]
        else:
            retMe = []
            i = 0
            while scoreList[0]['vitPValue']==scoreList[i]['vitPValue']:
                retMe.append(scoreList[i])
                i += 1
            return retMe

    # Get the scores for a PSSM.
    # Returns all miRNAs which conatin the miRNA name.
    def getmiRNAHit(self,pssmName,miRNAname):
        scoreList = self.getScoreList(pssmName)
        retMe = []
        for i in range(len(scoreList)):
            if not scoreList[i]['miRNA.name'].find(miRNAname)==-1:
                scoreList[i]['rank'] = i
                retMe.append(scoreList[i])
        return retMe

    # Strucuture of elements:
    # [pssm, getConsensus(pssms[pssm]), miRNA, miRNAs[miRNA], forVit]
    # <tr><td>miRNA Name</td><td>Alignment Start<sup>*</sup></td><td>Alignment Stop<sup>*</sup></td><td>Alignment Length</td><td>Motif (Length)</td><td>Alignment Type</td><td>Alignment</td><td>P(Alignment)</td></tr>
    def outHtml(self, outMe):
        writeMe = ''
        # miRNA name
        writeMe += '<tr align=\'center\' valign=\'center\'><td>'+str(outMe[2])+'</td><td>'+str(outMe[3])+'</td>'
        alignment = outMe[4][1] # Grab the alignment
        alignment.pop() # Get rid of the extra state which is added by the forwardViterbi function
        start = 1
        if alignment[0]=='NM1':
            for i in alignment:
                if i=='NM1':
                    start += 1
        # Alignment
        seedAlign = ''
        seed = list(outMe[3])
        motifAlign = ''
        motif = list(outMe[1])
        alignMe = alignment
        aligned = ''
        lenMatch = 0
        # First get them zero'd to the same point
        if start>1:
            for i in range(start-1):
                seedAlign += seed.pop(0)
                aligned += ' '
                motifAlign += '-'
                alignMe.pop(0)
        if not alignMe[0]=='PSSM0' and not alignMe[0]=='WOBBLE0':
            for i in range(int(alignMe[0][4])):
                seedAlign += '-'
                aligned += ' '
                motifAlign += motif.pop(0)
        # Then add the parts that align
        while 1:
            if len(alignMe)==0 or alignMe[0]=='NM2':
                break
            seedAlign += seed.pop(0)
            if alignMe[0][0]=='P':
                aligned += '|'
            elif alignMe[0][0]=='W':
                aligned += ':'
            lenMatch += 1
            motifAlign += motif.pop(0)
            alignMe.pop(0)
        # Then do the ending
        if len(alignMe)>0:
            for i in alignMe:
                seedAlign += seed[0]
                seed = seed[1:]
                aligned += ' '
                motifAlign += '-'
            alignMe = []
        if len(motif)>0 and len(alignMe)==0:
            for i in motif:
                seedAlign += '-'
                aligned += ' '
                motifAlign += i
        writeMe += '<td>'+str(start)+'</td><td>'+str(start+lenMatch-1)+'</td><td>'+str(lenMatch)+'</td><td><font face="Courier New"><pre>Seed  '+str(seedAlign)+'\n      '+str(aligned)+'\nMotif '+str(motifAlign)+'</pre></font></td>'
        # P(Alignment)
        writeMe += '<td>'+str(outMe[4][0])+'</td><td>'+str(float(len([i for i in self.totPs if i>=outMe[4][0]]))/float(len(self.totPs)))+'</td><td>'+str(outMe[4][2])+'</td><td>'+str(float(len([i for i in self.vitPs if i>=outMe[4][2]]))/float(len(self.vitPs)))+'</td></tr>'

        return writeMe

    # Decide whether to add and add at correct position
    # Strucuture of elements:
    # input = [pssm, getConsensus(pssms[pssm]), miRNA, miRNAs[miRNA], forVit], totPs, vitPs (Ps are permuted probabilities)
    # output = [miRNAname,miRNAseed,AlignStart,AlignStop,AlignLength,SeedAlign,Align,MotifAlign,P(Total),P-valueTotal,P(Viterbi),P-valueViterbi]
    def outData(self,outMe,totP,vitP,compModel,fullSeed):
        #print outMe[4][1]
        output = []
        # miRNA name
        output += [outMe[2],self.reverseComplement(fullSeed)]
        alignment = outMe[4][1] # Grab the alignment
        alignment.pop() # Get rid of the extra state which is added by the forwardViterbi function
        start = 1
        if alignment[0]=='NM1':
            for i in alignment:
                if i=='NM1':
                    start += 1
        # Alignment
        seedAlign = ''
        seed = list(self.reverseMe(outMe[3]))
        motifAlign = ''
        motif = list(outMe[1])
        alignMe = alignment
        aligned = ''
        lenMatch = 0
        # First get them zero'd to the same point
        if start>1:
            for i in range(start-1):
                seedAlign += seed.pop(0)
                aligned += '_'
                motifAlign += '-'
                alignMe.pop(0)
        if len(alignMe)>0 and not alignMe[0]=='PSSM0' and not alignMe[0]=='WOBBLE0':
            if alignMe[0][0]=='P':
                upTo = int(alignMe[0][4])
            elif alignMe[0][0]=='W':
                upTo = int(alignMe[0][6])
            for i in range(upTo):
                seedAlign += '-'
                aligned += '_'
                motifAlign += motif.pop(0)
        # Then add the parts that align
        while 1:
            if len(alignMe)==0 or alignMe[0]=='NM2':
                break
            seedAlign += seed.pop(0)
            if alignMe[0][0]=='P':
                aligned += '|'
            elif alignMe[0][0]=='W':
                aligned += ':'
            lenMatch += 1
            motifAlign += motif.pop(0)
            alignMe.pop(0)
        # Then do the ending
        if len(alignMe)>0:
            for i in alignMe:
                seedAlign += seed[0]
                seed = seed[1:]
                aligned += '_'
                motifAlign += '-'
            alignMe = []
        if len(motif)>0 and len(alignMe)==0:
            for i in motif:
                seedAlign += '-'
                aligned += '_'
                motifAlign += i
        output += [start, start + lenMatch - 1, lenMatch, compModel,
                   "'" + motifAlign, "'" + aligned, "'" + seedAlign]
        # P(Alignment)
        output += [outMe[4][0], totP, outMe[4][2], vitP]
        return [str(i) for i in output]
