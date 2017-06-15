#!/usr/bin/env python
#################################################################
# @Program: sygnal.py                                           #
# @Version: 2                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  6/14/2016                      #
#################################################################

import math
import cPickle, os, re, sys
from sys import stdout, exit
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
import subprocess
from shutil import rmtree
from copy import deepcopy
from collections import defaultdict
import gzip

import rpy2.robjects as robj
from rpy2.robjects import FloatVector, IntVector, StrVector
from rpy2 import rinterface

from scipy.stats import pearsonr
from scipy.stats import ttest_ind
import numpy as np

# Custom offYerBack libraries
from cMonkeyWrapper import cMonkeyWrapper
from pssm import PSSM
import pssm as pssm_mod
from miRvestigator import miRvestigator
import utils
import config
import json

#################################################################
## rpy2 integration                                            ##
#################################################################


def make_rint_vector(a):
    return IntVector(map(lambda i: rinterface.NA_Integer if i == 'NA' else int(i), a))


def make_rfloat_vector(a):
    return FloatVector(map(lambda f: rinterface.NA_Real if f == 'NA' else float(f), a))


#################################################################
## Parameters                                                  ##
#################################################################

# global state for concurrent operations
g_weeder_args        = None
g_weeder_results     = None
g_meme_args          = None
g_cluster_meme_runs  = None
g_pred_dict          = None
g_pred_total_targets = None
g_biclusters         = None
RATIOS               = None

#################################################################
## Functions                                                   ##
#################################################################


def meme(num, seqfile, bgfile, nmotifs, min_motif_width, max_motif_width, revcomp, seed=None):
    """ Run meme and get the output into PSSMs"""
    global g_cluster_meme_runs

    # Arguments for meme
    args = '%s -bfile %s -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw %d -maxw %d -nmotifs %d' % (seqfile, bgfile, min_motif_width, max_motif_width, nmotifs)

    if revcomp:
        args += ' -revcomp'

    if not seed is None:
        args += ' -cons ' + str(seed)

    print "MEME args: '%s'" % args
    meme_proc = Popen("meme %s" % args, shell=True, stdout=PIPE)
    output = meme_proc.communicate()[0].split('\n')

    PSSMs = []

    for i in range(len(output)):
        desc_comps = output[i].strip().split(' ')

        if desc_comps[0] == 'Motif' and desc_comps[2] == 'position-specific' and desc_comps[3] == 'probability':
            i += 2  # Skip the separator line, go to the summary line
            summary_comps = output[i].strip().split(' ')
            width = int(summary_comps[5])
            sites = int(summary_comps[7])
            evalue = float(summary_comps[9])
            matrix = []

            for j in range(width):
                i += 1
                matrix += [[float(let) for let in output[i].strip().split(' ') if let]]
            PSSMs.append(PSSM('%s_motif%s_meme' % (os.path.basename(seqfile).split('_')[1].split('.')[0], desc_comps[1]),
                              sites, evalue, matrix, [], 'meme'))
    g_cluster_meme_runs[num] = PSSMs


def run_meme(runarg):
    """Wrapper function to run the meme function using a multiprocessing pool
    """
    global g_meme_args

    run_num, filepath = runarg
    meme(run_num, filepath, g_meme_args['bgfile'], g_meme_args['nmotifs'], g_meme_args['min_motif_width'],
         g_meme_args['max_motif_width'], g_meme_args['revcomp'])


def weeder(bicluster, seqfile, bgfile, size, enriched, revcomp):
    """
    Run weeder and parse its output
    First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
    """
    global g_weeder_results
    print "run weeder on '%s'" % seqfile

    # First run weederTFBS
    weeder_args = "%s %s %s %s" % (seqfile, bgfile, size, enriched)
    if revcomp:
        weeder_args += ' S'
    errout = open(cfg['tmpdir']+'/weeder/stderr.out','w')
    weeder_proc = Popen("weederlauncher %s" % weeder_args, shell=True, stdout=PIPE, stderr=errout)
    output = weeder_proc.communicate()

    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqfile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads will be
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    if size=='medium':
        # Scroll to where the 10bp reads wll be
        while 1:
            if len(outLines)<=1:
                break
            outLine = outLines.pop(0)
            if not outLine.find('Searching for motifs of length 10') == -1:
                break

        # Get top hit of 10bp look for "1)"
        while 1:
            if len(outLines)<=1:
                break
            outLine = outLines.pop(0)
            if not outLine.find('1) ') == -1:
                break
        hitBp[10] = outLine.strip().split(' ')[1:]

    # Scroll to where the 10bp reads will be
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break

    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        if len(outLines)<=1:
            break
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    motif_id = 1
    bicluster_id = int(os.path.basename(seqfile).split('_')[1].split('.')[0])
    while 1:
        if len(outLines)<=1:
            break

        if revcomp:
            name = outLines.pop(0).strip() # Get match
            name += '_'+outLines.pop(0).strip()
        else:
            name = outLines.pop(0).strip() # Get match

        if not name.find('(not highest-ranking)') == -1:
            break

        # Get redundant motifs
        outLines.pop(0)
        red_motifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [float(i.strip()) for i in col.strip().split('\t')[1].split(' ') if i]
            colsum = sum(nums)
            matrix += [[ nums[0] / colsum, nums[1] / colsum, nums[2] / colsum, nums[3] / colsum]]
            col = outLines.pop(0)
        PSSMs += [PSSM('%d_motif%d_weeder' % (bicluster_id, motif_id),
                       len(instances), hitBp[len(matrix)][1], matrix, red_motifs, 'weeder')]
        motif_id += 1
    g_weeder_results[bicluster] = PSSMs


def run_weeder(run_arg):
    global g_weeder_args

    run_num, filepath = run_arg
    weeder(run_num, filepath, g_weeder_args['bgfile'], g_weeder_args['size'],
           g_weeder_args['enriched'], g_weeder_args['revcomp'])


def tomtom(num, dist_meth='ed', q_thresh=1, min_overlap=6):
    args = '-dist %s -o tmp/tomtom_out -text -thresh %d -min-overlap %d -verbosity 1 tmp/query%d.tomtom tmp/target%d.tomtom' % (dist_meth, q_thresh, min_overlap, num, num)
    print args

    with open(cfg['tmpdir']+'/stderr_%d.out' % num,'w') as errout:
        tomtom_proc = Popen("tomtom %s" % args, shell=True, stdout=PIPE, stderr=errout)
        with open(cfg['tmpdir']+'/tomtom_out/tomtom%d.out' % num, 'w') as outfile:
            output = tomtom_proc.communicate()[0]
            outfile.write(output)


def run_tomtom(i):
    tomtom(i, 'ed', 1, 6)


def phyper(q, m, n, k, lower_tail=False):
    """calls the R function phyper, input values are lists and returns a list"""
    r_phyper = robj.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return [f for f in
            r_phyper(FloatVector(q), FloatVector(m), FloatVector(n), FloatVector(k), **kwargs)]


def correlation(a1, a2):
    """
    Calculate the correlation coefficient and p-value between two variables.
    Input: Two arrays of float or integers.
    Returns: Corrleation coefficient and p-value.
    """
    """
    cor_test = robj.r['cor.test']
    result = cor_test(make_rfloat_vector(a1), make_rfloat_vector(a2))
    res = dict(zip(result.names, list(result)))
    tmp_e = 'NA'
    if not str(res['estimate'][0])=='NA':
        tmp_e = float(res['estimate'][0])
    tmp_p = 'NA'
    if not str(res['p.value'][0])=='NA':
        tmp_p = float(res['p.value'][0])
    return tmp_e, tmp_p
    """
    tmp = [i for i in range(len(a1)) if not (a1[i]=='NA' or a2[i]=='NA')]
    a1 = [float(a1[i]) for i in tmp]
    a2 = [float(a2[i]) for i in tmp]
    res = pearsonr(a1, a2)
    if not np.isnan(res[0]):
        return res[0], res[1]
    else:
        return 'NA', 'NA'


def ttest(a1, a2):
    """
    Calculate the T-test statistic and p-value testing for a significant
    difference between between two groups.
    https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.ttest_ind.html
    Input: Two arrays of float or integers.
    Returns: T statistic and p-value.
    """
    tmp1 = [i for i in range(len(a1)) if not a1[i]=='NA']
    a1 = [float(a1[i]) for i in tmp1]
    tmp2 = [i for i in range(len(a2)) if not a2[i]=='NA']
    a2 = [float(a2[i]) for i in tmp2]
    res = ttest_ind(a1, a2)
    if not np.isnan(res[0]):
        return res[0], res[1]
    else:
        return 'NA', 'NA'

def anova_oneway(a1, a2, a3, a4):
    """
    Perform one-way anova test, tests null hypothesis that 2+ groups have same population mean.
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.f_oneway.html
    Input: multiple arrays of float or integers
    Returns: F-value and pvalue associated with F-value
    """
    tmp1 = [i for i in range(len(a1)) if not a1[i] == 'NA']
    a1   = [float(a1[i]) for i in tmp1]
    tmp2 = [i for i in range(len(a2)) if not a2[i] == 'NA']
    a2   = [float(a2[i]) for i in tmp2]
    tmp3 = [i for i in range(len(a3)) if not a3[i] == 'NA']
    a3   = [float(a1[i]) for i in tmp1]
    tmp4 = [i for i in range(len(a4)) if not a4[i] == 'NA']
    a4   = [float(a4[i]) for i in tmp1]
    res  = f_oneway(a1, a2, a3, a4)
    if not np.isnan(res[0]):
        return res[0], res[1]
    else: 
	return 'NA', 'NA'


def mean_fold_change(a1, a2):
    """
    Calculate fold change between mean expression values of  two groups.
    Input: Two arrays of float or integers, assumes that values are log2-transformed
    Returns: fold change 2^(a1-a2) (assumes input is log2-transformed data).
    """
    tmp1 = [i for i in range(len(a1)) if not a1[i]=='NA']
    a1 = [float(a1[i]) for i in tmp1]
    tmp2 = [i for i in range(len(a2)) if not a2[i]=='NA']
    a2 = [float(a2[i]) for i in tmp2]
    a1_mean = np.mean(a1)
    a2_mean = np.mean(a2)    
    log_fc  = a1_mean - a2_mean
    res = np.power(2, log_fc) 
    if not np.isnan(res):
        return res
    else:
        return 'NA'


def survival(survival, dead, pc1, age, sex):
    """
    Calculate the survival correlation coefficient and p-value between two variables.
    Input: Four arrays of float or integers.
    Returns:
    """
    surv = robj.r("""
    library('survival')
    surv <- function(s, dead, pc1, age, sex) {
        scph1 = summary(coxph(Surv(s,dead == 1) ~ pc1))
        scph2 = summary(coxph(Surv(s,dead == 1) ~ pc1 + age))
        scph3 = summary(coxph(Surv(s,dead == 1) ~ pc1 + age + sex))
        c(scph1$coef[1,4], scph1$coef[1,5], scph2$coef[1,4], scph2$coef[1,5], scph3$coef[1,4], scph3$coef[1,5])
    }
    """)
    try:
        res = surv(make_rfloat_vector(survival), StrVector(dead), make_rfloat_vector(pc1), make_rfloat_vector(age), make_rint_vector(sex))
        for i in range(len(res)):
            if np.isnan(res[i]):
                res[i] = 'NA'
        return [[res[0], res[1]], [res[2], res[3]], [res[4], res[5]]]
    except:
        return [['NA','NA'],['NA','NA'],['NA','NA']]

def compareMiRNANames(a, b):
    """
    Match miRNA names.
    Input: Two miRNA names to be compared.
    Returns:
    """
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
            re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(b+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(a):
            return 1
    return 0


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_synonyms(cfg):
    """ Load synonym thesaurus to get a mapping from entrez id to a list
    of UCSC IDs.

    The synonym file is assumed to be in the format
    <UCSC ID>,<alt1>;<alt2>;...
    and the Entrez ID is assumed to be the first numerical value in the
    alternatives."""
    entrez2id = defaultdict(list)

    with gzip.open(cfg['synonyms-file'], 'r') as infile:
        for line in infile:
            row = line.strip().split(',')
            id = row[0]
            entrez = [i for i in row[1].split(';') if is_number(i)]
            if len(entrez) == 1:
                entrez2id[entrez[0]].append(id)

    return entrez2id


def miRNA_mappings(cfg):
    """Create a dictionary to convert the miRNAs to there respective ids"""
    mirna_ids = {}
    mirna_ids_rev = {}

    with open(cfg['mirna-fasta-file'], 'r') as infile:
        for line in infile:
            splitUp = line.split(' ')
            if not splitUp[1] in mirna_ids_rev:
                mirna_ids_rev[splitUp[1]] = splitUp[0].lower()

            if not splitUp[0].lower() in mirna_ids:
                mirna_ids[splitUp[0].lower()] = splitUp[1]
            else:
                print 'Uh oh!', splitUp
    return mirna_ids, mirna_ids_rev


def read_cmonkey_run(cfg, gene_conv):
    """ Load cMonkey Object - turns cMonkey data into objects"""
    output_path = cfg.outdir_path('c1.pkl')

    # If this is the first time then load from the RData file
    if not os.path.exists(output_path):
        c1 = cMonkeyWrapper(cfg['cmonkey-rundb'], meme_upstream=False, weeder_upstream=False,
                            weeder_3pUTR=False, tfbs_db=False,
                            pita_3pUTR=False,
                            targetscan_3pUTR=False,
                            promoterSeq=cfg['promoterSeq'],
                            p3utrSeq=cfg['p3utrSeq'],
                            geneConv=False)

        with open(output_path, 'wb') as pklFile:
            cPickle.dump(c1, pklFile)

    # Otherwise load the dumped pickle file if it exists
    else:
        print 'Loading precached cMonkey object...'
        with open(output_path, 'rb') as pklFile:
            c1 = cPickle.load(pklFile)
        print 'Done.\n'

    return c1


def compute_upstream_motifs_meme(cfg, c1):
    """
    A. Upstream motifs (MEME) ##
    If MEME hasn't been run on the biclusters upstream sequences then do so
    """
    global g_meme_args, g_cluster_meme_runs

    if not c1.meme_upstream:
        print 'Running MEME on Upstreams:'
        pkl_path = cfg.outdir_path('meme_upstream.pkl')
        # Use already run MEME results if they exist
        if not os.path.exists(pkl_path):
            # Make needed directories
            cfg.clear_tmp()
            cfg.create_tmpdir('meme/fasta')

            # Run MEME on all biclusters
            mgr = Manager()
            run_args = []

            # First make fasta files for all biclusters
            print 'Making Files for MEME Upstream run...'
            for cluster_num in c1.biclusters:
                seqs = c1.bicluster_seqs_upstream(cluster_num)
                if len(seqs) > 0:
                    cluster_filename = cfg.tmpdir_path('meme/fasta/bicluster_%d.fasta' % cluster_num)
                    run_args.append((cluster_num, cluster_filename))
                    with open(cluster_filename, 'w') as outfile:
                        outfile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

            # Where all the results will be stored
            g_cluster_meme_runs = mgr.dict()

            # Parameters to use for running MEME
            g_meme_args = mgr.dict({'bgfile': cfg['meme']['upstream']['bgfile'],
                                    'nmotifs': cfg['meme']['upstream']['nmotifs'],
                                    'min_motif_width': cfg['meme']['upstream']['motif-widths'][0],
                                    'max_motif_width': cfg['meme']['upstream']['motif-widths'][1],
                                    'revcomp': bool(cfg['meme']['upstream']['revcomp'])})


            print 'Running MEME on Upstream sequences...'
            print 'There are ',cfg['cores'],' CPUs available.'
            pool = Pool(processes=int(cfg['cores']))
            pool.map(run_meme, run_args)
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            with open(pkl_path, 'wb') as outfile:
                cPickle.dump(deepcopy(g_cluster_meme_runs), outfile)
        else:
            print 'Loading from precached object...'
            with open(pkl_path, 'rb') as outfile:
                g_cluster_meme_runs = cPickle.load(outfile)

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for cluster_num, pssms in g_cluster_meme_runs.items():
            for pssm1 in pssms:
                bicluster = c1.biclusters[cluster_num]
                pssm1.de_novo_method = 'meme'
                bicluster.add_pssm_upstream(pssm1)

        print 'Done with MEMEing.\n'

        # MEME upstream has been run on cMonkey run
        c1.meme_upstream = True


def __compute_motifs_weeder(cfg, pkl_path, biclusters, add_result, bicluster_seqs, args_dict):
    """Generic motif detection function with weeder
    """
    global g_weeder_args, g_weeder_results

    if not os.path.exists(pkl_path):
        cfg.clear_tmp()
        cfg.create_tmpdir('weeder/fasta')

        # Run WEEDER on all biclusters
        mgr = Manager()
        run_args = []

        # First make fasta files for all biclusters
        for cluster_num in biclusters:
            seqs = bicluster_seqs(cluster_num)
            if len(seqs) > 0:
                cluster_filename = cfg.tmpdir_path('weeder/fasta/bicluster_%d.fasta' % cluster_num)
                run_args.append((cluster_num, cluster_filename))
                with open(cluster_filename, 'w') as outfile:
                    outfile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

        # Where all the results will be stored
        g_weeder_results = mgr.dict()

        # Parameters to use for running Weeder
        # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
        g_weeder_args = mgr.dict(args_dict)

        print 'Running Weeder...'
        print 'There are ', cfg['cores'],' CPUs available.'
        pool = Pool(processes=int(cfg['cores']))
        pool.map(run_weeder, run_args)
        pool.close()
        pool.join()

        # Dump weeder results as a pickle file
        with open(pkl_path,'wb') as outfile:
            cPickle.dump(deepcopy(g_weeder_results), outfile)

    else:
        print 'Loading from precached object...'
        with open(pkl_path,'rb') as infile:
            g_weeder_results = cPickle.load(infile)

    print 'Storing output...'
    for i, pssms in g_weeder_results.items():
        for p in pssms:
            p.de_novo_method = 'weeder'
            add_result(i, p)

    print 'Done with WEEDERing.\n'


def compute_upstream_motifs_weeder(cfg, c1):
    print '\nRunning WEEDER on Upstreams:'
    if not c1.weeder_upstream:
        __compute_motifs_weeder(cfg, cfg.outdir_path('weeder_upstream.pkl'),
                                c1.biclusters,
                                lambda bi, p: c1.biclusters[bi].add_pssm_upstream(p),
                                lambda bi: c1.bicluster_seqs_upstream(bi),
                                {'bgfile': cfg['weeder']['upstream']['bgfile'],
                                 'size': cfg['weeder']['upstream']['size'],
                                 'enriched': cfg['weeder']['upstream']['enriched'],
                                 'revcomp': bool(cfg['weeder']['upstream']['revcomp'])})
        c1.weeder_upstream = True


def compute_3pUTR_weeder(cfg, c1):
    print '\nRunning WEEDR on 3\' UTRs:'
    if not c1.weeder_3pUTR:
        __compute_motifs_weeder(cfg, cfg.outdir_path('weeder_3pUTR.pkl'),
                                c1.biclusters,
                                lambda bi, p: c1.biclusters[bi].add_pssm_3putr(p),
                                lambda bi: c1.bicluster_seqs_3putr(bi),
                                {'bgfile': cfg['weeder']['3pUTR']['bgfile'],
                                 'size': cfg['weeder']['3pUTR']['size'],
                                 'enriched': cfg['weeder']['3pUTR']['enriched'],
                                 'revcomp': bool(cfg['weeder']['3pUTR']['revcomp'])})
        c1.weeder_3pUTR = True


def cluster_hypergeo(bicluster_id):
    """concurrent computation
    k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
    Take gene list and compute overlap with each miRNA
    """
    global g_pred_dict, g_pred_total_targets, g_biclusters
    db = g_pred_dict
    all_genes = g_pred_total_targets[0]

    genes = all_genes.intersection(g_biclusters[bicluster_id].genes)
    m1s = []
    q = []
    m = []
    n = []
    k = []
    for m1, m1_genes in db.items():
        m1s.append(m1)
        miRNAGenes = all_genes.intersection(m1_genes)
        q.append(len(set(miRNAGenes).intersection(genes)))
        m.append(len(miRNAGenes))
        n.append(len(all_genes) - len(miRNAGenes))
        k.append(len(genes))

    results = phyper(q,m,n,k)
    min_mirna = []
    perc_targets = []
    min_pvalue = 1.0

    for i in range(1, len(results)):
        if float(results[i]) <= 0.05 / 674.0 and q[i] != 0 and float(q[i]) / float(k[i]) >= 0.1:
            if min_mirna == [] or float(results[i]) < min_pvalue:
                min_mirna = [i]
                perc_targets = [ float(q[i]) / float(k[i]) ]
                min_pvalue = float(results[i])
            elif float(results[i]) == min_pvalue:
                min_mirna.append(i)
                perc_targets.append(float(q[i]) / float(k[i]))

    print 'Bicluster #%d' % bicluster_id, ' '.join([m1s[miRNA] for miRNA in min_mirna])
    return [bicluster_id, ' '.join([m1s[i] for i in min_mirna]),
            ' '.join(map(str, perc_targets)), min_pvalue]


def __bicluster_genes(c1):
    result = set()
    for key, bicluster in c1.biclusters.items():
        for gene in bicluster.genes:
            result.add(gene)
    return result


def __read_predictions(pred_path, pkl_path, genes_in_biclusters, manager):
    if not os.path.exists(pkl_path):
        print 'loading predictions...'
        tmp_set = set()
        tmp_dict = {}
        inFile = open(pred_path,'r')
        tmp_dict = json.load(inFile)
        inFile.close()
        tmp_set = set([item for sublist in tmp_dict.values() for item in sublist])
        with open(pkl_path,'wb') as pklFile:
            cPickle.dump(tmp_dict, pklFile)
            cPickle.dump(tmp_set, pklFile)

    # Otherwise load the dumped pickle file if it exists
    else:
        print 'Loading pickled predictions...'
        with open(pkl_path,'rb') as pklFile:
            tmp_dict = cPickle.load(pklFile)
            tmp_set = cPickle.load(pklFile)

    pred_dict = manager.dict(tmp_dict)
    pred_total_targets = manager.list()
    pred_total_targets.append(set(tmp_set))
    return pred_dict, pred_total_targets


def __compute_enrichment(c1, name,  pkl_path, pred_path, pred_pkl_path):
    """General enrichment analysis function"""
    global g_pred_dict, g_pred_total_targets, g_biclusters

    print 'Running %s enrichment on Biclusters:' % name

    if not os.path.exists(pkl_path):
        print 'Get a list of all genes in run...'
        mgr = Manager()
        g_biclusters = mgr.dict(c1.biclusters)
        g_pred_dict, g_pred_total_targets = __read_predictions(pred_path, pred_pkl_path,
                                                               __bicluster_genes(c1), mgr)
        print '%s prediction has %d TFs.' % (name, len(g_pred_dict))

        print 'Running %s enrichment analyses...' % name
        pool = Pool(processes=int(cfg['cores']))
        result = pool.map(cluster_hypergeo, c1.biclusters.keys())
        pool.close()
        pool.join()

        with open(pkl_path, 'wb') as pklFile:
            cPickle.dump(result, pklFile)
    else:
        print 'Loading precached analysis...'
        with open(pkl_path, 'rb') as pklFile:
            result = cPickle.load(pklFile)
    return result


def compute_tfbsdb_enrichment(cfg, c1):
    """D. Upstream TFBS DB enrichment Analysis"""
    if not c1.tfbs_db:
        res1 = __compute_enrichment(c1, 'TFBS_DB', cfg.outdir_path('tfbs_db.pkl'),
                                    cfg['tfbsdb'],
                                    'TF/tfbs_db.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, tf(s), Percent Targets, P-Value]
            bicluster = c1.biclusters[r1[0]]
            bicluster.add_attribute('tfbs_db', {'tf':r1[1], 'percentTargets': r1[2],
                                                'pValue':r1[3]})
        print 'Done.\n'
        c1.tfbs_db = True


def compute_3pUTR_pita_set_enrichment(cfg, c1, mirna_ids):
    """E. 3' UTR PITA"""
    if not c1.pita_3pUTR:
        res1 = __compute_enrichment(c1, 'PITA', cfg.outdir_path('pita_3pUTR.pkl'),
                                    cfg['pita'],
                                    'miRNA/pita.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            bicluster = c1.biclusters[r1[0]]
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += utils.mirna_in_dict(m1.lower(), mirna_ids)
            bicluster.add_attribute('pita_3pUTR', {'miRNA': r1[1],
                                                   'percentTargets': r1[2],
                                                   'pValue': r1[3],
                                                   'mature_seq_ids': miRNA_mature_seq_ids})
        print 'Done.\n'
        c1.pita_3pUTR = True


def compute_3pUTR_targetscan_set_enrichment(cfg, c1, mirna_ids):
    """F. 3' UTR TargetScan"""
    if not c1.targetscan_3pUTR:
        res1 = __compute_enrichment(c1, 'TargetScan', cfg.outdir_path('targetscan_3pUTR.pkl'),
                                    cfg['targetscan'],
                                    'miRNA/targetScan.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            bicluster = c1.biclusters[r1[0]]
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += utils.mirna_in_dict(m1.lower(), mirna_ids)
            bicluster.add_attribute('targetscan_3pUTR',
                                    {'miRNA': r1[1], 'percentTargets': r1[2],
                                     'pValue':r1[3], 'mature_seq_ids': miRNA_mature_seq_ids})
        print 'Done.\n'
        c1.targetscan_3pUTR = True


def compute_additional_info(cfg, mirna_ids, gene_conv):
    cm_pkl_path = cfg.outdir_path('c1_all.pkl')
    if not os.path.exists(cm_pkl_path):
        c1 = read_cmonkey_run(cfg, gene_conv)

        ################################################################
        ## Fill in the missing parts                                  ##
        ################################################################
        #  Check to see if all parts are there:                        #
        #   A. Upstream motifs (MEME)                                  #
        #   B. Upstream motif (Weeder)                                 #
        #   C. 3' UTR Weeder-miRvestigator (Weeder)                    #
        #   D. TFBS DB Enrichment                                      #
        #   E. 3' UTR PITA (Set Enrichment)                            #
        #   F. 3' UTR TargetScan (Set Enrichment)                      #
        ################################################################

        compute_upstream_motifs_meme(cfg, c1)
        compute_upstream_motifs_weeder(cfg, c1)
        compute_3pUTR_weeder(cfg, c1)
        compute_tfbsdb_enrichment(cfg, c1)
        compute_3pUTR_pita_set_enrichment(cfg, c1, mirna_ids)
        compute_3pUTR_targetscan_set_enrichment(cfg, c1, mirna_ids)

        with open(cm_pkl_path, 'wb') as outfile:
            cPickle.dump(c1, outfile)
    else:
        print 'Loading prechached cMonkey Object (c1_all.pkl):'
        with open(cm_pkl_path, 'rb') as infile:
            c1 = cPickle.load(infile)
    return c1

def convertSex(x):
    if x=='NA':
        return 'NA'
    if x=='M':
        return 0
    if x=='F':
        return 1

PHENOTYPES = ['preop_treatment','sex']

def post_process(arg):
    global RATIOS
    cluster_num, bicluster, phenotypes = arg
    attributes = {}
    print ' Postprocessing cluster:', cluster_num

    attributes['k'] = cluster_num
    attributes['k.rows'] = bicluster.num_genes()
    attributes['k.cols'] = bicluster.num_conditions()

    # Get matrix of expression for genes
    genes = bicluster.genes

    # Get first principal component variance explained
    key0 = RATIOS.keys()[0]
    gene0_conditions = set(RATIOS[key0].keys())
    dtd_conditions = set(phenotypes['survival'].keys())
    g0_dtd_conds = gene0_conditions.intersection(dtd_conditions)
    fpc = bicluster.attributes['pc1']
    fpc_conds = fpc.keys()
    result_conds = g0_dtd_conds.intersection(fpc_conds)
    pc1_1 = [fpc[i] for i in result_conds]

    # Phenotype tests with patient traits
    for phenotype in PHENOTYPES:
        if phenotype=='sex':
            case = [fpc[i] for i in result_conds if phenotypes[phenotype][i]=='M']
            control = [fpc[i] for i in result_conds if phenotypes[phenotype][i]=='F']
        else:
            case = [fpc[i] for i in result_conds if phenotypes[phenotype][i]=='1']
            control = [fpc[i] for i in result_conds if phenotypes[phenotype][i]=='0']

        attr_name = phenotype
        t_stat, t_pvalue = ttest(case, control)
        attributes[attr_name] = {'T': t_stat, 'pValue': t_pvalue }

    # Association of bicluster expression with patient survival
    surv = [phenotypes['survival'][c] for c in result_conds]
    dead = [phenotypes['status'][c] for c in result_conds]
    age = [phenotypes['age_at_surgery'][c] for c in result_conds]
    sex = [convertSex(phenotypes['sex'][c]) for c in result_conds]
    surv_0, surv_age, surv_ageSex = survival(surv, dead, pc1_1, age, sex)
    attributes['Survival'] = {'z': surv_0[0], 'pValue': surv_0[1] }
    attributes['Survival.AGE'] = {'z': surv_age[0], 'pValue': surv_age[1]}
    attributes['Survival.AGE_SEX'] = {'z': surv_ageSex[0], 'pValue': surv_ageSex[1]}
    return attributes

def __read_ratios(cfg):
    print "reading ratios matrix"
    ratios = {}
    with open(cfg['all-ratios-file'], 'r') as infile:
        conditions = [i.strip('"') for i in infile.readline().strip().split(',')]
        for line in infile:
            comps = line.strip().split(',')
            ratios[comps[0].strip('"')] = dict(zip(conditions, comps[1:]))
    return ratios

def dump_cluster_members(cfg, c1):
    print "dump cluster row members"

    with open(cfg.outdir_path('cluster.members.genes.txt'), 'w') as outfile:
        for cluster_num in c1.biclusters:
            outfile.write('%s %s\n' % (cluster_num, ' '.join(c1.biclusters[cluster_num].genes)))

    print "dump cluster condition members"
    with open(cfg.outdir_path('cluster.members.conditions.txt'), 'w') as outfile:
        for cluster_num in c1.biclusters:
            outfile.write('%s %s\n' % (cluster_num, ' '.join(c1.biclusters[cluster_num].conditions)))


def __get_cluster_eigengenes(cfg, c1):
    # Calculate bicluster eigengene (first principal components)
    print "Compute bicluster eigengenes"
    cluster_eigengenes_path = cfg.outdir_path('biclusterEigengenes.csv')
    if not os.path.exists(cluster_eigengenes_path):
        ret = subprocess.check_call(['./getEigengene.R',
                                     '-r', cfg['all-ratios-file'],
                                     '-o', cfg['outdir'],
                                     '-c', str(cfg['cores'])],
                                    stderr=subprocess.STDOUT)
        if ret == 1:
            print "could not create Eigengenes"
            exit(1)

    # Read in bicluster eigengene
    print "read bicluster eigengenes"
    with open(cluster_eigengenes_path, 'r') as infile:
        patients = [i.strip('"') for i in infile.readline().strip().split(',')]
        patients.pop(0) # Get rid of rowname placeholder

        for line in infile:
            eigengene = line.strip().split(',')
            cluster_num = int(eigengene.pop(0).strip('"'))
            bicluster = c1.biclusters[cluster_num]
            bicluster.add_attribute('pc1', dict(zip(patients, eigengene)))


def __get_cluster_variance_explained(cfg, c1):
    print "read bicluster variance explained"
    with open(cfg.outdir_path('biclusterVarianceExplained.csv'), 'r') as infile:
        infile.readline() # Get rid of header
        for line in infile:
            varExplained = line.strip().split(',')
            cluster_num = int(varExplained.pop(0).strip('"'))
            bicluster = c1.biclusters[cluster_num]
            bicluster.add_attribute('pc1.var.exp', varExplained[0])


def __get_phenotype_info(cfg, c1):
    print "read phenotype information"
    phenotypes = {}
    with open(cfg['phenotypes-file'], 'r') as infile:
        ids = infile.readline().strip().split(',')[1:]
        phenotypes = {}
        for i in ids:
            phenotypes[i] = {}

        for line in infile:
            splitUp = line.strip().split(',')
            for i in range(len(ids)):
                phenotypes[ids[i]][splitUp[0]] = splitUp[i+1]
    return phenotypes


def __do_postprocess(cfg, c1, phenotypes):
    global RATIOS

    multi = False
    postprocess_pkl_path = cfg.outdir_path('postProcessed.pkl')
    if not os.path.exists(postprocess_pkl_path):
        RATIOS = __read_ratios(cfg)
        args = [[bicluster, c1.biclusters[bicluster], phenotypes] for bicluster in c1.biclusters]
        if multi:
            pool = Pool()
            res1 = pool.map(post_process, args)
            pool.close()
            pool.join()
        else:
            res1 = [post_process(arg) for arg in args]

        # Dump res1 into a pkl
        with open(postprocess_pkl_path, 'wb') as outfile:
            cPickle.dump(res1, outfile)
    else:
        with open(postprocess_pkl_path, 'rb') as infile:
            res1 = cPickle.load(infile)

    # Put results in cMonkey object
    for entry in res1:
        bicluster = c1.biclusters[entry['k']]
        for attribute in entry:
            if not attribute == 'k':
                bicluster.add_attribute(attribute, entry[attribute])
    print 'Done.\n'


def __tomtom_upstream_motifs(cfg):
    #################################################################
    ## TomTom Upstream motifs versus Jaspar and Transfac           ##
    #################################################################
    print 'Running TOMTOM on Upstream Motifs:'

    # Make needed directories
    cfg.clear_tmp()
    cfg.create_tmpdir('tomtom_out')

    pssms = c1.pssms_upstream()
    upstreamMatches = {}

    comparison_pkl_path = cfg.outdir_path('upstreamJasparTransfacComparison.pkl')
    comparison_csv_path = cfg.outdir_path('upstreamComparison_jaspar_transfac.csv')

    target_pssms_in = []
    for motif_file in cfg['tomtom']['upstream']['motif-files']:
        tmpPssms = pssm_mod.load_pssms_json(motif_file)
        for pssm in tmpPssms.values():
            pssm.de_novo_method = 'meme'
        target_pssms_in.append(tmpPssms)

    if not os.path.exists(comparison_pkl_path):

        # Write out results
        with open(comparison_csv_path, 'w') as outFile:
            outFile.write('Motif Name,Original E-Value,Consensus,JASPAR Motif,JASPAR Consensus,TomTom.pValue,TomTom.qValue,Probe In Bicluster,Bicluster Residual')

            # Making MEME formatted files (make_files function in utils)
            print 'Making files...'
            for i, target_pssms in enumerate(target_pssms_in):
                utils.make_files(c1.nucFreqsUpstream, pssms.values(),
                                 target_pssms.values(), i)

            # Run TomTom
            print 'Comparing Upstream motifs against databases...'
            pool = Pool(processes=int(cfg['cores']))
            res1 = pool.map(run_tomtom, [i for i in range(len(target_pssms_in))])
            pool.close()
            pool.join()

            print 'Reading in Tomtom run...'
            output_lines = []
            for i in range(len(target_pssms_in)):
                with open(cfg.tmpdir_path('tomtom_out/tomtom%d.out' % i), 'r') as tomtom_outfile:
                    tomtom_outfile.readline()  # skip header
                    output_lines += [line.strip().split('\t') for line in tomtom_outfile
                                    if float(line.split('\t')[5]) <= 0.05]

            # Now iterate through output and save data
            print 'Now parsing output for %d matches...' % len(output_lines)
            for outputLine in output_lines:
                if len(outputLine) == 10 and float(outputLine[5]) <= 0.05:
                    tfName = outputLine[1].split('_')[0]
                    if not outputLine[0] in upstreamMatches:
                        upstreamMatches[outputLine[0]] = [{'factor':outputLine[1],'confidence':outputLine[5]}]
                    else:
                        upstreamMatches[outputLine[0]].append({'factor':outputLine[1],'confidence':outputLine[5]})

        with open(comparison_pkl_path,'wb') as pklFile:
            cPickle.dump(upstreamMatches, pklFile)
    else:
        print 'Loading precached upstream matches...'
        with open(comparison_pkl_path, 'rb') as pklFile:
            upstreamMatches = cPickle.load(pklFile)

    # Pump into pssms
    for pssmName in upstreamMatches:
        for match in upstreamMatches[pssmName]:
            pssms[pssmName].add_match(factor=match['factor'], confidence=match['confidence'])
    print 'We matched '+str(len(upstreamMatches))+' upstream motifs.\n'


def __expand_tf_factor_list(entrez2id):
    #################################################################
    ## Prepare to get expanded set of TFs from TFClass             ##
    ## (http://tfclass.bioinf.med.uni-goettingen.de/)              ##
    #################################################################
    print 'Expanding TF factor list with TFClass families...'
    # Read in humanTFs_All.csv with <TF Name>,<Entrez ID>
    tfName2entrezId = {}
    with open(cfg['tfExpansion']['tfs'],'r') as inFile:
        inFile.readline() # Get rid of header
        for line in inFile:
            splitUp = line.strip().split(',')
            """
            if splitUp[2] in entrez2id:
                for i in entrez2id[splitUp[2]]:
                    tfName2entrezId[splitUp[0]] = i
            """
            tfName2entrezId[splitUp[0]] = splitUp[2]

    # Read in tfFamilies.csv for expanded list of possible TFs
    tfFamilies = {}
    with open(cfg['tfExpansion']['tfFamilies'],'r') as inFile:
        inFile.readline() # Get rid of header
        for line in inFile:
            splitUp = line.strip().split(',')
            """
            tmp = []
            for i in splitUp[2].split(' '):
                if i in entrez2id:
                    tmp += entrez2id[i]
            """
            tmp = splitUp[2].split(' ')
            tfFamilies[splitUp[0]] = tmp

    # Add expanded TF regulators
    pssms = c1.pssms_upstream()
    for pssm in pssms.values():
        expanded_factors = {}

        # Collapse matches to a set of entrez IDs and add expanded factors
        if len(pssm.matches) > 0:
            for match in pssm.matches:
                if match['factor'] in tfName2entrezId:
                    factor = tfName2entrezId[match['factor']]
                    if not factor in expanded_factors:
                        expanded_factors[factor] = [factor]
                        for family in tfFamilies:
                            if factor in tfFamilies[family]:
                                expanded_factors[factor] += tfFamilies[family]
            # Push expanded TF factor list into the PSSM object
            for factor in expanded_factors:
                for expanded_factor in list(set(expanded_factors[factor])):
                    pssm.add_expanded_match(expanded_factor, factor)
    print 'Finished expanding TF factor list.\n'
    return tfName2entrezId, tfFamilies


def __correlate_tfs_with_cluster_eigengenes(cfg, c1):

    #######################################################################
    ## Filter expanded TFs through correlation with bicluster eigengenes ##
    #######################################################################
    print 'Correlate TFs with eigengenes...'

    # Get microarray data
    exp_data = {}
    with open(cfg['all-ratios-file'], 'r') as infile:
        all_names = [i for i in map(lambda n: n.strip('"').replace('-','.'), infile.readline().strip().lstrip(',').split(',')) if not i=='']
        for line in infile:
            row = line.strip().split(',')
            exp_data[row[0].strip('"')] = dict(zip(all_names, row[1:]))

    # [rho, pValue] = correlation(a1,a2)
    for cluster_num, bicluster in c1.biclusters.items():
        for pssm in bicluster.pssms_upstream:
            compared = []

            # Get maximum amount of correlation positive or negative
            if len(pssm.expanded_matches) > 0:
                print cluster_num, pssm.name, len(pssm.expanded_matches)
                for factor in pssm.expanded_matches:
                    corMax = []
                    if factor['factor'] in exp_data.keys():
                        eigengene = bicluster.attributes['pc1']
                        if (not factor['factor'] in compared) and (not sum([1 for i in eigengene.values() if i=='NA'])==len(eigengene)) and (float(sum([1 for i in exp_data[factor['factor']].values() if i=='NA']))/float(len(exp_data[factor['factor']])))<=0.8:
                            compared.append(factor['factor'])
                            tmp_names = set(eigengene.keys()).intersection(all_names)
                            #print [eigengene[i] for i in tmp_names]
                            #print [exp_data[factor['factor']][i] for i in tmp_names]
                            cor1 = correlation([eigengene[i] for i in tmp_names], [exp_data[factor['factor']][i] for i in tmp_names])
                            print factor['factor'], cor1
                            if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                                corMax = cor1
                    if not corMax==[]:
                        pssm.add_correlated_match(factor['factor'],corMax[0],corMax[1])
    print 'Done.\n'
    return exp_data, all_names


def __expand_and_correlate_tfbsdb_tfs(c1, tf_name2entrezid, tf_families, exp_data,
                                      all_names):
    print 'Expand and correlate TFBS_DB TFs...'
    for bicluster in c1.biclusters.values():
        # Get the tfbs_db attribute and for each TF get the list of expanded factors
        tfs = bicluster.attributes['tfbs_db']
        expanded_factors = {}
        if not tfs is None:
            for tf in tfs['tf'].split(' '):
                if tf[0:2]=='V_':
                    tf = 'V$'+tf[2:]
                # Get the list of expanded factors
                if tf in tf_name2entrezid:
                    factor = tf_name2entrezid[tf]
                    if not factor in expanded_factors:
                        expanded_factors[factor] = [factor]
                        for family_factors in tf_families.values():
                            if factor in family_factors:
                                expanded_factors[factor] += family_factors
                        expanded_factors[factor] = list(set(expanded_factors[factor]))

        # Push expanded TF factor list into the bicluster object
        if len(expanded_factors) > 0:
            print factor, expanded_factors
        bicluster.add_attribute('tfbs_db_expanded', expanded_factors)

    # [rho, pValue] = correlation(a1,a2)
    for bicluster in c1.biclusters.values():
        factors = bicluster.attributes['tfbs_db_expanded']
        compared = []
        correlatedFactor = []

        # Get maximum amount of correlation positive or negative
        for factor1 in factors:
            for factor2 in factors[factor1]:
                corMax = []
                if factor2 in exp_data:
                    eigengene = bicluster.attributes['pc1']
                    if (not factor2 in compared) and (not sum([1 for i in eigengene.values() if i=='NA'])==len(eigengene)) and (float(sum([1 for i in exp_data[factor2].values() if i=='NA']))/float(len(exp_data[factor2])))<=0.8:
                        compared.append(factor2)
                        tmp_names = set(eigengene.keys()).intersection(all_names)
                        cor1 = correlation([eigengene[i] for i in tmp_names], [exp_data[factor2][i] for i in tmp_names])
                        print 'Expanded:', factor2, cor1
                        if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                            corMax = cor1
                    if not corMax==[]:
                        correlatedFactor.append({'factor':factor2,'rho':corMax[0],'pValue':corMax[1]})
    	bicluster.add_attribute('tfbs_db_correlated', correlatedFactor)
    print 'Done.\n'


def __write_first_principal_components(cfg, c1):
    print 'Write biclusterFirstPrincComponents.csv...'
    # Get all first principal components for each bicluster
    fpcWrite = []
    conditions = c1.biclusters[1].attributes['pc1'].keys()
    for i in sorted(c1.biclusters.keys()):
        pc1 = c1.biclusters[i].attributes['pc1']
        fpcWrite.append(str(i)+','+','.join([str(pc1[j]) for j in conditions]))

    # Dump out file
    with open(cfg.outdir_path(cfg["first-principal-comps-result-file"]), 'w') as outfile:
        outfile.write('Bicluster,'+','.join([j.strip() for j in conditions])+'\n')
        outfile.write('\n'.join(fpcWrite))
    print 'Done.\n'


def __get_permuted_pvalues_for_upstream_meme_motifs(cfg, c1):
    cfg.clear_tmp()
    cfg.create_tmpdir('tomtom_out')

    # Compare the random motifs to the original motif in TOMTOM
    permPValues = {}
    pssms = c1.pssms_upstream(de_novo_method='meme')
    output_path = cfg.outdir_path('upstreamMotifPermutedPValues.csv')
    if not os.path.exists(output_path):
        outFile = open(output_path, 'w')
        outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
        pssmsNames = pssms.keys()
        print pssmsNames
        print 'Loading precached random PSSMs...'
        randPssmsDict = {}
        for i in range(5, 70, 5):  # [5,10,15,... ,65]:
            stdout.write(str(i)+' ')
            stdout.flush()
            filepath = os.path.join(cfg['rand_pssms_dir'], 'pssms_upstream_%d.json' % i)
            randPssmsDict[i] = pssm_mod.load_pssms_json(filepath)

            for pssm1 in randPssmsDict[i]:
                randPssmsDict[i][pssm1].de_novo_method = 'meme'
            delMes = []
            for randPssm in randPssmsDict[i]:
                if not float(randPssmsDict[i][randPssm].evalue) <= cfg['max-evalue']:
                    delMes.append(randPssm)
            for j in delMes:
                del randPssmsDict[i][j]

        print '\nMaking files...'
        for i in range(len(pssms)):
            clustSize = utils.rand_pssm_clust_size(c1.biclusters[int(pssmsNames[i].split('_')[0])].num_genes())
            utils.make_files(c1.nucFreqsUpstream, [pssms[pssmsNames[i]]],
                             randPssmsDict[clustSize].values(), i)

        # Run this using all cores available
        cpus = int(cfg['cores'])
        print 'There are', cpus,'CPUs available.'
        print 'Running TOMTOM to compare PSSMs...'
        pool = Pool(processes=cpus)
        pool.map(run_tomtom, range(len(pssms)))
        pool.close()
        pool.join()

        print 'Reading in Tomtom run...'
        for run in range(len(pssms)):
            tomtomPValues = {}
            with open(cfg.tmpdir_path('tomtom_out/tomtom%d.out' % run), 'r') as infile:
                output = infile.readlines()

            # Now iterate through output and save data
            output.pop(0) # Get rid of header
            while len(output) > 0:
                outputLine = output.pop(0).strip().split('\t')
                if len(outputLine)==10:
                    tomtomPValues[outputLine[1]] = float(outputLine[3])
            pValues = tomtomPValues.values()
            similar = 0
            for pValue in pValues:
                if float(pValue) <= float(0.05):
                    similar += 1

            # Write out the results
            mot = outputLine[0].split('_')[1]
            permPValues[outputLine[0]] = { mot + '.consensus': pssm_mod.consensus_motif(pssms[outputLine[0]]),
                                           mot + '.permutedEV<=10': str(len(pValues)),
                                           mot + '.similar':str(similar),
                                           mot+'.permPV':str(float(similar) / 1000.0) }
            pssms[outputLine[0]].permuted_pvalue = float(similar) / 1000.0
            outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].evalue)+',' + pssm_mod.consensus_motif(pssms[outputLine[0]]) + ',' + str(len(pValues)) + ',' + str(similar) + ',' + str(1000) + ',' + str(float(similar)/float(1000)))
        outFile.close()
    else:
        print 'Using precalculated upstream p-values...'
        with open(output_path, 'r') as infile:
            infile.readline()
            upPValues = [line.strip().split(',') for line in infile]

        for line in upPValues:
            pssms[line[0]].permuted_pvalue = float(line[7]) / 1000.0
    print 'Done.\n'


def __run_mirvestigator_3putr(cfg, c1):
    """Compare 3' UTR Weeder Motifs to miRBase using miRvestigator"""
    print 'Running miRvestigator on 3\' UTR Motifs:'
    pkl_path = cfg.outdir_path('m2m.pkl')
    if not os.path.exists(pkl_path):
        print 'Computing miRNA matches...'
        pssms = c1.pssms_3putr()
        seqs3pUTR = c1.seqs3pUTR.values()

        m2m = miRvestigator(pssms.values(),
                            seqs3pUTR,
                            seedModel=cfg['mirvestigator']['seedModel'],
                            minor=bool(cfg['mirvestigator']['minor']),
                            p5=bool(cfg['mirvestigator']['p5']),
                            p3=bool(cfg['mirvestigator']['p3']),
                            wobble=bool(cfg['mirvestigator']['wobble']),
                            wobbleCut=cfg['mirvestigator']['wobbleCut'],
                            baseDir='output',
                            species=cfg['mirvestigator']['species'])
        with open(pkl_path, 'wb') as outfile:
            cPickle.dump(m2m, outfile)
    else:
        print 'Loading precached miRNA matches...'
        with open(pkl_path, 'rb') as infile:
            m2m = cPickle.load(infile)
    print 'Done.\n'


def __convert_mirvestigator_3putr_results(cfg, c1, mirna_ids):
    """Convert miRNAs and Get permuted p-values for 3' UTR motifs"""

    pssms = c1.pssms_3putr()
    print 'Loading miRvestigator results...'
    # Convert miRvestigator results
    pkl_path = cfg.outdir_path('miRvestigatorResults.pkl')
    if not os.path.exists(pkl_path):
        with open(cfg.outdir_path('miRNA/scores.csv'), 'r') as infile:
            infile.readline()  # get rid of header
            lines = [i.strip().split(',') for i in infile.readlines()]

        miRNA_matches = {}
        for line in lines:
            if not line[1] == 'NA':
                miRNA_mature_seq_ids = []
                miRNA_mature_seq_ids = line[1].split('_')
                miRNA_matches[line[0]] = {'miRNA':line[1],'model':line[2],'mature_seq_ids':miRNA_mature_seq_ids}
                for m1 in miRNA_mature_seq_ids:
                    pssms[line[0]].add_match(factor=m1, confidence=line[2])

        with open(pkl_path, 'wb') as outfile:
            cPickle.dump(miRNA_matches, outfile)
    else:
        with open(pkl_path, 'rb') as infile:
            miRNA_matches = cPickle.load(infile)
            for m1 in miRNA_matches:
                for m2 in miRNA_matches[m1]['mature_seq_ids']:
                    pssms[m1].add_match(factor=m2, confidence=miRNA_matches[m1]['model'])

    # Compile results to put them into the postProcessed
    print 'Get perumted p-values for 3\' UTR motifs...'
    with open('randPSSMs/weederRand.pkl', 'rb') as infile:
        weederRand = cPickle.load(infile)

    with open('randPSSMs/weederRand_all.pkl', 'rb') as infile:
        weederRand_all = cPickle.load(infile)

    clustSizes = sorted(weederRand['8bp'].keys())
    for pssm1 in pssms:
        seqNum = pssms[pssm1].num_genes()
        consensus = pssm_mod.consensus_motif(pssms[pssm1])
        width = len(consensus)
        splitUp = pssm1.split('_')
        clustInd = 5
        for c in clustSizes:
            if seqNum > c:
                clustInd = c
        pValue = float(sum(1 for i in weederRand[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].evalue)))/float(len(weederRand[str(width)+'bp'][clustInd]))
        pValue_all = float(sum(1 for i in weederRand_all[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].evalue)))/float(len(weederRand_all[str(width)+'bp'][clustInd]))

        # note the inconsistent usage of the permuted_pvalue attribute
        pssms[pssm1].permuted_pvalue = {'pValue':pValue,'pValue_all':pValue_all}
    print 'Done.\n'


def run_replication(rep_script):
    ## TODO:
    ## these R scripts are external to the project and therefore
    ## we can not generalize them safely. The script has a "loc1" variable which
    ## hardcodes the output path, we need to make both
    ## What actually needs to be done is to move the script into this project and
    ## make it configurable
    print '  Replication running for %s...' % rep_script
    ret = subprocess.check_call('Rscript replication_'+rep_script+'/replicationDatasetPermutation.R -d '+rep_script+' -o output',
                                stderr=subprocess.STDOUT, shell=True)
    if ret == 1:
        raise Exception('could not run replication - check dependencies')


def __make_replication_pvalues(cfg, c1):
    # Dump a file containing all the genes for each cluster
    with open(cfg.outdir_path('cluster.members.genes.txt'), 'w') as cmgFile:
        writeMe = []
        for cluster_num, bicluster in c1.biclusters.items():
            writeMe.append(str(cluster_num) + ' ' + ' '.join(bicluster.genes))
        cmgFile.write('\n'.join(writeMe))

    # Dump a file containing all the genes for each cluster
    with open(cfg.outdir_path('cluster.members.conditions.txt'), 'w') as cmcFile:
        writeMe = []
        for cluster_num, bicluster in c1.biclusters.items():
            writeMe.append(str(cluster_num) + ' ' + ' '.join(bicluster.conditions))
        cmcFile.write('\n'.join(writeMe))

    # Run replication on all datasets
    run_sets = [name for name in cfg["replication-dataset-names"]
                if not os.path.exists(cfg.outdir_path('replicationPvalues_%s.csv' % name))]

    if len(run_sets) > 0:
        print 'Run replication..'
        # Run this using all cores available
        cpus = int(cfg['cores'])
        print 'There are', cpus,'CPUs avialable.'
        pool = Pool(processes=cpus)
        pool.map(run_replication, run_sets)
        pool.close()
        pool.join()


def __read_replication_pvalues(cfg, c1):
    """TODO: This function still contains the largest refactoring potential
    It is hardcodes the dataset names and is inherently redundant. The way
    to go is to ensure uniform input and output formats and then use the same
    code to handle all data set types."""
    # Read in replication p-values - mesoTCGA Dataset
    print 'Loading replication p-values...'
    with open(cfg.outdir_path('replicationPvalues_mesoTCGA.csv'), 'r') as infile:
        infile.readline()
        for line in infile:
            splitUp = line.strip().split(',')
            bicluster = c1.biclusters[int(splitUp[0].replace('"',''))]
            bicluster.add_attribute(key='replication_mesoTCGA',value=dict(zip(['mesoTCGA_var.exp','mesoTCGA_avg.pc1.var.exp','mesoTCGA_pc1.perm.p','mesoTCGA_OS','mesoTCGA_OS.p','mesoTCGA_OS.age','mesoTCGA_OS.age.p', 'mesoTCGA_OS.age.sex','mesoTCGA_OS.age.sex.p'],splitUp[3:])))
    print 'Done.\n'


def __make_permuted_pvalues(cfg, c1):
    ###########################################################################
    ## Run permuted p-value for variance epxlained first principal component ##
    ###########################################################################
    pvalues_path = cfg.outdir_path('residualPermutedPvalues_permAll.csv')
    if not os.path.exists(pvalues_path):
        print 'Calculating FPC permuted p-values...'
        ret = subprocess.check_call(['./permutedResidualPvalues_permAll_mc.R',
                                     '-o', cfg['outdir'],
                                     '-r', cfg['ratios-file'],
                                     '-c', str(cfg['cores'])],
                                    stderr=subprocess.STDOUT)
        if ret == 1:
            print "error in calling R script."
            exit(1)

        print 'Done.\n'

    #################################################################
    ## Read in residual permutations to use for filtering          ##
    #################################################################
    print 'Load residual permuted p-values...'
    with open(pvalues_path, 'r') as infile:
        # "","bicluster","n.rows","n.cols","orig.resid","avg.perm.resid","perm.p","orig.resid.norm","avg.norm.perm.resid","norm.perm.p","pc1.var.exp","avg.pc1.var.exp","pc1.perm.p"
        infile.readline()
        for line in infile:
            splitUp = line.strip().split(',')
            bicluster = c1.biclusters[int(splitUp[0].strip('"'))]
            bicluster.add_attribute(key='pc1.perm.p',value=str(splitUp[7]))
    print 'Done.\n'


def __make_functional_enrichment(cfg, c1):
    #################################################################
    ## Run functional enrichment                                   ##
    #################################################################
    # Note that these are external to the project and have hard-coded paths !!!
    if not os.path.exists(cfg.outdir_path('biclusterEnrichment_GOBP.csv')):
        print 'Run functional enrichment...'
        ret = subprocess.check_call("./enrichment.R -o "+cfg['outdir'], # +" -c "+cfg['gene_conv']
                                    stderr=subprocess.STDOUT, shell=True)
        if ret == 1:
            raise Exception('could not run functional enrichment')

        print 'Done.\n'

    #################################################################
    ## Read in functional enrichment                               ##
    #################################################################
    print 'Load GO Biological Process functional enrichment...'
    with open(cfg.outdir_path('biclusterEnrichment_GOBP.csv'), 'r') as infile:
        infile.readline()  # Get rid of header
        lines = [line.strip().split(',') for line in infile]

    for line in lines:
        bicluster = c1.biclusters[int(line[0].strip('"'))]
        bicluster.add_attribute(key='goTermBP',value=line[2].strip('"').split(';'))
    print 'Done.\n'


def __make_go_term_semantic_similarity(cfg, c1):
    #################################################################
    ## Run GO term semantic similarity                             ##
    #################################################################
    if not os.path.exists(cfg.outdir_path('jiangConrath_hallmarks.csv')):
        print 'Run semantic similarity...'
        ret = subprocess.check_call("./goSimHallmarksOfCancer.R -o %s" % cfg.outdir,
                                    stderr=subprocess.STDOUT, shell=True)
        if ret == 1:
            raise Exception('could not run semantic similarity')

        print 'Done.\n'

    print 'Load Jiang-Conrath semantic similarity to Hallmarks of Cancer...'
    with open(cfg.outdir_path('jiangConrath_hallmarks.csv'), 'r') as infile:
        hallmarks = [i for i in infile.readline().split(',') if not i.strip('"')=='']
        lines = [line.strip().split(',') for line in infile]

        for line in lines:
            bicluster = c1.biclusters[int(line[0].strip('"'))]
            bicluster.add_attribute(key='hallmarksOfCancer',value=dict(zip(hallmarks,line[1:])))
        print 'Done.\n'


def perform_postprocessing(cfg, c1, entrez2id, mirna_ids):
    pkl_path = cfg.outdir_path('c1_postProc.pkl')
    if not os.path.exists(pkl_path):
        interim_path = cfg.outdir_path('c1_postInterim.pkl')
        if not os.path.exists(interim_path):
            dump_cluster_members(cfg, c1)
            __get_cluster_eigengenes(cfg, c1)
            __get_cluster_variance_explained(cfg, c1)
            phenotypes = __get_phenotype_info(cfg, c1)

            __do_postprocess(cfg, c1, phenotypes)

            __tomtom_upstream_motifs(cfg)
            tf_name2entrezid, tf_families = __expand_tf_factor_list(entrez2id)
            exp_data, all_names = __correlate_tfs_with_cluster_eigengenes(cfg, c1)
            __expand_and_correlate_tfbsdb_tfs(c1, tf_name2entrezid, tf_families, exp_data, all_names)
            print 'Dumping interim post-processing Object:'
            with open(interim_path, 'wb') as outfile:
                cPickle.dump(c1, outfile)
            print 'Done.\n'
        else:
            print 'Loading from precached cMonkey Object:'
            with open(interim_path, 'rb') as infile:
                c1 = cPickle.load(infile)
            print 'Done.\n'
        __write_first_principal_components(cfg, c1)
        __get_permuted_pvalues_for_upstream_meme_motifs(cfg, c1)
        __run_mirvestigator_3putr(cfg, c1)
        __convert_mirvestigator_3putr_results(cfg, c1, mirna_ids)
        __make_replication_pvalues(cfg, c1)
        __read_replication_pvalues(cfg, c1)
        __make_permuted_pvalues(cfg, c1)
        __make_functional_enrichment(cfg, c1)
        __make_go_term_semantic_similarity(cfg, c1)

        print 'Dumping Final cMonkey Object:'
        with open(pkl_path, 'wb') as outfile:
            cPickle.dump(c1, outfile)
        print 'Done.\n'

    else:
        print 'Loading from precached cMonkey Object:'
        with open(pkl_path, 'rb') as infile:
            c1 = cPickle.load(infile)
        print 'Done.\n'
    return c1


def run_neo(cfg):
    """Run NEO and integrate some form of results """
    if not os.path.exists(cfg.outdir_path('causal')):
        ## Run the runNEO.R script and do the causality analyses
        print '  Network edge orienting (NEO)...'
        print ['./NEO/runNEO_V2.R', '-o', cfg['outdir'], '-r', cfg['all-ratios-file'], '-s', cfg['som-muts-file'], '-t', 'meso', '-e', 'output/biclusterEigengenes.csv', '-c', str(cfg['cores'])]
        ret = subprocess.check_call(['./NEO/runNEO_VnoMir.R',
                                         '-o', cfg['outdir'],
                                         '-r', cfg['all-ratios-file'],
                                         #'-m', cfg['mirna-file'],
                                         '-s', cfg['som-muts-file'],
                                         '-t', 'meso',
                                         '-e', 'output/biclusterEigengenes.csv',
                                         '-c', str(cfg['cores'])],
                                        stderr=subprocess.STDOUT)
        if ret == 1:
            raise Exception('could not run causality analyses')


def write_neo_summary(cfg):
    ## Pull together analysis into cohesive output
    causalSummary = []
    # For each mutation
    for file1 in os.listdir(cfg.outdir_path('causal')):
        if file1[0:7]=='causal_':
            with open(cfg.outdir_path('causal/%s' % (file1)), 'r') as inFile:
                inLine = inFile.readline() # Get rid of header
                while 1:
                    inLine = inFile.readline()
                    if not inLine:
                        break
                    splitUp = inLine.strip().split(',')
                    if (float(splitUp[6]) > cfg['leo-nb-atob'] and
                        float(splitUp[12]) <= cfg['mlogp-m-atob']):
                        # Somatic Mutation(1), Regulator(3), Biclster(5), leo.nb.AtoB(6), mlogp.M.AtoB(12), PathAB(17), SEPathAB(18), ZPathAB(19), PPathAB(20), BLV.AtoB(25), RMSEA.AtoB(28)
                        causalSummary.append({'Mutation': splitUp[1].strip('"').lstrip('M:'), 'Regulator': splitUp[3].strip('"').lstrip('A:X'), 'Bicluster': splitUp[5].strip('"').lstrip('B:bic_'), 'leo.nb.AtoB': splitUp[6], 'mlogp.M.AtoB': splitUp[12], 'PathAB': splitUp[17], 'SEPathAB': splitUp[18], 'ZPathAB': splitUp[19], 'PPathAB': splitUp[20], 'BLV.AtoB': splitUp[25], 'RMSEA.AtoB': splitUp[28]})

    ## Output:  Somatic Mutation(1), Regulator(3), Biclster(5), leo.nb.AtoB(6), mlogp.M.AtoB(12), PathAB(17), SEPathAB(18), ZPathAB(19), PPathAB(20), BLV.AtoB(25), RMSEA.AtoB(28)
    header = ['Mutation', 'Regulator', 'Bicluster', 'leo.nb.AtoB', 'mlogp.M.AtoB', 'PathAB', 'SEPathAB', 'ZPathAB', 'PPathAB', 'BLV.AtoB', 'RMSEA.AtoB']
    with open(cfg.outdir_path('causalitySummary.csv'), 'w') as outfile:
        outfile.write(','.join(header)+'\n')
        outfile.write('\n'.join([','.join([i[j] for j in header]) for i in causalSummary]))
    return causalSummary


def add_correspondent_regulators(cfg, c1, causal_summary, mirna_ids_rev):
    """Dump out correspondent regulators (both mechanistically and causally predicted)"""
    correspondentRegulators = {}
    for causalFlow in causal_summary:
        bicluster = c1.biclusters[int(causalFlow['Bicluster'])]
        ## Upstream (TFs)
        tfs = []
        # 1. MEME and WEEDER Upstream motifs
        for pssm in bicluster.pssms_upstream:
            matches = pssm.correlated_matches()
            if not matches is None:
                for corTf in matches:
                    if (corTf['pValue'] <= cfg['pvalue-cut'] and
                        abs(corTf['rho']) >= cfg['rho-cut']):
                        tfs.append(corTf['factor'])
            # 2. TFBS_DB
            for corTf in bicluster.attributes['tfbs_db_correlated']:
                if corTf['pValue'] <= cfg['pvalue-cut'] and abs(corTf['rho']) >= cfg['rho-cut']:
                    tfs.append(corTf['factor'])
        # 3. Find Correspondent TF regulators
        if causalFlow['Regulator'] in tfs:
            if not int(causalFlow['Bicluster']) in correspondentRegulators:
                correspondentRegulators[int(causalFlow['Bicluster'])] = {'tf':[],'miRNA':[]}
            correspondentRegulators[int(causalFlow['Bicluster'])]['tf'].append(causalFlow['Regulator'])

        ## 3' UTR (miRNA)
        miRNAs = []
        # 1. WEEDER 3'UTR
        for pssm in bicluster.pssms_3putr:
            for miR in pssm.matches:
                if miR['confidence'] in ['8mer','7mer_a1','7mer_m8']:
                    print miR['factor'], miR['factor'] in mirna_ids_rev
                    #miRNAs += mirna_ids_rev[miR['factor']]
                    miRNAs += miR['factor']

        # 2. PITA (not correlated)
        if (float(bicluster.attributes['pita_3pUTR']['pValue']) <= cfg['pvalue-cut'] and
            float(bicluster.attributes['pita_3pUTR']['percentTargets'].split(' ')[0]) >= cfg['perc-targets']):
            miRNAs += bicluster.attributes['pita_3pUTR']['miRNA'].split(' ')
        # 3. TargetScan (not correlated)
        if (float(bicluster.attributes['targetscan_3pUTR']['pValue']) <= cfg['pvalue-cut'] and
            float(bicluster.attributes['targetscan_3pUTR']['percentTargets'].split(' ')[0]) >= cfg['perc-targets']):
            miRNAs += bicluster.attributes['targetscan_3pUTR']['miRNA'].split(' ')
        # 4. Find Correspondent miRNA regulators
        for miR in miRNAs:
            if compareMiRNANames(causalFlow['Regulator'].lower(), miR.lower()):
                if not int(causalFlow['Bicluster']) in correspondentRegulators:
                    correspondentRegulators[int(causalFlow['Bicluster'])] = {'tf':[],'miRNA':[]}
                correspondentRegulators[int(causalFlow['Bicluster'])]['miRNA'].append(causalFlow['Regulator'])

    ## Put correspondent regulators into cMonkeyWrapper object
    for cluster_num in correspondentRegulators.keys():
        bicluster = c1.biclusters[cluster_num]
        bicluster.add_attribute(key='correspondentRegulators', value=correspondentRegulators[cluster_num])


def write_final_result(cfg, c1, mirna_ids_rev):
    #################################################################
    ## Write out the final post-processed file                     ##
    #################################################################
    print 'Write postProcessedVFinal.csv...'
    postOut = []
    hallmarksOfCancer = c1.biclusters[1].attributes['hallmarksOfCancer'].keys()

    for cluster_num in sorted(c1.biclusters.keys()):
        writeMe = []
        bicluster = c1.biclusters[cluster_num]
        # Write file line by line
        #   a. Bicluster basics:  id, genes, conditions, resid, resid.norm, resid.norm.perm.p
        print cluster_num #, bicluster.attributes
        writeMe += [str(cluster_num),
                    str(bicluster.num_genes()),  # genes
                    str(bicluster.num_conditions())]  # conditions
        writeMe += [str(bicluster.attributes['pc1.var.exp']), # Variance explained by first principal component
                    str(bicluster.attributes['pc1.perm.p'])] # Variance explained by first principal component
        #   b. Upstream motifs:  meme.motif1.E, meme.motif1.consensus, meme.motif1.matches, meme.motif1.permPV
        motifNames = bicluster.pssm_names_upstream()
        upstreamMotifs = { 'meme_motif1':None, 'meme_motif2':None, 'weeder_motif1':None, 'weeder_motif2':None }
        for m1 in motifNames:
            splitUp = m1.split('_')
            if splitUp[1]=='motif1' and splitUp[2]=='meme':
                upstreamMotifs['meme_motif1'] = m1
            if splitUp[1]=='motif2' and splitUp[2]=='meme':
                upstreamMotifs['meme_motif2'] = m1
            if splitUp[1]=='motif1' and splitUp[2]=='weeder':
                upstreamMotifs['weeder_motif1'] = m1
            if splitUp[1]=='motif2' and splitUp[2]=='weeder':
                upstreamMotifs['weeder_motif2'] = m1
        #   - MEME motifs
        for meme1 in ['meme_motif1','meme_motif2']:
            if not upstreamMotifs[meme1]==None:
                pssm1 = bicluster.pssm_upstream(upstreamMotifs[meme1])
                original = []
                matches = 'NA'
                expandedMatches = 'NA'
                correlatedMatches = {}
                originalExpanded = {}
                minCorrelated = {}
                correlatedMatches = 'NA'
                originalExpanded = 'NA'
                minCorrelated = 'NA'
                if len(pssm1.matches) > 0:
                    matches = ' '.join([match1['factor'] for match1 in pssm1.matches])
                    if len(pssm1.expanded_matches) > 0:
                        expandedMatches = {}
                        for i in pssm1.expanded_matches:
                            if not i['seedFactor'] in original:
                                original.append(i['seedFactor'])
                            if not i['seedFactor'] in expandedMatches:
                                expandedMatches[i['seedFactor']] = []
                            expandedMatches[i['seedFactor']].append(i['factor'])
                        expandedMatches = ' '.join([seedFactor+':'+';'.join(expandedMatches[seedFactor]) for seedFactor in expandedMatches])
                    tmp = pssm1.correlated_matches()
                    if not tmp is None:
                        for match1 in tmp:
                            if (match1['pValue'] <= cfg['pvalue-cut'] and
                                abs(match1['rho']) >= cfg['rho-cut']):
                                if correlatedMatches=='NA':
                                    correlatedMatches = []
                                correlatedMatches.append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                                if minCorrelated=='NA' or match1['pValue']<minCorrelated['pValue']:
                                    minCorrelated = match1
                        if not correlatedMatches=='NA':
                            correlatedMatches = ' '.join(correlatedMatches)
                        if not minCorrelated=='NA':
                            if minCorrelated['factor'] in original:
                                originalExpanded = 'Original'
                            else:
                                originalExpanded = 'Expanded'
                            minCorrelated = minCorrelated['factor']+':'+str(minCorrelated['rho'])+':'+str(minCorrelated['pValue'])
                writeMe += ([str(pssm1.evalue),
                             str(pssm1.permuted_pvalue),
                             pssm_mod.consensus_motif(pssm1),
                             matches,  # motif natches from TransFac and Jaspar
                             expandedMatches]  # Expanded matches (TFClass)
                            + [correlatedMatches+','+originalExpanded+','+minCorrelated])
            else:
                writeMe += (['NA', # E-value
                            'NA', # Permuted p-value for motif
                            'NA', # Motif consensus sequence
                            'NA', # Matches to the motif from TransFac and Jaspar
                            'NA'] # Expanded matches using TFClass
                            + ['NA,NA,NA'])
        #   - WEEDER motifs
        for weeder1 in ['weeder_motif1','weeder_motif2']:
            if not upstreamMotifs[weeder1]==None:
                pssm1 = bicluster.pssm_upstream(upstreamMotifs[weeder1])
                original = []
                matches = 'NA'
                expandedMatches = 'NA'
                correlatedMatches = {}
                originalExpanded = {}
                minCorrelated = {}
                correlatedMatches = 'NA'
                originalExpanded = 'NA'
                minCorrelated = 'NA'
                if len(pssm1.matches) > 0:
                    matches = ' '.join([match1['factor'] for match1 in pssm1.matches])
                    if len(pssm1.expanded_matches) > 0:
                        expandedMatches = {}
                        for i in pssm1.expanded_matches:
                            if not i['seedFactor'] in original:
                                original.append(i['seedFactor'])
                            if not i['seedFactor'] in expandedMatches:
                                expandedMatches[i['seedFactor']] = []
                            expandedMatches[i['seedFactor']].append(i['factor'])
                        expandedMatches = ' '.join([seedFactor+':'+';'.join(expandedMatches[seedFactor]) for seedFactor in expandedMatches])
                    tmp = pssm1.correlated_matches()
                    if not tmp is None:
                        for match1 in tmp:
                            if (match1['pValue'] <= cfg['pvalue-cut'] and
                                abs(match1['rho']) >= cfg['rho-cut']):
                                if correlatedMatches=='NA':
                                    correlatedMatches = []
                                correlatedMatches.append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                                if minCorrelated=='NA' or match1['pValue']<minCorrelated['pValue']:
                                    minCorrelated = match1
                        if not correlatedMatches=='NA':
                            correlatedMatches = ' '.join(correlatedMatches)
                        if not minCorrelated=='NA':
                            if minCorrelated['factor'] in original:
                                originalExpanded = 'Original'
                            else:
                                originalExpanded = 'Expanded'
                            minCorrelated = minCorrelated['factor']+':'+str(minCorrelated['rho'])+':'+str(minCorrelated['pValue'])
                writeMe += ([str(pssm1.evalue),
                             pssm_mod.consensus_motif(pssm1),
                             matches,  # motif matches (TransFac and Jaspar)
                             expandedMatches]  # Expanded matches (TFClass)
                            + [correlatedMatches+','+originalExpanded+','+minCorrelated])
            else:
                writeMe += (['NA', # E-value
                             'NA', # Motif consensus sequence
                             'NA', # Matches to the motif from TransFac and Jaspar
                             'NA'] # Expanded matches using TFClass
                            + ['NA,NA,NA'])
        #   c. Enriched TFs:  TFBS_DB.TFs,TFBS_DB.percTargets,TFBS_DB.pValue
        for association in ['tfbs_db']:
            a1 = bicluster.attributes[association]
            if a1['tf'] != '':
                expandedMatches = 'NA'
                correlatedMatches = {}
                originalExpanded = {}
                minCorrelated = {}
                correlatedMatches = 'NA'
                originalExpanded = 'NA'
                minCorrelated = 'NA'
                tmp = bicluster.attributes['tfbs_db_expanded']
                if not tmp is None:
                    expandedMatches = ' '.join([seedFactor+':'+';'.join(tmp[seedFactor]) for seedFactor in tmp])
                tmp = bicluster.attributes['tfbs_db_correlated']
                if not tmp is None:
                    for match1 in tmp:
                        if (match1['pValue'] <= cfg['pvalue-cut'] and
                            abs(match1['rho']) >= cfg['rho-cut']):
                            if correlatedMatches=='NA':
                                correlatedMatches = []
                            correlatedMatches.append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                            if minCorrelated=='NA' or match1['pValue']<minCorrelated['pValue']:
                                minCorrelated = match1
                    if not correlatedMatches=='NA':
                        correlatedMatches = ' '.join(correlatedMatches)
                    if not minCorrelated=='NA':
                        if minCorrelated['factor'] in original:
                            originalExpanded = 'Original'
                        else:
                            originalExpanded = 'Expanded'
                        minCorrelated = minCorrelated['factor']+':'+str(minCorrelated['rho'])+':'+str(minCorrelated['pValue'])
                if len(correlatedMatches)==0:
                    correlatedMatches = 'NA'
                writeMe += [str(a1['tf']).replace(';',' '),
                            str(a1['percentTargets']).replace(';',' '),
                            str(a1['pValue']),
                            expandedMatches] + [correlatedMatches + ',' + originalExpanded + ',' + minCorrelated]
            else:
                writeMe += (['NA','NA','NA','NA'] + ['NA,NA,NA'])

        #   d. 3' UTR motifs:  weeder.motif1.E, weeder.motif1.permPV, weeder.motif1.permPV_all, weeder.motif1.consensus, weeder.motif1.matches, weeder.motif1.model
        motifNames = bicluster.pssm_names_3putr()
        p3utrMotifs = { 'weeder_motif1':None, 'weeder_motif2':None }
        for m1 in motifNames:
            splitUp = m1.split('_')
            if splitUp[1]=='motif1' and splitUp[2]=='weeder':
                p3utrMotifs['weeder_motif1'] = m1
            if splitUp[1]=='motif2' and splitUp[2]=='weeder':
                p3utrMotifs['weeder_motif2'] = m1
        #   - WEEDER motifs
        for weeder1 in ['weeder_motif1','weeder_motif2']:
            if not p3utrMotifs[weeder1]==None:
                pssm1 = bicluster.pssm_3putr(p3utrMotifs[weeder1])
                matches = 'NA'
                model = 'NA'
                if len(pssm1.matches) > 0:
                    #matches = ' '.join([mirna_ids_rev[match1['factor']] for match1 in pssm1.matches])
                    matches = ' '.join([match1['factor'] for match1 in pssm1.matches])
                    model = pssm1.matches[0]['confidence']
                permutedPValue = 'NA' if pssm1.permuted_pvalue is None else pssm1.permuted_pvalue
                writeMe += [str(pssm1.evalue),
                            str(permutedPValue['pValue']), # Permuted p-value for motif
                            str(permutedPValue['pValue_all']), # Permuted p-value for motif (all)
                            pssm_mod.consensus_motif(pssm1),
                            matches, # Matches to the motif to miRBase
                            model] # Model fit by miRvestigator
            else:
                writeMe += ['NA', # E-value
                            'NA', # Permuted p-value for motif
                            'NA', # Permuted p-value for motif (all)
                            'NA', # Motif consensus sequence
                            'NA', # Matches to the motif to miRBase
                            'NA'] # Model fit by miRvestigator

        #   e. Enriched miRNAs:  3pUTR_pita.miRNAs,3pUTR_pita.percTargets,3pUTR_pita.pValue,3pUTR_targetScan.miRNAs,3pUTR_targetScan.percTargets,3pUTR_targetScan.pValue
        for association in ['pita_3pUTR', 'targetscan_3pUTR']:
            a1 = bicluster.attributes[association]
            if a1['miRNA'] != '':
                writeMe += [str(a1['miRNA']).replace(';',' '), str(a1['percentTargets']).replace(';',' '), str(a1['pValue'])]
            else:
                writeMe += ['NA','NA','NA']

        #   f. Correpsondent regulators
        if cfg['run-neo']:
            corrRegs = bicluster.attributes['correspondentRegulators'] if 'correspondentRegulators' in bicluster.attributes else None
            if not corrRegs is None:
                # Split up by TF and miRNA
                if not len(corrRegs['tf'])==0:
                    writeMe += [' '.join(sorted(set(corrRegs['tf'])))]
                else:
                    writeMe += ['NA']
                """if not len(corrRegs['miRNA'])==0:
                    writeMe += [' '.join(sorted(set(corrRegs['miRNA'])))]
                else:
                    writeMe += ['NA']
                """
            else:
                #writeMe += ['NA','NA']
                writeMe += ['NA']

        #   g. Associations with traits:  age, sex.bi, chemo_therapy, radiation_therapy
        for association in PHENOTYPES:
            if association in bicluster.attributes:
                print 'ass1', bicluster.attributes[association]
                ass1 = bicluster.attributes[association]
                writeMe += [str(ass1['T']), str(ass1['pValue'])]
            else:
                writeMe += ['NA', 'NA']
        surv1 = bicluster.attributes['Survival']
        survAge1 = bicluster.attributes['Survival.AGE']
        survAgeSex1 = bicluster.attributes['Survival.AGE_SEX']
        writeMe += [str(surv1['z']), str(surv1['pValue']), str(survAge1['z']), str(survAge1['pValue']), str(survAgeSex1['z']), str(survAgeSex1['pValue'])]

        #   h. Independent replication:
        replications_mesoTCGA = bicluster.attributes['replication_mesoTCGA']
        for replication in ['mesoTCGA_var.exp','mesoTCGA_avg.pc1.var.exp','mesoTCGA_pc1.perm.p','mesoTCGA_OS','mesoTCGA_OS.p','mesoTCGA_OS.age','mesoTCGA_OS.age.p', 'mesoTCGA_OS.age.sex','mesoTCGA_OS.age.sex.p']:
            writeMe.append(str(replications_mesoTCGA[replication]))

        #   i. Functional enrichment of biclusters using GO term Biological Processes
        bfe1 = bicluster.attributes['goTermBP']
        if bfe1 == ['']:
            writeMe.append('NA')
        else:
            writeMe.append(';'.join(bfe1))

        #   j. Hallmarks of Cancer:  Hanahan and Weinberg, 2011
        bhc1 = bicluster.attributes['hallmarksOfCancer']
        for hallmark in hallmarksOfCancer:
            writeMe.append(str(bhc1[hallmark]))

        #   k. Glioma sub-type enrichment: 'NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM'
        #for overlap in ['NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']:
        #    writeMe.append(str(bicluster.attributes[overlap]))

        #   l. Add genes
        writeMe.append(' '.join(bicluster.genes))

        # Add to the final output file
        postOut.append(deepcopy(writeMe))

    with open(cfg.outdir_path(cfg['postprocessing-result-file']), 'w') as postFinal:
        header = ['Bicluster', 'Genes', 'Patients']
        header += ['Var. Exp. First PC', 'Var. Exp. First PC Perm. P-Value']
        header += ['MEME Motif1 E-Value', 'MEME Motif1 Perm. P-Value', 'Up.MEME Motif1 Consensus', 'Up.MEME Motif1 Matches','Up.MEME Motif1 Expanded Matches'] + ['Up.MEME Motif1 Correlated Matches', 'Up.MEME Motif1 Original/Expanded', 'Up.MEME Motif1 Minimum Correlated'] + \
         ['Up.MEME Motif2 E-Value', 'Up.MEME Motif2 Perm. P-Value', 'Up.MEME Motif2 Consensus', 'Up.MEME Motif2 Matches', 'Up.MEME Motif2 Expanded Matches'] + ['Up.MEME Motif2 Correlated Matches', 'Up.MEME Motif2 Original/Expanded', 'Up.MEME Motif2 Minimum Correlated'] + \
         ['Up.WEEDER Motif1 Score', 'Up.WEEDER Motif1 Consensus', 'Up.WEEDER Motif1 Matches', 'Up.WEEDER Motif1 Expanded Matches'] + ['Up.WEEDER Motif1 Correlated Matches', 'Up.WEEDER Motif1 Original/Expanded', 'Up.WEEDER Motif1 Minimum Correlated'] + \
         ['Up.WEEDER Motif2 Score', 'Up.WEEDER Motif2 Consensus', 'Up.WEEDER Motif2 Matches', 'Up.WEEDER Motif2 Expanded Matches'] + ['Up.WEEDER Motif2 Correlated Matches', 'Up.WEEDER Motif2 Original/Expanded', 'Up.WEEDER Motif2 Minimum Correlated'] + \
         ['TFBS_DB.TFs', 'TFBS_DB.percTargets', 'TFBS_DB.pValue', 'TFBS_DB.Exapnded Matches'] + ['TFBS_DB.Correlated Matches', 'TFBS_DB.Original/Expanded','TFBS_DB.Minimum Correlated'] + \
         ['3pUTR.WEEDER Motif1 E-Value', '3pUTR.WEEDER Motif1 Perm. P-Value', '3pUTR.WEEDER Motif1 Perm. P-Value (All)', '3pUTR.WEEDER Motif1 Consensus', '3pUTR.WEEDER Motif1 Matches', '3pUTR.WEEDER Motif1 Model'] + \
         ['3pUTR.WEEDER Motif2 E-Value', '3pUTR.WEEDER Motif2 Perm. P-Value', '3pUTR.WEEDER Motif2 Perm. P-Value (All)', '3pUTR.WEEDER Motif2 Consensus', '3pUTR.WEEDER Motif2 Matches', '3pUTR.WEEDER Motif2 Model'] + \
         ['3pUTR_pita.miRNAs', '3pUTR_pita.percTargets', '3pUTR_pita.pValue'] + \
         ['3pUTR_targetScan.miRNAs', '3pUTR_targetScan.percTargets', '3pUTR_targetScan.pValue']
        #if cfg['run-neo']:
        #    header += ['Correspondent.TFs', 'Correspondent.miRNAs']
        if cfg['run-neo']:
            header += ['Correspondent.TFs']
        header += sum([[association, association+'.p'] for association in PHENOTYPES+['OS','OS.covAge','OS.covAgeSex']],[]) + \
         ['mesoTCGA_var.exp','mesoTCGA_avg.pc1.var.exp','mesoTCGA_pc1.perm.p','mesoTCGA_OS','mesoTCGA_OS.p','mesoTCGA_OS.age','mesoTCGA_OS.age.p', 'mesoTCGA_OS.age.sex','mesoTCGA_OS.age.sex.p'] + \
         ['GO_Term_BP'] + \
         [i.strip() for i in hallmarksOfCancer] + \
         ['Genes']
        postFinal.write(','.join(header)+'\n'+'\n'.join([','.join(i) for i in postOut]))

    print 'Done.\n'


if __name__ == '__main__':
    cfg = config.sygnal_init()

    entrez2id = read_synonyms(cfg)
    if bool(cfg['gene_conv']):
        gene_conv = entrez2id
    else:
        gene_conv = False
    mirna_ids, mirna_ids_rev = miRNA_mappings(cfg)
    c1 = compute_additional_info(cfg, mirna_ids, gene_conv)
    c1 = perform_postprocessing(cfg, c1, entrez2id, mirna_ids)
    # Run neo only if True
    if cfg['run-neo']:
        run_neo(cfg)
        causal_summary = write_neo_summary(cfg)
        add_correspondent_regulators(cfg, c1, causal_summary, mirna_ids_rev)
    write_final_result(cfg, c1, mirna_ids_rev)
