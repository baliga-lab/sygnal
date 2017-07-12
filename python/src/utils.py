#################################################################
# @Program: utils.py                                            #
# @Version: 1                                                   #
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
# Copyrighted by Chris Plaisier  5/17/2012                      #
#################################################################

import re
import pssm


def mirna_in_dict(mirna, mirna_ids):
    """convert miRNA names to miRBase IDs"""
    return [mirna_id for name, mirna_id in mirna_ids.items()
            if __compare_mirna_names(mirna, name)]


def __compare_mirna_names(a, b):
    if a == b:
        return True

    if len(a) < len(b):
        re1 = re.compile(a+'[a-z]$')
        if re1.match(b):
            return True
    else:
        re1 = re.compile(b+'[a-z]$')
        if re1.match(a):
            return True
    return False


def rand_pssm_clust_size(clust_size):
    if clust_size <= 5:
        return 5
    elif clust_size <= 10:
        return 10
    elif clust_size <= 15:
        return 15
    elif clust_size <= 20:
        return 20
    elif clust_size <= 25:
        return 25
    elif clust_size <= 30:
        return 30
    elif clust_size <= 35:
        return 35
    elif clust_size <= 40:
        return 40
    elif clust_size <= 45:
        return 45
    elif clust_size <= 50:
        return 50
    elif clust_size <= 55:
        return 55
    elif clust_size <= 60:
        return 60
    else:
        return 65


def make_files(nuc_freqs, query_pssms, target_pssms, num, strands='+ -'):
    """Make the files for a TomTom run"""
    meme_header = ''
    meme_header += 'MEME version 3.0\n\n'
    meme_header += 'ALPHABET= ACGT\n\n'

    # Here is where we tell it what strand: for miRNAs this would just be '+'
    meme_header += 'strands: '+strands+'\n\n'
    meme_header += 'Background letter frequencies (from genome):\n'
    meme_header += 'A '+str(round(float(nuc_freqs['A']),3))+' C '+str(round(float(nuc_freqs['C']),3))+' G '+str(round(float(nuc_freqs['G']),3))+' T '+str(round(float(nuc_freqs['T']),3))+'\n\n'

    # Make query PSSM file
    with open('tmp/query%d.tomtom' % num, 'w') as query_file:
        query_file.write(meme_header)
        query_file.write('\n\n'.join([pssm.to_meme_str(p) for p in query_pssms]))

    # Make target PSSM file
    with open('tmp/target%d.tomtom' % num,'w') as target_file:
        target_file.write(meme_header)
        target_file.write('\n\n'.join([pssm.to_meme_str(p) for p in target_pssms]))
