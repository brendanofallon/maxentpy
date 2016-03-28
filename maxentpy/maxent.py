'''
maxent.py
Calculate splice site strength
Modified from MaxEntScan perl scripts developed by Gene Yeo and Christopher
Burge
Yeo G and Burge C. Maximum entropy modeling of short sequence motifs with
applications to RNA splicing signals. Journal of Computational Biology,
2004; 11(2-3):377-94.
'''

import math
from collections import defaultdict
from string import maketrans
import os.path


DONOR_BASES = 9 #Total length of sequence to be supplied to score5 (donor seq)
DONOR_JUNCTION_OFFSET=3 #Number of bases into donor seq that splice junction begins

ACCEPTOR_BASES = 23 #Total length of sequence to be given to score3 (acceptor)
ACCEPTOR_JUNCTION_OFFSET=20

dir_path = os.path.dirname(os.path.abspath(__file__))

class SpliceScorer(object):

    _DATA_DIR = "data"
    _DONOR_DATA = "score5_matrix.txt"
    _ACCEPTOR_DATA = "score3_matrix.txt"

    def __init__(self):
        self.donor_lookup = self._init_donor()
        self.acceptor_lookup = self._init_acceptor()

    def _init_donor(self):
        """
        Read in scoring matrix data for splice donor prediction
        """
        dir_path = os.path.dirname(os.path.abspath(__file__))
        matrix_f = dir_path + os.path.join(dir_path, SpliceScorer._DATA_DIR, SpliceScorer._DONOR_DATA)
        data = {}
        with open(matrix_f, 'r') as f:
            for line in f:
                entry = line.split()
                data[entry[0]] = float(entry[1])
        return data

    def _init_acceptor(self):
        """
        Read in matrix for splice acceptor prediction
        :return:
        """
        dir_path = os.path.dirname(os.path.abspath(__file__))
        matrix_f = dir_path + os.path.join(dir_path, SpliceScorer._DATA_DIR, SpliceScorer._ACCEPTOR_DATA)
        matrix = defaultdict(dict)
        with open(matrix_f, 'r') as f:
            for line in f:
                n, m, s = line.split()
                matrix[int(n)][int(m)] = float(s)
        return matrix

    def score5(self, fa):
        '''
        Calculate 5' (donor) splice site strength, requires 9 bp of sequence
        (exon)XXX|XXXXXX(intron)
                  **
        '''
        # for key elements
        if len(fa) != DONOR_BASES:
            raise ValueError('Need at least ' + str(DONOR_BASES) + ' nucleotides for computation')
        key = fa[3:5].upper()
        bgd = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
        cons1 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
        cons2 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}
        score = cons1[key[0]] * cons2[key[1]] / (bgd[key[0]] * bgd[key[1]])
        # for rest elements
        rest = (fa[:3] + fa[5:]).upper()
        rest_score = self.donor_lookup[rest]
        return math.log(score * rest_score, 2)


    def score3(self, fa):
        '''
        Calculate 3' (acceptor) splice site strength
        (intron)XXXXXXXXXXXXXXXXXXXX|XXX(exon)
                                  **
        '''
        # for key elements
        if len(fa) != ACCEPTOR_BASES:
            raise ValueError('Need ' + str(ACCEPTOR_BASES) + ' nucleotides for computation')
        key = fa[18:20].upper()
        bgd = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
        cons1 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
        cons2 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}
        score = cons1[key[0]] * cons2[key[1]] / (bgd[key[0]] * bgd[key[1]])
        # for rest elements
        rest = (fa[:18] + fa[20:]).upper()
        matrix = self.acceptor_lookup
        rest_score = 1
        rest_score *= matrix[0][hashseq(rest[:7])]
        rest_score *= matrix[1][hashseq(rest[7:14])]
        rest_score *= matrix[2][hashseq(rest[14:])]
        rest_score *= matrix[3][hashseq(rest[4:11])]
        rest_score *= matrix[4][hashseq(rest[11:18])]
        rest_score /= matrix[5][hashseq(rest[4:7])]
        rest_score /= matrix[6][hashseq(rest[7:11])]
        rest_score /= matrix[7][hashseq(rest[11:14])]
        rest_score /= matrix[8][hashseq(rest[14:18])]
        # final score
        return math.log(score * rest_score, 2)


def hashseq(fa):
    table = maketrans('ACGT', '0123')
    seq = fa.translate(table)
    return sum(int(j) * 4**(len(seq) - i - 1) for i, j in enumerate(seq))
