import sys
import itertools
import numpy as np
from scipy.stats import itemfreq
from collections import defaultdict


def read_aln(aln_filename):
    """ Read multiple sequence alngnment in jones format,
    i.e. each row represents one alngned protein sequence.

    @param  aln_filename    path to alngnment file
    @return aln             numpy array of sequences in numerical representation.
    """

    with open(aln_filename) as aln_file:
        for i, seq in enumerate(aln_file):
            pass
    
    # number of sequences
    B = i+1
    
    with open(aln_filename) as aln_file:
        seq = aln_file.readline().strip()

    # sequence length
    N = len(seq)

    aln = np.zeros((B,N))
    aa_dict = {'R':1, 'H':2, 'K':3, 'D':4, 'E':5, 'S':6, 'T':7,\
            'N':8, 'Q':9, 'C':10, 'G':11, 'P':12, 'A':13, 'I':14,\
            'L':15, 'M':16, 'F':17, 'W':18, 'Y':19, 'V':20, '-':21}

    with open(aln_filename) as aln_file:
        for i, seq in enumerate(aln_file):
            seq = seq.strip()
            aln[i,:] = [aa_dict[aa] for aa in seq]
        
    return aln


def seq_identity(seq1, seq2):
    """ Calculate sequence identity between two sequences.

    @param  seq1    numpy array of first sequence
    @param  seq2    numpy array of second sequence
    @return ident   fraction of identical positions
    """
    
    # ignore gaps
    seq1_m = np.ma.masked_where(seq1 == 21, seq1)
    seq2_m = np.ma.masked_where(seq2 == 21, seq2)

    # calc identity
    diff = seq2_m - seq1_m
    N = np.ma.count(diff)
    n_diff = np.size(np.ma.nonzero(diff)[0], 0)
    ident = 1 - (n_diff / float(N))

    return ident



def get_per_residue_numseq(aln, threshold):
    """ Calculate number of similar sequences for each column in the
    alngnment. Similarity is defined as percentage of identical
    residues. Sequences less similar than threshold are counted. 
    
    @param  aln             numpy array of the alignment
    @return numseq_lst      per residue sequence counts
    """
    numseq_lst = []

    B,N = aln.shape
    seq0 = aln[0,:]
    ident_lst = np.empty(B)

    for i in xrange(1, B):
        ident_lst[i] = seq_identity(seq0, aln[i,:])

    mask = np.empty((B,N), dtype=bool)
    mask[:,:] = (ident_lst > threshold)[:,np.newaxis]
    
    #TODO: apply mask to get sequence counts
    print mask


if __name__ == "__main__":

    aln_filename = sys.argv[1]
    threshold = float(sys.argv[2])
    aln = read_aln(aln_filename)
    get_per_residue_numseq(aln, threshold)
