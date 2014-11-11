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

    # count number of sequences, ignore gap-only sequences
    empty_seqs = []

    with open(aln_filename) as aln_file:
        for i, seq in enumerate(aln_file):
            if seq.count('-') == len(seq.strip()):
                empty_seqs.append(i)
    
    print empty_seqs
    # number of sequences
    B = i+1 - len(empty_seqs)
    
    with open(aln_filename) as aln_file:
        seq = aln_file.readline().strip()

    # sequence length
    N = len(seq)

    aln = np.zeros((B,N))
    aa_dict = {'R':1, 'H':2, 'K':3, 'D':4, 'E':5, 'S':6, 'T':7,\
            'N':8, 'Q':9, 'C':10, 'G':11, 'P':12, 'A':13, 'I':14,\
            'L':15, 'M':16, 'F':17, 'W':18, 'Y':19, 'V':20, '-':21}

    i_nonempty = 0
    with open(aln_filename) as aln_file:
        for i, seq in enumerate(aln_file):
            if i not in empty_seqs:
                seq = seq.strip()
                aln[i_nonempty,:] = [aa_dict[aa] for aa in seq]
                i_nonempty += 1
        
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

    # no identity to sequence with only gaps
    if N == 0:
        return 0.0

    n_diff = np.size(np.ma.nonzero(diff)[0], 0)
    ident = 1 - (n_diff / float(N))

    return ident



def get_per_residue_numseq(aln, threshold):
    """ Calculate number of similar sequences for each column in the
    alignment. Similarity is defined as percentage of identical
    residues. Sequences less similar than threshold are counted. 
    
    @param  aln             numpy array of the alignment
    @param  threshold       similarity threshold
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
    #print mask


def get_weights(aln, threshold):
    """ Calculate sequence weight based on similarity. Similarity is defined 
    as percentage of identical residues. Sequences less similar than 
    threshold are counted. 
    
    @param  aln             numpy array of the alignment
    @param  threshold       similarity threshold
    @return weights         sequence weights
    """

    B,N = aln.shape
    sim_mat = np.empty((B,B))
    aln_m = np.ma.masked_where(aln == 21, aln)

    for i in xrange(B):
        diff = aln_m - aln_m[i,:]
        #n_diff = np.zeros(B, dtype=float)
        n_diff = np.bincount(np.ma.nonzero(diff)[0])
        n_diff.resize(B)
        N = np.ma.count(diff, 1).astype(np.float)
        #print n_diff
        #print N
        #print 1 - (n_diff / N)
        sim_mat[i,:] = 1 - (n_diff / N)

    num_sim = np.sum(sim_mat >= threshold, 1)
    #print sim_mat[np.where(num_sim == 0),:]
    weights = 1./num_sim
    #print weights
    return weights



def get_beff(aln, threshold):
    """ Calculate number of effective sequences.     

    @param  aln             numpy array of the alignment
    @param  threshold       similarity threshold
    @return b_eff           number of effective sequences
    """

    weights = get_weights(aln, threshold)
    b_eff = sum(weights)
    return b_eff



if __name__ == "__main__":

    aln_filename = sys.argv[1]
    threshold = float(sys.argv[2])
    aln = read_aln(aln_filename)
    numseq_lst = get_per_residue_numseq(aln, threshold)

    b_eff = get_beff(aln, threshold)
    print b_eff
