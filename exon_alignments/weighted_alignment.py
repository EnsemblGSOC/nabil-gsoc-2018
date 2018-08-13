"""
    Weighted Global Alignment

    Implementation of our proposed global alignment algorithm
    A modified Needleman Wunsch Algorithm
"""

import math
import sys
from subprocess import Popen, PIPE
import os


def weight_value(index1, index2, len1, len2, mode='gaussian'):
    """
    This function takes the coordinates and returns the weight score

    The weights are calculated from a statistical distribution motivated from 
    Gaussian distribution, the difference being it is the highest at the terminals
    and the lowest at the middle

    Please refer to the manuscript for more idea
    
    Arguments:
        index1 {int} -- current index of the 1st sequence
        index2 {int} -- current index of the 1st sequence
        len1 {int} -- length of 1st sequence
        len2 {int} -- length of 2nd sequence
    
    Keyword Arguments:
        mode {str} -- weight mode (default: {'gaussian'})
    
    Returns:
        float -- weight at the particular position
    """


    if(mode=='gaussian'):
        # nearer to the left terminal ?

        x1_left = index1
        x2_left = index2

        x_left = min(x1_left, x2_left)


        # nearer to the right terminal ?

        x1_right = len1 - index1 -1 
        x2_right = len2 - index2 -1 

        x_right = min(x1_right, x2_right)

        # pick the smaller one
        x = min(x1_left,x1_right)

        # normalize
        x = x / (max(len1,len2)/2)

        # obtain weight from a modified gaussian function

        miu = 0.0
        sigma = 0.4

        weight = math.e**(-(x-miu)**2/(2*sigma*sigma))


    elif(mode=='uniform'):

        return 1

    return weight


def dp(index1, index2, continuing_gap, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo,weight_mode='gaussian'):
    """
    Weighted Global Alignment 

    Simply put, we have implementd the Needleman Wunsch Algorithm, with putting higher weight
    along the splice site

    dynamic programming states : 3 (index1, index2, continuing_gap)
    
    Arguments:
        index1 {int} -- index of seq1
        index2 {int} -- index of seq2
        continuing_gap {int} -- are we following a gap ? 1-> True, 0 -> False
        seq1 {str} -- exon sequence 1
        seq2 {str} -- exon sequence 2
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_start {float} -- gap open penalty
        gap_extend {float} -- gap extend penalty
        memo {2d list} -- dp memo
        backtrack_memo {2d list} -- memo used for alignment construction
    
    Keyword Arguments:
        weight_mode {str} -- [description] (default: {'gaussian'})
    
    Returns:
        float -- optimal score
    """

                                    # avoiding recomputation
    if(memo[index1][index2][continuing_gap]!=None):

        return memo[index1][index2][continuing_gap]

    if(index1 == len(seq1)) and (index2 == len(seq2)):      # terminating condition, both seq1 and seq2 ended
        
        backtrack_memo[len(seq1)][len(seq2)][continuing_gap] = 'X'
        return 0

    elif (index1 == len(seq1)) and (index2 != len(seq2)):       # special case, seq1 ended
                                                                                            
        maxx = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * (gap_extend if continuing_gap==1 else gap_start) + dp(index1, index2+1, 1, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)
        backtrack_memo[index1][index2][continuing_gap] = 'V'

    elif (index1 != len(seq1)) and (index2 == len(seq2)):           # special case, seq2 ended

        maxx = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * (gap_extend if continuing_gap==1 else gap_start) + dp(index1+1, index2, 1, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)
        backtrack_memo[index1][index2][continuing_gap] = 'H'

    else:           

        maxx = None         #initialization

        # match 

        score1 = None
        if(seq1[index1] == seq2[index2]):
            score1 = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * match_score + dp(index1+1, index2+1, 0, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)

        # mismatch
        
        score2 = None
        if(seq1[index1] != seq2[index2]):
            score2 = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * mismatch_penalty + dp(index1+1, index2+1, 0, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)

        # gap

        # add gap in seq1
        score3 = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * (gap_extend if continuing_gap==1 else gap_start) + dp(index1, index2+1, 1, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)

        # add gap in seq2
        score4 = weight_value(index1, index2, len(seq1), len(seq2), weight_mode) * (gap_extend if continuing_gap==1 else gap_start) + dp(index1+1, index2, 1, seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)


        if(maxx == None and score1 != None):        
            maxx = score1                       
            backtrack_memo[index1][index2][continuing_gap] = 'D'        # traverse diagonally

        if(maxx == None and score2 != None) or (score2!= None and score2>maxx) :
            maxx = score2   
            backtrack_memo[index1][index2][continuing_gap] = 'D'        # traverse diagonally

        if(maxx == None and score3 != None) or (score3!= None and score3>maxx) :
            maxx = score3
            backtrack_memo[index1][index2][continuing_gap] = 'V'        # traverse vertically

        if(maxx == None and score4 != None) or (score4!= None and score4>maxx) :
            maxx = score4   
            backtrack_memo[index1][index2][continuing_gap] = 'H'        # traverse horizontally

    memo[index1][index2][continuing_gap] = maxx         # updating memo

    return maxx


def get_alignment(backtrack_memo, seq1, seq2):
    """
    Construct alignment from backtrack memo
    
    Arguments:
        backtrack_memo {2d list} -- backtracking memo
        seq1 {str} -- exon seq 1
        seq2 {str} -- exon seq 2
    
    Returns:
        tuple -- a tuple containing the aligned sequences
    """

                # initialization
    line1 = ""
    line2 = ""

    index1 = 0
    index2 = 0
    gap_extend = 0

    while(backtrack_memo[index1][index2][gap_extend] != 'X'):   # until termination met

        if(backtrack_memo[index1][index2][gap_extend] == 'D'):      # move diagonally

            line1 += seq1[index1]
            line2 += seq2[index2]

            index1 += 1 
            index2 += 1

            gap_extend = 0

        elif(backtrack_memo[index1][index2][gap_extend] == 'H'):    # move horizontally, skip along seq2, gap introduced

            line1 += seq1[index1]
            line2 += '-'

            index1 += 1 

            gap_extend = 1
        
        elif(backtrack_memo[index1][index2][gap_extend] == 'V'):    # move vertically, skip along seq1, gap introduced

            line1 += '-'
            line2 += seq2[index2]

            index2 += 1 

            gap_extend = 1
        
    return (line1,line2)


def weighted_needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, weight_mode, score_only):
    """
    Python backend for Weighted Needleman Wunsch
    
    Arguments:
        seq1 {str} -- exon sequence 1
        seq2 {str} -- exon sequence 2 
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_start {float} -- gap open penalty
        gap_extend {float} -- gap extend penalty
        weight_mode {str} -- gaussian or uniform
        score_only {bool} -- return only score ? otherwise return alignment with score
    
    Returns:
        float or tuple -- returns either score or both the score and the alignment
    """

                    # initializing the memo
    memo = []    

    for i in range(len(seq1)+1):
        li = []
        for j in range(len(seq2)+1):
            li.append([None]*2)
        
        memo.append(li)
                                    # increasing recursion limit
    sys.setrecursionlimit(len(memo)*len(memo[0])*len(memo[0][0]))

    backtrack_memo = []     # initializing the backtrack memo

    for i in range(len(seq1)+1):
        li = []
        for j in range(len(seq2)+1):
            li.append([None]*2)
        
        backtrack_memo.append(li)

                            # obtaining score from dp
    score = dp(0,0,0,seq1, seq2, match_score, mismatch_penalty, gap_start, gap_extend, memo, backtrack_memo, weight_mode)

    if(score_only==True):       # return only score

        return score
    
    else:                       # return both score and alignment
        
        alignment = get_alignment(backtrack_memo, seq1, seq2)

        return (score, alignment)


def weighted_alignment_wrapper(exon1_seq, exon2_seq, match_score, mismatch_penalty, gap_start, gap_extend, score_only, backend='C++'):
    """
    Wrapper for weighted global alignment
    
    Arguments:
        exon1_seq {str} -- exon sequence 1
        exon2_seq {str} -- exon sequence 1
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_start {float} -- gap open penalty
        gap_extend {float} -- gap extend penalty        
        score_only {bool} -- return only score ? otherwise return alignment with score
    
    Keyword Arguments:
        backend {str} -- which backend to use C++ or Python (default: {'C++'})
    
    Returns:
        float or tuple -- returns either score or both the score and the alignment
    """


    if(backend=='C++'):     # use C++
                                                # start sub process
        process = Popen([os.path.join('..','weighted_alignment','exe_file'), exon1_seq, exon2_seq, str(match_score),str(mismatch_penalty), str(gap_start), str(gap_extend), str(score_only)], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()      # obtain the response
        stdout = stdout.decode("utf-8")             # convert the byte to string
        stdout = stdout.split('\n')                 # following the format of returned data from C++ code
        
        if(score_only):         # return only score
        
            score = float(stdout[0])
            return score

        else:                   # return both score and alignment

            score = float(stdout[0])
            alignment = [stdout[1],stdout[2]]

            return (score, alignment)

    elif(backend=='Python'):            # use Python

        out = weighted_needleman_wunsch(exon1_seq, exon2_seq, match_score, mismatch_penalty, gap_start, gap_extend, 'gaussian', score_only)

        return out


def main():
    # example code

    seq1 = 'CATGACTGTTATCTTTTTGATGGAAAAGATGAAGTTCAGCGAAGAACAAATCTACAACTTGGCCAGAAGACTGTAGACTTTGCCTTCAAAGAG'
    seq2 = 'GAAAATGGACAAAGGAAGTATGGTGGCCCTCCACCAGGCTGGGATTCTACACCCCCAGAAAGGGGCTGCGAGATTTTCATTGGGAAACTTCCCCGGGACCTCTTTGAGGATGAACTCATACCATTGTGTGAAAAA'
    
    print(weighted_needleman_wunsch(seq1, seq2, 1, 0, 0, 0, 'gaussian',False))
    

if __name__ == '__main__':
    main()

