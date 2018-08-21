"""
    Exon pairing using dynamic programming
"""


from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from weighted_alignment import weighted_alignment_wrapper



def compute_similarity(seq1, seq2, index1, index2, params, similarity_memo, weight_mode):
    """
    Compute similarity score of the two exons
    
    Arguments:
        seq1 {str} -- exon 1 sequence
        seq2 {str} -- exon 2 sequence
        index1 {int} -- index1
        index2 {int} -- index2
        params {dict} -- alignment algorithm parametrs
        similarity_memo {2d list} -- lookup table for storing the similarity values
        weight_mode {str} -- 'gaussian' or 'uniform
    
    Returns:
        float -- similarity score
    """

    # avoid recomputation
    if(similarity_memo[index1][index2] != None ):

        return similarity_memo[index1][index2]

    # uniform alignment
    if(weight_mode=='uniform'):        

        if(len(seq1)==0) and (len(seq2)==0):            # special case

            similarity_memo[index1][index2] = 0.0

        elif(len(seq1)==0):                             # special case

            similarity_memo[index1][index2] = round(params["gap_open"] + params["gap_extend"]*(len(seq2)-1),2)
            
        elif(len(seq2)==0):                             # special case

            similarity_memo[index1][index2] = round(params["gap_open"] + params["gap_extend"]*(len(seq1)-1),2)
        
        
        else:
            
            similarity_memo[index1][index2] = round(pairwise2.align.globalms(seq1, seq2, params["match_score"], params["mismatch_penalty"], params["gap_open"], params["gap_extend"], score_only=True, one_alignment_only=True) / max(len(seq1), len(seq2),1) , 3 )
    
    # weighted alignment
    elif(weight_mode=='gaussian'):
                                                
        similarity_memo[index1][index2] = round(weighted_alignment_wrapper(seq1, seq2, params["match_score"], params["mismatch_penalty"], params["gap_open"], params["gap_extend"], score_only=True)  , 3 )
        
    
    return similarity_memo[index1][index2]


def dp(index1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode):
    """
    Exon pairing dynamic programming
    
    Arguments:
        index1 {int} -- index 1
        index2 {int} -- index 2
        query_transcript {Transcript} -- [description]
        ortholog_transcript {Transcript} -- [description]
        params {dict} -- alignment algorithm parameters
        memo {2d list} -- dp memo
        backtrack_memo {2d list} -- memo for reconstruction
        similarity_memo {2d list} -- memo for similarity scores
        weight_mode {str} -- 'gaussian' or 'uniform'
    
    Returns:
        [type] -- [description]
    """



    if(memo[index1][index2]!=None):

        return memo[index1][index2]

    elif(index1==len(query_transcript.exon_sequence) and (index2==len(ortholog_transcript.exon_sequence))):

        memo[index1][index2] = 0
        backtrack_memo[index1][index2] = (None, None)

        return memo[index1][index2]

    elif(index1==len(query_transcript.exon_sequence)):

        memo[index1][index2] = params["skip_penalty"] + dp(index1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)
        backtrack_memo[index1][index2] = (index1, index2+1)

        return memo[index1][index2]

    elif(index2==len(ortholog_transcript.exon_sequence)):

        memo[index1][index2] = params["skip_penalty"] + dp(index1+1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)
        backtrack_memo[index1][index2] = (index1+1, index2)

        return memo[index1][index2]
        
    # no skip        
    best = compute_similarity(query_transcript.exon_sequence[index1], ortholog_transcript.exon_sequence[index2], index1, index2, params, similarity_memo, weight_mode) + dp(index1+1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)
    
    backtrack_memo[index1][index2] = (index1+1, index2+1)

    #skip seq 1
    score = params["skip_penalty"] + dp(index1+1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)

    if(score>best):

        best = score 
        backtrack_memo[index1][index2] = (index1+1, index2)

    #skip seq 2
    score = params["skip_penalty"] + dp(index1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)

    if(score>best):

        best = score
        backtrack_memo[index1][index2] = (index1, index2+1)
    
    memo[index1][index2] = best
    return best


def dp_wrapper(query_transcript, ortholog_transcript, params, weight_mode):
    """
    Wrapper for the exon pairing dynamic programming
    
    Arguments:
        query_transcript {Transcript} -- [description]
        ortholog_transcript {Transcript} -- [description]
        params {dict} -- alignment algorithm parameters
        weight_mode {str} -- 'gaussian' or 'uniform'
    
    Returns:
        dict -- {"transcript1_total_exons",
                 "transcript2_total_exons",
                 "best_score",
                 "best_pairing"}
    """


    memo = []
    backtrack_memo = []
    similarity_memo = []




    for i in range(len(query_transcript.exon_sequence)+2):

        li = [None] * (len(ortholog_transcript.exon_sequence)+2)

        memo.append(li)

        li2 = [None] * (len(ortholog_transcript.exon_sequence)+2)

        backtrack_memo.append(li2)

        li3 = [None] * (len(ortholog_transcript.exon_sequence)+2)

        similarity_memo.append(li3)

    score = dp(0, 0, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo, weight_mode)
    
    normalization = len(query_transcript.exon_sequence)


    for seq in query_transcript.exon_sequence:

        if(len(seq)==0):

            normalization -= 1

    if(normalization!=0):
        score /= normalization
    else:
        score = 0

    pairing = backtrack(backtrack_memo)


    return {"transcript1_total_exons" : len(query_transcript.exon_sequence), 
            "transcript2_total_exons" : len(ortholog_transcript.exon_sequence),
            "best_score": round(score,3), 
            "best_pairing":pairing}


def backtrack(backtrack_memo):
    """
    Construct the exon pairings
    
    Arguments:
        backtrack_memo {2d list} -- backtrack memo for reconstruction
    
    Returns:
        list -- resulting alignment
    """

    pairing = []


    oldIndex1 = 0
    oldIndex2 = 0

    while(True):

        newIndex1, newIndex2 = backtrack_memo[oldIndex1][oldIndex2]
        

        if(newIndex1 ==None and newIndex2==None):

            break

        if(newIndex1!=oldIndex1 and newIndex2!=oldIndex2):

            pairing.append({ 'transcript1_exon_index' : oldIndex1, 
                             'transcript2_exon_index' : oldIndex2})

        oldIndex1 = newIndex1
        oldIndex2 = newIndex2

    return pairing

    