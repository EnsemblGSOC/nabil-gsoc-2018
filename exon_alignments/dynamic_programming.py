from Bio.pairwise2 import format_alignment
from Bio import pairwise2


def compute_similarity(seq1, seq2, index1, index2, params, similarity_memo):

    if(similarity_memo[index1][index2] != None ):

        return similarity_memo[index1][index2]

    similarity_memo[index1][index2] = round(pairwise2.align.globalms(seq1, seq2, params["match_score"], params["mismatch_penalty"], params["gap_open"], params["gap_extend"], score_only=True, one_alignment_only=True) / max(len(seq1), len(seq2)) , 2 )

    return similarity_memo[index1][index2]

def dp(index1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo):

    if(memo[index1][index2]!=None):

        return memo[index1][index2]

    elif(index1==len(query_transcript.exon_sequence) and (index2==len(ortholog_transcript.exon_sequence))):

        memo[index1][index2] = 0
        backtrack_memo[index1][index2] = (None, None)

        return memo[index1][index2]

    elif(index1==len(query_transcript.exon_sequence)):

        memo[index1][index2] = params["skip_penalty"] + dp(index1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)
        backtrack_memo[index1][index2] = (index1, index2+1)

        return memo[index1][index2]

    elif(index2==len(ortholog_transcript.exon_sequence)):

        memo[index1][index2] = params["skip_penalty"] + dp(index1+1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)
        backtrack_memo[index1][index2] = (index1, index2+1)

        return memo[index1][index2]
        
    # no skip

    best = compute_similarity(query_transcript.exon_sequence[index1], ortholog_transcript.exon_sequence[index2], index1, index2, params, similarity_memo) + dp(index1+1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)
    backtrack_memo[index1][index2] = (index1+1, index2+1)

    #skip seq 1
    score = params["skip_penalty"] + dp(index1+1, index2, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)

    if(score>best):

        best = score 
        backtrack_memo[index1][index2] = (index1+1, index2)

    #skip seq 2
    score = params["skip_penalty"] + dp(index1, index2+1, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)

    if(score>best):

        best = score
        backtrack_memo[index1][index2] = (index1, index2+1)
    
    memo[index1][index2] = best
    return best


def dp_wrapper(query_transcript, ortholog_transcript, params):

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

    score = dp(0, 0, query_transcript, ortholog_transcript, params, memo, backtrack_memo, similarity_memo)

    score /= len(query_transcript.exon_sequence)

    pairing = backtrack(backtrack_memo)



    return {"transcript1_total_exons" : len(query_transcript.exon_sequence), 
            "transcript2_total_exons" : len(ortholog_transcript.exon_sequence),
            "best_score": score, 
            "best_pairing":pairing}

def backtrack(backtrack_memo):

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

    