"""
    

"""


from exon_similarity import exon_similarity

global pairing
global pairing_index
global query_transcript
global ortholog_transcript
global match_score
global mismatch_penalty
global gap_open
global gap_extend
global best_score
global best_pairing


def process_result():

    global best_score
    global best_pairing

    try:
        if(best_score * len(query_transcript.exon_sequence) > pairing_index):
            return
    except:
        pass 

    total = 0

    print(pairing[:pairing_index])

    for i in range(pairing_index):

        if(pairing[i][0]==None or pairing[i][1]==None):

            continue

        else:

            similarity_score = exon_similarity(query_transcript.exon_sequence[pairing[i][0]], ortholog_transcript.exon_sequence[pairing[i][1]], match_score, mismatch_penalty, gap_open, gap_extend)            
            total += similarity_score

    total /= len(query_transcript.exon_sequence)

    

    if((best_score == None) or (best_score < total)):

        best_score = total

        best_pairing = []

        for i in range(pairing_index):

            if(pairing[i][0]==None or pairing[i][1]==None):

                continue

            else:
                
                best_pairing.append({"transcript1_exon_index" : pairing[i][0], "transcript2_exon_index" : pairing[i][1]})

    print(total,best_score)

def recursive_backtracking(index1, index2, skips_left):

    global pairing
    global pairing_index

    if( (index1==len(query_transcript.exon_sequence)) or (index2==len(ortholog_transcript.exon_sequence)) ):

        process_result()

    else:
        
        # no skip
        pairing[pairing_index][0] = index1
        pairing[pairing_index][1] = index2
        pairing_index += 1
        
        recursive_backtracking(index1+1, index2+1, skips_left)
        
        pairing_index -= 1

        if(skips_left>0):

            # skip from query transcript

            recursive_backtracking(index1+1, index2, skips_left-1)

        if(skips_left>0):
            
            # skip from ortholog transcript

            recursive_backtracking(index1, index2+1, skips_left-1)

        


def process_transcript_recursive(transcript1, transcript2, custom_match_score = 1, custom_mismatch_penalty = 0, custom_gap_open = 0, custom_gap_extend = 0, skips_allowed= 2):

    global query_transcript
    global ortholog_transcript
    global match_score
    global mismatch_penalty
    global gap_open
    global gap_extend
    global best_score
    global pairing
    global pairing_index

    query_transcript = transcript1
    ortholog_transcript = transcript2
    match_score = custom_match_score
    mismatch_penalty = custom_mismatch_penalty
    gap_open = custom_gap_open
    gap_extend = custom_gap_extend
    
    best_score = None 
    pairing = []
    pairing_index = 0

    for i in range( max(len(query_transcript.exon_sequence), len(ortholog_transcript.exon_sequence) ) + skips_allowed):
        
        pairing.append([None,None])

    recursive_backtracking(0,0,skips_allowed)

    return {"transcript1_total_exons" : len(query_transcript.exon_sequence), 
            "transcript2_total_exons" : len(ortholog_transcript.exon_sequence),
            "best_score": best_score, 
            "best_pairing":best_pairing}

