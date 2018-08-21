"""

    Backend code for the Flask server

"""


import sqlite3
import os 
import time 
import sys

sys.path.append(os.path.join('..','exon_alignments'))          # adding the exon_alignment codes to path

from db_interface import get_transcript, process_ortholog_pairs
from muscle_interface import protein_similarity_wrapper, transcript_similarity_wrapper
from biopython_alignment_interface import protein_similarity_wrapper_biopython, transcript_similarity_wrapper_biopython
from blosum_matrices import get_blosum_scores
from weighted_alignment import weighted_alignment_wrapper

global db_path
db_path = os.path.join("..","data_acquisition","data","db.db")


def load_from_db(data_type, species, gene_id, protein_coding, gencode_basic):
    """
    Load data from database
    
    Arguments:
        data_type {str} -- data type consists of gene_name, gene_id, transcript_name, transcript_id
        species {str} -- name of the species 'human' or 'mouse'
        gene_id {str} -- gene id of transcript ( if we're interested with extracting transcript_name or transcript_id)
        protein_coding {str} -- 'true' or 'false'
        gencode_basic {str} -- 'true' or 'false'
    
    Returns:
        dict -- {"data" : list of extracted data}
    """
    global db_path
    conn = sqlite3.connect(db_path)         # connecting to the database
    c = conn.cursor()

    out = []        # list to store the data

    if( (data_type=='gene_name') or (data_type=='gene_id') ) :      # extracting gene_id or gene_name

        c.execute("SELECT DISTINCT "+data_type+" FROM Genes WHERE Species='"+species+"' ORDER BY "+data_type)

        all_data = c.fetchall()

        for i in range(len(all_data)):
            
            out.append(all_data[i][0])


    else:                   # extracting transcript_id or transcrip_name

        if(gene_id=='all'):

            if(protein_coding=='true'):
                if(gencode_basic=='false'):
                    c.execute("SELECT DISTINCT "+ data_type +" FROM geneTranscripts WHERE biotype='protein_coding' AND Species ='"+species+"' ORDER BY "+data_type)

                else:
                    c.execute("SELECT DISTINCT "+ data_type +" FROM geneTranscripts WHERE biotype='protein_coding' AND gencode_basic='1' AND Species ='"+species+"' ORDER BY "+data_type)

        else:

            if(protein_coding=='true'):
                if(gencode_basic=='false'):
                    c.execute("SELECT DISTINCT "+ data_type +" FROM geneTranscripts WHERE biotype='protein_coding' AND gene_id='"+gene_id+"'"+ " ORDER BY "+data_type)

                else:
                    c.execute("SELECT DISTINCT "+ data_type +" FROM geneTranscripts WHERE biotype='protein_coding' AND gencode_basic='1'AND gene_id='"+gene_id+"' ORDER BY "+data_type)

        all_data = c.fetchall()

        for i in range(len(all_data)):

            out.append(all_data[i][0])
            

    conn.close()

    return {"data":out}


def extract_from_gene(gene_name,gene_id,species, protein_coding, gencode_basic):
    """
    Checks validity of a gene and collect the transcripts
    
    Arguments:
        gene_name {str} -- name of the gene
        gene_id {str} -- gene id
        species {str} -- 'human' or 'mouse'
        protein_coding {str} -- 'true' or 'false'
        gencode_basic {str} -- 'true' or 'false'
    
    Returns:
        dict -- {"valid",
                 "gene_id",
                 "gene_name",
                 "transcript_ids",
                 "transcript_names"}
    """

    
    global db_path
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    if(gene_id==None):
        # get gene_id
        c.execute("SELECT Gene_ID FROM Genes WHERE Gene_Name='"+gene_name+"' AND Species='"+species+"' ORDER BY Gene_ID")
        all_data = c.fetchall()

        if(len(all_data)==0):
            conn.close()
            return {"valid":False}
        
        gene_id = all_data[0][0]

    else:
        # get gene_name
        c.execute("SELECT Gene_Name FROM Genes WHERE Gene_ID='"+gene_id+"' AND Species='"+species+"'")
        all_data = c.fetchall()

        if(len(all_data)==0):
            conn.close()
            return {"valid":False}
        
        gene_name = all_data[0][0]

    transcript_ids = []
    transcript_names = []

    if(protein_coding=='true'):
        
        if(gencode_basic=='false'):
                    
            c.execute("SELECT Transcript_ID, Transcript_Name FROM geneTranscripts WHERE biotype='protein_coding' AND gene_id='"+gene_id+"'"+ " ORDER BY Transcript_ID")

        else:
            c.execute("SELECT Transcript_ID, Transcript_Name FROM geneTranscripts WHERE biotype='protein_coding' AND gencode_basic='1'AND gene_id='"+gene_id+"'"+ " ORDER BY Transcript_ID")
   

        all_data = c.fetchall()

        for i in range(len(all_data)):

            transcript_ids.append(all_data[i][0])
            transcript_names.append(all_data[i][1])

    conn.close()

    return {
            "valid" : True,
            "gene_id" : gene_id,
            "gene_name" : gene_name,
            "transcript_ids" : transcript_ids,
            "transcript_names" : transcript_names
            }


def extract_from_transcript(transcript_name, transcript_id):
    """
    Checks validity of a transcript and also extract the related information
    
    Arguments:
        transcript_name {str} -- name of the transcript
        transcript_id {str} -- transcript id
    
    Returns:
        dict -- {"valid",
                 "gene_id",
                 "gene_name",
                 "transcript_id",
                 "transcript_name"}
    """
        
    global db_path
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    gene_id = ''
    gene_name = ''

    if(transcript_id==None):
        # get transcript_id
        c.execute("SELECT Transcript_ID, Gene_ID FROM GeneTranscripts WHERE Transcript_Name='"+transcript_name+"'")
        all_data = c.fetchall()

        if(len(all_data)==0):
            conn.close()
            return {"valid":False}
        
        transcript_id = all_data[0][0]
        gene_id = all_data[0][1]

    else:
        # get transcript_name
        c.execute("SELECT Transcript_Name, Gene_ID FROM GeneTranscripts WHERE Transcript_ID='"+transcript_id+"'")
        all_data = c.fetchall()

        if(len(all_data)==0):
            conn.close()
            return {"valid":False}
        
        transcript_name = all_data[0][0]
        gene_id = all_data[0][1]

    
    c.execute("SELECT Gene_Name FROM Genes WHERE Gene_ID='"+gene_id+"'")
    all_data = c.fetchall()
    gene_name = all_data[0][0]

    conn.close()

    return {
            "valid" : True,
            "gene_id" : gene_id,
            "gene_name" : gene_name,
            "transcript_id" : transcript_id,
            "transcript_name" : transcript_name
            }


def get_scoring_metrics(match_score, mismatch_penalty, gap_open, blosum_mat):
    """
    Obtain the scoring metrics
    
    Arguments:
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_open {float} -- gap open penalty
        blosum_mat {str} -- blosum matrix name
    
    Returns:
        dict -- {   "transcript_scores" : {dict} -> transcript scores,
                    "protein_scores" : {dict} -> protein_scores,
                    "transcript_scores_max" : {float} transcript_scores_max,
                    "transcript_scores_min" : {float} transcript_scores_min,
                    "protein_scores_max" : {float} protein_scores_max,
                    "protein_scores_min" : {float} protein_scores_min
                }
    """

                                # initializing the transcript scores
    transcript_scores = {}
    transcript_scores['-']={}

    for char1 in ['A','T','C','G']:    

        transcript_scores[char1]={}

        for char2 in ['A','T','C','G','-']:

            if(char1==char2):       # match
                
                transcript_scores[char1][char2] = match_score
            
            elif(char2=='-'):       # gap penalty

                transcript_scores[char1][char2] = gap_open
                transcript_scores[char2][char1] = gap_open

            else:                   # mismatch

                transcript_scores[char1][char2] = mismatch_penalty
    
    transcript_scores_max = max(match_score, mismatch_penalty, gap_open)        # used for colormap
    transcript_scores_min = min(match_score, mismatch_penalty, gap_open)

    protein_scores, protein_scores_max, protein_scores_min = get_blosum_scores(blosum_mat, gap_open)    # obtain blosum scores

    return { "transcript_scores" : transcript_scores,
             "protein_scores" : protein_scores,
             "transcript_scores_max" : transcript_scores_max,
             "transcript_scores_min" : transcript_scores_min,
             "protein_scores_max" : protein_scores_max,
             "protein_scores_min" : protein_scores_min
            }
    

def load_orthologs(transcript_id, protein_coding, gencode_basic):
    """
    
    Load all orthologs of the transcript
    
    Arguments:
        transcript_id {str} -- id of the transcript
        protein_coding {str} -- 'true' or 'false'
        gencode_basic {str} -- 'true' or 'false'
    
    Returns:
        dict -- { "data":ortholog_transcripts -> {dict} {"transcript_id", 
                                                         "biotype", 
                                                         "transcript_name", 
                                                         "gencode_basic", 
                                                         "species"}

    """
                                            # connecting to the database
    global db_path
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute("SELECT gene_id FROM geneTranscripts WHERE transcript_id='"+transcript_id+"'")    # load the gene id

    row = c.fetchall()

    gene_id = row[0][0]

    c.execute("SELECT human_gene_id, mouse_gene_id FROM orthologPairs WHERE human_gene_id='"+gene_id+"' OR mouse_gene_id='"+gene_id+"'")      # load ortholog gene id

    row = c.fetchall()

    ortholog_gene_id = row[0][0] if row[0][1]==gene_id else row[0][1]

    # filtering 
    if(protein_coding=='true'):
    
        if(gencode_basic=='false'):
            
            c.execute("SELECT Transcript_ID, BioType, transcript_name, gencode_basic, species FROM geneTranscripts WHERE gene_id='"+ortholog_gene_id+"' AND biotype='protein_coding' ORDER BY Transcript_ID")

        else:

            c.execute("SELECT Transcript_ID, BioType, transcript_name, gencode_basic, species FROM geneTranscripts WHERE gene_id='"+ortholog_gene_id+"' AND  biotype='protein_coding' AND gencode_basic='1' ORDER BY Transcript_ID")


    rows = c.fetchall()

    ortholog_transcripts = []   # list to store orthologs

    for row in rows:

        ortholog_transcripts.append({
                                        "transcript_id" : row[0],
                                        "biotype" : row[1],
                                        "transcript_name":row[2],                                        
                                        "gencode_basic": 'GENCODE basic' if row[3]=='1' else '',
                                        "species" : row[4]
                                    })

    conn.close()        # closing the connection

    return {"data":ortholog_transcripts}


def load_exon_sequence(transcript_id, exon_number):
    """
    Load the sequence of an Exon specified by its number and transcript
    
    Arguments:
        transcript_id {str} -- id of the transcript
        exon_number {str} -- number of the exon
    
    Returns:
        dict -- {   "transcript_id": {str} transcript_id,
                    "exon_id": {str} exon id
                    "exon_sequence": {str} coding sequence
                    "utr_sequence": {str} UTR sequence
                    "utr_pos": {int} -1 means UTR to the left, 1 means UTR to the right
                }
    """


    transcript = get_transcript(transcript_id)

    out = {"transcript_id":transcript_id,
           "exon_id":transcript.exon_id[int(exon_number)],
           "exon_sequence":transcript.exon_sequence[int(exon_number)],
           "utr_sequence":transcript.utr_sequence[int(exon_number)],
           "utr_pos":transcript.utr_pos[int(exon_number)]
            }

    return out


def make_pairs(transcript1_id,transcript2_id, weight_mode, match_score, mismatch_penalty, gap_open, gap_extend, skip_penalty):
    """
    Pair exons
    
    Arguments:
        transcript1_id {str} -- id of 1st transcript
        transcript2_id {str} -- id of 2nd transcript
        weight_mode {str} -- weight mode 'gaussian' or 'uniform'
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_open {float} -- gap open penalty
        gap_extend {float} -- gap extend penalty
        skip_penalty {float} -- exon skip penalty
    
    Returns:
        [type] -- [description]
    """

    
    return process_ortholog_pairs(transcript1_id, transcript2_id, weight_mode, match_score, mismatch_penalty, gap_open, gap_extend, skip_penalty)


def get_splice_site_info(transcript_id, exon_id):
    """
    Obtain splice cite coordinates
    
    Arguments:
        transcript_id {str} -- transcript id
        exon_id {str} -- exon id
    
    Returns:
        dict -- {"absolute_start",
                 "absolute_end",
                 "relative_start",
                 "relative_end",
                 "strand"}
    """

    global db_path
    conn = sqlite3.connect(db_path)         # connecting to database
    c = conn.cursor()

    table = 'Transcripts' + transcript_id[-3:]      # table name

    c.execute("SELECT start, end, strand FROM "+table+" WHERE Transcript_ID='"+transcript_id+"' AND Exon_ID='"+exon_id+"'")

    row = c.fetchall()
                            # extract all the information
    start = row[0][0]
    end = row[0][1]
    strand = row[0][2]

    c.execute("SELECT start FROM "+table+" WHERE Transcript_ID='"+transcript_id+"' AND Exon_No='1'")

    row = c.fetchall()

    exon1_start = row[0][0]

    absolute_start = start
    absolute_end = end
    relative_start = int(start) - int(exon1_start)
    relative_end = int(end) - int(exon1_start)

    conn.close()

    return { "absolute_start" : absolute_start,
             "absolute_end" : absolute_end,
             "relative_start" : relative_start,
             "relative_end" : relative_end,
             "strand" : strand }

                
def get_protein_similarity(transcript1_id, transcript2_id, exon1_id, exon2_id):
    """
    Obtains the protein alignment using MUSCLE
    
    Arguments:
        transcript1_id {str} -- id of the 1st transcript
        transcript2_id {str} -- id of the 2nd transcript
        exon1_id {str} -- id of the 1st exon
        exon2_id {str} -- id of the 2nd exon
    Returns:
        list -- aligned sequences
    """

    # retrieve the transcripts
    transcript_1 = get_transcript(transcript1_id)
    transcript_2 = get_transcript(transcript2_id)
    
    exon1_seq = ''
    exon2_seq = ''

    # O(N) is good enough for small values of N
    for i in range(len(transcript_1.exon_id)):

        if(transcript_1.exon_id[i]==exon1_id):

            exon1_seq = transcript_1.exon_sequence[i]

    for i in range(len(transcript_2.exon_id)):

        if(transcript_2.exon_id[i]==exon2_id):

            exon2_seq = transcript_2.exon_sequence[i]

    global db_path
    conn = sqlite3.connect(db_path)     # connecting to the database
    c = conn.cursor()

    table = 'Transcripts' + transcript1_id[-3:]

    c.execute("SELECT  start_phase, end_phase FROM "+table+" WHERE Transcript_ID='"+transcript1_id+"' AND Exon_ID='"+exon1_id+"'")

    rows = c.fetchall()

    # integrating the phase information
    exon1_start_phase = int(rows[0][0])
    exon1_end_phase = int(rows[0][1])

    if(exon1_start_phase==-1 and exon1_end_phase==-1):
        return {'error':'Exon '+str(exon1_id)+' is non-coding'}     # error, non coding
    
    if(exon1_start_phase==-1):              # please see the ensembl api doc for better understanding
        exon1_start_phase = 0
    else:
        exon1_start_phase = (3-exon1_start_phase) % 3
    
    if(exon1_end_phase==-1):              # please see the ensembl api doc for better understanding
        exon1_end_phase = 0
    else:
        exon1_end_phase = (3-exon1_end_phase) % 3
    

    exon1_seq = exon1_seq[exon1_start_phase:len(exon1_seq)-exon1_end_phase]

    table = 'Transcripts' + transcript2_id[-3:]

    c.execute("SELECT start_phase, end_phase FROM "+table+" WHERE Transcript_ID='"+transcript2_id+"' AND Exon_ID='"+exon2_id+"'")

    rows = c.fetchall()

    exon2_start_phase = int(rows[0][0])
    exon2_end_phase = int(rows[0][1])

    if(exon2_start_phase==-1 and exon2_end_phase==-1):
        return {'error':'Exon '+str(exon2_id)+' is non-coding'}     # error, non coding


    if(exon2_start_phase==-1):              # please see the ensembl api doc for better understanding
        exon2_start_phase = 0
    else:
        exon2_start_phase = (3-exon2_start_phase) % 3
    
    if(exon2_end_phase==-1):              # please see the ensembl api doc for better understanding
        exon2_end_phase = 0
    else:
        exon2_end_phase = (3-exon2_end_phase) % 3

    exon2_seq = exon2_seq[exon2_start_phase:len(exon2_seq)-exon2_end_phase] 

    conn.close()

    muscle_err = 'None'

    try:
        alignment = protein_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id)
    
    except:

        alignment = protein_similarity_wrapper_biopython(exon1_seq, exon1_id, exon2_seq, exon2_id)
        muscle_err = 'error'

    return {'alignment' : alignment,
            'muscle_error' : muscle_err }


def get_transcript_similarity(transcript1_id, transcript2_id, exon1_id, exon2_id, match_score, mismatch_penalty, gap_start, gap_extend, weight_mode):
    """
    Obtains the transcript alignment and score using MUSCLE or weighted alignment algorithm
    
    Arguments:
        transcript1_id {str} -- id of the 1st transcript
        transcript2_id {str} -- id of the 2nd transcript
        exon1_id {str} -- id of the 1st exon
        exon2_id {str} -- id of the 2nd exon
        match_score {float} -- match score
        mismatch_penalty {float} -- mismatch penalty
        gap_start {float} -- gap open penalty
        gap_extend {float} -- gap continue penalty
        weight_mode {[type]} -- 'gaussian' or 'uniform'
    
    Returns:
        dict -- {"score" {float} -> alignment score, 
                 "alignment" {list of str} -> aligned sequences }
    """

    # obtain the transcripts
    transcript_1 = get_transcript(transcript1_id)
    transcript_2 = get_transcript(transcript2_id)
    
    exon1_seq = ''
    exon2_seq = ''

    # O(n) is good enough since n is too small
    for i in range(len(transcript_1.exon_id)):

        if(transcript_1.exon_id[i]==exon1_id):

            exon1_seq = transcript_1.exon_sequence[i]

    for i in range(len(transcript_2.exon_id)):

        if(transcript_2.exon_id[i]==exon2_id):

            exon2_seq = transcript_2.exon_sequence[i]

    if(weight_mode=='uniform'):         # using MUSCLE

        muscle_err = 'None'

        alignment = ''
        score = 0
        gap_continue = 0

        if(len(exon1_seq)==0) and (len(exon2_seq)==0):      # special case

            alignment = ['','']

        elif (len(exon1_seq)==0):                             # special case

            alignment = ['-'*len(exon2_seq),exon2_seq]

        
        elif (len(exon2_seq)==0):                            # special case   

            alignment = [exon1_seq,'-'*len(exon1_seq)]

        else:                           # MUSCLE interface
            
            try:
                alignment = transcript_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id)

            except:
                alignment = transcript_similarity_wrapper_biopython(exon1_seq, exon1_id, exon2_seq, exon2_id, match_score, mismatch_penalty, gap_start, gap_extend)
                muscle_err = 'error'

        # calculating the score
        for i in range(len(alignment[0])):

            if(alignment[0][i]=='-' or alignment[1][i]=='-'):

                score += gap_extend if gap_continue==1 else gap_start
                gap_continue = 1

            elif(alignment[0][i]==alignment[1][i]):

                score += match_score
                gap_continue = 0

            else:

                score += mismatch_penalty
                gap_continue = 0
        
        return {'score' : score/(len(alignment[0])*max(match_score, mismatch_penalty, gap_extend, gap_start)),
                'alignment' : alignment,
                'muscle_error' : muscle_err }
    
    else:               # use Weight Alignment Algorithm
        
        (score,alignment) = weighted_alignment_wrapper(exon1_seq, exon2_seq, match_score, mismatch_penalty, gap_start, gap_extend, score_only=False, backend='C++')

        return {'score' : score,
                'alignment' : alignment}
                

def extract_species(transcript_id):
    """
    Check the species
    
    Arguments:
        transcript_id {str} -- Transcript ID
    
    Returns:
        dict -- {"species"}
    """ 
    global db_path
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute("SELECT species FROM GeneTranscripts WHERE Transcript_ID='"+transcript_id+"'")
    all_data = c.fetchall()

    species = all_data[0][0]

    if(species=='human'):

        species = 'Homo_sapiens'

    elif(species=='mouse'):

        species = 'Mus_musculus'
    
    conn.close()

    return {"species" : species}