'''
    This code interfaces the python script for exon alignment
    with the local database constructed using Ensembl perl api
'''


from data_structures import Transcript_From_DB
import sqlite3
from tqdm import tqdm
from exon_similarity import process_transcript
import os
from recursive_backtracking import process_transcript_recursive
from dynamic_programming import dp_wrapper

db_path = os.path.join("..","data_acquistion","data","db.db")       # path to the local database
gene_list_path = "gene_ids_from_publication.txt"  # path to the 


def load_from_db():

    conn = sqlite3.connect(db_path)    # connecting to the local database
    c = conn.cursor()

    # extracting information from the database
    c.execute('SELECT human_gene_id, human_gene_name, mouse_gene_id FROM OrthologPairs')
    all_rows = c.fetchall()

    id_mapper = {}        # map gene name to stable ids
    ortholog_mapper = {}  # map human gene id to ortholog mouse gene id

    for i in range(len(all_rows)):

        id_mapper[all_rows[i][1]] = all_rows[i][0]
        ortholog_mapper[all_rows[i][0]] = all_rows[i][2]

    conn.close()

    return id_mapper,ortholog_mapper

def process_all_gene(id_mapper, gene_list, ortholog_mapper):

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    try:
        os.makedirs('results')   # create a directory to store the results
    except:
        pass 

    for gene in tqdm(gene_list,total=len(gene_list)):

        if(gene not in id_mapper):
            tqdm.write('Not Found')
            continue

        try:
            os.makedirs(os.path.join('results',gene+':'+id_mapper[gene]))
        except:
            pass


        c.execute("SELECT transcript_id FROM GeneTranscripts WHERE gene_id='" + id_mapper[gene] + "' AND Biotype='protein_coding'")
        all_human_transcipts = c.fetchall()

        c.execute("SELECT transcript_id FROM GeneTranscripts WHERE gene_id='" + ortholog_mapper[id_mapper[gene]] + "' AND Biotype='protein_coding'")
        all_mouse_transcipts = c.fetchall()

        for i in range(len(all_human_transcipts)):

            all_human_transcipts[i] = all_human_transcipts[i][0]
 
        for i in range(len(all_mouse_transcipts)):

            all_mouse_transcipts[i] = all_mouse_transcipts[i][0]
 
        
        if( (len(all_human_transcipts)==0) or (len(all_mouse_transcipts)==0) ):
            continue

        ortholog_transcripts = [ Transcript_From_DB(db_path,mouse_transcript) for mouse_transcript in all_mouse_transcipts ]

        for transcript in tqdm(all_human_transcipts, total=len(all_human_transcipts)):

            query_transcript = Transcript_From_DB(db_path, transcript)

            out_str = process_transcript(query_transcript, ortholog_transcripts)

            fp = open(os.path.join('results',gene+':'+id_mapper[gene],transcript+'.txt'),'w')
            fp.write(out_str)
            fp.close()

def create_report():
    
    id_mapper,ortholog_mapper = load_from_db()

    fp = open(gene_list_path,'r')
    gene_list = fp.read().split('\n')[:-1]
    fp.close()

    tot = 0

    for i in gene_list:

        if i not in id_mapper:
            
            tot += 1
    

    process_all_gene(id_mapper, gene_list, ortholog_mapper)

def process_ortholog_pairs(transcript_id_1, transcript_id_2):

    query_transcript = Transcript_From_DB(db_path, transcript_id_1)
    ortholog_transcript = Transcript_From_DB(db_path, transcript_id_2)

    params = {}
    params["match_score"] = 1
    params["mismatch_penalty"] = 0
    params["gap_open"] = 0
    params["gap_extend"] = 0
    params["skip_penalty"]= -1


    return(dp_wrapper(query_transcript,ortholog_transcript, params))

def get_transcript(transcript_id):

    return Transcript_From_DB(db_path, transcript_id)



def main():

    process_ortholog_pairs('ENST00000301894','ENSMUST00000113458')
    

    

if __name__ == '__main__':
    main()
