import sqlite3
import os 
import time 
import sys

sys.path.append(os.path.join('..','exon_alignments'))

from db_interface import get_transcript, process_ortholog_pairs
from protein_similarity import protein_similarity_wrapper

db_path = os.path.join("..","data_acquistion","data","db.db")

def load_from_db(data_type, species, gene_id, protein_coding, gencode_basic):


    t1 = time.time()

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    out = []

    if( (data_type=='gene_name') or (data_type=='gene_id') ) :

        c.execute("SELECT DISTINCT "+data_type+" FROM Genes WHERE Species='"+species+"' ORDER BY "+data_type)

        all_data = c.fetchall()

        for i in range(len(all_data)):
            
            out.append(all_data[i][0])


    else:

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

    t2 = time.time()

    return {"data":out}


def load_orthologs(transcript_id, protein_coding, gencode_basic):

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute("SELECT gene_id FROM geneTranscripts WHERE transcript_id='"+transcript_id+"'")

    row = c.fetchall()

    gene_id = row[0][0]

    c.execute("SELECT human_gene_id, mouse_gene_id FROM orthologPairs WHERE human_gene_id='"+gene_id+"' OR mouse_gene_id='"+gene_id+"'")

    row = c.fetchall()

    ortholog_gene_id = row[0][0] if row[0][1]==gene_id else row[0][1]

    if(protein_coding=='true'):
    
        if(gencode_basic=='false'):
            
            c.execute("SELECT Transcript_ID, BioType, transcript_name, gencode_basic, species FROM geneTranscripts WHERE gene_id='"+ortholog_gene_id+"' AND biotype='protein_coding' ORDER BY Transcript_ID")

        else:

            c.execute("SELECT Transcript_ID, BioType, transcript_name, gencode_basic, species FROM geneTranscripts WHERE gene_id='"+ortholog_gene_id+"' AND  biotype='protein_coding' AND gencode_basic='1' ORDER BY Transcript_ID")


    rows = c.fetchall()

    ortholog_transcripts = []

    for row in rows:

        ortholog_transcripts.append({
                                        "transcript_id" : row[0],
                                        "biotype" : row[1],
                                        "transcript_name":row[2],                                        
                                        "gencode_basic": 'GENCODE basic' if row[3]=='1' else '',
                                        "species" : row[4]
                                    })

    conn.close()

    print(ortholog_transcripts)

    return {"data":ortholog_transcripts}


def load_exon_sequence(transcript_id, exon_number):

    transcript = get_transcript(transcript_id)

    out = {"transcript_id":transcript_id,
           "exon_id":transcript.exon_id[int(exon_number)],
           "exon_sequence":transcript.exon_sequence[int(exon_number)]
            }

    return out


def make_pairs(transcript1_id,transcript2_id):

    return process_ortholog_pairs(transcript1_id,transcript2_id)


def get_splice_site_info(transcript_id, exon_id):

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    table = 'Transcripts' + transcript_id[-3:]

    c.execute("SELECT start, end, strand FROM "+table+" WHERE Transcript_ID='"+transcript_id+"' AND Exon_ID='"+exon_id+"'")

    row = c.fetchall()

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

    return { "absolute_start" : absolute_start,
             "absolute_end" : absolute_end,
             "relative_start" : relative_start,
             "relative_end" : relative_end,
             "strand" : strand }

                
def get_protein_similarity(transcript1_id, transcript2_id, exon1_id, exon2_id):

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    table = 'Transcripts' + transcript1_id[-3:]

    c.execute("SELECT seq FROM "+table+" WHERE Transcript_ID='"+transcript1_id+"' AND Exon_ID='"+exon1_id+"'")

    rows = c.fetchall()

    exon1_seq = rows[0][0]

    table = 'Transcripts' + transcript2_id[-3:]

    c.execute("SELECT seq FROM "+table+" WHERE Transcript_ID='"+transcript2_id+"' AND Exon_ID='"+exon2_id+"'")

    rows = c.fetchall()

    exon2_seq = rows[0][0]

    return protein_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id)