'''
    The data structures used to model the transcripts.
    We have support for both fasta and local sqlite db.
'''


import sqlite3
from Bio import SeqIO



class Transcript_From_Fasta(object):

    '''
        Data structure for a transcript

        For now it's just a simple list 
        of exon sequences

        Works with fasta
    '''

    def __init__(self, fasta_file, speices_name=None, gene_id=None, gene_name=None, transcript_id=None, transcript_name=None):

        '''
            Loads the fasta_file
            and extracts the exon sequences

            Also contains the metadata (optional)

            Args:
                fasta_file : string = path to the fasta file                
                speices_name : string (optional) = Name of the species
                gene_id : string (optional) = Ensembl stable id of the gene
                gene_name : string (optional) = Name of the gene
                transcript_id : string = Ensembl stable id of the transcript
                transcript_name : string = Name of the transcript
        '''

        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_name = fasta_file.split('_')[-3] + fasta_file.split('_')[-2]  # Extract the transcript name 
                                                                                      # from name of the fasta file

        self.exon_sequence = []
        self.exon_id = []
        
        exons = SeqIO.parse(fasta_file, "fasta")  # parses all the exons from the fasta file

        for exon in exons:
            self.exon_id.append(exon.description.split(' ')[1])
            self.exon_sequence.append(exon.seq)
            


class Transcript_From_DB(object):

    '''
        Data structure for a transcript

        For now it's just a simple list 
        of exon sequences

        Works with sqlite database
    '''

    def __init__(self, path_to_db, transcript_id, speices_name=None, gene_id=None, gene_name=None):

        '''
            connects to the database
            and extracts the exon sequences

            Also contains the metadata (optional)

            Args:
                path_to_db : string = path to the database
                transcript_id : string = Ensembl stable id of the transcript
                speices_name : string (optional) = Name of the species
                gene_id : string (optional) = Ensembl stable id of the gene
                gene_name : string (optional) = Name of the gene
        '''

        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id

        self.exon_sequence = []
        self.exon_id = []

        conn = sqlite3.connect(path_to_db)  # connecting to the database
        c = conn.cursor()

        # Extract all the information related to the transcript from the database
        c.execute('SELECT Exon_ID, seq, Exon_No FROM Transcripts'+transcript_id[-3:]+" WHERE transcript_id='"+ transcript_id +"' ORDER BY Exon_No")
        all_rows = c.fetchall()

        for row in all_rows:

            self.exon_id.append(row[0])
            self.exon_sequence.append(row[1])

        conn.close()