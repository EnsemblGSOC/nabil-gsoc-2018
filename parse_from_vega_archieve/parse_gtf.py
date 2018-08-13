'''
    This script parses the .gtf files 
    to extract informations related to the 
    transcripts, and stores them in a local
    sqlite database
'''


import sqlite3
from tqdm import tqdm 

def process_gtf_file(path,species):

    fp = open(path,'r')   # loading the file

    conn = sqlite3.connect('data/'+species+'.db')   # connecting the database
    c = conn.cursor()

    lines = fp.read().split('\n')[:-1]   # reading all the lines of the .gtf file

    for line in tqdm(lines,total=len(lines)):

        process_line(line,species,c)

    conn.commit()       # commiting to the database
    conn.close()        # closing the database

def process_line(line,species,cursor):


    elements = line.split('\t')
    
    seqname = elements[0]
    biotype = elements[1]
    category = elements[2]  # category => exon / CDS / stop_codon / start_codon
    start = elements[3]
    end = elements[4]
    strand = elements[6]
    attributes = elements[8].split(';')[:-1] # avoiding the last ';'

    dictionary = {
                    'gene_id' : '',
                    'transcript_id' : '',
                    'exon_number' : '',
                    'gene_name' : '',
                    'transcript_name' : '',
                    'exon_id' : '',
                    'protein_id' : ''
                 }
    


    for attribute in attributes:

        attribute = attribute.split(' ')[1:]  # avoiding the first space
        field = attribute[0]  
        value = attribute[1][1:-1]  # avoiding the quotation marks " "

        dictionary[field] = value

    # Add new gene

    try:

        cursor.execute("INSERT INTO " + species + " (gene_id, gene_name) VALUES  \
                                    ('"+dictionary['gene_id']+"','"+dictionary['gene_name']+"')")

    except Exception as e:  # species table hasn't been created

        try:
            
            cursor.execute('CREATE TABLE ' + species + ' ( gene_id VARCHAR(50), gene_name VARCHAR(50),  ensembl_id VARCHAR(50), UNIQUE (gene_id) )')

            cursor.execute("INSERT INTO " + species + " (gene_id, gene_name) VALUES  \
                                        ('"+dictionary['gene_id']+"','"+dictionary['gene_name']+"')")

        except:   # this is done for resume support
            pass

        

    # Add new transcript to gene table
    try:

        cursor.execute("INSERT INTO " + dictionary['gene_id'] + " (biotype, transcript_id, transcript_name) VALUES  \
                                        ('"+biotype+"','"+dictionary['transcript_id']+"','"+dictionary['transcript_name']+"')")

    except Exception as e:  # gene table hasn't been created

        try:
            cursor.execute('CREATE TABLE ' + dictionary['gene_id'] + ' ( biotype VARCHAR(50), transcript_id VARCHAR(50),  transcript_name VARCHAR(50), UNIQUE (transcript_id) )')

            cursor.execute("INSERT INTO " + dictionary['gene_id'] + " (biotype, transcript_id, transcript_name) VALUES  \
                                            ('"+biotype+"','"+dictionary['transcript_id']+"','"+dictionary['transcript_name']+"')")

        except:   # this is done for resume support
            pass


    # Add new data to transcript
    try:
        cursor.execute("INSERT INTO " + dictionary['transcript_id'] + " (category, start, end, strand, exon_number, exon_id, protein_id ) VALUES  \
                                        ('"+category+"','"+start+"','"+end+"','"+strand+"','"+dictionary['exon_number']+"','"+dictionary['exon_id']+"','"+dictionary['protein_id']+"')")

    except Exception as e:  # transcript table hasn't been created

        try:
            cursor.execute('CREATE TABLE ' + dictionary['transcript_id'] + ' ( seqname VARCHAR(50), category VARCHAR(30), start VARCHAR(10), end  VARCHAR(10), strand VARCHAR(5), exon_number VARCHAR(10), exon_id VARCHAR(50) , protein_id VARCHAR(50) , sequence VARCHAR(60000), UNIQUE (exon_id, protein_id,category,exon_number) )')

            cursor.execute("INSERT INTO " + dictionary['transcript_id'] + " (seqname, category, start, end, strand, exon_number, exon_id, protein_id ) VALUES  \
                                            ('"+seqname+"','"+category+"','"+start+"','"+end+"','"+strand+"','"+dictionary['exon_number']+"','"+dictionary['exon_id']+"','"+dictionary['protein_id']+"')")
        except:   # this is done for resume support
            pass


def main():

    # parse the gtf file for human   
    process_gtf_file('vega_files/Homo_sapiens.GRCh38.68.gtf','human')
    
    # parse the gtf file for mouse
    process_gtf_file('vega_files/Mus_musculus.GRCm38.68.gtf','mouse')

if __name__ == '__main__':
    main()
