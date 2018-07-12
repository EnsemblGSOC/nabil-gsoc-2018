"""
    This python script creates a latex file from the 
    exon similarity results. 

    From the latex file we can create a pdf file, which
    more convenient to read than raw text files
"""

import os 
from tqdm import tqdm

def create_latex_document():

    """
        This function encapsulates all the steps in creating the 
        latex file.

        This function stores the file as document.text

        After the file is generated we can create a pdf file
    """

    # latex document
    document = """
\\documentclass[a4paper]{article}

%% Language and font encodings
\\usepackage[english]{babel}
\\usepackage[utf8x]{inputenc}
\\usepackage[T1]{fontenc}

%% Sets page size and margins
\\usepackage[a4paper,top=3cm,bottom=2cm,left=3cm,right=3cm,marginparwidth=1.75cm]{geometry}

%% Useful packages
\\usepackage{longtable}
\\usepackage{amsmath}
\\usepackage{graphicx}
\\usepackage[colorinlistoftodos]{todonotes}
\\usepackage[colorlinks=true, allcolors=blue]{hyperref}

\\usepackage{array}
\\newcolumntype{C}[1]{>{\centering\let\\newline\\\\\\arraybackslash\hspace{0pt}}m{#1}}

\\title{Transcript Comparisons}
\\author{Exon Alignments (Number of skips allowed = 0)}

\\begin{document}

\\maketitle

\\tableofcontents

"""
    # browse all the results
    genes = next(os.walk('results'))[1]

    for gene in tqdm(genes,total=len(genes)):

        # adding the gene information
        document += "\\newpage\\section{"+gene.split(':')[0]+" ("+gene.split(':')[1]+ ")}\n\n"

        # collect all transcripts
        transcripts = next(os.walk(os.path.join('results',gene)))[2]

        # read all the similarity scores
        for transcript in tqdm(transcripts,total=len(transcripts)):
            
            fp = open(os.path.join('results',gene,transcript),'r')
            data = fp.read().split('\n')[:-1]
            fp.close()

            transcript_name = data[0].split('\t')[0]
            exon_number = data[0].split('\t')[1]

            # adding the transcript information
            document += "\\subsection{"+ transcript_name +"}\nTranscript \\href{https://www.ensembl.org/Homo_sapiens/Transcript/Exons?t=" + transcript_name + "}{" + transcript_name + "} (" + exon_number + ")\n"

            # creating a table
            document += """
\\begin{center}
\\begin{longtable}{| C{9cm} | C{1.4cm} | C{3.8cm} | C{1.4cm} |}
\\hline
Exon Similarity Scores & Transcript Score & Mouse Transcript & Number of exons \\\\ 
\\hline
"""
            for i in range(2,len(data)):
                
                row_data = data[i].split('\t')

                # adding data to the table
                document += row_data[0] + " & " + row_data[1] + " & " + "\\href{https://www.ensembl.org/Mus_musculus/Transcript/Exons?t=" + row_data[2] + "}{" + row_data[2] + "}" + " & " + row_data[3] + "\\\\\n\\hline\n"

            # closing the table
            document += "\\end{longtable}\n\\end{center}\n\n"


             

        
    # completing the latex file
    document += "\\end{document}"

    #saving the latex file
    fp = open('document.tex','w')
    fp.write(document)
    fp.close()


def main():
    
    create_latex_document()

if __name__ == '__main__':
    main()