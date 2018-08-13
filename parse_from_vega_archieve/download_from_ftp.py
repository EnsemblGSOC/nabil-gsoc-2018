'''
    This script downloads the .gtf and .fa files 
    from vega archieve using the ftp 
'''


import wget
import gzip

def download_from_archieve(url,file_name):

    # download the file in compressed format
    output_directory = 'vega_files'
    downloaded_file = wget.download(url, out=output_directory)
    
    # read the file content
    f = gzip.open(downloaded_file, 'rb')
    file_content = f.read()
    f.close()

    # save the file
    fp = open('vega_files/'+file_name,'wb')
    fp.write(file_content)
    fp.close()




def main():


    # download the gtf file for human
    download_from_archieve("ftp://ftp.ensembl.org/pub/vega/human/Homo_sapiens.GRCh38.68.gtf.gz","Homo_sapiens.GRCh38.68.gtf")
    
    # download the gtf file for mouse
    download_from_archieve("ftp://ftp.ensembl.org/pub/vega/mouse/Mus_musculus.GRCm38.68.gtf.gz","Mus_musculus.GRCm38.68.gtf")
    
    # download the cdna fasta file for human
    download_from_archieve("ftp://ftp.ensembl.org/pub/vega/human/cdna/Homo_sapiens.VEGA68.cdna.all.fa.gz","Homo_sapiens.VEGA68.cdna.all.fa")

    # download the cdna fasta file for mouse
    download_from_archieve("ftp://ftp.ensembl.org/pub/vega/mouse/cdna/Mus_musculus.VEGA68.cdna.all.fa.gz","Mus_musculus.VEGA68.cdna.all.fa")


if __name__ == '__main__':
    main()

