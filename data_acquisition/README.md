# Data Acquisition

Perl scripts for obtaining data from Ensembl using the Ensembl apis.

## Instructions

Download the scripts and create a 'data' folder

Run the following script to collect all the one2one orthologs between mouse and human and store them in a sqlite database, situated in the data directory.

```
$ perl downloadOrthologGenes.pl
```
Next running the following script will download all the transcript ids associated with the orthologous genes, and store them in the same database on a separate table called 'GeneTranscripts'

```
$ perl downloadGeneTranscripts.pl
```
