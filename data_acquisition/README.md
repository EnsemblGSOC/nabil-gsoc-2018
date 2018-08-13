# Data Acquisition

Perl scripts for obtaining data from Ensembl using the Ensembl apis.

## Instructions

Download the scripts and create a 'data' folder

Run the following script to collect all the one2one orthologs between mouse and human and store them in a sqlite database, situated in the data directory.

```
$ perl download_ortholog_genes.pl
```
Next, running the following script will download all the transcript ids associated with the orthologous genes, and store them in the same database on a separate table called 'GeneTranscripts'

```
$ perl download_gene_transcripts.pl
```
Next, download the exon sequences of all the transcripts.
```
$ perl download_exon_sequences.pl
```

Next, download the name of the ortholog genes, so that we can make a bridge between the current Ensembl database and the archieved vega database
```
$ perl download_gene_names.pl
```

Next, we can download the name of the transcripts, for better readability and represenation instead of using long transcript stable ids
```
$ perl download_transcript_names.pl
```

Next, we need to download the phase information of the exons for translating them to protein correctly
```
$ perl download_phase_information.pl
```

Next, we need to download the coding region start and end information to identify the UTRs
```
$ perl refine_coding_region.pl
```

## Shortcut

Instead of following the whole pipeline of data acquisition, you can easily download the processed data from this [link](https://drive.google.com/open?id=1X51bFX-WVbsp346byt36VW3kCUb6mPtb)


## Troubleshooting

If you face difficultiese while connecting to Ensembl, please use a host that is most suitable based on your location. Since I'm in Asia, I'm using the Asia mirror.

In all the codes the connection is made in the following way

```
# connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );
```
Just change the host here based on your location.
