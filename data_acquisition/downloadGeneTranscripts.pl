=pod
This script downloads all the 
transcripts of a gene specified by
its stable id and saves them in the
database 
=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Compara::Homology;
use Class::Inspector;
use DBI;
use Try::Tiny;

sub get_gene_transcripts {
   
   # Takes the speices name and gene stable_id as input
   # collects the transcripts of the gene
   # saves the transcript id and biotype in database

    # connecting to Ensembl db
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org',  # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # extracting the values from function arguments
    my $species = $_[0];
    my $gene_stable_id = $_[1];

    # initialize the gene adaptor
    my $ga = $reg->get_adaptor($species, 'core', 'gene');

    # collect the gene
    my $gene = $ga->fetch_by_stable_id($gene_stable_id);

    # collect the transcripts associated with the gene
    my $transcripts = $gene->get_all_Transcripts();


    # Saving data to the database

    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=db.db","","");

    try {
        $dbh->do("CREATE TABLE " . $gene_stable_id. "Transcripts (Transcript_ID VARCHAR(1000),
                                                       Biotype VARCHAR(1000),
                                                       PRIMARY KEY (Transcript_ID))");
    } catch {
        warn "caught error: $_";
    };

    # store the data in the database table
    for( my $i = 0; $i < @$transcripts; $i++ ) {

        my $transcript_id = $transcripts->[$i]->stable_id();
        my $biotype = $transcripts->[$i]->biotype();
        try {
                $dbh->do("INSERT INTO " . $gene_stable_id. "Transcripts 
                (Transcript_ID, Biotype) VALUES 
                ('$transcript_id', '$biotype')");
        } catch {
        warn "caught error: $_";
        };
    }

    $dbh->disconnect();

}

get_gene_transcripts('human', 'ENSG00000276626');

