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
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # store the data in the database table
    for( my $i = 0; $i < @$transcripts; $i++ ) {

        my $transcript_id = $transcripts->[$i]->stable_id();
        my $biotype = $transcripts->[$i]->biotype();
        try {
                $dbh->do("INSERT INTO GeneTranscripts 
                (Transcript_ID, Gene_ID, Biotype) VALUES 
                ('$transcript_id', '$gene_stable_id' , '$biotype')");
        } catch {
        warn "caught error: $_";
        };
    }

    $dbh->disconnect();

}

sub get_all_transcripts {

    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create table to store gene transcripts
    try {
        $dbh->do("CREATE TABLE GeneTranscripts (Transcript_ID VARCHAR(1000),
                                                Gene_ID VARCHAR(1000),                                                
                                                Biotype VARCHAR(1000),
                                                PRIMARY KEY (Transcript_ID))");
    } catch {
        warn "caught error: $_";
    };

    my @human_genes = ();
    my @mouse_genes = ();

    my $sth = $dbh->prepare("SELECT human_gene_id, mouse_gene_id from OrthologPairs");
    $sth->execute();
    
    while(my @row = $sth->fetchrow_array()){
        push @human_genes, $row[0];
        push @mouse_genes, $row[1];
    }

    $dbh->disconnect();

    for(my $i = 0; $i < scalar @human_genes; $i++ ){
        get_gene_transcripts('human',$human_genes[$i]);
    }

    for(my $i = 0; $i < scalar @mouse_genes; $i++ ){
        get_gene_transcripts('mouse',$mouse_genes[$i]);
    }

}


get_all_transcripts()

