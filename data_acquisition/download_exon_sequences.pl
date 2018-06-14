=pod
/* This script downloads the exon sequence 
 * of a transcript specified by
 * its stable id and saves it in the 
 * database
*/ 
=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Compara::Homology;
use DBI;
use Try::Tiny;
use Data::Dumper;


sub get_exons_from_transcript {

    # Takes the gene_id and transcript_id as input
    # collects the exon sequences with data like start, end, strand
    # saves those information in the local database

    my $gene_stable_id = $_[0];
    my $transcript_stable_id = $_[1];


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

   
    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );


    my $species = 'human';  # default value
    
    if(substr($gene_stable_id, 3, 1) ne 'G'){        
        my $species = 'mouse';                 # test for the species name
    } 

    # To avoid millions of table or millions of rows, we stored transcripts in 1000 tables
    my $table_name = "Transcripts" . substr($transcript_stable_id, -3);

    try {
        $dbh->do("CREATE TABLE $table_name (Transcript_ID VARCHAR(1000),
                                                      Exon_No VARCHAR(100),
                                                      Exon_ID VARCHAR(1000),
                                                      start VARCHAR(100),                                                
                                                      end VARCHAR(100),
                                                      strand VARCHAR(5),
                                                      seq VARCHAR(5000),                                                    
                                                      PRIMARY KEY (Exon_ID))");
    } catch {
        warn "caught error: $_";
    };

    my $ta = $reg->get_adaptor($species, 'core', 'Transcript'); # transcript adaptor

    my $transcript = $ta->fetch_by_stable_id($transcript_stable_id);

    my $exons = $transcript->get_all_Exons(); # get all exons from this transcript

    my $total_exons =scalar @{$exons};


    for (my $index=0; $index< $total_exons; $index++){

        my $exon_no = "".($index+1);
        my $exon_id = $exons->[$index]->stable_id();
        my $start = $exons->[$index]->start();
        my $end = $exons->[$index]->end();
        my $strand = $exons->[$index]->strand();
        my $seq = $exons->[$index]->seq()->primary_seq->seq();

        # save the exon sequence in database

        try {
                $dbh->do("INSERT INTO $table_name 
                (Transcript_ID, Exon_No, Exon_ID, start, end, strand, seq) VALUES 
                ('$transcript_stable_id', '$exon_no', '$exon_id', '$start', '$end', '$strand', '$seq')");
        } catch {
        warn "caught error: $_";
        };
    }

    # close the database connection
    $dbh->disconnect();

}

sub get_exons_from_all_transcripts{

    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # selecting only the protein_coding transcripts
    my $sth = $dbh->prepare("SELECT Gene_ID, Transcript_ID FROM GeneTranscripts WHERE Biotype='protein_coding' ORDER BY Gene_ID");
    $sth->execute();

    my @genes = ();
    my @transcripts = ();

    while(my @row = $sth->fetchrow_array()){
        push @genes, $row[0];
        push @transcripts, $row[1];
    }

    $dbh->disconnect();

    # collecting exon sequence for all the transcripts
    for(my $i = 0; $i < scalar @genes; $i++ ){
        get_exons_from_transcript($genes[$i],$transcripts[$i]);
    }
}


get_exons_from_all_transcripts();
