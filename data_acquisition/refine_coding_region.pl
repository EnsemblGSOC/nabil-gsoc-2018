=pod
/* This script downloads the coding region start and
 * coding region end information of  all the human-mouse transcripts
 * and saves it in the database
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


sub get_coding_region {

    # This collects all the coding region start and coding region end of the 
    # transcripts and saves those information in the local database


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create column for storing coding region start and end
    try{
        $dbh->do("ALTER TABLE GeneTranscripts ADD coding_region_start VARCHAR(20)");
        $dbh->do("ALTER TABLE GeneTranscripts ADD coding_region_end VARCHAR(20)");
    } catch {
        warn "caught error: $_";
    };

    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # collect all the transcripts
    my $sth = $dbh->prepare("SELECT transcript_id, species FROM GeneTranscripts WHERE biotype='protein_coding'");
    $sth->execute();

    my @transcripts = (); # array to store transcript ids
    my @species_names = (); # array to store species name

    while(my @row = $sth->fetchrow_array()){
        push @transcripts, $row[0];
        push @species_names, $row[1];
    }

    my $total_transcripts = scalar @transcripts;

    for(my $i = 0; $i < $total_transcripts; $i++ ){
    
        my $ta = $reg->get_adaptor($species_names[$i], 'core', 'Transcript'); # transcript adaptor

        my $transcript = $ta->fetch_by_stable_id($transcripts[$i]);

        my $coding_start = $transcript -> cdna_coding_start();
        my $coding_end = $transcript -> cdna_coding_end();
        
            # save the coding region start and end in database

            try {
                    $dbh->do("UPDATE GeneTranscripts
                    SET coding_region_start = '$coding_start', coding_region_end = '$coding_end' 
                    WHERE Transcript_ID='$transcripts[$i]'");
            } catch {
            warn "caught error: $_";
            };
    }

    

   
    # close the database connection
    $dbh->disconnect();

}



get_coding_region();