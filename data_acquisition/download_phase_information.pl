=pod
/* This script downloads the start phase and
 * end phase of  all the human-mouse exons
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


sub get_phase {

    # This collects all the start phase and end phase of the 
    # exons and saves those information in the local database


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # iterate over all transcripts tables

    for( my $i= 0; $i<1000; $i++){

        my $table_name = "Transcripts";

        if($i<10){
            $table_name = $table_name . "00$i";
        }

        elsif($i<100){
            $table_name = $table_name . "0$i";
        }

        else{
            $table_name = $table_name . "$i";
        }

        # create column for storing starting and ending phase 
        try{
            $dbh->do("ALTER TABLE $table_name ADD start_phase VARCHAR(10)");
            $dbh->do("ALTER TABLE $table_name ADD end_phase VARCHAR(10)");
        } catch {
            warn "caught error: $_";
        };

    }


    # collect all the transcripts
    my $sth = $dbh->prepare("SELECT transcript_id, species FROM GeneTranscripts WHERE biotype='protein_coding' AND gencode_basic='1'");
    $sth->execute();

    my @transcripts = (); # array to store transcript ids
    my @species_names = (); # array to store species name

    while(my @row = $sth->fetchrow_array()){
        push @transcripts, $row[0];
        push @species_names, $row[1];
    }

    my $total_transcripts = scalar @transcripts;

    for(my $i = 0; $i < $total_transcripts; $i++ ){
    
        my $table_name = "Transcripts" . substr($transcripts[$i], -3);

        my $ta = $reg->get_adaptor($species_names[$i], 'core', 'Transcript'); # transcript adaptor

        my $transcript = $ta->fetch_by_stable_id($transcripts[$i]);

        my $exons = $transcript->get_all_Exons(); # get all exons from this transcript

        my $total_exons =scalar @{$exons};

        for (my $index=0; $index< $total_exons; $index++){
            
            my $exon_id = $exons->[$index]->stable_id();
            my $start_phase = $exons->[$index]->phase();
            my $end_phase = $exons->[$index]->end_phase();

            # save the exon phase in database

            try {
                    $dbh->do("UPDATE $table_name 
                    SET start_phase = '$start_phase', end_phase = '$end_phase' 
                    WHERE Transcript_ID='$transcripts[$i]' AND Exon_ID='$exon_id'");
            } catch {
            warn "caught error: $_";
            };
        }

    }

   
    # close the database connection
    $dbh->disconnect();

}



get_phase();