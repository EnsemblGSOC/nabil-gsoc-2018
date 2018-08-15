=pod
/* This script checks whether a transcript
 * is gencode_basic or not
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


sub filter_gencode_basic {

    # This checks all the transcripts for gencode_basic 
    # saves those information in the local database


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create column for storing gencode_basic
    try{
        $dbh->do("ALTER TABLE GeneTranscripts ADD gencode_basic VARCHAR(100)");
    } catch {
        warn "caught error: $_";
    };

    # since we're only working with the protein coding transcripts, we're only checking for those
    my $sth = $dbh->prepare("SELECT transcript_id, gene_id FROM GeneTranscripts WHERE Biotype='protein_coding'");
    $sth->execute();

    my @transcripts = (); # array to store transcript ids
    my @genes = ();       # array to store gene ids

    while(my @row = $sth->fetchrow_array()){
        push @transcripts, $row[0];
        push @genes, $row[1];
    }


    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    for(my $i = 0; $i < scalar @transcripts; $i++ ){
        
        try{
            my $species = 'human';  # default value
        
            if(substr($genes[$i], 3, 1) ne 'G'){        
                $species = 'mouse';                 # test for the species name
            }

            my $ta = $reg->get_adaptor($species, 'core', 'transcript'); # transcript adaptor

            my $transcript = $ta->fetch_by_stable_id($transcripts[$i]);
            
            my $gencode_basic = $transcript->gencode_basic();


            try { # saving the gencode name in the local database
                    $dbh->do("UPDATE GeneTranscripts 
                    SET gencode_basic = '$gencode_basic' 
                    WHERE transcript_id='$transcripts[$i]'");
            } catch {
                warn "caught error: $_";
            }

        } catch {  # in case there are no name saved for that transcript
            warn "caught error: $_";

        }
    }
   
    # close the database connection
    $dbh->disconnect();

}


filter_gencode_basic();
