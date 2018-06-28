=pod
/* This script downloads the gene names
 * of the human-mouse ortholog pairs
 * and saves it in the database
 * 
 * so that we can retrieve the genes mentioned
 * in the publication
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


sub get_human_gene_names {

    # This collects the human gene names and
    # saves those information in the local database


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create column for storing human gene names
    try{
        $dbh->do("ALTER TABLE OrthologPairs ADD human_gene_name VARCHAR(100)");
    } catch {
        warn "caught error: $_";
    };

    my $sth = $dbh->prepare("SELECT human_gene_id FROM OrthologPairs"); # load all the gene ids
    $sth->execute();

    my @genes = ();    # array to store gene ids

    while(my @row = $sth->fetchrow_array()){
        push @genes, $row[0];
    }

    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    my $ga = $reg->get_adaptor('human', 'core', 'gene'); # gene adaptor
 
    for(my $i = 0; $i < scalar @genes; $i++ ){
        
        my $gene = $ga->fetch_by_stable_id($genes[$i]);
        
        my $gene_name = $gene->external_name(); # retrieve the gene name

        try { # save the gene name in the local database
                $dbh->do("UPDATE OrthologPairs 
                 SET human_gene_name = '$gene_name' 
                 WHERE human_gene_id='$genes[$i]'");
        } catch {
        warn "caught error: $_";
        };
    }
   
    # close the database connection
    $dbh->disconnect();

}

sub get_mouse_gene_names {

    # This collects the mouse gene names and
    # saves those information in the local database


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create column for storing mouse gene names
    try{
        $dbh->do("ALTER TABLE OrthologPairs ADD mouse_gene_name VARCHAR(100)");
    } catch {
        warn "caught error: $_";
    };

    my $sth = $dbh->prepare("SELECT mouse_gene_id FROM OrthologPairs"); # load all the gene ids
    $sth->execute();

    my @genes = ();    # array to store gene ids

    while(my @row = $sth->fetchrow_array()){
        push @genes, $row[0];
    }

    # connecting to Ensembl
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org', # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    my $ga = $reg->get_adaptor('mouse', 'core', 'gene'); # gene adaptor
 
    for(my $i = 0; $i < scalar @genes; $i++ ){
        
        my $gene = $ga->fetch_by_stable_id($genes[$i]);
        
        my $gene_name = $gene->external_name(); # retrieve the gene name

       

        try { # save the gene name in the local database
                $dbh->do("UPDATE OrthologPairs 
                 SET mouse_gene_name = '$gene_name' 
                 WHERE mouse_gene_id='$genes[$i]'");
        } catch {
        warn "caught error: $_";
        };
    }
   
    # close the database connection
    $dbh->disconnect();

}


get_human_gene_names();
get_mouse_gene_names();