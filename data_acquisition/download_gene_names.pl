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

sub create_gene_table {

    # This function creates a new table for genes where
    # the gene_name, gene_id and species information are stored
    # for better handling of the data


    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create table for storing gene information
    try{
        $dbh->do("CREATE TABLE Genes (Gene_ID VARCHAR(100), Gene_Name VARCHAR(100), Species VARCHAR(100), PRIMARY KEY (gene_id) )");
    } catch {
        warn "caught error: $_";
    };


    my $sth = $dbh->prepare("SELECT mouse_gene_id, mouse_gene_name FROM OrthologPairs"); # load all the mouse genes
    $sth->execute();

    my @gene_ids = ();    # array to store gene ids
    my @gene_names = ();    # array to store gene names

    while(my @row = $sth->fetchrow_array()){
        push @gene_ids, $row[0];
        push @gene_names, $row[1];
    }

    for(my $i = 0; $i < scalar @gene_ids; $i++ ){
        
        try { # save the gene information in the database
                $dbh->do("INSERT INTO genes 
                (gene_name,gene_id,species) 
                 VALUES ('$gene_names[$i]', '$gene_ids[$i]', 'mouse')");
        } catch {
        warn "caught error: $_";
        };
    }

    $sth = $dbh->prepare("SELECT human_gene_id, human_gene_name FROM OrthologPairs"); # load all the human genes
    $sth->execute();

    @gene_ids = ();    # array to store gene ids
    @gene_names = ();    # array to store gene names

    while(my @row = $sth->fetchrow_array()){
        push @gene_ids, $row[0];
        push @gene_names, $row[1];
    }

    for(my $i = 0; $i < scalar @gene_ids; $i++ ){
        
        try { # save the gene information in the database
                $dbh->do("INSERT INTO genes 
                 (gene_name,gene_id,species) 
                 VALUES ('$gene_names[$i]', '$gene_ids[$i]', 'human')");
        } catch {
        warn "caught error: $_";
        };
    }
   
    # close the database connection
    $dbh->disconnect();

}


get_human_gene_names();
get_mouse_gene_names();
create_gene_table();
