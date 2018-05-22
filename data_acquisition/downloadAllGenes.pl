=pod
This script downloads all the 
genes of a species
and saves them in a sqlite database
=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Compara::Homology;
use Class::Inspector;
use DBI;
use Try::Tiny;

sub download_all_genes{

    my $reg = 'Bio::EnsEMBL::Registry';    

    $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org',   # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # extracting the values from function arguments
    my $species = $_[0];

    my $ga = $reg->get_adaptor($species, 'core', 'gene');  # creating the gene adaptor for human

    my $genes = $ga->fetch_all();

    # connecting to the database
    my $dbh = DBI->connect("dbi:SQLite:dbname=db.db","","");

    try {
        $dbh->do("CREATE TABLE ".$species."Genes (gene_stable_ID VARCHAR(1000), PRIMARY KEY (gene_stable_ID))");
    } catch {
        warn "caught error: $_";
    };

    # store the data in the database table
    for( my $i = 0; $i < @{$genes}; $i++) {
        
        my $gene_id = $genes->[$i]->stable_id();
        
        try {
                $dbh->do("INSERT INTO ".$species."Genes 
                (gene_stable_ID) VALUES 
                ('$gene_id')");
        } catch {
        warn "caught error: $_";
        };
    }
}

download_all_genes('human');
download_all_genes('mouse');