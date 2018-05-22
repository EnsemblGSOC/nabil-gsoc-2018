=pod
This script downloads all the 
orthologous mouse genes of human 
genes. Stores them as a pair in database.
Also stores the goc_score, wga_coverage
and is_high_confidence
=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Compara::Homology;
use Class::Inspector;
use DBI;
use Try::Tiny;

sub get_gene_orthologs {
   
    # connecting to Ensembl db
    my $reg = 'Bio::EnsEMBL::Registry';

     $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org',  # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # connecting to the local database
    my $dbh = DBI->connect("dbi:SQLite:dbname=db.db","","");

    # extracting the values from function arguments
    my $gene_stable_id = $_[0];
    
    # initialize the gene adaptor
    my $ga = $reg->get_adaptor('human', 'core', 'gene');

    try{
        # collect the gene
        my $gene = $ga->fetch_by_stable_id($gene_stable_id);
    
        # collect the orthologs of the gene
        my $orthologs = $gene->get_all_homologous_Genes();

        for(my $i=0; $i<@$orthologs; $i++){
            
            my $ortholog_gene_id = $orthologs->[$i]->[0]->stable_id();
            
            # testing whether this gene is a mouse gene or not
            my $sth = $dbh->prepare("SELECT gene_stable_ID from mouseGenes WHERE gene_stable_ID='$ortholog_gene_id'");
            $sth->execute();

            if(my @row = $sth->fetchrow_array()){
                
                my $goc_score = "".$orthologs->[$i]->[1]->goc_score();
                my $wga_coverage = "".$orthologs->[$i]->[1]->wga_coverage();
                my $is_high_confidence = "".$orthologs->[$i]->[1]->is_high_confidence();
                try {
                    $dbh->do("INSERT INTO OrthologPairs  
                    (human_gene_id, mouse_gene_id, goc_score, wga_coverage, is_high_confidence) VALUES 
                    ('$gene_stable_id', '$ortholog_gene_id', '$goc_score', '$wga_coverage', '$is_high_confidence' )");
                } catch {
                warn "caught error: $_";
                };
            }
        }

        $dbh->do("UPDATE humanGenes
                 SET done = '1'
                 WHERE gene_stable_ID='$gene_stable_id'");

    } catch{
        warn "caught error: $_";

        $dbh->do("UPDATE humanGenes
                 SET done = '-1'
                 WHERE gene_stable_ID='$gene_stable_id'");

    };

    $dbh->disconnect();
}

sub batch_download
{
    my @genes = @{$_[0]};

    # connecting to Ensembl db
    my $reg = 'Bio::EnsEMBL::Registry';

        $ reg->load_registry_from_db(
        -host => 'asiadb.ensembl.org',  # Please use the suitable host according to your location
        -user => 'anonymous'
    );

    # initialize the gene adaptor
    my $ga = $reg->get_adaptor('human', 'core', 'gene');

    # connecting to the local database
    my $dbh = DBI->connect("dbi:SQLite:dbname=db.db","","");

    foreach (@genes)
    {
        my $gene_stable_id = $_;

        try{
            # collect the gene
            my $gene = $ga->fetch_by_stable_id($gene_stable_id);

            # collect the orthologs of the gene
            my $orthologs = $gene->get_all_homologous_Genes();
            
            for(my $i=0; $i<@$orthologs; $i++){
                
                my $ortholog_gene_id = $orthologs->[$i]->[0]->stable_id();
                
                # testing whether this gene is a mouse gene or not
                my $sth = $dbh->prepare("SELECT gene_stable_ID from mouseGenes WHERE gene_stable_ID='$ortholog_gene_id'");
                $sth->execute();

                if(my @row = $sth->fetchrow_array()){
                    
                    my $goc_score = "".$orthologs->[$i]->[1]->goc_score();
                    my $wga_coverage = "".$orthologs->[$i]->[1]->wga_coverage();
                    my $is_high_confidence = "".$orthologs->[$i]->[1]->is_high_confidence();
                    try {
                        $dbh->do("INSERT INTO OrthologPairs  
                        (human_gene_id, mouse_gene_id, goc_score, wga_coverage, is_high_confidence) VALUES 
                        ('$gene_stable_id', '$ortholog_gene_id', '$goc_score', '$wga_coverage', '$is_high_confidence' )");
                    } catch {
                    warn "caught error: $_";
                    };
                }
            }

            $dbh->do("UPDATE humanGenes
                    SET done = '1'
                    WHERE gene_stable_ID='$gene_stable_id'");

        } catch{
            warn "caught error: $_";

            $dbh->do("UPDATE humanGenes
                    SET done = '-1'
                    WHERE gene_stable_ID='$gene_stable_id'");

        };
    }

    $dbh->disconnect();
}


sub download_all_orthologs{

    my $dbh = DBI->connect("dbi:SQLite:dbname=db.db","","");

    try {
        $dbh->do("CREATE TABLE OrthologPairs (human_gene_id VARCHAR(1000),
                                              mouse_gene_id VARCHAR(1000),
                                              goc_score VARCHAR(10),
                                              wga_coverage VARCHAR(10),
                                              is_high_confidence VARCHAR(10),
                                              PRIMARY KEY (human_gene_id,mouse_gene_id))");
    } catch {
        warn "caught error: $_";
    };

    try {
        $dbh->do("ALTER TABLE humanGenes 
                ADD done VARCHAR(5) DEFAULT(0)");
    } catch {
        warn "caught error: $_";
    };

    my $sth = $dbh->prepare("SELECT gene_stable_ID from humanGenes WHERE done='0'");
    $sth->execute();

    my @genes = ();

    while(my @row = $sth->fetchrow_array()){
        push @genes, $row[0];
        
    }

    $dbh->disconnect();

    # Download all orthologes
    batch_download(\@genes);

}

download_all_orthologs();
