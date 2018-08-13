=pod
This script downloads all the 
orthologous mouse genes of human 
genes. Stores them as a pair in database.
Also stores the goc_score, wga_coverage
and is_high_confidence
=cut


use strict;
use warnings;
use DBI;
use Try::Tiny;
use Bio::EnsEMBL::Registry;
use Data::Dumper;


sub download_all_orthologs{

    my $reg = "Bio::EnsEMBL::Registry";

    $reg->load_registry_from_db(
    -host => "asiadb.ensembl.org", # Please use the suitable host according to your location
    -user => "anonymous"
    );

    # connecting to the local database
    my $dbh = DBI->connect("dbi:SQLite:dbname=data/db.db","","");

    # create table for ortholog pairs
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

    # Getting the GenomeDB adaptor
    my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
        'Multi', 'compara', 'GenomeDB');

    # Fetching GenomeDB objects for human and mouse
    my $human_genome_db = $genome_db_adaptor->fetch_by_name_assembly('homo_sapiens');
    my $mouse_genome_db = $genome_db_adaptor->fetch_by_name_assembly('mus_musculus');

    # Getting the MethodLinkSpeciesSet adaptor
    my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
        'Multi', 'compara', 'MethodLinkSpeciesSet');

    # Fetching the MethodLinkSpeciesSet object corresponding to orthologs between human and mouse genomic sequences:
    my $human_mouse_orthologs =
        $method_link_species_set_adaptor->fetch_by_method_link_type_genome_db_ids(
            'ENSEMBL_ORTHOLOGUES',
            [$human_genome_db, $mouse_genome_db]
        );

    # Getting the Homology adaptor
    my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');

    # Fetch all the human-mouse one2one orthologs
    my $homologies = $homology_adaptor->fetch_all_by_MethodLinkSpeciesSet($human_mouse_orthologs,  -ORTHOLOGY_TYPE => 'ortholog_one2one');

    for(my $i = 0; $i < scalar @$homologies ; $i++){
    
        my $gene_members = $homologies->[$i]->gene_list();
        my $goc_score = $homologies->[$i]->goc_score();
        my $wga_coverage = $homologies->[$i]->wga_coverage();
        my $is_high_confidence = $homologies->[$i]->is_high_confidence();
        my $human_gene = $gene_members->[0]->stable_id();
        my $mouse_gene = $gene_members->[1]->stable_id();

        if(substr($mouse_gene, 3, 1) ne 'M')
        {
            $human_gene = $gene_members->[1]->stable_id();
            $mouse_gene = $gene_members->[0]->stable_id();
        }

        try {
                $dbh->do("INSERT INTO OrthologPairs  
                (human_gene_id, mouse_gene_id, goc_score, wga_coverage, is_high_confidence) VALUES 
                ('$human_gene', '$mouse_gene', '$goc_score', '$wga_coverage', '$is_high_confidence' )");
            } catch {
                warn "caught error: $_";
            };
    }

}


download_all_orthologs();


