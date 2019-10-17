****************
Querying Compara
****************

The Ensembl compara database is represented by ``ensembldb3.compara.Compara``. This object provides a means for querying for relationships among genomes and obtaining multiple alignments. 

Creating a compara instance
===========================

Instantiating ``Compara`` requires the ensembl release, the series of species of interest and optionally an account (we also use our local account for speed). For the purpose of illustration we'll use the human, chimpanzee and macaque genomes. The resulting object has a ``Genome`` instance added as an attribute with name corresponding to the capitalised common name for each species.

.. doctest::

    >>> import os
    >>> from ensembldb3 import HostAccount, Compara
    >>> account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    >>> compara = Compara(['human', 'chimp', 'macaque'], release=85, account=account)
    >>> compara.Human
    Genome(species='Homo sapiens'; release='85')
    >>> compara.Chimp
    Genome(species='Pan troglodytes'; release='85')
    >>> compara.Macaque
    Genome(species='Macaca mulatta'; release='85')

.. note::
    Use ``Species.get_compara_name(species_name)`` to see what the attribute name will be.

Get the species tree
====================

This is accessed from the ``Compara`` instance.

.. doctest::
    
    >>> tree = compara.get_species_tree(just_members=True)
    >>> print(tree.ascii_art())
                        /-Pan troglodytes
              /Homininae
    -root----|          \-Homo sapiens
             |
              \-Macaca mulatta

What alignment types are available
==================================

What alignments are available for the species chosen can be displayed printing the ``method_species_links`` attribute of ``Compara``.

.. doctest::
    
    >>> print(compara.method_species_links)
    Align Methods/Clades
    ===================================================================================================================
    method_link_species_set_id  method_link_id  species_set_id      align_method                            align_clade
    -------------------------------------------------------------------------------------------------------------------
                           756              13           35886               EPO                         8 primates EPO
                           780              13           36102               EPO               17 eutherian mammals EPO
                           781              14           36103  EPO_LOW_COVERAGE  39 eutherian mammals EPO_LOW_COVERAGE
                           788              10           36176             PECAN           23 amniota vertebrates Pecan
    -------------------------------------------------------------------------------------------------------------------

.. note::
    Any queries on this instance of compara will only return results for the indicated species. If you want to query about other species, create another instance.

Get a syntenic region
=====================

Get genomic alignment for the *BRCA2* gene region. We can specify the alignment set we want the data from using ``align_method`` and ``align_clade``. We then use the ``get_alignment()`` method. We can further identify feature types we want the sequences in the alignment to be annotated with.

.. doctest::
    
    >>> human_brca2 = compara.Human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> regions = compara.get_syntenic_regions(region=human_brca2, align_method='EPO', align_clade='primate')
    >>> for region in regions:
    ...     print(region)
    SyntenicRegions:
      Coordinate(Human,chro...,13,32315473-32400266,1)
      Coordinate(Chimp,chro...,13,31957346-32041418,-1)
      Coordinate(Macaque,chro...,17,11686607-11779396,-1)
    >>> aln = region.get_alignment(feature_types=["gene", "repeat"])
    >>> aln
    3 x 99492 dna alignment: Homo sapiens:chromosome:13:32315473-32400266:1[GGGCTTGTGGC...], Pan troglodytes:chromosome:13:31957346-32041418:1[GGGCTTGTGGC...], Macaca mulatta:chromosome:17:11686607-11779396:1[GGGCTTGTGGC...]

A Cogent3 annotated alignment object is returned. This can be queried to get annotations corresponding to specific features, or for masking those features, etc.. See the Cogent3 documentation for more information on `using annotations <http://cogent3.readthedocs.io/en/latest/examples/complete_seq_features.html>`_.

.. doctest::
    
    >>> print(aln.get_annotations_from_any_seq('CDS'))
    [CDS "ENST00000380152" at [1000:1067, 3676:3925, 10390:10499...

You can specify an equivalent query using the corresponding ``method_link_species_set_id`` value.

.. doctest::
    
    >>> regions = compara.get_syntenic_regions(region=human_brca2, method_clade_id=756)
    >>> for region in regions:
    ...     print(region)
    ...     # tree = region.get_species_tree()
    ...     # print(tree.ascii_art())
    SyntenicRegions:
      Coordinate(Human,chro...,13,32315473-32400266,1)
      Coordinate(Chimp,chro...,13,31957346-32041418,-1)
      Coordinate(Macaque,chro...,17,11686607-11779396,-1)

Get related genes
=================

Types of relationships
----------------------

.. doctest::
    
    >>> compara.get_distinct('relationship') # doctest: +SKIP
    ['ortholog_many2many', 'within_species_paralog', 'gene_split', 'ortholog_one2one', 'ortholog_one2many', 'alt_allele', 'other_paralog']

One-to-one orthologs
--------------------

We get the one-to-one orthologs for *BRCA2*.

.. doctest::

    >>> orthologs = compara.get_related_genes(stableid='ENSG00000139618',
    ...                  relationship='ortholog_one2one')

We iterate over the related members, which are gene instances (see :ref:`gene`)

.. doctest::
    
    >>> for ortholog in orthologs.members:
    ...     print(ortholog)
    Gene(species='Pan troglodytes'; biotype='protein_coding'; description='BRCA2, DNA repair...'; location=Coordinate(Chimp,chro...,13,31957352-32040817,1); stableid='ENSPTRG00000005766'; status='KNOWN'; symbol='BRCA2')
    Gene(species='Macaca mulatta'; biotype='protein_coding'; description='BRCA2, DNA repair...'; location=Coordinate(Macaque,chro...,17,11687583-11777925,1); stableid='ENSMMUG00000007197'; status='KNOWN_BY_PROJECTION'; symbol='BRCA2')
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='BRCA2, DNA repair...'; location=Coordinate(Human,chro...,13,32315473-32400266,1); stableid='ENSG00000139618'; status='KNOWN'; symbol='BRCA2')

We can identify which species are in the set (it may not always be the same as those the ``Compara`` instance was created with).

.. doctest::
    
    >>> orthologs.get_species_set() # doctest: +SKIP
    {'Pan troglodytes', 'Macaca mulatta', 'Homo sapiens'}

Within species paralogs
-----------------------

I'm using the haemoglobin B locus identifier ``ENSG00000244734``.

.. doctest::
    
    >>> paras = compara.get_related_genes(stableid='ENSG00000244734',
    ...                 relationship="within_species_paralog")
    >>> print(paras)
    RelatedGenes:
     relationships=within_species_paralog
      Gene(species='Homo sapiens'; biotype='protein_coding'; description='hemoglobin subunit beta...'; location=Coordinate(Human,chro...,11,5225463-5229395,-1); stableid='ENSG00000244734'; status='KNOWN'; symbol='HBB')
      Gene(species='Homo sapiens'; biotype='protein_coding'; description='hemoglobin subunit delta...'; location=Coordinate(Human,chro...,11,5232677-5235370,-1); stableid='ENSG00000223609'; status='KNOWN'; symbol='HBD')...
    >>> tree = paras.get_tree()
    >>> print(tree.ascii_art())
                                  /-ENSG00000223609
                        /edge.0--|
                       |          \-ENSG00000244734
              /edge.3--|
             |         |                    /-ENSG00000196565
             |         |          /edge.1--|
             |          \edge.2--|          \-ENSG00000213934
             |                   |
             |                    \-ENSG00000213931
    -root----|
             |                                        /-ENSG00000188536
             |                              /edge.4--|
             |                    /edge.5--|          \-ENSG00000206172
             |                   |         |
             |          /edge.6--|          \-ENSG00000086506
             |         |         |
              \edge.7--|          \-ENSG00000206177
                       |
                        \-ENSG00000130656

Getting the CDS from related genes
----------------------------------

Note the name for the sequence is derived from the coordinate.

.. doctest::
    
    >>> cds_seqs = []
    >>> for gene in paras.members:
    ...     cds = gene.canonical_transcript.cds
    ...     cds_seqs.append([gene.stableid, cds])
    ...     
    >>> cds_seqs
    [['ENSG00000244734', DnaSequence(ATGGTGC... 444)], ['ENSG00000223609', DnaSequence(ATGGTGC... 444)], ['ENSG00000213934', DnaSequence(ATGGGTC... 444)], ['ENSG00000196565', DnaSequence(ATGGGTC... 444)], ['ENSG00000213931', DnaSequence(ATGGTGC... 444)], ['ENSG00000130656', DnaSequence(ATGTCTC... 429)], ['ENSG00000206177', DnaSequence(ATGCTCA... 426)], ['ENSG00000188536', DnaSequence(ATGGTGC... 429)], ['ENSG00000206172', DnaSequence(ATGGTGC... 429)], ['ENSG00000086506', DnaSequence(ATGGCGC... 429)]]
