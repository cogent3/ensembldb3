Note that much more extensive documentation is available in :ref:`query-ensembl`.

Connecting
----------

.. Gavin Huttley

`Ensembl <http://www.ensembl.org>`_ provides access to their MySQL databases directly or users can download and run those databases on a local machine. To use the Ensembl's UK servers for running queries, nothing special needs to be done as this is the default setting for PyCogent's ``ensembl`` module. To use a different Ensembl installation, you create an account instance:

.. doctest::

    >>> from ensembldb import HostAccount
    >>> account = HostAccount('fastcomputer.topuni.edu', 'username',
    ...                       'somepass')

To specify a specific port to connect to MySQL on:

.. doctest::

    >>> from ensembldb import HostAccount
    >>> account = HostAccount('anensembl.server.edu', 'someuser',
    ...                       'somepass', port=3306)

.. we create valid account now to work on my local machines here at ANU

.. doctest::
    :hide:

    >>> import os
    >>> hotsname, uname, passwd = os.environ['ENSEMBL_ACCOUNT'].split()
    >>> account = HostAccount(hotsname, uname, passwd)

Species to be queried
---------------------

To see what existing species are available

.. doctest::

    >>> from ensembldb import Species
    >>> print Species
    ================================================================================
           Common Name                   Species Name              Ensembl Db Prefix
    --------------------------------------------------------------------------------
             A.aegypti                  Aedes aegypti                  aedes_aegypti
            A.clavatus           Aspergillus clavatus           aspergillus_clavatus...

If Ensembl has added a new species which is not yet included in ``Species``, you can add it yourself.

.. doctest::

    >>> Species.amend_species('A latinname', 'a common name')

You can get the common name for a species

.. doctest::

    >>> Species.get_common_name('Procavia capensis')
    'Rock hyrax'

and the Ensembl database name prefix which will be used for all databases for this species.

.. doctest::

    >>> Species.get_ensembl_db_prefix('Procavia capensis')
    'procavia_capensis'

Species common names are used to construct attributes on PyCogent ``Compara`` instances). You can get the name that will be using the ``get_compara_name`` method. For species with a real common name

.. doctest::
    
    >>> Species.get_compara_name('Procavia capensis')
    'RockHyrax'

or with a shortened species name

.. doctest::
    
    >>> Species.get_compara_name('Caenorhabditis remanei')
    'Cremanei'

Get genomic features
--------------------

Find a gene by gene symbol
^^^^^^^^^^^^^^^^^^^^^^^^^^

We query for the *BRCA2* gene for humans.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> print human
    Genome(species='Homo sapiens'; release='76')
    >>> genes = human.get_genes_matching(symbol='BRCA2')
    >>> for gene in genes:
    ...     if gene.symbol == 'BRCA2':
    ...         print gene
    ...         break
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast cancer 2,...'; stableid='ENSG00000139618'; Status='KNOWN'; symbol='BRCA2')

Find a gene by Ensembl Stable ID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the stable ID for *BRCA2*.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> gene = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> print gene
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast cancer 2,...'; stableid='ENSG00000139618'; Status='KNOWN'; symbol='BRCA2')

Find genes matching a description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We look for breast cancer related genes that are estrogen induced.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> genes = human.get_genes_matching(description='breast cancer anti-estrogen')
    >>> for gene in genes:
    ...     print gene
    Gene(species='Homo sapiens'; biotype='lincRNA'; description='breast cancer anti-estrogen...'; stableid='ENSG00000262117'; Status='NOVEL'; symbol='BCAR4')...

We can also require that an exact (case insensitive) match to the word(s) occurs within the description by setting ``like=False``.

.. doctest::
    
    >>> genes = human.get_genes_matching(description='breast cancer anti-estrogen',
    ...                                  like=False)
    >>> for gene in genes:
    ...     print gene
    Gene(species='Homo sapiens'; biotype='lincRNA'; description='breast cancer anti-estrogen...'; stableid='ENSG00000262117'; Status='NOVEL'; symbol='BCAR4')...

Get canonical transcript for a gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We get the canonical transcripts for *BRCA2*.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> print transcript
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315473; end=32400266; length=84793; strand='+')

Get the CDS for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> cds = transcript.Cds
    >>> print type(cds)
    <class 'cogent.core.sequence.DnaSequence'>
    >>> print cds
    ATGCCTATTGGATCCAAAGAGAGGCCA...

Look at all transcripts for a gene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> for transcript in brca2.Transcripts:
    ...     print transcript
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315473; end=32400266; length=84793; strand='+')
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315504; end=32333291; length=17787; strand='+')...

Get the first exon for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We show just for the canonical transcript.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> print brca2.CanonicalTranscript.Exons[0]
    Exon(stableid=ENSE00001184784, rank=1)

Get the introns for a transcript
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We show just for the canonical transcript.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> for intron in brca2.CanonicalTranscript.Introns:
    ...     print intron
    Intron(TranscriptId=ENST00000380152, rank=1)
    Intron(TranscriptId=ENST00000380152, rank=2)
    Intron(TranscriptId=ENST00000380152, rank=3)...


Inspect the genomic coordinate for a feature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> print brca2.location.coord_name
    13
    >>> print brca2.location.start
    32315473
    >>> print brca2.location.strand
    1

Get repeat elements in a genomic interval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query the genome for repeats within a specific coordinate range on chromosome 13.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> repeats = human.get_features(coord_name='13', start=32305473, end=32315473, feature_types='repeat')
    >>> for repeat in repeats:
    ...     print repeat.RepeatClass
    ...     print repeat
    ...     break
    SINE/Alu
    Repeat(coord_name='13'; start=32305225; end=32305525; length=300; strand='-', Score=2770.0)

Get CpG island elements in a genomic interval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query the genome for CpG islands within a specific coordinate range on chromosome 11.

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> islands = human.get_features(coord_name='11', start=2129111, end=2149604, feature_types='cpg')
    >>> for island in islands:
    ...     print island
    ...     break
    CpGisland(coord_name='11'; start=2137721; end=2141254; length=3533; strand='-', Score=3254.0)

Get SNPs
--------

For a gene
^^^^^^^^^^

We find the genetic variants for the canonical transcript of *BRCA2*.

.. note:: The output is significantly truncated!

.. doctest::

    >>> from ensembldb import Genome
    >>> human = Genome('human', release=76, account=account)
    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> transcript = brca2.CanonicalTranscript
    >>> print transcript.Variants
    (<cogent.db.ensembl.region.Variation object at ...
    >>> for variant in transcript.Variants:
    ...     print variant
    ...     break
    Variation(symbol='rs370721506'; effect=['non_coding_exon_variant', 'nc_transcript_variant', '5_prime...

Get a single SNP
^^^^^^^^^^^^^^^^

We get a single SNP and print it's allele frequencies.

.. doctest::
    
    >>> snp = list(human.get_variation(symbol='rs34213141'))[0]
    >>> print snp.AlleleFreqs
    =================================
    allele      freq    population_id
    ---------------------------------
         A    0.0303              933
         G    0.9697              933
         G    1.0000            11208
         G    1.0000            11519
         A                      11961
         G                      11961...

What alignment types available
------------------------------

We create a ``Compara`` instance for human, chimpanzee and macaque.

.. doctest::

    >>> from ensembldb import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], release=76,
    ...                  account=account)
    >>> print compara.method_species_links
    Align Methods/Clades
    ===================================================================================================================
    method_link_species_set_id  method_link_id  species_set_id      align_method                            align_clade
    -------------------------------------------------------------------------------------------------------------------
                           753              10           35883             PECAN           22 amniota vertebrates Pecan
                           741              13           35734               EPO               16 eutherian mammals EPO
                           742              13           35735               EPO                         7 primates EPO
                           743              14           35736  EPO_LOW_COVERAGE  38 eutherian mammals EPO_LOW_COVERAGE
    -------------------------------------------------------------------------------------------------------------------

Get genomic alignment for a gene region
---------------------------------------

We first get the syntenic region corresponding to human gene *BRCA2*.

.. doctest::

    >>> from ensembldb import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], release=76,
    ...                  account=account)
    >>> human_brca2 = compara.Human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> regions = compara.get_syntenic_regions(region=human_brca2, align_method='EPO', align_clade='primates')
    >>> for region in regions:
    ...     print region
    SyntenicRegions:
      Coordinate(Human,chro...,13,32315473-32400266,1)
      Coordinate(Chimp,chro...,13,31957346-32041418,-1)
      Coordinate(Macaque,chro...,17,11686607-11779396,-1)...

We then get a cogent ``Alignment`` object, requesting that sequences be annotated for gene spans.

.. doctest::

    >>> aln = region.getAlignment(feature_types='gene')
    >>> print repr(aln)
    3 x 99457 dna alignment: Homo sapiens:chromosome:13:3231...

Parsing syntenic regions
------------------------

Not all regions in a given genome have a syntenic alignment, and some have more than one alignment.
To illustrate these cases, we can consider an alignment between mouse and human, using the ``PECAN`` 
alignment method in the vertebrates clade:

.. doctest::

    >>> species = ["mouse", "human"]
    >>> compara = Compara(species, release=66)
    >>> clade = "vertebrates"
    >>> chrom, start, end, strand = "X", 155754928, 155755079, "-"
    >>> regions = compara.get_syntenic_regions(species="mouse", coord_name=chrom, 
    ...                                      start=start, end=end, align_method="PECAN", 
    ...                                      align_clade=clade, strand=strand)     
    >>> aligned_pairs = [r for r in regions]
    >>> alignment = aligned_pairs[0]                                                            
    >>> aligned_regions = [m for m in alignment.members
    ...                    if m.Region is not None]
    >>> source_region, target_region = aligned_regions
    >>> print source_region.location.coord_name, source_region.location.start, source_region.location.end
    X 155754928 155755079
    >>> print target_region.location.coord_name, target_region.location.start, target_region.location.end
    X 20222659 20223163

.. note:: We took the aligned regions from the ``regions`` generator and put them in a list for convenience.

If there are no regions returned (i.e. ``num_pairs`` is zero), then no alignment could be found. In the case of 
the above region, an exon in the *Hccs* gene, there is only one alignment. We then accessed the coordinates of the 
alignment using the ``members`` attribute of the region. Each element of ``aligned_regions`` is a ``SyntenicRegion``
instance, whose coordinates can be pulled from the ``location`` attribute.

This example shows that mouse region ``X:155754928-155755079`` aligns only to human region ``X:20222659-20223163``.

.. note:: Sometimes, the genomic coordinates given to ``get_syntenic_regions`` will contain multiple alignments between the pair of genomes, in which case two or more regions will be returned in ``aligned_pairs``.

Getting related genes
---------------------

What gene relationships are available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::

    >>> from ensembldb import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], release=76,
    ...                  account=account)
    >>> print compara.get_distinct('relationship')
    [u'gene_split', u'alt_allele', u'ortholog_one2many', u'ortholog_one2one'...

Get one-to-one orthologs
^^^^^^^^^^^^^^^^^^^^^^^^

We get the one-to-one orthologs for *BRCA2*.

.. doctest::

    >>> from ensembldb import Compara
    >>> compara = Compara(['human', 'chimp', 'macaque'], release=76,
    ...                  account=account)
    >>> orthologs = compara.get_related_genes(stableid='ENSG00000139618',
    ...                  Relationship='ortholog_one2one')
    >>> print orthologs
    RelatedGenes:
     relationships=ortholog_one2one
      Gene(species='Macaca mulatta'; biotype='protein_coding'; description=...

We iterate over the related members.

.. doctest::
    
    >>> for ortholog in orthologs.members:
    ...     print ortholog
    Gene(species='Macaca mulatta'; biotype='protein_coding'; description=...

We get statistics on the ortholog CDS lengths.

.. doctest::
    
    >>> print orthologs.get_max_cds_lengths()
    [10008, 10257, 10257]

We get the sequences as a sequence collection, with annotations for gene.

.. doctest::
    
    >>> seqs = orthologs.get_seqCollection(feature_types='gene')

Get CDS for all one-to-one orthologs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We sample all one-to-one orthologs for a group of species, generating a FASTA formatted string that can be written to file. We check all species have an ortholog and that all are translatable.

.. doctest::
    
    >>> from cogent3.core.alphabet import AlphabetError
    >>> common_names = ["mouse", "rat", "human", "opossum"]
    >>> latin_names = set([Species.get_species_name(n) for n in common_names])
    >>> latin_to_common = dict(zip(latin_names, common_names))
    >>> compara = Compara(common_names, release=76, account=account)
    >>> for gene in compara.Human.get_genes_matching(biotype='protein_coding'):
    ...     orthologs = compara.get_related_genes(gene,
    ...                                  Relationship='ortholog_one2one')
    ...     # make sure all species represented
    ...     if orthologs is None or orthologs.get_species_set() != latin_names:
    ...         continue
    ...     seqs = []
    ...     for m in orthologs.members:
    ...         try: # if sequence can't be translated, we ignore it
    ...             # get the CDS without the ending stop
    ...             seq = m.CanonicalTranscript.Cds.trim_stop_codon()
    ...             # make the sequence name
    ...             seq.name = '%s:%s:%s' % \
    ...         (latin_to_common[m.genome.species], m.stableid, m.location)
    ...             aa = seq.getTranslation()
    ...             seqs += [seq]
    ...         except (AlphabetError, AssertionError):
    ...             seqs = [] # exclude this gene
    ...             break
    ...     if len(seqs) == len(common_names):
    ...         fasta = '\n'.join(s.to_fasta() for s in seqs)
    ...         break

Get within species paralogs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> paralogs = compara.get_related_genes(stableid='ENSG00000164032',
    ...             Relationship='within_species_paralog')
    >>> print paralogs
    RelatedGenes:
     relationships=within_species_paralog
      Gene(species='Homo sapiens'; biotype='protein_coding'; description='H2A histone...

