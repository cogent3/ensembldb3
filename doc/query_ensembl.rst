.. _query-ensembl:

Querying Ensembl
================

.. sectionauthor:: Gavin Huttley, Hua Ying

NEEDS TO BE COMPLETELY REWRITTEN


We begin this documentation with a note on dependencies, performance and code status. Ensembl_ makes their data available via MySQL servers, so the ``cogent.db.ensembl`` module has additional dependencies of `mysql-connector-python`_ or `MySQL-python`_ and SQLAlchemy_. You can use pip_ to install both `mysql-connector-python` and SQLAlchemy_, but MySQL-python is more involved. If you experience trouble, please create a new issue on the `issue tracker <https://github.com/pycogent/pycogent/issues>`_ and tag it with the *question* tag. Regarding performance, significant queries to the UK servers from Australia are very slow. The examples in this documentation, for instance, take ~15 minutes to run when pointed at the UK servers. Running against a local installation, however, is ~50 times faster. On status, the ``cogent.db.ensembl`` module should be considered beta level code. We still strongly advise users to check results for a subset of their analyses against the data from the UK Ensembl web site.

.. _`MySQL-python`: http://sourceforge.net/projects/mysql-python
.. _`mysql-connector-python`: https://pypi.python.org/pypi/mysql-connector-python
.. _SQLAlchemy: http://www.sqlalchemy.org/
.. _pip: https://pypi.python.org/pypi/pip

The major components of Ensembl_ are compara and individual genomes. In all cases extracting data requires connecting to MySQL databases on a server and the server may be located remotely or locally. For convenience, the critical objects you'll need to query a database are provided in the top-level module import, ie immediately under ``cogent.db.ensembl``.

.. _Ensembl: http://www.ensembl.org

Specifying a Host and Account
-----------------------------

So the first step is to specify what host and account are to be used. On my lab's machines, I have set an environment variable with the username and password for our local installation of the Ensembl_ MySQL databases, e.g. ``ENSEMBL_ACCOUNT="host username password"``. So I'll check for that (since the documentation runs much quicker when this is true) and if it's absent, we just set ``account=None`` and the account used defaults to the UK Ensembl_ service. I also define which release of ensembl we'll use in one place to allow easier updating of this documentation.

.. doctest::
    
    >>> import os
    >>> release = 76
    >>> from ensembldb import HostAccount
    >>> if 'ENSEMBL_ACCOUNT' in os.environ:
    ...     host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
    ...     account = HostAccount(host, username, password)
    ... else:
    ...     account = None

What Species Are Available?
---------------------------

Another key element, of course, is the species available. Included as part of ``cogent.db.ensembl`` is the module ``species``. This module contains a class that translates between latin names, common names and ensembl database prefixes. The species listed are guaranteed to be incomplete, given Ensembl's update schedule, so it's possible to dynamically add to this listing, or even change the common name for a given latin name.

.. doctest::

    >>> from ensembldb import Species
    >>> print Species
    ================================================================================
           Common Name                   Species Name              Ensembl Db Prefix
    --------------------------------------------------------------------------------
             A.aegypti                  Aedes aegypti                  aedes_aegypti
            A.clavatus           Aspergillus clavatus           aspergillus_clavatus...

In Australia, the common name for *Gallus gallus* is chook, so I'll modify that.

.. doctest::

    >>> Species.amend_species('Gallus gallus', 'chook')
    >>> assert Species.get_common_name('Gallus gallus') == 'chook'

You can also add new species for when they become available using ``Species.amend_species``.

Species common names are used to construct attributes on PyCogent ``Compara`` instances). You can get the name that will be using the ``get_compara_name`` method. For species with a real common name

.. doctest::
    
    >>> Species.get_compara_name('Procavia capensis')
    'RockHyrax'

or with a shortened species name

.. doctest::
    
    >>> Species.get_compara_name('Caenorhabditis remanei')
    'Cremanei'

The ``Species`` class is basically used to translate between latin names and ensembl's database naming scheme. It also serves to allow the user to simply enter the common name for a species in order to reference it's genome databases. The queries are case-insensitive.

Interrogating a Genome
----------------------

As implied above, Ensembl databases are versioned, hence you must explicitly state what release you want. Aside from that, getting an object for querying a genome is simply a matter of importing the ``HostAccount`` and ``Genome`` classes. Here I'm going to use the ``cogent.db.ensembl`` level imports.

.. doctest::

    >>> from ensembldb import HostAccount, Genome
    >>> human = Genome(species='human', release=release, account=account)
    >>> print human
    Genome(species='Homo sapiens'; release='76')

Notice I used the common name rather than full name. The ``Genome`` provides an interface to obtaining different attributes. It's primary role is to allow selection of genomic regions according to some search criteria. The type of region is presently limited to ``Gene``, ``Est``, ``CpGisland``, ``Repeat`` and ``Variation``. There's also a ``GenericRegion``. The specific types are also capable of identifying information related to themselves, as we will demonstrate below.

A Note on Coordinate Systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The positions employed on Ensembl's web-site, and in their MySQL database differ from those used internally by ``cogent.db.ensembl``. In all cases where you are querying ``cogent.db.ensembl`` objects directly inputting nucleotide positions you can indicate you are using Ensembl coordinates by setting ``ensembl_coord=True``. If you are explicitly passing in a ``cogent.db.ensembl`` region, that argument has no effect.

Selecting Gene's
^^^^^^^^^^^^^^^^

The genome can be queried for gene's in a number of ways. You can search for genes using the ``Genome.get_gene_by_stableid`` method which requires you know the Ensembl stable id.

.. doctest::
    
    >>> brca1 = human.get_gene_by_stableid(stableid='ENSG00000012048')
    >>> print brca1.description
    breast cancer 1, early onset...

Alternatively, you can query using the ``Genome.get_genes_matching`` method. This method allows querying for gene(s) by the following identifiers: HGNC symbol; Ensembl ``stable_id``; description; or coding type.

.. note:: When querying by description, you can specify that the exact words in the query must be present in the description by setting the argument ``like=True``. The default is ``like=False``.

In general for such queries, case shouldn't matter. For instance, find the *BRCA2* gene by it's HGNC symbol.

.. doctest::

    >>> genes = human.get_genes_matching(symbol='brca2')

Because there can be multiple hits from a ``get_genes_matching`` query, and because we wish to not spend time doing things (like talking to the database) unnecessarily, the result of the query is a python generator. This acts like a series and allows you to iterate over the database hits until you find the one you want and then terminate the record collection.

.. doctest::

    >>> for gene in genes:
    ...     if gene.symbol.lower() == 'brca2':
    ...         break
    ...
    >>> brca2 = gene # so we keep track of this reference for later on
    >>> print brca2.symbol
    BRCA2
    >>> print brca2.description
    breast cancer 2...
    >>> print brca2
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast...

This code serves to illustrate a few things. First, the sorts of properties that exist on the object. These can be directly accessed as illustrated above. Secondly, that the argument names to ``get_genes_matching`` match the properties.

Gene's also have a location. The length of a gene is the difference between its start and end location.

.. note:: Unfortunately all gene coordinates can vary between genome builds. So start, end and length can all differ between Ensembl releases for the same gene.

.. doctest::

    >>> print brca2.location
    Homo sapiens:chromosome:13:32315473...
    >>> print len(brca2)
    84793
    
Each location is directly tied to the parent genome and the coordinate above also shows the coordinates' *type* (chromosome in this case), name (13), start, end and strand. The start and end positions are python indices and will differ from the Ensembl indices in that start will be the Ensembl index - 1. This is because python counts from 0, not 1. In querying for regions using a specific set of coordinates, it is possible to put in the Ensembl coordinates (demonstrated below).

``Gene`` has several useful properties, including the ability to directly get their own DNA sequence and their ``canonical_transcript`` and ``Transcripts``. ``canonical_transcript`` is the characteristic transcript for a gene, as defined by Ensembl. ``Transcripts`` is a ``tuple`` attribute containing individual region instances of type ``Transcript``. A ``Transcript`` has ``Exons``, ``Introns``, a ``Cds`` and, if the ``biotype`` is protein coding, a protein sequence. In the following we grab the cannonical transcript from ``brca2``

.. doctest::

    >>> print brca2.biotype
    protein_coding
    >>> print brca2.seq
    GGGCTTGTGGCGC...
    >>> print brca2.canonical_transcript.Cds
    ATGCCTATTGGATC...
    >>> print brca2.canonical_transcript.ProteinSeq
    MPIGSKERPTF...

It is also possible to iterate over a transcript's exons, over their translated exons, or to obtain their coding DNA sequence. We grab the second transcript for this.

.. doctest::
    
    >>> transcript = brca2.Transcripts[0]
    >>> for exon in transcript.Exons:
    ...     print exon, exon.location
    Exon(stableid=ENSE00001184784, rank=1) Homo sapiens:chromosome:13:...
    >>> for exon in transcript.TranslatedExons:
    ...     print exon, exon.location
    Exon(stableid=ENSE00001484009, rank=2) Homo sapiens:chromosome:13:...
    >>> print transcript.Cds
    ATGCCTATTGGATCCAAA...

The ``Cds`` sequence includes the stop-codon, if present. The reason for this is there are many annotated transcripts in the Ensembl database the length of whose transcribed exons are not divisible by 3. Hence we leave it to the user to decide how to deal with that, but mention here that determining the number of complete codons is trivial and you can slice the ``Cds`` so that it's length is divisible by 3.

The ``Exons`` and ``TranslatedExons`` properties are tuples that are evaluated on demand and can be sliced. Each ``Exon/TranslatedExon`` is itself a region, with all of the properties of generic regions (like having a ``seq`` attribute). Similar descriptions apply to the ``Introns`` property and ``Intron`` class. We show just for the canonical transcript.

.. doctest::

    >>> for intron in brca2.canonical_transcript.Introns:
    ...     print intron
    Intron(TranscriptId=ENST00000380152, rank=1)
    Intron(TranscriptId=ENST00000380152, rank=2)
    Intron(TranscriptId=ENST00000380152, rank=3)...


The ``Gene`` region also has convenience methods for examining properties of it's transcripts, in presenting the ``Cds`` lengths and getting the ``Transcript`` encoding the longest ``Cds``.

.. doctest::

    >>> print brca2.get_cds_lengths()
    [10257, 10257]
    >>> longest = brca2.get_longest_cds_transcript()
    >>> print longest.Cds
    ATGCCTATTGGATCCAAA...

All Regions have a ``get_features`` method which differs from that on genome only in that the genomic coordinates are automatically entered for you. Regions also have the ability to return their sequence as an annotated ``cogent`` sequence. The method on ``Gene`` simply queries the parent genome using the gene's own location as the coordinate for the currently supported region types. We will query ``brca2`` asking for gene features, the end-result will be a ``cogent`` sequence that can be used to obtain the CDS, for instance, using the standard ``cogent`` annotation capabilities.

.. doctest::

    >>> annot_brca2 = brca2.get_annotated_seq(feature_types='gene')
    >>> cds = annot_brca2.get_annotations_matching('CDS')[0].get_slice()
    >>> print cds
    ATGCCTATTGGATCCAAA...

Those are the properties of a ``Gene``, at present, of direct interest to end-users.

There are obviously different types of genes, and the ``Genome`` object provides an ability to establish exactly what distinct types are defined in Ensembl.

.. doctest::

    >>> print human.get_distinct('biotype')
    [u'unitary_pseudogene', u'rRNA', u'lincRNA'...

The genome can be queried for any of these types, for instance we'll query for ``rRNA``. We'll get the first few records and then exit.

.. doctest::

    >>> rRNA_genes = human.get_genes_matching(biotype='rRNA')
    >>> count = 0
    >>> for gene in rRNA_genes:
    ...     print gene
    ...     count += 1
    ...     if count == 1:
    ...         break
    ...
    Gene(species='Homo sapiens'; biotype='rRNA'; description='RNA, 5S...

This has the effect of returning any gene whose ``biotype`` includes the phrase ``rRNA``. If a gene is not a protein coding gene, as in the current case, then it's ``Transcripts`` will have ``ProteinSeq==None`` and ``TranslatedExons==None``, but it will have ``Exons`` and a ``Cds``.

.. doctest::

    >>> transcript = gene.Transcripts[0]
    >>> assert transcript.ProteinSeq == None
    >>> assert transcript.TranslatedExons == None
    >>> assert transcript.Cds != None

Getting ESTs
^^^^^^^^^^^^

Ensembl's ``otherfeatures`` database mirrors the structure of the ``core`` database and contains EST information. Hence, the ``Est`` region inherits directly from ``Gene`` (ie has many of the same properties). ``est`` is a supported ``feature_types`` for the ``get_features`` method. You can also directly query for an EST using Ensembl's ``StableID``. Here, however, we'll just query for ``Est`` that map to the ``brca2`` region.

.. doctest::

    >>> ests = human.get_features(feature_types='est', region=brca2)
    >>> for est in ests:
    ...     print est
    Est(species='Homo sapiens'; biotype='protein_coding'; description='None';...

Getting Variation
^^^^^^^^^^^^^^^^^

``Variation`` regions also have distinctive properties worthy of additional mention. As for genes, there are distinct types stored in Ensembl that may be of interest. Those types can likewise be discovered from the genome,

.. doctest::

    >>> print human.get_distinct('effect')
    ['3_prime_UTR_variant', 'splice_acceptor_variant', 'intergenic_variant'...

and that information can be used to query the genome for all variation of that effect. 

.. note:: What we term ``effect``, Ensembl terms consequence. We use ``effect`` because it's shorter.

We allow the query to be an inexact match by setting ``like=True``. Again we'll just iterate over the first few.

.. doctest::

    >>> nsyn_variants = human.get_variation(effect='missense_variant')
    >>> for i, nsyn_variant in enumerate(nsyn_variants):
    ...     if nsyn_variant.effect == 'missense_variant' and\
    ...                          nsyn_variant.AlleleFreqs:
    ...         break
    ...     
    >>> print nsyn_variant
    Variation(symbol='rs1638319'; effect='missense_variant'; Alleles='A/G')
    >>> print nsyn_variant.AlleleFreqs
    =================================
    allele      freq    population_id
    ---------------------------------
         A    0.5000             9088
         G    0.5000             9088
    ---------------------------------

``Variation`` objects also have other useful properties, such as a location, the number of alleles and the allele frequencies. The length of a ``Variation`` instance is the length of it's longest allele.

.. doctest::

    >>> assert len(nsyn_variant) == 1
    >>> print nsyn_variant.location
    Homo sapiens:chromosome:1:69967-69968:1
    >>> assert nsyn_variant.NumAlleles == 2

``Variation`` objects have ``FlankingSeq`` and ``seq`` attributes which in the case of a SNP is a single nucleotide long and should correspond to one of the alleles. In the latter case, this property is a tuple with the 0th entry being the 5'- 300 nucleotides and the 1st entry being the 3' nucleotides.

.. note:: The flanking sequence is only returned when the SNPs flank matches reference (according to Ensembl).

.. doctest::

    >>> print nsyn_variant.FlankingSeq[0]
    TTGCTAACAGT...
    >>> print nsyn_variant.FlankingSeq[1]
    GCTGAGAAAAT...
    >>> assert str(nsyn_variant.seq) in nsyn_variant.Alleles, str(nsyn_variant.seq)

As a standard feature, ``Variation`` within a specific interval can also be obtained. Using the ``brca2`` gene region instance created above, we can find all the genetic variants using the ``variants`` property of genome regions. We use this example to also demonstrate the ``PeptideAlleles`` and ``TranslationLocation`` attributes. ``PeptideAlleles`` is the amino-acid variation resulting from the nucleotide variation while ``TranslationLocation`` is the position in the translated peptide of the variant. If a variant does not affect protein coding sequence (either it's not exonic or it's a synonymous variant) then these properties have the value ``None``.
We illustrate their use.

.. doctest::

    >>> for variant in brca2.variants:
    ...     if variant.PeptideAlleles is None:
    ...         continue
    ...     print variant.PeptideAlleles, variant.TranslationLocation
    P/L 1...

.. note:: These are Python coordinates, add 1 to get the Ensembl value.

We can also use a slightly more involved query to find all variants within the gene of a specific type. (Of course, you could also simply iterate over the ``variants`` attribute to grab these out too.)

.. doctest::

    >>> brca2_snps = human.get_features(feature_types='variation',
    ...                      region=brca2)
    >>> for snp in brca2_snps:
    ...     if 'missense_variant' in snp.effect:
    ...         break
    >>> print snp
    Variation(symbol='rs80358836'; effect=['non_coding_exon_variant', 'nc_transcript_variant', 'upstream_gene_variant', 'missense_variant'...
    >>> print snp.location
    Homo sapiens:chromosome:13:32316464-32316465:1


Other Region Types
^^^^^^^^^^^^^^^^^^

These can be obtained from the genome instance using the genomes ``get_features`` method. At present, only repeats, CpG islands, variation, EST's and genes can be obtained through this method. There's also ``GenericRegion``, which is precisely that.

In Ensembl's databases, each type of feature may be recorded at multiple coordinate levels. Accordingly, each level is checked to obtain full information of that feature. 

.. doctest::

   >>> chicken = Genome(species='chook', release=release, account=account)
   >>> print chicken.feature_coord_levels
   Gallus gallus
   ============================================
        type                             levels
   --------------------------------------------
        gene                         chromosome
      repeat                             contig
         est                         chromosome
   variation                         chromosome
         cpg    chromosome, supercontig, contig
   --------------------------------------------

Comparative Analyses
--------------------

The Ensembl compara database is represented by ``cogent.db.ensembl.compara.Compara``. This object provides a means for querying for relationships among genomes and obtaining multiple alignments. For convenience the class is made available through the top-level module for importing  (i.e. ``cogent.db.ensembl.Compara``). Instantiating ``Compara`` requires, as before, the ensembl release, the series of species of interest and optionally an account (we also use our local account for speed). For the purpose of illustration we'll use the human, mouse and rat genomes.

.. note:: Any queries on this instance of compara will only return results for the indicated species. If you want to query about other species, create another instance.

.. doctest::

    >>> from ensembldb import Compara
    >>> compara = Compara(['human', 'mouse', 'rat'], account=account,
    ...                  release=release)
    >>> print compara
    Compara(species=('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'); release=76...

The ``Compara`` object loads the corresponding ``Genome``'s and attaches them to itself as named attributes (use ``Species.get_compara_name`` to find out what the attribute will be). The genome instances are named according to their common name in CamelCase, or Scase. For instance, if we had created a ``Compara`` instance with the American pika species included, then that genome would be accessed as ``compara.AmericanPika``. Common names containing a '.' are treated differently. For instance, the common name for *Caenorhabditis remanei* is ``C.remanei`` which becomes ``compara.Cremanei``. We access the human genome in this ``Compara`` instance and conduct a gene search.

.. doctest::

    >>> brca2 = compara.Human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> print brca2
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast...

We can now use this result to search compara for related genes. We note here that like ``Genome``, ``Compara`` has the ``get_distinct`` method to assist in identifying appropriate search criteria. What are the distinct types of gene relationships recorded in Ensembl, for instance?

.. doctest::

    >>> relationships = compara.get_distinct('relationship')
    >>> print relationships
    [u'gene_split', u'alt_allele', u'ortholog_one2many', u'ortholog_one2one'...

So we use the ``brca2`` instance above and search for orthologs among the human, mouse, rat genomes.

.. doctest::

    >>> orthologs = compara.get_related_genes(gene_region=brca2,
    ...                 Relationship='ortholog_one2one')
    >>> print orthologs
    RelatedGenes:
     relationships=ortholog_one2one
      Gene(species='Rattus norvegicus'; biotype='protein_coding'; description='breast cancer ...

I could also have done that query using a ``stableid``, which I now do using the Ensembl mouse identifier for *Brca2*.

.. doctest::

    >>> orthologs = compara.get_related_genes(stableid='ENSMUSG00000041147',
    ...                 Relationship='ortholog_one2one')
    >>> print orthologs
    RelatedGenes:
     relationships=ortholog_one2one
      Gene(species='Rattus norvegicus'; biotype='protein_coding'; description='breast cancer...

The ``RelatedGenes`` object has a number of properties allowing you to get access to data. A ``members`` attribute holds each of the ``Gene`` instances displayed above. The length of this attribute tells you how many hits there were, while each member has all of the capabilities described for ``Gene`` above, eg. a ``Cds`` property. There is also a ``get_seqLengths`` method which returns the vector of sequence lengths for the members. This method returns just the lengths of the individual genes.

.. doctest::

    >>> print orthologs.members
    (Gene(species='Rattus norvegicus'; biotype='protein_coding'; Descr...
    >>> print orthologs.get_seqLengths()
    [40748, 84793, 47117]

In addition there's a ``get_max_cds_lengths`` method for returning the lengths of the longest ``Cds`` from each member.

.. doctest::

    >>> print orthologs.get_max_cds_lengths()
    [10032, 10257, 9990]

You can also obtain the sequences as a ``cogent`` ``SequenceCollection`` (unaligned), with the ability to have those sequences annotated as described above. The sequences are named in accordance with their genomic coordinates.

.. doctest::

    >>> seqs = orthologs.get_seqCollection(feature_types='gene')
    >>> print seqs.names
    ['Rattus norvegicus:chromosome:12:491...

We can also search for other relationship types, which we do here for a histone.

.. doctest::

    >>> paralogs = compara.get_related_genes(stableid='ENSG00000164032',
    ...             Relationship='within_species_paralog')
    >>> print paralogs
    RelatedGenes:
     relationships=within_species_paralog
      Gene(species='Homo sapiens'; biotype='protein_coding'; description='H2A...

Getting Comparative Alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ensembl stores multiple sequence alignments for selected species. For a given group of species, you can examine what alignments are available by printing the ``method_species_links`` attribute of ``Compara``. This will return something like

    >>> print compara.method_species_links
    Align Methods/Clades
    ============================================================================...
    method_link_species_set_id  method_link_id  species_set_id      align_method...
    ----------------------------------------------------------------------------...
                           753              10           35883             PECAN...
                           741              13           35734               EPO...
                           743              14           35736  EPO_LOW_COVERAGE...
    ----------------------------------------------------------------------------...

The ``align_method`` and ``align_clade`` columns can be used as arguments to ``get_syntenic_regions``. This method is responsible for returning ``SyntenicRegions`` instances for a given coordinate from a species. As it's possible that multiple records may be found from the multiple alignment for a given set of coordinates, the result of calling this method is a python generator. The returned regions have a length, defined by the full set of aligned sequences. If the ``omit_redundant`` argument is used, then positions with gaps in all sampled species will be removed in the alignment to be returned. The length of the syntenic region, however, is the length of the unfiltered alignment.

.. note:: It's important to realise that multiple alignments are from these clades. Hence, sequence regions that you might expect would result in a contiguous alignment in the species subset of interest may be returned as separate ``SyntenicRegions`` due to the influence on the alignment of the other species.

.. doctest::

    >>> syntenic_regions = compara.get_syntenic_regions(region=brca2,
    ...                      align_method='EPO', align_clade='eutherian')
    >>> for syntenic_region in syntenic_regions:
    ...     print syntenic_region
    ...     print len(syntenic_region)
    ...     print repr(syntenic_region.getAlignment(omit_redundant=False))
    SyntenicRegions:
      Coordinate(Human,chro...,13,32316515-32356518,1)
      Coordinate(Mouse,chro...,5,150523267-150550922,1)
      Coordinate(Rat,chro...,12,491349-514449,1)
    118454
    3 x 118454 dna alignment: Homo sapiens:chromosome:13:32316515-32356518...

We consider a species for which pairwise alignments are available -- the bush baby.

.. doctest::

    >>> compara_pair = Compara(['Human', 'Bushbaby'], release=release,
    ...                        account=account)
    >>> print compara_pair
    Compara(species=('Homo sapiens', 'Otolemur garnettii'); release=76; connected=True)


Printing the ``method_species_links`` table provides all the necessary information for specifying selection conditions.

    >>> print compara_pair.method_species_links
    Align Methods/Clades
    ============================================================================...
    method_link_species_set_id  method_link_id  species_set_id      align_method...
    ----------------------------------------------------------------------------...
                           743              14           35736  EPO_LOW_COVERAGE...
                           682              16           35678         LASTZ_NET...
    ----------------------------------------------------------------------------...

.. doctest::
    
    >>> gene = compara_pair.Bushbaby.get_gene_by_stableid(
    ...                             stableid='ENSOGAG00000003166')
    ...
    >>> print gene
    Gene(species='Otolemur garnettii'; biotype='protein_coding'...
    >>> syntenic = compara_pair.get_syntenic_regions(region=gene,
    ...          align_method='LASTZ_NET', align_clade='H.sap-O.gar')
    ...
    >>> for region in syntenic:
    ...     print region
    ...     break
    SyntenicRegions:
      Coordinate(Bushbaby,scaf...,GL87...,8624867-8626121,1)
      Coordinate(Human,chro...,7,135726146-135727496,1)