****************
Querying Genomes
****************

.. _genome:

Getting a ``Genome``
====================

We use the species common name to get the genome instance we want.

.. doctest::

    >>> import os
    >>> from ensembldb3 import HostAccount, Genome
    >>> account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    >>> human = Genome(species='human', release=85, account=account)
    >>> print(human)
    Genome(species='Homo sapiens'; release='85')

.. note::

    The positions employed on Ensembl's web-site, and in their MySQL database differ from those used internally by ``ensembldb3`` -- Ensembl starts counting from 1, ``ensembldb3`` starts counting from 0. In all cases where you are querying ``ensembldb3`` objects directly inputting nucleotide positions you can indicate you are using Ensembl coordinates by setting ``ensembl_coord=True``. If you are explicitly passing in a ``ensembldb3`` region, that argument has no effect.

.. _gene:

Getting a ``Gene``
==================

.. note::
     A ``Gene`` is a type of region. All region types have ``.location`` and ``.seq`` attributes.

By biotype
----------

We can search for all genes of a specific biotype. We will use the ``limit`` argument to reduce the amount of results to be processed.

.. doctest::
    
    >>> genes = human.get_genes_matching(biotype="protein_coding", limit=2)
    >>> for gene in genes:
    ...     print(gene)
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='RING1 and YY1...'; stableid='ENSG00000281766'; status='KNOWN'; symbol='RYBP')
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='forkhead box O6...'; stableid='ENSG00000281518'; status='KNOWN'; symbol='FOXO6')

By gene symbol
--------------

We query for the *BRCA2* gene for humans.

.. doctest::

    >>> genes = human.get_genes_matching(symbol='BRCA2')
    >>> for gene in genes:
    ...     if gene.symbol == 'BRCA2':
    ...         print(gene)
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='BRCA2, DNA repair...'; stableid='ENSG00000139618'; status='KNOWN'; symbol='BRCA2')
    Gene(species='Homo sapiens'; biotype='LRG_gene'; description='BRCA2, DNA repair...'; stableid='LRG_293'; status='KNOWN'; symbol='BRCA2')

By Ensembl Stable ID
--------------------

We use the stable ID for *BRCA2*.

.. doctest::

    >>> gene = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> print(gene)
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='BRCA2, DNA repair...'; stableid='ENSG00000139618'; status='KNOWN'; symbol='BRCA2')

Matching a description
----------------------

We look for breast cancer related genes that are estrogen induced.

.. doctest::

    >>> genes = human.get_genes_matching(description='breast cancer anti-estrogen')
    >>> for gene in genes:
    ...     print(gene)
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast cancer anti-estrogen...'; stableid='ENSG00000137936';...

We can also require that an exact (case insensitive) match to the word(s) occurs within the description by setting ``like=False``.

.. doctest::

    >>> genes = human.get_genes_matching(description='breast cancer anti-estrogen',
    ...                                  like=False)
    >>> for gene in genes:
    ...     print(gene)
    Gene(species='Homo sapiens'; biotype='protein_coding'; description='breast cancer anti-estrogen...'; stableid='ENSG00000137936'; status='KNOWN'; symbol='BCAR3')...

``Gene`` attributes
-------------------

Stable ID, Symbol, Biotype, etc..
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doctest::
    
    >>> print(gene.stableid, gene.symbol, gene.biotype, gene.status)
    ENSG00000262117 BCAR4 lincRNA KNOWN

Location
^^^^^^^^

.. doctest::

    >>> gene.location
    Coordinate(Human,chro...,16,11819828-11828845,-1)
    >>> print(gene.location)
    Homo sapiens:chromosome:16:11819828-11828845:-1
    >>> print(gene.location.coord_name)
    16
    >>> print(gene.location.start)
    11819828
    >>> print(gene.location.strand)
    -1

The gene sequence
-----------------

This is an attribute of the gene.

.. doctest::
    
    >>> gene.seq
    DnaSequence(GATTCTT... 9017)

Getting a ``Transcript``
========================

This is done via the ``Gene``.

Canonical transcript
--------------------

We get the canonical transcripts for *BRCA2*.

.. doctest::

    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> transcript = brca2.canonical_transcript
    >>> print(transcript)
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315473; end=32400266; length=84793; strand='+')

Get the CDS for a transcript
----------------------------

.. doctest::

    >>> transcript = brca2.canonical_transcript
    >>> cds = transcript.cds
    >>> cds
    DnaSequence(ATGCCTA... 10257)
    >>> print(cds)
    ATGCCTATTGGATCCAAAGAGAGGCCA...

Look at all transcripts for a gene
----------------------------------

Done via the ``Gene.transcripts`` attribute

.. doctest::

    >>> for transcript in brca2.transcripts:
    ...     print(transcript)
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315473; end=32400266; length=84793; strand='+')
    Transcript(species='Homo sapiens'; coord_name='13'; start=32315504; end=32333291; length=17787; strand='+')...

Transcript exons
----------------

We show just for the canonical transcript.

.. doctest::

    >>> print(brca2.canonical_transcript.exons[0])
    Exon(stableid=ENSE00001184784, rank=1)

Get the introns for a transcript
--------------------------------

We show just for the canonical transcript.

.. doctest::

    >>> for intron in brca2.canonical_transcript.introns:
    ...     print(intron)
    Intron(TranscriptId=ENST00000380152, rank=1)
    Intron(TranscriptId=ENST00000380152, rank=2)
    Intron(TranscriptId=ENST00000380152, rank=3)...

Other region types
==================

Getting a generic genomic region
--------------------------------

Genomic regions can be obtained just using by coordinates. They can also be used to get features that lay within themselves.

.. doctest::
    
    >>> region = human.get_region(coord_name='1', start=1000000, end=1010000)
    >>> print(region)
    GenericRegion(species='Homo sapiens'; coord_name='1'; start=1000000; end=1010000; length=10000; strand='+')
    >>> region.seq
    DnaSequence(GTGGAGC... 10000)
    >>> repeats = region.get_features('repeat', limit=2)
    >>> for repeat in repeats:
    ...     print(repeat)
    Repeat(coord_name='1'; start=1000268; end=1000292; length=24; strand='+', Score=24.0)
    Repeat(coord_name='1'; start=1002182; end=1002202; length=20; strand='-', Score=0.0)

Get repeat elements in a genomic interval
-----------------------------------------

We query the genome for repeats within a specific coordinate range on chromosome 13.

.. doctest::

    >>> repeats = human.get_features(coord_name='13', start=32305473, end=32315473, feature_types='repeat', limit=2)
    >>> for repeat in repeats:
    ...     print(repeat.repeat_class)
    ...     print(repeat)
    SINE/Alu
    Repeat(coord_name='13'; start=32305225; end=32305525; length=300; strand='-', Score=2770.0)
    SINE/Alu
    Repeat(coord_name='13'; start=32305225; end=32305525; length=300; strand='-', Score=2770.0)

Get CpG island elements in a genomic interval
---------------------------------------------

We query the genome for CpG islands within a specific coordinate range on chromosome 11.

.. doctest::

    >>> islands = human.get_features(coord_name='11', start=2129111, end=2149604, feature_types='cpg', limit=2)
    >>> for island in islands:
    ...     print(island)
    CpGisland(coord_name='11'; start=2137721; end=2141254; length=3533; strand='-', Score=3254.0)
    CpGisland(coord_name='11'; start=2143905; end=2144442; length=537; strand='-', Score=652.0)

