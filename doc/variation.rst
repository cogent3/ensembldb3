******************
Querying Variation
******************

Variation databases are attributes of a species ``Genome``, so we create an instance of the human genome for our examples.

.. doctest::

    >>> import os
    >>> from ensembldb3 import HostAccount, Genome
    >>> account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    >>> human = Genome(species='human', release=85, account=account)

What variation effects are there?
=================================

There are distinct effect types stored in Ensembl which can be discovered from the genome,

.. doctest::

    >>> print(human.get_distinct('effect')) # doctest: +SKIP
    ['3_prime_UTR_variant', 'missense_variant', 'non_coding_transcript_exon_variant',...

.. note::
    What we term ``effect``, Ensembl terms consequence. We use ``effect`` because it's shorter.

Get variants by effect
======================

Genetic variants are represented by the region type class ``Variation``. We query the genome for a specific effect type, using tyhe limit argument in the interests of speed.

.. doctest::

    >>> snps = human.get_variation(effect='missense_variant', limit=3)
    >>> for snp in snps:
    ...     print(snp)
    ...
    Variation(symbol='rs2853516'; effect=['downstream_gene_variant', 'missense_variant', 'upstream_gene_variant']; alleles='G/A')
    Variation(symbol='rs28416101'; effect=['downstream_gene_variant', 'missense_variant', 'upstream_gene_variant']; alleles='T/G')
    Variation(symbol='rs201969351'; effect=['downstream_gene_variant', 'missense_variant', 'upstream_gene_variant']; alleles='T/C')

Get variantion from a genomic region
====================================

We can also use a slightly more involved query to find all variants within the gene of a specific type. (Of course, you could also simply iterate over the ``variants`` attribute to grab these out too.)

.. doctest::

    >>> region = human.get_region(coord_name='16', start=1000000, end=2000000)
    >>> snps = human.get_features(feature_types='variation',
    ...                      region=region, limit=2)
    >>> for snp in snps:
    ...     print(snp)
    ...     
    Variation(symbol='rs560025441'; effect=['upstream_gene_variant', 'non_coding_transcript_variant', 'intron_variant']; alleles='T/C')
    Variation(symbol='rs773514182'; effect=['upstream_gene_variant', 'non_coding_transcript_variant', 'intron_variant']; alleles='A/T')

Get variants from a gene
========================

Gene's have a ``variation`` attribute.

.. doctest::

    >>> brca2 = human.get_gene_by_stableid(stableid='ENSG00000139618')
    >>> for snp in brca2.variants:
    ...     if snp.peptide_alleles is not None:
    ...         break # so we can demonstrate stuff below

Attributes of ``Variation``
===========================

.. doctest::
    
    >>> snp.location, snp.num_alleles, snp.alleles
    (Coordinate(Human,chro...,13,32316464-32316465,1), 2, 'C/T')

They can also have peptide alleles and position of the polymorphism in the translated peptide.

.. doctest::
    
    >>> print(snp.peptide_alleles, snp.translation_location)
    P/L 1...

.. note::
    If a variant does not affect protein coding sequence (either it's not exonic or it's a synonymous variant) then these properties have the value ``None``.

They have length, so single nucleotide polymorphisms have alength of 1. (The length is the length of it's longest allele.)

.. doctest::

    >>> len(snp)
    1

Flanking sequence of a variant
------------------------------

``Variation`` objects have a ``seq`` attribute which in the case of a SNP is a single nucleotide long and should correspond to one of the alleles. They also have a ``flanking_seq`` attribute. This property is a tuple with the 0th entry being the 5'- 300 nucleotides and the 1st entry being the 3' nucleotides.

.. doctest::

    >>> print(snp.flanking_seq[0]) # 5' flank
    TGATTGAAACT...
    >>> print(snp.flanking_seq[1]) # 3' flank
    TATTGGATCCA...

.. note::
    The flanking sequence is only returned when the SNPs flank matches reference.

Getting allele frequencies
--------------------------

``Variation`` objects have an allele frequency table.

.. doctest::
    
    >>> snp = list(human.get_variation(symbol="rs55880202"))[0]
    >>> print(snp.allele_freqs)
    =================================
    allele      freq    population_id
    ---------------------------------
         C                       5474
         T                       5474
         C    0.9661            11888
         T    0.0339            11888
    ---------------------------------

Validation status
-----------------

.. doctest::
    
    >>> snp.validation
    {'Phenotype_or_Disease', 'Frequency', '1000Genomes'}
