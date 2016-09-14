****************
Advanced Queries
****************

Direct SQL based queries
========================

You can use ``ensembldb3`` to perform your own custom SQL queries of the Ensembl databases. ``ensembldb3`` relies on `SQLAlchemy <http://www.sqlalchemy.org/>`_ for SQL. Hence, you need to learn how to use the SQLAlchemy `expression language <http://docs.sqlalchemy.org/en/latest/core/tutorial.html>`_.

.. note::
    The `Ensembl MySQL schema <http://asia.ensembl.org/info/docs/api/core/core_schema.html>`_ is the essential reference.

To start with, we create a ``Genome`` instance.

.. doctest::

    >>> import os
    >>> from pprint import pprint
    >>> import sqlalchemy as sql
    >>> from ensembldb3 import Genome, HostAccount
    >>> account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    >>> human = Genome('human', release=85, account=account, pool_recycle=10000)

We perform a custom query of the human variation database, sampling non-somatic SNPs whose flanks match the reference genome. We construct our ``where`` clause and use that in a ``select`` statement that will be performed on the Ensembl's ``variation_feature`` table. In the interests of computational speed we will limit the query to 1 record and display it using ``pprint`` to make it more legible.

.. doctest::

    >>> variation_feature = human.VarDb.get_table("variation_feature")
    >>> whereclause = sql.and_(variation_feature.c.somatic==0, variation_feature.c.alignment_quality==1)
    >>> query = sql.select([variation_feature], whereclause)
    >>> query = query.limit(1)
    >>> for record in query.execute():
    ...     pprint(record.items())
    [('variation_feature_id', 19004134),
     ('seq_region_id', 131537),
     ('seq_region_start', 60360),
     ('seq_region_end', 60360),
     ('seq_region_strand', 1),
     ('variation_id', 20065981),
     ('allele_string', 'C/G'),
     ('variation_name', 'rs111660247'),
     ('map_weight', 1),
     ('flags', None),
     ('source_id', 1),
     ('consequence_types', {'downstream_gene_variant'}),
     ('variation_set_id', set()),
     ('class_attrib_id', 2),
     ('somatic', 0),
     ('minor_allele', None),
     ('minor_allele_freq', None),
     ('minor_allele_count', None),
     ('alignment_quality', Decimal('1.0000000000')),
     ('evidence_attribs', None),
     ('clinical_significance', None),
     ('display', 1)]

.. note::
    Values in the SQLAlchemy result objects can be accessed in multiple ways, including like a dictionary.
