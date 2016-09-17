********************
Accounts and Species
********************

.. _accounts:

Accounts
========

You need to specify the host name (the domain name where there's a running MySQL server hosting the Ensembl databases), username and password. Use the ``HostAccount`` class for this

.. doctest::
    
    >>> from ensembldb3 import HostAccount
    >>> account = HostAccount("mysqlhostname.anu.edu.au", "username", "password")

You can also specify the port number, e.g. ``HostAccount(..., port=5306)`` if it differs from the default 3306.

I find it convenient to specify the MySQL server account details as an environment variable called ``ENSEMBL_ACCOUNT`` by adding the following to my ``.bashrc``::

    export ENSEMBL_ACCOUNT="mysqlhostname.anu.edu.au username password"
    

In my scripts I then create the ``HostAccount`` instance as

.. doctest::
    
    >>> import os
    >>> from ensembldb3 import HostAccount
    >>> account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())

.. note::
    ``ensembldb3`` defaults to using the Ensembl UK MySQL servers if you don't specify an account. 

.. _species:

Species
=======

The ``Species`` class is a top level import that is used to translate between latin names and Ensembl's database naming scheme. It also serves to allow the user to use just a species common name to reference it's genome databases. The queries are case-insensitive.

.. doctest::

    >>> from ensembldb3 import Species
    >>> print(Species)
    =========================================================================================================
                 Common name                Species name           Ensembl Db Prefix                Synonymns
    ---------------------------------------------------------------------------------------------------------
                      Alpaca               Vicugna pacos               vicugna_pacos                         
                Amazon molly            Poecilia formosa            poecilia_formosa                         
                Anole Lizard         Anolis carolinensis         anolis_carolinensis                         ...
 
You can directly extend the list of species, or modify an existing entry, using ``Species.amend_species``. If you wish to edit the species list on a larger scale or just do it once so all your scripts can rely on that change, you can directly modify the reference species data used by ``ensembldb3``.

#. See  :ref:`exportrc` to obtain the species data distributed with ``ensembldb3`` plus other configuration files and edit the ``species.tsv`` file
#. Add an environment variable ``ENSEMBLDBRC`` to your ``.bashrc`` as follows::
    
    export ENSEMBLDBRC="~/path/to/ensembldbrc/"

.. note::
    The ``species.tsv`` file has at least 2, and up to 3, fields per line: species name, common name, species name synonymn.


Look up a species common name
-----------------------------

This can be done using the species name. Note that some synonymns are supported.

.. doctest::
    
    >>> Species.get_common_name("Felis catus")
    'Cat'
    >>> Species.get_common_name("Canis familiaris")
    'Dog'
    >>> Species.get_common_name("Canis lupus familiaris")
    'Dog'

