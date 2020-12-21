**********************
Database Optimisations
**********************

Basic Info
----------
SQL performance is influenced by the amount of data being requested and how efficiently that data is accessed. So it is important to make sure we aren't fetching data unnecessarily, and our indexes are setup correctly, matching the 'join' or 'where' clauses used in our queries.


Checking what indexes we have:
------------------------------
::

    sudo mysql
    show databases;
    use (database)
    show tables;
    show indexes from (table);
    # In particular, check the cardinalities to see how many unique values there are in each index.

Getting some sample results:
----------------------------
::

    select * from (table) limit 1;

Testing query speed:
--------------------
::

    pager cat > /dev/null	# removes the output display for more consistent times. To undo, use: "pager"
    explain select * from (table) where (column)=(value)
    # In particular, the 'rows' attribute should be low if using an index effectively, but high (close to the cardinality) if not.

Removing indexes:
-----------------
::

    drop index (index_name) on (table);  
    # Should see a drop in performance if you test the query speed now

Updating indexes:
-----------------
::

    analyze table (table_name);     # For large tables, may want to do: nohup sudo mysql (database) -e "analyze table (table_name);"

Creating indexes:
-----------------
::

    CREATE INDEX (index_name) ON (table) (column(s));
    # can specify multiple columns separated by commas if making a composite index
