import sqlalchemy as sql
from cogent3.core.location import Map

from .species import Species as _Species
from .util import asserted_one, convert_strand, DisplayString
from .host import DbConnection

__author__ = "Hua Ying"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "BSD"
__version__ = "1.0.a"
__maintainer__ = "Hua Ying"
__email__ = "Hua.Ying@anu.edu.au"
__status__ = "alpha"


def location_query(table, query_start, query_end,
                   start_col='seq_region_start', end_col='seq_region_end', query=None,
                   where='overlap'):
    # TODO should we allow for spans, overlaps, within?
    # the union result is a complex query and has be appended to any other queries
    # in which it's being employed
    # should we be setting default values here regarding the columns that start/end
    # are pulled from, or explicitly state which columns
    if query is None:
        query = sql.select([table])

    if where == 'within':
        query.append_whereclause(sql.and_(table.c[start_col] < query_start,
                                          table.c[end_col] > query_end))
    else:
        query.append_whereclause(
            sql.or_(sql.and_(table.c[start_col] < query_start,
                             table.c[end_col] > query_end),
                    sql.and_(table.c[start_col] >= query_start,
                             table.c[start_col] <= query_end),
                    sql.and_(table.c[end_col] >= query_start,
                             table.c[end_col] <= query_end)))
    # the union is only being used here to order the results
    # that usage imposes the limitation this function must be appended to
    # other queries components being built into a fuller SQL query
    # makes me think it shouldn't be here?
    query = query.order_by(table.c[start_col])
    return query


def _get_coord_type_and_seq_region_id(coord_name, core_db):
    seq_region_table = core_db.get_table('seq_region')
    rows = sql.select([seq_region_table]).\
        where(seq_region_table.c.name == str(coord_name)).execute().fetchall()
    species_coord_sys = CoordSystem(species=core_db.db_name.species,
                                    core_db=core_db)
    try:
        selected_row = asserted_one(rows)
    except ValueError:
        selected_row = None
        for row in rows:
            # not a default_version
            if not row['coord_system_id'] in species_coord_sys:
                continue
            elif not selected_row:
                selected_row = row
                break
        if selected_row is None:
            raise ValueError("Ambigous coordinate name: %s" % coord_name)
    coord_type = species_coord_sys[selected_row['coord_system_id']].name
    return selected_row, coord_type


class Coordinate(object):

    def __init__(self, genome, coord_name, start, end, strand=1,
                 CoordType=None, seq_region_id=None, ensembl_coord=False):
        if not CoordType or not (seq_region_id or start or end):
            seq_region_data, CoordType = \
                _get_coord_type_and_seq_region_id(coord_name, genome.CoreDb)
            seq_region_id = seq_region_data['seq_region_id']
            start = start or 0
            end = end or seq_region_data['length']
        # TODO allow creation with just seq_region_id
        self.species = genome.species
        self.CoordType = DisplayString(CoordType, repr_length=4,
                                       with_quotes=False)
        self.coord_name = DisplayString(coord_name, repr_length=4,
                                       with_quotes=False)
        # if start == end, we +1 to end, unless these are ensembl_coord's
        if ensembl_coord:
            start -= 1
        elif start == end:
            end += 1

        if start > end:
            assert strand == -1,\
                "strand incorrect for start[%s] > end[%s]" % (start, end)
            start, end = end, start

        self.start = start
        self.end = end
        self.strand = convert_strand(strand)
        self.seq_region_id = seq_region_id
        self.genome = genome

    def __len__(self):
        return self.end - self.start

    def __lt__(self, other):
        return (self.coord_name, self.start) < (other.coord_name, other.start)

    def __eq__(self, other):
        return (self.coord_name, self.start) == (other.coord_name, other.start)

    def _get_ensembl_start(self):
        # ensembl counting starts from 1
        return self.start + 1

    ensembl_start = property(_get_ensembl_start)

    def _get_ensembl_end(self):
        return self.end

    ensembl_end = property(_get_ensembl_end)

    def __str__(self):
        return '%s:%s:%s:%d-%d:%d' % (self.species, self.CoordType,
                                      self.coord_name, self.start, self.end, self.strand)

    def __repr__(self):
        my_type = self.__class__.__name__
        name = _Species.get_common_name(self.species)
        coord_type = self.CoordType
        c = '%s(%r,%r,%r,%d-%d,%d)' % (my_type, name, coord_type,
                                       self.coord_name, self.start, self.end, self.strand)
        return c.replace("'", "")

    def adopted(self, other, shift=False):
        """adopts the seq_region_id (including coord_name and CoordType) of
        another coordinate.

        Arguments:
            - shift: an int or True/False. If int, it's added to start/end.
              If bool, other.start is added to start/end"""
        if type(shift) == bool:
            shift = [0, other.start][shift]
        return self.__class__(other.genome, coord_name=other.coord_name,
                              start=self.start + shift, end=self.end + shift,
                              strand=other.strand,
                              seq_region_id=other.seq_region_id)

    def shifted(self, value):
        """adds value to start/end coords, returning a new instance."""
        new = self.copy()
        new.start += value
        new.end += value
        assert len(new) > 0, 'shift generated a negative length'
        return new

    def copy(self):
        """returns a copy"""
        return self.__class__(genome=self.genome, coord_name=self.coord_name,
                              start=self.start, end=self.end, strand=self.strand,
                              CoordType=self.CoordType, seq_region_id=self.seq_region_id)

    def resized(self, from_start, from_end):
        """returns a new resized Coordinate with the
        start=self.start+from_start and end = self.end+from_end.

        If you want to shift start upstream, add a -ve number"""
        new = self.copy()
        new.start += from_start
        new.end += from_end
        try:
            assert len(
                new) >= 0, 'resized generated a negative length: %s' % new
        except (ValueError, AssertionError):
            raise ValueError
        return new

    def make_relative_to(self, other, make_relative=True):
        """returns a new coordinate with attributes adopted from other, and
        positioned relative to other."""

        if other.strand != self.strand:
            start = other.end - self.end
        elif make_relative:
            start = self.start - other.start
        else:
            start = self.start + other.start

        end = start + len(self)

        return self.__class__(other.genome, coord_name=other.coord_name,
                              start=start, end=end, strand=other.strand,
                              seq_region_id=other.seq_region_id)


class _CoordRecord(object):
    """store one record of the coord"""

    def __init__(self, attrib, rank, name=None, coord_system_id=None):
        self.coord_system_id = coord_system_id
        self.name = name
        self.rank = rank
        self.attr = attrib

    def __str__(self):
        return "coord_system_id = %d; name = %s; rank = %d; attr = %s "\
            % (self.coord_system_id, self.name, self.rank, self.attr)


class CoordSystemCache(object):
    """store coord_system table from core database.
    (only read default_version as stated in attrib column)
    There are two ways to get information about coordinate system:
    (1) use coord_type (e.g contig) which is at coord_system.c.name,
        and which are keys of _species_coord_systems[species]
    (2) use coord_system_id (e.g 17 refers to chromosome) which are also keys
        of _species_coord_systems[species]
    (3) to get which level of system is used for storing dna table, check
        'attrib' column of coord_system as default_version, sequence_level.
    """
    # Problem: multiple species (for compara) --> organized as {species: coordsystem}
    # TODO: simplify _species_coord_systems?
    # we place each species coord-system in _species_coord_systems, once, so
    # this attribute is a very _public_ attribute, and serves as a cache to
    # reduce unecessary lookups
    _species_coord_systems = {}
    # columns needed from coord_system table
    columns = ['coord_system_id', 'name', 'rank', 'attrib']
    # the attrib property has sequence_level, which means this the coordinate
    # system employed for sequence

    def _set_species_system(self, core_db, species):
        if species in self._species_coord_systems:
            return
        self._species_coord_systems[species] = {}
        coord_table = core_db.get_table('coord_system')
        records = sql.select([coord_table]).where(coord_table.c.attrib.like('default%')).\
            execute().fetchall()    # only select default version
        for record in records:
            attr = self._species_coord_systems[species]
            for key in ['coord_system_id', 'name']:
                key_val = record[key]
                vals = {}
                for column in self.columns:
                    val = record[column]
                    if isinstance(val, set):  # join items in set to one string
                        try:
                            val = ", ".join(sorted(val))
                        except TypeError:
                            pass
                    vals[column] = val
                attr[key_val] = _CoordRecord(**vals)

    def _get_seq_level_system(self, species):
        """returns the sequence level system for species"""
        sp_sys = self._species_coord_systems[species]
        for key, val in list(sp_sys.items()):
            if 'sequence_level' in val.attr:
                return val.name

        raise RuntimeError('no coord system for %s' % species)

    def __call__(self, coord_type=None, core_db=None, species=None,
                 seq_level=False):
        """coord_type can be coord_type or coord_system_id"""
        # TODO should only pass in core_db here, not that and species, or just
        # the genome - what if someone wants to compare different ensembl
        # releases? keying by species is then a bad idea! better to key by
        # id(object)
        # change identifier to coord_system, handle either string val or int
        # (see MySQL table) as is this shouldn't be a __call__, see line 168
        # for reason why we should have a method to set data: setSpeciesCoord
        # call then then just returns the coords for the named species
        species = _Species.get_species_name(species or core_db.db_name.species)
        self._set_species_system(core_db, species)
        if seq_level:
            result = self._get_seq_level_system(species)
        elif coord_type:
            result = self._species_coord_systems[species][coord_type]
        else:
            result = self._species_coord_systems[species]
        return result

CoordSystem = CoordSystemCache()


def _rank_checking(query_coord_type, target_coord_type, core_db, species):
    # assiting in constructingthe query language for assembly

    # in order to convert between coordinate systems, we need to establish the
    # ranking for coordinate types
    # rank defines the order of conversion between coord system 'levels'
    # chromosome has rank 1
    # super contig has rank 2
    # contig has rank 4
    # clone has rank 3
    # converting requires changing columns between 'asm' and 'cmp'
    # converting from clone -> contig, use 'asm' column
    # converting from contig -> clone, use 'cmp' column
    query_rank = CoordSystem(core_db=core_db, species=species,
                             coord_type=query_coord_type).rank
    target_rank = CoordSystem(core_db=core_db, species=species,
                              coord_type=target_coord_type).rank

    if query_rank < target_rank:
        query_prefix, target_prefix = 'asm', 'cmp'
    elif query_rank > target_rank:
        query_prefix, target_prefix = 'cmp', 'asm'
    else:
        query_prefix, target_prefix = '', ''
    return query_prefix, target_prefix


def _get_equivalent_coords(query_coord, assembly_row, query_prefix,
                           target_prefix, target_coord_type):
    # TODO better function name
    start = query_coord.ensembl_start
    end = query_coord.ensembl_end
    strand = query_coord.strand

    ori = assembly_row['ori']
    q_strand, t_strand = strand, strand * ori
    if 'seq_region' not in query_prefix:
        q_seq_region_id = assembly_row['%s_seq_region_id' % query_prefix]
        t_seq_region_id = assembly_row['%s_seq_region_id' % target_prefix]
    else:
        q_seq_region_id = assembly_row['_'.join([query_prefix, 'id'])]
        t_seq_region_id = assembly_row['_'.join([target_prefix, 'id'])]

    # d -- distance
    d_start = max(0, start - int(assembly_row['%s_start' % query_prefix]))
    d_end = max(0, int(assembly_row['%s_end' % query_prefix]) - end)
    # q -- query (to differ from the origin query block)
    q_start = int(assembly_row['%s_start' % query_prefix]) + d_start
    q_end = int(assembly_row['%s_end' % query_prefix]) - d_end

    if int(assembly_row['ori']) == -1:
        d_start, d_end = d_end, d_start
    # t -- target
    t_start = int(assembly_row['%s_start' % target_prefix]) + d_start
    t_end = int(assembly_row['%s_end' % target_prefix]) - d_end

    q_location = Coordinate(coord_name=query_coord.coord_name, start=q_start,
                            end=q_end, strand=q_strand,
                            CoordType=query_coord.CoordType,
                            seq_region_id=q_seq_region_id,
                            genome=query_coord.genome, ensembl_coord=True)
    t_location = Coordinate(coord_name=assembly_row['name'], start=t_start,
                            end=t_end, strand=t_strand, CoordType=target_coord_type,
                            seq_region_id=t_seq_region_id,
                            genome=query_coord.genome,
                            ensembl_coord=True)
    return [q_location, t_location]


def assembly_exception_coordinate(loc):
    """returns a coordinate conversion for one with an assembly exception"""
    genome = loc.genome
    assemb_except_table = genome.CoreDb.get_table('assembly_exception')
    seq_region_table = genome.CoreDb.get_table('seq_region')

    query = sql.select([assemb_except_table, seq_region_table.c.name],
                       sql.and_(
        assemb_except_table.c.seq_region_id ==
        loc.seq_region_id,
        assemb_except_table.c.exc_seq_region_id ==
        seq_region_table.c.seq_region_id))
    query = location_query(assemb_except_table,
                           loc.start, loc.end, query=query)
    record = asserted_one(query.execute().fetchall())
    s, conv_loc = _get_equivalent_coords(loc, record, "seq_region",
                                         "exc_seq_region", loc.CoordType)
    return conv_loc


def get_coord_conversion(query_location, target_coord_type, core_db, where=None):
    """returns the ???"""
    where = where or 'overlap'
    # TODO better function name
    species = core_db.db_name.species
    assert query_location.species == species
    assembly = core_db.get_table('assembly')
    seq_region = core_db.get_table('seq_region')
    target_coord_system_id = CoordSystem(target_coord_type, core_db=core_db,
                                         species=species).coord_system_id

    query_prefix, target_prefix = _rank_checking(query_location.CoordType,
                                                 target_coord_type, core_db, species)
    if query_prefix == target_prefix:
        return [[query_location, query_location]]
    # TODO: deal with query_prefix == target_prefix == '' --> could happen
    # when query features.
    query = sql.select([assembly, seq_region.c.name], sql.and_(assembly.c
                                                               ['%s_seq_region_id' %
                                                                   target_prefix] == seq_region.c.seq_region_id,
                                                               seq_region.c.coord_system_id == target_coord_system_id,
                                                               assembly.c['%s_seq_region_id' % query_prefix] ==
                                                               query_location.seq_region_id))
    query = location_query(assembly, query_location.ensembl_start,
                           query_location.ensembl_end,
                           start_col="%s_start" % query_prefix,
                           end_col="%s_end" % query_prefix, query=query,
                           where=where)
    assembly_rows = query.execute().fetchall()
    results = []
    for assembly_row in assembly_rows:
        results.append(_get_equivalent_coords(query_location, assembly_row,
                                              query_prefix, target_prefix, target_coord_type))
    return results
