import sys

import sqlalchemy as sql
from cogent3 import DNA
from cogent3.core.annotation import Feature
from cogent3.core.location import Map
from cogent3.util.table import Table

from .util import LazyRecord, asserted_one, DisplayString, \
    NoItemError
from .assembly import Coordinate, CoordSystem, \
    location_query, assembly_exception_coordinate
from .sequence import get_sequence
from .database import cached_attribs

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"

DEFAULT_PARENT_LENGTH = 2 ** 30

# some common string display formatting
_quoted = lambda x: DisplayString(x, with_quotes=True)
_limit_words = lambda x: DisplayString(x, with_quotes=True, num_words=3)


class _Region(LazyRecord):
    """a simple genomic region object"""
    type = None

    def __init__(self):
        super(_Region, self).__init__()
        self._attr_ensembl_table_map = None
        self._location_column_prefix = 'seq_region_'

    def __len__(self):
        return len(self.location)

    def __lt__(self, other):
        try:
            return self.location < other.location
        except AttributeError:
            raise TypeError("cannot compare %s to %s" % (self.__class__,
                                                         other.__class__))

    def __eq__(self, other):
        try:
            val = self.location == other.location
        except AttributeError:
            val = False
        return val

    def __ne__(self, other):
        try:
            val = self.location != other.location
        except AttributeError:
            val = True
        return val

    def _make_location(self):
        row = self._table_rows[self._attr_ensembl_table_map['location']]
        if row is None:
            return
        seq_region_id = row['%sid' % self._location_column_prefix]
        start = row['%sstart' % self._location_column_prefix]
        end = row['%send' % self._location_column_prefix]
        strand = row['%sstrand' % self._location_column_prefix]
        seq_region_table = self.db.get_table('seq_region')
        query = sql.select([seq_region_table.c.name],
                           seq_region_table.c.seq_region_id == seq_region_id)
        result = asserted_one(query.execute().fetchall())
        coord_name = result['name']

        coord = Coordinate(genome=self.genome, coord_name=coord_name,
                           start=start, end=end, strand=strand,
                           seq_region_id=seq_region_id,
                           ensembl_coord=True)
        self._cached['location'] = coord

    def _get_location_record(self):
        """makes the location data"""
        if not self._attr_ensembl_table_map['location'] in self._table_rows:
            # we use a bit of magic to figure out what method will be required
            # this magic assumes the method for obtaining a record from a table
            # are named _get_tablename_record
            dep_record_func = getattr(self, '_get_%s_record' %
                                      self._attr_ensembl_table_map['location'])
            dep_record_func()
        self._make_location()

    def _get_location(self):
        return self._get_cached_value('location', self._get_location_record)

    location = property(_get_location)

    def _get_sequence(self):
        if 'seq' not in self._cached:
            try:
                seq = get_sequence(self.location)
            except NoItemError:
                try:
                    alt_loc = assembly_exception_coordinate(self.location)
                    seq = get_sequence(alt_loc)
                except NoItemError:
                    seq = DNA.make_seq("N" * len(self))
            seq.name = str(self.location)
            self._cached['seq'] = seq
        return self._cached['seq']

    seq = property(_get_sequence)

    def _get_symbol(self):
        # override in subclasses
        return None

    symbol = property(_get_symbol)

    def get_features(self, feature_types, where_feature=None, limit=None):
        """queries the parent genome for feature types corresponding to this
        region
        where_feature: the returned region can either lie 'within' this region,
        'overlap' this region, or 'span' this region"""
        return self.genome.get_features(self.location,
                                       feature_types=feature_types,
                                       where_feature=where_feature,
                                       limit=limit)

    def _get_variants(self):
        """constructs the variants attribute"""
        if 'variants' not in self._cached:
            variants = self.genome.get_features(
                feature_types='variation', region=self)
            self._cached['variants'] = tuple(variants)

        return self._cached['variants']

    variants = property(_get_variants)

    def _feature_data(self, parent_map):
        symbol = self.symbol or getattr(self, 'stableid', '')
        assert not parent_map.reverse
        feat_map = parent_map[self.location.start:self.location.end]
        if feat_map.useful:
            if self.location.strand == -1:
                # this map is relative to + strand
                feat_map = feat_map.reversed()
            data = (self.type, str(symbol), feat_map)
        else:
            data = None
        return data

    def get_annotated_seq(self, feature_types=None, where_feature=None):
        regions = list(self.get_features(feature_types=feature_types,
                                        where_feature=where_feature))
        # seq_map is on the + strand, regardless the actual strand of sequence
        seq_map = Map(locations=[(self.location.start, self.location.end)],
                      parent_length=DEFAULT_PARENT_LENGTH)
        seq_map = seq_map.inverse()

        for region in regions:
            data = region._feature_data(seq_map)
            if data is None:
                continue
            # this will consider the strand information of actual sequence
            feature_map = [data[-1],
                           data[-1].nucleic_reversed()][self.location.strand == -1]
            self.seq.add_annotation(Feature, data[0], data[1], feature_map)

            if region.type == 'gene':  # TODO: SHOULD be much simplified
                sub_data = region._sub_feature_data(seq_map)
                for feature_type, feature_name, feature_map in sub_data:
                    if self.location.strand == -1:
                        # again, change feature map to -1 strand sequence if
                        # needed.
                        feature_map = feature_map.nucleic_reversed()
                    self.seq.add_annotation(Feature, feature_type,
                                           feature_name, feature_map)

        return self.seq


class GenericRegion(_Region):
    """a generic genomic region"""

    type = 'generic_region'

    def __init__(self, genome, db, location=None, coord_name=None, start=None,
                 end=None, strand=1, ensembl_coord=False):
        super(GenericRegion, self).__init__()
        self.genome = genome
        self.db = db

        if location is None and coord_name:
            self._get_seq_region_record(str(coord_name))
            if end is not None:
                assert self._table_rows['seq_region']['length'] > end, \
                    'Requested end[%s] too large' % end
            seq_region_id = self._table_rows['seq_region']['seq_region_id']
            location = Coordinate(genome=genome, coord_name=str(coord_name),
                                  start=start, end=end, strand=strand,
                                  seq_region_id=seq_region_id,
                                  ensembl_coord=ensembl_coord)

        if location is not None:
            self._cached['location'] = location

    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(species='%s'; coord_name='%s'; start=%s; end=%s;"\
               " length=%s; strand='%s')" % (my_type,
                                             self.genome.species,
                                             self.location.coord_name, self.location.start,
                                             self.location.end, len(self), '-+'[self.location.strand > 0])

    def _get_seq_region_record(self, coord_name):
        # override the _Region class method, since, we take the provided start
        # etc .. attributes
        # coord_name comes from seq_region_table.c.name
        # matched, by coord_system_id, to default coord system
        seq_region_table = self.genome.db.get_table('seq_region')
        coord_systems = CoordSystem(core_db=self.genome.CoreDb)
        coord_system_ids = [k for k in coord_systems if not isinstance(k, str)]
        record = sql.select([seq_region_table],
                            sql.and_(seq_region_table.c.name == coord_name,
                                     seq_region_table.c.coord_system_id.in_(coord_system_ids)))
        record = asserted_one(record.execute().fetchall())
        self._table_rows['seq_region'] = record


class _StableRegion(GenericRegion):
    """region with a stable_id"""

    _member_types = None

    def __init__(self, genome, db, **kwargs):
        super(_StableRegion, self).__init__(genome, db, **kwargs)

    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(%s; %s)' % (my_type, self.genome.species, self.stableid)

    def _get_record_for_stable_id(self):
        # subclasses need to provide a function for loading the correct
        # record for obtaining a stable_id
        table_name = self._attr_ensembl_table_map['stableid']

        if self.genome.general_release >= 65:
            func_name = '_get_%s_record' % (table_name + '_stable_id')
        else:
            func_name = '_get_%s_record' % table_name

        func = getattr(self, func_name)
        func()
        attr_column_map = [('stableid', 'stable_id', _quoted)]

        self._populate_cache_from_record(attr_column_map, table_name)

    def _get_stable_id(self):
        return self._get_cached_value('stableid',
                                      self._get_record_for_stable_id)

    stableid = property(_get_stable_id)

    def get_member(self, stableid, member_types=None):
        """returns the associated member with matching stableid or None if not
        found.

        Arguments:
            - member_types: the property to be searched, depends on self.type.
              transcripts for genes, exons/translated_exons for transcripts."""

        member_types = member_types or self._member_types
        if type(member_types) == str:
            member_types = [member_types]

        for member_type in member_types:
            member = getattr(self, member_type, None)
            if member is None:
                raise AttributeError(
                    "%s doesn't have property %s" % (self.type, member_type))
            for element in member:
                if element.stableid == stableid:
                    return element
        return None


class Gene(_StableRegion):
    """a gene region"""
    type = 'gene'
    _member_types = ['transcripts']

    def __init__(self, genome, db, stableid=None, symbol=None, location=None, data=None):
        """constructed by a genome instance"""
        super(Gene, self).__init__(genome, db, location=location)

        self._attr_ensembl_table_map = dict(stableid=['gene_stable_id',
                                                      'gene'][genome.general_release >= 65],
                                            symbol='xref',
                                            description='gene', biotype='gene', location='gene',
                                            canonical_transcript='gene',
                                            transcripts='transcript',
                                            exons='transcript')

        if data is None:
            args = [dict(stableid=stableid), dict(
                symbol=symbol)][stableid is None]
            assert args
            data = asserted_one(
                list(self.genome._get_gene_query(db, **args).execute()))

        for name, func in \
            [('stableid', self._get_gene_stable_id_record),
             ('biotype', self._get_gene_record),
             ('description', self._get_gene_record),
             ('symbol', self._get_xref_record),
             ('location', self._get_gene_record)]:
            # For EST
            if name == 'symbol' and 'display_label' not in list(data.keys()):
                continue
            self._table_rows[self._attr_ensembl_table_map[name]] = data
            func()  # this populates the attributes

    def __str__(self):
        my_type = self.__class__.__name__
        vals = ['%s=%r' % (key, val) for key, val in list(self._cached.items())
                if val is not None]
        vals.sort()
        vals.insert(0, "species='%s'" % self.genome.species)
        return '%s(%s)' % (my_type, '; '.join(vals))

    def __repr__(self):
        my_type = self.__class__.__name__
        vals = ['%s=%r' % (key, val) for key, val in list(self._cached.items())
                if val is not None]
        vals.sort()
        vals.insert(0, 'species=%r' % self.genome.species)
        return '%s(%s)' % (my_type, '; '.join(vals))

    def _get_gene_record(self):
        """adds the gene data to self._table_rows"""
        attr_column_map = [('biotype', 'biotype', _quoted),
                           ('status', 'status', _quoted),
                           ('description', 'description', _limit_words)]
        # we set all the attributes that derive from this
        self._populate_cache_from_record(attr_column_map, 'gene')
        return

    def _get_gene_stable_id_record(self):
        """adds the gene_stable_id data to self._table_rows"""
        attr_column_map = [('stableid', 'stable_id', _quoted)]
        self._populate_cache_from_record(attr_column_map,
                                         self._attr_ensembl_table_map['stableid'])
        return

    def _get_xref_record(self):
        attr_column_map = [('symbol', 'display_label', _quoted)]
        self._populate_cache_from_record(attr_column_map, 'xref')
        return

    def _get_biotype(self):
        return self._get_cached_value('biotype', self._get_gene_record)

    biotype = property(_get_biotype)

    def _get_symbol(self):
        if 'xref' in self._table_rows:
            return self._get_cached_value('symbol', self._get_xref_record)
        self._set_null_values(['symbol'])
        return self._cached['symbol']

    symbol = property(_get_symbol)

    def _get_description(self):
        return self._get_cached_value('description', self._get_gene_record)

    description = property(_get_description)

    def _get_status(self):
        return self._get_cached_value('status', self._get_gene_record)

    status = property(_get_status)

    def _make_canonical_transcript(self):
        if 'gene' not in self._table_rows:
            self._get_gene_record()
        canonical_id = self._table_rows['gene']['canonical_transcript_id']
        transcript_table = self.db.get_table('transcript')
        query = sql.select([transcript_table],
                           transcript_table.c.transcript_id == canonical_id)
        records = query.execute().fetchall()
        assert len(records) == 1,\
            "wrong number of records from canonical_transcript"
        record = records[0]
        transcript = Transcript(self.genome, self.db, canonical_id,
                                data=record)
        self._cached['canonical_transcript'] = transcript

    def _get_canonical_transcript(self):
        return self._get_cached_value('canonical_transcript',
                                      self._make_canonical_transcript)

    canonical_transcript = property(_get_canonical_transcript)

    def _make_transcripts(self):
        if 'gene' not in self._table_rows:
            self._get_gene_record()
        gene_id = self._table_rows['gene']['gene_id']
        transcript_table = self.db.get_table('transcript')
        query = sql.select([transcript_table],
                           transcript_table.c.gene_id == gene_id)
        records = query.execute().fetchall()
        if not records:
            self._set_null_values(['transcripts'], 'transcript')
            return
        transcripts = []
        for record in records:
            transcript_id = record['transcript_id']
            transcripts.append(Transcript(self.genome, self.db, transcript_id,
                                          data=record))
        self._cached['transcripts'] = tuple(transcripts)

    def _get_transcripts(self):
        return self._get_cached_value('transcripts', self._make_transcripts)

    transcripts = property(_get_transcripts)

    def _sub_feature_data(self, parent_map):
        """returns data for making a cogent Feature. These can be
        automatically applied to the seq by the get_annotated_seq method.
        Returns None if self lies outside parent's span.
        """
        features = []
        for transcript in self.transcripts:
            transcript_data = transcript._feature_data(parent_map)
            if transcript_data:
                features.append(transcript_data)
                data = transcript._sub_feature_data(parent_map)
                features.extend(data)
        return features

    def get_cds_lengths(self):
        """returns the cds lengths from transcripts with the same biotype.
        returns None if no transcripts."""
        if self.transcripts is self.NULL_VALUE:
            return None
        l = [ts.get_cds_length() for ts in self.transcripts
             if ts.biotype == self.biotype]
        return l

    def get_longest_cds_transcript(self):
        """returns the Transcript with the longest cds and the same biotype"""
        result = sorted([(ts.get_cds_length(), ts) for ts in self.transcripts
                         if ts.biotype == self.biotype])

        if result:  # last result is longest
            result = result[-1][1]

        return result


class Transcript(_StableRegion):
    type = 'transcript'
    _member_types = ['exons', 'translated_exons']

    def __init__(self, genome, db, transcript_id, data, location=None):
        """created by Gene"""
        super(Transcript, self).__init__(genome, db, location=location)

        self._attr_ensembl_table_map = dict(stableid=['transcript_stable_id',
                                                      'transcript'][genome.general_release >= 65],
                                            location='transcript',
                                            status='transcript',
                                            translated_exons='translation')

        self._am_prot_coding = None
        self.transcript_id = transcript_id
        self._table_rows['transcript'] = data
        self.gene_id = self._table_rows['transcript']['gene_id']
        self._set_transcript_record()

    def _set_transcript_record(self):
        attr_column_map = [('biotype', 'biotype', _quoted),
                           ('status', 'status', _quoted)]
        self._populate_cache_from_record(attr_column_map, 'transcript')
        self._am_prot_coding = self._cached[
            'biotype'].lower() == 'protein_coding'

    def _get_status(self):
        return self._cached['status']

    status = property(_get_status)

    def _get_biotype(self):
        return self._cached['biotype']

    biotype = property(_get_biotype)

    def _get_gene(self):
        gene_id = self.gene_id
        gene_table = self.db.get_table('gene')
        query = sql.select([gene_table.c.stable_id],
                           gene_table.c.gene_id == gene_id)
        record = asserted_one(query.execute())
        gene = self.genome.get_gene_by_stableid(record[0])
        return gene

    gene = property(_get_gene)

    def _get_transcript_stable_id_record(self):
        table_name = self._attr_ensembl_table_map['stableid']
        if table_name in self._table_rows:
            return
        transcript_id = self.transcript_id
        table = self.db.get_table(table_name)
        query = sql.select([table], table.c.transcript_id == transcript_id)
        record = asserted_one(query.execute())
        self._table_rows[table_name] = record

    def _get_exon_transcript_records(self):
        transcript_id = self.transcript_id
        exon_transcript_table = self.db.get_table('exon_transcript')
        query = sql.select([exon_transcript_table],
                           exon_transcript_table.c.transcript_id == transcript_id)
        records = query.execute()
        exons = []
        for record in records:
            exons.append(Exon(self.genome, self.db, record['exon_id'],
                              record['rank']))
        exons.sort()
        self._cached['exons'] = tuple(exons)

    def _get_exons(self):
        return self._get_cached_value('exons',
                                      self._get_exon_transcript_records)

    exons = property(_get_exons)

    def _get_intron_transcript_records(self):
        if len(self.exons) < 2:
            self._set_null_values(["introns"])
            return

        exon_positions = [(exon.location.start, exon.location.end)
                          for exon in self.exons]
        exon_positions.sort()
        end = exon_positions[-1][-1]
        exon_map = Map(locations=exon_positions, parent_length=end)
        intron_map = exon_map.shadow()

        intron_positions = [(span.start, span.end)
                            for span in intron_map.spans if span.start != 0]

        chrom = self.location.coord_name
        strand = self.location.strand
        introns = []
        rank = 1
        if strand == -1:
            intron_positions.reverse()
        for s, e in intron_positions:
            coord = self.genome.make_location(coord_name=chrom, start=s, end=e,
                                             strand=strand, ensembl_coord=False)
            introns.append(Intron(self.genome, self.db, rank, self.stableid,
                                  coord))
            rank += 1

        self._cached['introns'] = tuple(introns)

    def _get_introns(self):
        return self._get_cached_value('introns',
                                      self._get_intron_transcript_records)

    introns = property(_get_introns)

    def _get_translation_record(self):
        transcript_id = self.transcript_id
        translation_table = self.db.get_table('translation')
        query = sql.select([translation_table],
                           translation_table.c.transcript_id == transcript_id)
        try:
            record = asserted_one(query.execute())
        except NoItemError:
            self._set_null_values(['translated_exons'], 'translation')
            return

        self._table_rows['translation'] = record

    def _get_transcript(self):
        self._get_translation_record()
        record = self._table_rows['translation']
        if record == self.NULL_VALUE:
            return

        start_exon_id = record['start_exon_id']
        end_exon_id = record['end_exon_id']
        # because this will be used to shift a coord. Note: these are relative
        # to the exon start but ignore strand, so we have to decide whether
        # the coord shifts need to be flipped
        seq_start = record['seq_start'] - 1
        seq_end = record['seq_end']
        flip_coords = self.exons[0].location.strand == -1

        start_index = None
        end_index = None
        for index, exon in enumerate(self.exons):
            if exon.exon_id == start_exon_id:
                start_index = index
            if exon.exon_id == end_exon_id:
                end_index = index
        assert None not in (start_index, end_index), \
            'Error in matching transcript and exons'

        start_exon = self.exons[start_index]

        if start_index == end_index:
            shift_start = [seq_start, len(start_exon) - seq_end][flip_coords]
            shift_end = [seq_end - len(start_exon), -
                         1 * seq_start][flip_coords]
        else:
            shift_start = [seq_start, 0][flip_coords]
            shift_end = [0, -1 * seq_start][flip_coords]

        coord = start_exon.location.resized(shift_start, shift_end)

        DEBUG = False
        if DEBUG:
            out = ['\nseq_start=%d; seq_end=%d' % (seq_start, seq_end),
                   'shift_start=%d; shift_end=%d' % (shift_start, shift_end),
                   'len=%s' % len(coord)]
            sys.stderr.write('\n'.join(map(str, out)) + '\n')

        new_start_exon = Exon(self.genome, self.db, start_exon.exon_id,
                              start_exon.rank, location=coord)
        translated_exons = (new_start_exon,) +\
            self.exons[start_index + 1:end_index]
        if start_index != end_index:
            end_exon = self.exons[end_index]
            shift_start = [0, len(end_exon) - seq_end][flip_coords]
            shift_end = [seq_end - len(end_exon), 0][flip_coords]
            coord = end_exon.location.resized(shift_start, shift_end)
            new_end_exon = Exon(self.genome, self.db, end_exon.exon_id,
                                end_exon.rank, location=coord)
            translated_exons += (new_end_exon,)
        self._cached['translated_exons'] = translated_exons

    def _get_translated_exons(self):
        return self._get_cached_value('translated_exons', self._get_transcript)

    translated_exons = property(_get_translated_exons)

    def _calculate_Utr_exons(self):
        # TODO clean up this code
        exons = self.exons
        translated_exons = self.translated_exons
        num_exons = len(self.exons)
        if not translated_exons:
            self._set_null_values(["untranslated_exons_5", "untranslated_exons_3"])
            return
        untranslated_5exons, untranslated_3exons = [], []
        start_exon, end_exon = translated_exons[0], translated_exons[-1]
        flip_coords = start_exon.location.strand == -1

        for exon in exons[0:start_exon.rank]:   # get 5'UTR
            coord = exon.location.copy()
            if exon.stableid == start_exon.stableid:
                coord.start = [coord.start,
                               start_exon.location.end][flip_coords]
                coord.end = [start_exon.location.start, coord.end][flip_coords]
            if len(coord) != 0:
                untranslated_5exons.append(Exon(self.genome, self.db,
                                                exon.exon_id, exon.rank, location=coord))
        for exon in exons[end_exon.rank - 1: num_exons]:  # get 3'UTR
            coord = exon.location.copy()
            if exon.stableid == end_exon.stableid:
                coord.start = [end_exon.location.end, coord.start][flip_coords]
                coord.end = [coord.end, end_exon.location.start][flip_coords]
            if len(coord) != 0:
                untranslated_3exons.append(Exon(self.genome, self.db,
                                                exon.exon_id, exon.rank, location=coord))

        self._cached["untranslated_exons_5"] = tuple(untranslated_5exons)
        self._cached["untranslated_exons_3"] = tuple(untranslated_3exons)

    def _get_5prime_untranslated_exons(self):
        return self._get_cached_value("untranslated_exons_5",
                                      self._calculate_Utr_exons)

    untranslated_exons_5 = property(_get_5prime_untranslated_exons)

    def _get_3prime_untranslated_exons(self):
        return self._get_cached_value("untranslated_exons_3",
                                      self._calculate_Utr_exons)

    untranslated_exons_3 = property(_get_3prime_untranslated_exons)

    def _make_utr_seq(self):
        if self.untranslated_exons_5 is None and self.untranslated_exons_3 is None:
            self._cached["utr5"] = self.NULL_VALUE
            self._cached["utr3"] = self.NULL_VALUE
            return
        Utr5_seq, Utr3_seq = DNA.make_seq(""), DNA.make_seq("")
        for exon in self.untranslated_exons_5:
            Utr5_seq += exon.seq
        for exon in self.untranslated_exons_3:
            Utr3_seq += exon.seq
        self._cached["utr5"] = Utr5_seq
        self._cached["utr3"] = Utr3_seq

    def _get_utr5_seq(self):
        return self._get_cached_value("utr5", self._make_utr_seq)

    utr5 = property(_get_utr5_seq)

    def _get_utr3_seq(self):
        return self._get_cached_value("utr3", self._make_utr_seq)

    utr3 = property(_get_utr3_seq)

    def _make_cds_seq(self):
        if self.exons is self.NULL_VALUE:
            self._cached['cds'] = self.NULL_VALUE
            return

        exons = [self.exons, self.translated_exons][self._am_prot_coding]
        full_seq = None
        for exon in exons:
            if full_seq is None:
                full_seq = exon.seq
                continue
            full_seq += exon.seq

        # check first exon phase_start is 0 and last exon phase_end
        if exons[0].phase_start > 0:
            fill = DNA.make_seq(
                'N' * exons[0].phase_start, name=full_seq.name)
            full_seq = fill + full_seq

        if exons[-1].phase_end > 0:
            fill = DNA.make_seq(
                'N' * exons[-1].phase_end, name=full_seq.name)
            full_seq += fill

        self._cached['cds'] = full_seq

    def _get_cds(self):
        return self._get_cached_value('cds', self._make_cds_seq)

    cds = property(_get_cds)

    def get_cds_length(self):
        """returns the length of the cds. If this property is not available,
        returns None."""
        if self.cds is self.NULL_VALUE:
            return None
        exons = [self.exons, self.translated_exons][self._am_prot_coding]
        return sum(map(len, exons))

    def _make_protein_seq(self):
        if not self._am_prot_coding or self.cds is self.NULL_VALUE:
            self._cached['protein_seq'] = self.NULL_VALUE
            return

        DEBUG = False
        # enforce multiple of 3
        cds = self.cds
        length = len(cds)
        cds = cds[: length - (length % 3)]
        try:
            cds = cds.trim_stop_codon()
        except AssertionError:
            if not DEBUG:
                raise
            out = ['\n****\nFAILED=%s' % self.stableid]
            for exon in self.translated_exons:
                out += ['TranslatedExon[rank=%d]\n' % exon.rank, exon,
                        exon.location,
                        '%s ... %s' % (exon.seq[:20], exon.seq[-20:])]
                sys.stderr.write('\n'.join(map(str, out)) + '\n')
            raise

        self._cached['protein_seq'] = cds.get_translation()

    def _get_protein_seq(self):
        return self._get_cached_value('protein_seq', self._make_protein_seq)

    protein_seq = property(_get_protein_seq)

    def _get_exon_feature_data(self, parent_map):
        """returns the exon feature data"""
        features = []
        if self.exons is self.NULL_VALUE:
            return features
        for exon in self.exons:
            feature_data = exon._feature_data(parent_map)
            if feature_data is None:
                continue
            features.append(feature_data)
        return features

    def _get_intron_feature_data(self, parent_map):
        """return the intron feature data"""
        features = []
        if self.introns is self.NULL_VALUE:
            return features
        for intron in self.introns:
            feature_data = intron._feature_data(parent_map)
            if feature_data is None:
                continue
            features.append(feature_data)
        return features

    def _get_translated_exon_feature_data(self, parent_map):
        """returns featureD data for translated exons"""
        features = []
        if self.translated_exons is self.NULL_VALUE:
            return features
        cds_spans = []
        for exon in self.translated_exons:
            feature_data = exon._feature_data(parent_map)
            if feature_data is None:
                continue
            cds_spans.extend(feature_data[-1].spans)
        if cds_spans:
            # TODO: check strand
            cds_map = Map(spans=cds_spans,
                          parent_length=parent_map.parent_length)
            features.append(('CDS', str(self.stableid), cds_map))
        return features

    def _get_Utr_feature_data(self, parent_map):
        # TODO: Simplify this part
        features = []
        utr5_spans, utr3_spans = [], []
        for exon in self.untranslated_exons_5:
            feature_data = exon._feature_data(parent_map)
            if feature_data is None:
                continue
            utr5_spans.extend(feature_data[-1].spans)
        for exon in self.untranslated_exons_3:
            feature_data = exon._feature_data(parent_map)
            if feature_data is None:
                continue
            utr3_spans.extend(feature_data[-1].spans)
        if utr5_spans:
            utr5_map = Map(spans=utr5_spans,
                           parent_length=parent_map.parent_length)
            features.append(("5'UTR", str(self.stableid), utr5_map))
        if utr3_spans:
            utr3_map = Map(spans=utr3_spans,
                           parent_length=parent_map.parent_length)
            features.append(("3'UTR", str(self.stableid), utr3_map))
        return features

    def _sub_feature_data(self, parent_map):
        """returns data for making a cogent Feature. This can be automatically
        applied to the seq by the get_annotated_seq method. Returns None if
        self lies outside parent's span.
        """
        features = self._get_exon_feature_data(parent_map)
        features += self._get_intron_feature_data(parent_map)
        features += self._get_translated_exon_feature_data(parent_map)
        if self.translated_exons:
            features += self._get_Utr_feature_data(parent_map)
        return features


class Exon(_StableRegion):
    type = 'exon'

    def __init__(self, genome, db, exon_id, rank, location=None):
        """created by a Gene"""
        _StableRegion.__init__(self, genome, db, location=location)

        self._attr_ensembl_table_map = dict(stableid=['exon_stable_id',
                                                      'exon'][genome.general_release >= 65],
                                            location='exon')

        self.exon_id = exon_id
        self.rank = rank

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(stableid=%s, rank=%s)' % (my_type, self.stableid, self.rank)

    def __lt__(self, other):
        return self.rank < other.rank

    def __eq__(self, other):
        return self.rank == other.rank

    def __ne__(self, other):
        return self.rank != other.rank

    def _get_exon_stable_id_record(self):
        if self.genome.general_release >= 65:
            # release >= 65, data is just in the exon table
            self._get_exon_record()
            return

        table_name = self._attr_ensembl_table_map['stableid']
        exon_stable_id_table = self.db.get_table(table_name)
        query = sql.select([exon_stable_id_table.c.stable_id],
                           exon_stable_id_table.c.exon_id == self.exon_id)
        records = query.execute()
        record = asserted_one(records.fetchall())
        self._table_rows[table_name] = record

    def _get_exon_record(self):
        # this will be called by _Region parent class to make the location
        exon_table = self.db.get_table('exon')
        query = sql.select([exon_table], exon_table.c.exon_id == self.exon_id)
        records = query.execute()
        record = asserted_one(records.fetchall())
        self._table_rows['exon'] = record

    def _make_symbol(self):
        self._cached['symbol'] = '%s-%s' % (self.stableid, self.rank)

    def _get_symbol(self):
        return self._get_cached_value('symbol', self._make_symbol)

    symbol = property(_get_symbol)

    def _make_phase(self):
        """creates the exon phase attributes"""
        if 'exon' not in self._table_rows:
            self._get_exon_record()

        exon = self._table_rows['exon']
        self._cached['phase_start'] = exon['phase']
        self._cached['phase_end'] = exon['end_phase']

    @property
    def phase_start(self):
        """reading frame start for this exon"""
        return self._get_cached_value('phase_start', self._make_phase)

    @property
    def phase_end(self):
        """reading frame end for this exon"""
        return self._get_cached_value('phase_end', self._make_phase)


class Intron(GenericRegion):
    type = 'intron'

    def __init__(self, genome, db, rank, transcript_stable_id, location=None):
        GenericRegion.__init__(self, genome, db, location=location)
        self.TranscriptStableId = transcript_stable_id
        self.rank = rank

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        my_type = self.__class__.__name__
        return '%s(TranscriptId=%s, rank=%s)' % (my_type,
                                                 self.TranscriptStableId, self.rank)

    def _make_symbol(self):
        self._cached['symbol'] = '%s-%s' % (self.TranscriptStableId, self.rank)

    def _get_symbol(self):
        return self._get_cached_value('symbol', self._make_symbol)

    symbol = property(_get_symbol)


class Est(Gene):
    """an EST region"""
    type = 'est'


def _set_to_string(val):
    if type(val) in (str, type(None)):
        return val
    val = list(val)
    while len(val) == 1 and type(val) in (tuple, list):
        val = val[0]
    return val


class Variation(_Region):
    """genomic variation"""
    type = 'variation'

    def __init__(self, genome, db=None, effect=None, symbol=None, data=None):
        self.genome = genome

        get_table = genome.VarDb.get_table
        self.variation_feature_table = get_table('variation_feature')
        self.transcript_variation_table = get_table('transcript_variation')
        self.allele_table = get_table('allele')
        self.variation_table = get_table('variation')
        self.attr_table = get_table("attrib")
        self.attrib_type_tab = get_table("attrib_type")
        try:
            self.allele_code_table = get_table('allele_code')
        except sql.exceptions.ProgrammingError:
            self.allele_code_table = None

        super(Variation, self).__init__()

        if genome.general_release < 70:
            self._get_flanking_seq_data = self._get_flanking_seq_data_lt_70
        else:
            self._get_flanking_seq_data = self._get_flanking_seq_data_ge_70

        self._attr_ensembl_table_map = dict(effect='variation_feature',
                                            symbol='variation_feature',
                                            validation='variation_feature',
                                            map_weight='variation_feature',
                                            flanking_seq='flanking_sequence',
                                            peptide_alleles='transcript_variation',
                                            translation_location='transcript_variation',
                                            location='variation_feature',
                                            allele_freqs='allele',
                                            ancestral='variation')

        assert data is not None, 'Variation record created in an unusual way'
        for name, value, func in \
            [('effect', effect, self._get_variation_table_record),
             ('symbol', symbol, self._get_variation_table_record)]:
            if value is not None:
                self._table_rows[self._attr_ensembl_table_map[name]] = data
                if func is not None:
                    func()  # this populates the attributes

        self.db = db or self.genome.CoreDb

    def __len__(self):
        """return the length of the longest allelic variant"""
        return max(list(map(len, self._split_alleles())))

    def __str__(self):
        my_type = self.__class__.__name__

        return "%s(symbol=%r; effect=%r; alleles=%r)" % \
               (my_type, self.symbol, self.effect, self.alleles)

    def _get_variation_table_record(self):
        # this is actually the variation_feature table
        consequence_type = 'consequence_type'
        if self.genome.general_release > 67:
            consequence_type += 's'  # change to plural column name

        attr_name_map = [('effect', consequence_type, _set_to_string),
                         ('alleles', 'allele_string', _quoted),
                         ('symbol', 'variation_name', _quoted),
                         ('map_weight', 'map_weight', int),
                         ('somatic', 'somatic', bool)]
        
        if self.genome.general_release < 83:
            attr_name_map.append(
                         ('validation', 'validation_status', _set_to_string))
        else:
            attr_name_map.append(('validation', 'evidence_attribs', lambda x: x))            
        
        self._populate_cache_from_record(attr_name_map, 'variation_feature')
        # TODO handle obtaining the variation_feature if we were created in
        # any way other than through the symbol or effect

    def _get_ancestral_data(self):
        # actually the variation table
        def str_or_none(val):
            if val is not None:
                val = _quoted(val)
            return val

        variation_id = self._table_rows['variation_feature']['variation_id']
        var_table = self.variation_table
        query = sql.select([var_table],
                           var_table.c.variation_id == variation_id)
        record = asserted_one(query.execute())
        self._table_rows['variation'] = record
        attr_name_map = [('ancestral', 'ancestral_allele', str_or_none)]
        self._populate_cache_from_record(attr_name_map, 'variation')

    def _get_seq_region_record(self, seq_region_id):
        # should this be on a parent class? or a generic function in assembly?
        seq_region_table = self.db.get_table('seq_region')
        query = sql.select([seq_region_table],
                           seq_region_table.c.seq_region_id == seq_region_id)
        record = asserted_one(query.execute())
        return record

    def _get_flanking_seq_data_ge_70(self):
        """return the flanking sequence data if release >= 70"""
        # variation_feature.alignment_quality == 1, means flanks match reference
        # genome, 0 means they don't
        aligned_ref = self._table_rows[
            'variation_feature']['alignment_quality'] == 1
        if not aligned_ref:
            self._cached['flanking_seq'] = self.NULL_VALUE
            return

        seqs = dict(up=self.NULL_VALUE, down=self.NULL_VALUE)
        for name, seq in list(seqs.items()):
            resized = [(-301, -1), (1, 301)][name == 'down']
            if self.location.strand == -1:
                resized = [(1, 301), (-301, -1)][name == 'down']
            flank = self.location.resized(*resized)
            flanking = self.genome.get_region(region=flank)
            seq = flanking.seq
            seqs[name] = seq

        self._cached[('flanking_seq')] = (seqs['up'][-300:], seqs['down'][:300])

    def _get_flanking_seq_data_lt_70(self):
        # maps to flanking_sequence through variation_feature_id
        # if this fails, we grab from genomic sequence
        variation_id = self._table_rows['variation_feature']['variation_id']
        flanking_seq_table = self.flanking_sequence_table
        query = sql.select([flanking_seq_table],
                           flanking_seq_table.c.variation_id == variation_id)
        record = asserted_one(query.execute())
        self._table_rows['flanking_sequence'] = record
        up_seq = record['up_seq']
        down_seq = record['down_seq']
        # the following two lines are because -- wait for it -- someone has
        # entered the string 'NULL' instead of NULL in the MySQL tables!!!
        up_seq = [up_seq, None][up_seq == 'NULL']
        down_seq = [down_seq, None][down_seq == 'NULL']
        seqs = dict(up=up_seq, down=down_seq)
        for name, seq in list(seqs.items()):
            if seq is not None:
                seq = DNA.make_seq(seq)
            else:
                resized = [(-301, -1), (1, 301)][name == 'down']
                if self.location.strand == -1:
                    resized = [(1, 301), (-301, -1)][name == 'down']
                flank = self.location.resized(*resized)
                flanking = self.genome.get_region(region=flank)
                seq = flanking.seq
            seqs[name] = seq

        self._cached[('flanking_seq')] = (seqs['up'][-300:], seqs['down'][:300])

    def _get_flanking_seq(self):
        return self._get_cached_value('flanking_seq',
                                      self._get_flanking_seq_data)

    flanking_seq = property(_get_flanking_seq)

    def _get_effect(self):
        return self._get_cached_value('effect',
                                      self._get_variation_table_record)

    effect = property(_get_effect)

    def _get_somatic(self):
        return self._get_cached_value('somatic',
                                      self._get_variation_table_record)

    somatic = property(_get_somatic)

    def _get_alleles(self):
        return self._get_cached_value('alleles',
                                      self._get_variation_table_record)

    alleles = property(_get_alleles)

    def _get_ancestral(self):
        return self._get_cached_value('ancestral', self._get_ancestral_data)

    ancestral = property(_get_ancestral)

    def _get_allele_table_record(self):
        variation_id = self._table_rows['variation_feature']['variation_id']
        allele_table = self.allele_table
        query = sql.select([allele_table],
                           allele_table.c.variation_id == variation_id)
        records = [r for r in query.execute()]

        if len(records) == 0:
            self._cached[('allele_freqs')] = self.NULL_VALUE
            return

        # property change from >= 65, allele ids need to be looked up in
        # the allele_code table
        allele_code = self.allele_code_table

        self._table_rows['allele_table'] = records
        data = []

        if self.genome.general_release > 71:
            sample_id = 'population_id'
        else:
            sample_id = 'sample_id'

        for rec in records:
            if not rec[sample_id]:
                continue

            if allele_code is None:
                allele = rec['allele']
            else:
                allele_query = sql.select([allele_code.c.allele],
                                          allele_code.c.allele_code_id == rec['allele_code_id'])
                allele = list(allele_query.execute())[0][0]

            data.append((allele, rec['frequency'], rec[sample_id]))

        if not data:
            self._cached[('allele_freqs')] = self.NULL_VALUE
            return

        table = Table(header=['allele', 'freq', sample_id], rows=data)
        self._cached[('allele_freqs')] = table.sorted([sample_id, 'allele'])

    def _get_allele_freqs(self):
        return self._get_cached_value('allele_freqs',
                                      self._get_allele_table_record)

    allele_freqs = property(_get_allele_freqs)

    def _get_symbol(self):
        return self._get_cached_value('symbol',
                                      self._get_variation_table_record)

    symbol = property(_get_symbol)

    @property
    def validation(self):
        result = self._get_cached_value('validation',
                                        self._get_variation_table_record)        
        if self.genome.general_release < 83:
            # we need to access cached descriptions
            return result
        
        mapping = cached_attribs[(self.genome, "validation")]
        if mapping is None:
            attr = self.attr_table
            attr_ty = self.attrib_type_tab
            query = sql.select([attr.c.attrib_id, attr.c.value],
                               sql.and_(attr.c.attrib_type_id==attr_ty.c.attrib_type_id,
                                        attr_ty.c.name == 'Variant evidence status'))
            mapping = dict((str(i), v) for i,v in query.execute().fetchall())
            cached_attribs.add_to_cache(self.genome, "validation", mapping)
        
        out = [mapping[k] for k in sorted(result)]
        result = set(out)
        
        return result
            

    def _get_map_weight(self):
        return self._get_cached_value('map_weight',
                                      self._get_variation_table_record)

    map_weight = property(_get_map_weight)

    def _get_transcript_record(self):
        if 'variation_feature' not in self._table_rows:
            raise NotImplementedError

        try:
            effects = [self.effect.lower()]
        except AttributeError:
            effects = [v.lower() for v in self.effect]

        effects = set(effects)
        # TODO swap what we use for nysn by Ensembl version, thanks Ensembl!
        nsyn = set(('coding_sequence_variant', 'missense_variant'))
        if not effects & nsyn:
            self._cached['peptide_alleles'] = self.NULL_VALUE
            self._cached['translation_location'] = self.NULL_VALUE
            return

        table_name = self._attr_ensembl_table_map['peptide_alleles']
        loc = lambda x: int(x) - 1

        # column name changed between releases, so we check to see which
        # one is being used for this instance and set the column strings

        # TODO can we modify the table on loading? This would give better
        # performance.

        if self.genome.VarDb.table_has_column(table_name, 'pep_allele_string'):
            pep_allele_string = 'pep_allele_string'
            consequence_type = 'consequence_types'
        else:
            pep_allele_string = 'peptide_allele_string'
            consequence_type = 'consequence_type'

        attr_column_map = [
            ('peptide_alleles', pep_allele_string, _quoted),
            ('translation_location', 'translation_start', loc)]

        if table_name in self._table_rows:
            self._populate_cache_from_record(attr_column_map, table_name)
            return

        var_feature_record = self._table_rows['variation_feature']
        var_feature_id = var_feature_record['variation_feature_id']
        table = self.genome.VarDb.get_table(table_name)
        self_effect = set([self.effect, [self.effect]]
                          [type(self.effect) == str])
        query = sql.select([table.c.variation_feature_id,
                            table.columns[pep_allele_string],
                            table.c.translation_start,
                            table.columns[consequence_type]],
                           sql.and_(table.c.variation_feature_id == var_feature_id,
                                    table.columns[pep_allele_string] is not None))
        records = query.execute().fetchall()
        pep_alleles = []
        translation_location = []
        for record in records:
            if not record[consequence_type] & self_effect:
                continue
            
            allele = record[pep_allele_string]
            if not allele:
                continue
            
            pep_alleles += [allele]
            translation_location += [record['translation_start']]

        if not pep_alleles:
            self._cached['peptide_alleles'] = self.NULL_VALUE
            self._cached['translation_location'] = self.NULL_VALUE
            return
        
        # we only want unique allele strings
        allele_location = dict(list(zip(pep_alleles, translation_location)))
        pep_alleles = list(set(pep_alleles))
        pep_alleles = [pep_alleles, pep_alleles[0]][len(pep_alleles) == 1]
        if type(pep_alleles) not in (str, str):
            for pep_allele in pep_alleles:
                translation_location = allele_location[pep_allele]
        else:
            translation_location = allele_location[pep_alleles]

        self._table_rows[table_name] = dict(pep_allele_string=pep_alleles,
                                            translation_start=translation_location)
        self._populate_cache_from_record(attr_column_map, table_name)

    def _get_peptide_variation(self):
        return self._get_cached_value('peptide_alleles',
                                      self._get_transcript_record)

    peptide_alleles = property(_get_peptide_variation)

    def _get_translation_location(self):
        return self._get_cached_value('translation_location',
                                      self._get_transcript_record)

    translation_location = property(_get_translation_location)

    def _split_alleles(self):
        return self.alleles.split('/')

    def _get_number_alleles(self):
        result = self._split_alleles()
        return len(result)

    num_alleles = property(_get_number_alleles)


class CpGisland(GenericRegion):
    type = 'CpGisland'

    def __init__(self, genome, db, location, Score):
        super(CpGisland, self).__init__(genome=genome, db=db,
                                        location=location)
        self.Score = Score

    def __str__(self):
        my_type = self.__class__.__name__

        return "%s(coord_name='%s'; start=%s; end=%s; length=%s;"\
               " strand='%s', Score=%.1f)" % (my_type,
                                              self.location.coord_name,
                                              self.location.start,
                                              self.location.end,
                                              len(self),
                                              '-+'[self.location.strand > 0], self.Score)


class Repeat(GenericRegion):
    type = 'repeat'

    def __init__(self, genome, db, location, Score, data):
        super(Repeat, self).__init__(genome=genome, db=db, location=location)
        self._attr_ensembl_table_map = dict(symbol='repeat_consensus',
                                            repeat_type='repeat_consensus',
                                            repeat_class='repeat_consensus',
                                            consensus='repeat_consensus')

        self.Score = Score
        # assume always created from repeat_feature table
        self._table_rows['repeat_feature'] = data

    def __str__(self):
        my_type = self.__class__.__name__

        return "%s(coord_name='%s'; start=%s; end=%s; length=%s;"\
               " strand='%s', Score=%.1f)" % (my_type,
                                              self.location.coord_name,
                                              self.location.start, self.location.end, len(
                                                  self),
                                              '-+'[self.location.strand > 0], self.Score)

    def _get_repeat_consensus_record(self):
        repeat_consensus_table = self.db.get_table('repeat_consensus')
        repeat_consensus_id = self._table_rows[
            'repeat_feature']['repeat_consensus_id']
        record = sql.select([repeat_consensus_table],
                            repeat_consensus_table.c.repeat_consensus_id == repeat_consensus_id)
        record = asserted_one(record.execute().fetchall())
        self._table_rows['repeat_consensus'] = record
        limit_length = lambda x: DisplayString(x, repr_length=10)
        attr_column_map = [('symbol', 'repeat_name', _quoted),
                           ('repeat_class', 'repeat_class', _quoted),
                           ('repeat_type', 'repeat_type', _quoted),
                           ('consensus', 'repeat_consensus', limit_length)]
        self._populate_cache_from_record(attr_column_map, 'repeat_consensus')

    def _get_symbol(self):
        return self._get_cached_value('symbol',
                                      self._get_repeat_consensus_record)

    symbol = property(_get_symbol)

    def _get_repeat_class(self):
        return self._get_cached_value('repeat_class',
                                      self._get_repeat_consensus_record)

    repeat_class = property(_get_repeat_class)

    def _get_repeat_type(self):
        return self._get_cached_value('repeat_type',
                                      self._get_repeat_consensus_record)

    repeat_type = property(_get_repeat_type)

    def _get_consensus(self):
        return self._get_cached_value('consensus',
                                      self._get_repeat_consensus_record)

    consensus = property(_get_consensus)
