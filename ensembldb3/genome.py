import re

import sqlalchemy as sql

from .assembly import (Coordinate, CoordSystem, get_coord_conversion,
                       location_query)
from .database import Database
from .feature_level import FeatureCoordLevels
from .host import get_ensembl_account, get_latest_release
from .region import (CpGisland, Est, Gene, GenericRegion, Repeat, Transcript,
                     Variation)
from .species import Species as _Species
from .util import (DisplayString, LazyRecord, asserted_one, convert_strand,
                   flatten)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "3.0a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


class FeatureTypeCache(LazyRecord):
    """stores critical indices for different feature types"""

    def __init__(self, genome):
        super(FeatureTypeCache, self).__init__()
        self.genome = genome
        self._type_func_map = dict(
            CpGisland=self._get_cpg_island_analysis_id, Repeat=self._get_repeat_id
        )

    def _get_cpg_island_analysis_id(self):
        analysis_description_table = self.genome.CoreDb.get_table(
            "analysis_description"
        )
        query = sql.select(
            [analysis_description_table.c.analysis_id],
            analysis_description_table.c.display_label.like("%CpG%"),
        )
        record = asserted_one(query.execute())
        self._table_rows["analysis_description"] = record
        self._populate_cache_from_record(
            [
                (
                    "CpGisland",
                    "analysis_id",
                    lambda x: DisplayString(x, with_quotes=True, num_words=2),
                )
            ],
            "analysis_description",
        )

    def _get_cpg_island_id(self):
        return self._get_cached_value("CpGisland", self._get_cpg_island_analysis_id)

    CpGisland = property(_get_cpg_island_id)

    def _get_repeat_id(self):
        raise NotImplementedError

    Repeat = property(_get_repeat_id)

    def get(self, feature_type):
        """returns the analysis_id for feature_type"""
        try:
            func = self._type_func_map[feature_type]
        except KeyError:
            raise RuntimeError("Unknown feature type: %s" % feature_type)
        return self._get_cached_value(feature_type, func)


class Genome(object):
    """An Ensembl Genome"""

    def __init__(self, species, release, account=None, pool_recycle=None):
        super(Genome, self).__init__()

        assert release, "invalid release specified"
        if account is None:
            account = get_ensembl_account(release=release)

        self._account = account
        self._pool_recycle = pool_recycle

        # TODO: check release may not be necessary because: assert release
        # above
        if release is None:
            release = get_latest_release(account=account)

        self._gen_release = None

        # TODO make name and release immutable properties
        self.species = _Species.get_species_name(species)
        self.release = str(release)

        # the db connections
        self._core_db = None
        self._var_db = None
        self._other_db = None
        self._feature_type_ids = FeatureTypeCache(self)
        self._feature_coord_levels = FeatureCoordLevels(self.species)

    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(species='%s'; release='%s')" % (my_type, self.species, self.release)

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.CoreDb < other.CoreDb

    def __eq__(self, other):
        return self.CoreDb == other.CoreDb

    def __ne__(self, other):
        return self.CoreDb != other.CoreDb

    def __hash__(self):
        return hash((self.species, self.release, self._account))

    def _connect_db(self, db_type):
        connection = dict(
            account=self._account,
            release=self.release,
            species=self.species,
            pool_recycle=self._pool_recycle,
        )
        if self._core_db is None and db_type == "core":
            self._core_db = Database(db_type="core", **connection)
            gen_rel = self.CoreDb.db_name.general_release
            gen_rel = int(re.findall(r"^\d+", str(gen_rel))[0])
            self._gen_release = gen_rel
        elif self._var_db is None and db_type == "variation":
            self._var_db = Database(db_type="variation", **connection)
        elif self._other_db is None and db_type == "otherfeatures":
            self._other_db = Database(db_type="otherfeatures", **connection)

    def _get_core_db(self):
        self._connect_db("core")
        return self._core_db

    CoreDb = property(_get_core_db)

    def _get_var_db(self):
        self._connect_db("variation")
        return self._var_db

    VarDb = property(_get_var_db)

    def _get_other_db(self):
        self._connect_db("otherfeatures")
        return self._other_db

    OtherFeaturesDb = property(_get_other_db)

    @property
    def general_release(self):
        """returns True if the general Ensembl release is >= 65"""
        # General release is used here as to support Ensembl genomes
        if self._gen_release is None:
            self.CoreDb

        return self._gen_release

    def _get_biotype_description_condition(
        self, gene_table, description=None, biotype=None, like=True
    ):
        assert description or biotype, "no valid argument provided"
        btype, descr = None, None

        if biotype:
            if like:
                btype = gene_table.c.biotype.like("%" + biotype + "%")
            else:
                btype = gene_table.c.biotype == biotype
        if description:
            if like:
                descr = gene_table.c.description.like("%" + description + "%")
            else:
                descr = gene_table.c.description.op("regexp")(
                    "[[:<:]]%s[[:>:]]" % description
                )

        if btype is not None and descr is not None:
            condition = sql.and_(btype, descr)
        elif btype is not None:
            condition = btype
        elif descr is not None:
            condition = descr

        return condition

    def _build_gene_query(
        self, db, condition, gene_table, gene_id_table, xref_table=None
    ):
        if gene_id_table is None:  # Ensembl releases later than >= 65
            join_obj = gene_table
            select_obj = [gene_table]
        else:
            join_obj = gene_id_table.join(
                gene_table, gene_id_table.c.gene_id == gene_table.c.gene_id
            )
            select_obj = [gene_id_table.c.stable_id, gene_table]

        if db.type == "core":
            join_obj = join_obj.outerjoin(
                xref_table, gene_table.c.display_xref_id == xref_table.c.xref_id
            )
            select_obj.append(xref_table.c.display_label)
        query = sql.select(select_obj, from_obj=[join_obj], whereclause=condition)
        return query

    def _get_symbol_from_synonym(self, db, synonym):
        """returns the gene symbol for a synonym"""
        synonym_table = db.get_table("external_synonym")
        xref_table = db.get_table("xref")
        joinclause = xref_table.join(
            synonym_table, xref_table.c.xref_id == synonym_table.c.xref_id
        )
        whereclause = synonym_table.c.synonym == synonym
        query = sql.select(
            [xref_table.c.display_label], from_obj=[joinclause], whereclause=whereclause
        ).distinct()
        result = query.execute().fetchall()
        if result:
            try:
                symbol = flatten(result)[0]
            except IndexError:
                symbol = None
        else:
            symbol = None
        return symbol

    def _get_gene_query(
        self,
        db,
        symbol=None,
        description=None,
        stableid=None,
        biotype=None,
        synonym=None,
        like=True,
    ):
        xref_table = [None, db.get_table("xref")][db.type == "core"]
        gene_table = db.get_table("gene")

        # after release 65, the gene_id_table is removed. The following is to
        # maintain support for earlier releases
        release_ge_65 = self.general_release >= 65
        if release_ge_65:
            gene_id_table = None
        else:
            gene_id_table = db.get_table("gene_stable_id")

        assert (
            symbol or description or stableid or biotype
        ), "no valid argument provided"
        if symbol:
            condition = xref_table.c.display_label == symbol
        elif stableid and release_ge_65:
            condition = gene_table.c.stable_id == stableid
        elif stableid:
            condition = gene_id_table.c.stable_id == stableid
        else:
            condition = self._get_biotype_description_condition(
                gene_table, description, biotype, like
            )

        query = self._build_gene_query(
            db, condition, gene_table, gene_id_table, xref_table
        )

        return query

    def make_location(
        self, coord_name, start=None, end=None, strand=1, ensembl_coord=False
    ):
        """returns a location in the genome"""
        return Coordinate(
            self,
            coord_name=coord_name,
            start=start,
            end=end,
            strand=strand,
            ensembl_coord=ensembl_coord,
        )

    def get_gene_by_stableid(self, stableid):
        """returns the gene matching stableid, or None if no record found"""
        query = self._get_gene_query(self.CoreDb, stableid=stableid)
        try:
            record = list(query.execute())[0]
            gene = Gene(self, self.CoreDb, data=record)
        except IndexError:
            gene = None
        return gene

    def get_genes_matching(
        self,
        symbol=None,
        description=None,
        stableid=None,
        biotype=None,
        like=True,
        limit=None,
    ):
        """returns a generator of Gene instances

        Arguments:
            - symbol: HGC gene symbol, case doesn't matter
            - description: a functional description
            - stableid: the ensebl identifier
            - biotype: the biological encoding type
            - like: allow incomplete matches
            - limit: only return this number of hits"""
        # TODO additional arguments to satisfy: external_ref, go_terms
        if symbol is not None:
            symbol = symbol.lower()
        # biotype -> gene
        # description -> gene
        # Symbols -> xref
        # stableid -> gene_stable_id
        # XREF table calls
        # for gene symbols, these need to be matched against the display_label
        # attribute of core.xref table
        # for description, these need to be matched against the description
        # field of the xref table

        # TODO catch conditions where user passes in both a symbol and a
        # biotype
        args = dict(
            symbol=symbol,
            description=description,
            stableid=stableid,
            biotype=biotype,
            like=like,
        )
        query = self._get_gene_query(self.CoreDb, **args)
        if limit is not None:
            query = query.limit(limit)

        records = query.execute()
        if records.rowcount == 0 and symbol is not None:
            # see if the symbol has a synonym
            symbol = self._get_symbol_from_synonym(self.CoreDb, symbol)
            if symbol is not None:
                args["symbol"] = symbol
                records = self._get_gene_query(self.CoreDb, **args).execute()
            else:
                records = []

        for record in records:
            gene = Gene(self, self.CoreDb, data=record)
            yield gene

    def get_transcript_by_stableid(self, stableid):
        """returns the transcript matching stableid,
        or None if no record found"""
        query = self._get_transcript_query(self.CoreDb, stableid=stableid)
        try:
            record = list(query.execute())[0]
            transcript_id = record[0]
            transcript = Transcript(self, self.CoreDb, transcript_id, data=record)
        except IndexError:
            transcript = None
        return transcript

    def _get_transcript_query(
        self,
        db,
        symbol=None,
        description=None,
        stableid=None,
        biotype=None,
        synonym=None,
        like=True,
    ):
        xref_table = [None, db.get_table("xref")][db.type == "core"]
        transcript_table = db.get_table("transcript")

        # after release 65, the transcript_id_table is removed. The following
        # is to maintain support for earlier releases
        release_ge_65 = self.general_release >= 65
        if release_ge_65:
            transcript_id_table = None
        else:
            transcript_id_table = db.get_table("transcript_stable_id")

        assert (
            symbol or description or stableid or biotype
        ), "no valid argument provided"
        if symbol:
            condition = xref_table.c.display_label == symbol
        elif stableid and release_ge_65:
            condition = transcript_table.c.stable_id == stableid
        elif stableid:
            condition = transcript_id_table.c.stable_id == stableid
        else:
            condition = self._get_biotype_description_condition(
                transcript_table, description, biotype, like
            )

        query = self._build_transcript_query(
            db, condition, transcript_table, transcript_id_table, xref_table
        )

        return query

    def _build_transcript_query(
        self, db, condition, transcript_table, transcript_id_table, xref_table=None
    ):
        if transcript_id_table is None:  # Ensembl releases later than >= 65
            join_obj = transcript_table
            select_obj = [transcript_table]
        else:
            join_obj = transcript_id_table.join(
                transcript_table,
                transcript_id_table.c.gene_id == transcript_table.c.transcript_id,
            )
            select_obj = [transcript_id_table.c.stable_id, transcript_table]

        if db.type == "core":
            join_obj = join_obj.outerjoin(
                xref_table, transcript_table.c.display_xref_id == xref_table.c.xref_id
            )
            select_obj.append(xref_table.c.display_label)
        query = sql.select(select_obj, from_obj=[join_obj], whereclause=condition)
        return query

    def get_est_matching(self, stableid):
        """returns an Est object from the otherfeatures db with the stableid"""
        query = self._get_gene_query(self.OtherFeaturesDb, stableid=stableid)
        records = query.execute()
        for record in records:
            yield Est(self, self.OtherFeaturesDb, stableid=stableid, data=record)

    def _get_seq_region_id(self, coord_name):
        """returns the seq_region_id for the provided coord_name"""
        seq_region_table = self.CoreDb.get_table("seq_region")
        coord_systems = CoordSystem(core_db=self.CoreDb)
        coord_system_ids = [k for k in coord_systems if type(k) not in (str, str)]
        record = sql.select(
            [seq_region_table.c.seq_region_id],
            sql.and_(
                seq_region_table.c.name == coord_name,
                seq_region_table.c.coord_system_id.in_(coord_system_ids),
            ),
        )
        record = asserted_one(record.execute().fetchall())
        return record["seq_region_id"]

    def _get_simple_features(
        self, db, klass, target_coord, query_coord, where_feature, limit=None
    ):
        """returns feature_type records for the query_coord from the
        simple_feature table. The returned coord is referenced to
        target_coord. At present, only CpG islands being queried."""
        simple_feature_table = db.get_table("simple_feature")
        feature_types = ["CpGisland"]
        feature_type_ids = [str(self._feature_type_ids.get(f)) for f in feature_types]
        # fix the following
        query = sql.select(
            [simple_feature_table],
            sql.and_(
                simple_feature_table.c.analysis_id.in_(feature_type_ids),
                simple_feature_table.c.seq_region_id == query_coord.seq_region_id,
            ),
        )
        query = location_query(
            simple_feature_table,
            query_coord.ensembl_start,
            query_coord.ensembl_end,
            query=query,
            where=where_feature,
        )
        if limit is not None:
            query = query.limit(limit)

        records = query.execute()
        for record in records:
            coord = Coordinate(
                self,
                coord_name=query_coord.coord_name,
                start=record["seq_region_start"],
                end=record["seq_region_end"],
                seq_region_id=record["seq_region_id"],
                strand=record["seq_region_strand"],
                ensembl_coord=True,
            )
            if query_coord.coord_name != target_coord.coord_name:
                coord = asserted_one(
                    get_coord_conversion(coord, target_coord.coord_type, self.CoreDb)
                )[1]

            yield klass(self, db, location=coord, Score=record["score"])

    def _get_repeat_features(
        self, db, klass, target_coord, query_coord, where_feature, limit=None
    ):
        """returns Repeat region instances"""
        # we build repeats using coordinates from repeat_feature table
        # the repeat_consensus_id is required to get the repeat name, class
        # and type
        repeat_feature_table = db.get_table("repeat_feature")
        query = sql.select(
            [repeat_feature_table],
            repeat_feature_table.c.seq_region_id == query_coord.seq_region_id,
        )
        query = location_query(
            repeat_feature_table,
            query_coord.ensembl_start,
            query_coord.ensembl_end,
            query=query,
            where=where_feature,
        )
        if limit is not None:
            query = query.limit(limit)

        for record in query.execute():
            coord = Coordinate(
                self,
                coord_name=query_coord.coord_name,
                start=record["seq_region_start"],
                end=record["seq_region_end"],
                seq_region_id=record["seq_region_id"],
                strand=record["seq_region_strand"],
                ensembl_coord=True,
            )
            if query_coord.coord_name != target_coord.coord_name:
                coord = asserted_one(
                    get_coord_conversion(coord, target_coord.coord_type, self.CoreDb)
                )[1]
            yield klass(self, db, location=coord, Score=record["score"], data=record)

    def _get_gene_features(
        self, db, klass, target_coord, query_coord, where_feature, limit=None
    ):
        """returns all genes"""
        xref_table = [None, db.get_table("xref")][db.type == "core"]
        gene_table = db.get_table("gene")

        # after release 65, the gene_id_table is removed. The following is
        # to maintain support for earlier releases.
        if self.general_release >= 65:
            gene_id_table = None
        else:
            gene_id_table = db.get_table("gene_stable_id")

        # note gene records are at chromosome, not contig, level
        condition = gene_table.c.seq_region_id == query_coord.seq_region_id
        query = self._build_gene_query(
            db, condition, gene_table, gene_id_table, xref_table
        )
        query = location_query(
            gene_table,
            query_coord.ensembl_start,
            query_coord.ensembl_end,
            query=query,
            where=where_feature,
        )
        if limit is not None:
            query = query.limit(limit)

        for record in query.execute():
            new = Coordinate(
                self,
                coord_name=query_coord.coord_name,
                start=record["seq_region_start"],
                end=record["seq_region_end"],
                strand=record["seq_region_strand"],
                seq_region_id=record["seq_region_id"],
                ensembl_coord=True,
            )

            gene = klass(self, db, location=new, data=record)
            yield gene

    def _get_variation_features(
        self, db, klass, target_coord, query_coord, where_feature, limit=None
    ):
        """returns variation instances within the specified region"""
        # variation features at supercontig level
        var_feature_table = self.VarDb.get_table("variation_feature")
        # note gene records are at chromosome, not contig, level
        query = sql.select(
            [var_feature_table],
            var_feature_table.c.seq_region_id == query_coord.seq_region_id,
        )
        query = location_query(
            var_feature_table,
            query_coord.ensembl_start,
            query_coord.ensembl_end,
            query=query,
            where=where_feature,
        )
        if limit is not None:
            query = query.limit(limit)

        for record in query.execute():
            yield klass(self, self.CoreDb, symbol=record["variation_name"], data=record)

    def _get_feature_coord_levels(self, feature_types):
        dbs = dict(core_db=self.CoreDb)
        if "variation" in feature_types:
            dbs["var_db"] = self.VarDb
        if "est" in feature_types:
            dbs["otherfeature_db"] = self.OtherFeaturesDb
        feature_coord_levels = self._feature_coord_levels(
            self.species, feature_types=feature_types, **dbs
        )
        return feature_coord_levels

    def _feature_coord_levels(self):
        if str(self._feature_coord_levels):
            return self._feature_coord_levels
        feature_types = ["gene", "est", "variation", "cpg", "repeat"]
        feature_coord_levels = self._get_feature_coord_levels(feature_types)
        return self._feature_coord_levels

    feature_coord_levels = property(_feature_coord_levels)

    def get_features(
        self,
        region=None,
        feature_types=None,
        where_feature=None,
        coord_name=None,
        start=None,
        end=None,
        strand=None,
        ensembl_coord=False,
        limit=None,
    ):
        """returns region instances for the specified location"""
        if isinstance(feature_types, str):
            feature_types = [feature_types]
        feature_types = [ft.lower() for ft in feature_types]
        feature_coord_levels = self._get_feature_coord_levels(feature_types)

        if region is None:
            seq_region_id = self._get_seq_region_id(coord_name)
            region = Coordinate(
                self,
                coord_name=coord_name,
                start=start,
                end=end,
                strand=convert_strand(strand),
                seq_region_id=seq_region_id,
                ensembl_coord=ensembl_coord,
            )
        elif hasattr(region, "location"):
            region = region.location

        coord = region
        # the coordinate system at which locations are to be referenced, and
        # the processing function
        target_coords_funcs = dict(
            cpg=(self._get_simple_features, CpGisland),
            repeat=(self._get_repeat_features, Repeat),
            gene=(self._get_gene_features, Gene),
            est=(self._get_gene_features, Est),
            variation=(self._get_variation_features, Variation),
        )

        known_types = set(target_coords_funcs.keys())
        if not set(feature_types) <= known_types:
            raise RuntimeError(
                "Unknown feature[%s], valid feature_types \
                are: %s"
                % (set(feature_types) ^ known_types, known_types)
            )

        for feature_type in feature_types:
            target_func, target_class = target_coords_funcs[feature_type]
            db = self.CoreDb
            if feature_type == "est":
                db = self.OtherFeaturesDb

            feature_coords = feature_coord_levels[feature_type].levels
            for feature_coord in feature_coords:
                chrom_other_coords = get_coord_conversion(
                    coord, feature_coord, db, where=where_feature
                )
                for chrom_coord, other_coord in chrom_other_coords:
                    for region in target_func(
                        db,
                        target_class,
                        chrom_coord,
                        other_coord,
                        where_feature,
                        limit=limit,
                    ):
                        yield region

    def get_variation(
        self,
        effect=None,
        symbol=None,
        like=True,
        validated=False,
        somatic=False,
        flanks_match_ref=False,
        limit=None,
    ):
        """returns a generator of Variation instances

        Arguments:
            - effect: the coding impact, eg. nonsynonymous
            - like: effect is exactly matched against records like that
              provided
            - symbol: the external or ensembl identifier - returns the exact
              match
            - validated: variant has validation != None
            - somatic: exclude somatic mutations
            - flanks_match_ref: flanking sequence matches the reference
            - limit: only return this number of hits"""
        var_feature_table = self.VarDb.get_table("variation_feature")

        assert effect or symbol, "No arguments provided"
        #  if we don't have symbol, then we deal with effect

        consequence_type = "consequence_type"
        if self.general_release > 67:
            consequence_type += "s"  # change to plural column name

        if effect is not None:
            if like:
                query = var_feature_table.columns[consequence_type].like(
                    "%" + effect + "%"
                )
            else:
                query = var_feature_table.columns[consequence_type] == effect
        else:
            query = var_feature_table.c.variation_name == symbol

        if validated:
            validated_col = "evidence_attribs"
            if self.general_release < 83:
                validated_col = "validation_status"

            null = None
            if int(self.release) >= 65:
                null = ""

            query = sql.and_(query, var_feature_table.columns[validated_col] != null)

        if not somatic:
            query = sql.and_(query, var_feature_table.c.somatic != 1)

        if flanks_match_ref:
            query = sql.and_(query, var_feature_table.c.alignment_quality == 1)

        query = sql.select([var_feature_table], query).order_by(
            var_feature_table.c.seq_region_start
        )

        if limit:
            query = query.limit(limit)

        for record in query.execute():
            yield Variation(
                self, self.CoreDb, effect=effect, symbol=symbol, data=record
            )

    def get_region(
        self,
        region=None,
        coord_name=None,
        start=None,
        end=None,
        strand=None,
        ensembl_coord=False,
    ):
        """returns a single generic region for the specified coordinates
        Arguments:
            - region: a genomic region or a Coordinate instance
            - ensembl_coords: if True, follows indexing system of Ensembl
              where indexing starts at 1"""
        if region is None:
            seq_region_id = self._get_seq_region_id(coord_name)
            region = Coordinate(
                self,
                coord_name=coord_name,
                start=start,
                end=end,
                strand=convert_strand(strand),
                seq_region_id=seq_region_id,
                ensembl_coord=ensembl_coord,
            )
        elif hasattr(region, "location"):
            region = region.location

        return GenericRegion(
            self,
            self.CoreDb,
            coord_name=coord_name,
            start=start,
            end=end,
            strand=strand,
            location=region,
            ensembl_coord=ensembl_coord,
        )

    def get_distinct(self, property_type):
        """returns the Ensembl data-bases distinct values for the named
        property_type.

        Arguments:
            - property_type: valid values are biotype, status (pre release 90)
              effect"""
        property_type = property_type.lower()
        if property_type == "effect":
            db = self.VarDb
        else:
            db = self.CoreDb

        consequence_type = "consequence_type"
        if self.general_release > 67:
            consequence_type += "s"  # change to plural column name

        property_map = {
            "effect": ("variation_feature", consequence_type),
            "biotype": ("gene", "biotype"),
        }

        if self.general_release < 90:
            property_map["status"] = ("gene", "status")

        if property_type not in property_map:
            raise RuntimeError("ERROR: Unknown property type: %s" % property_type)

        table_name, column = property_map[property_type]
        return list(db.get_distinct(table_name, column))
