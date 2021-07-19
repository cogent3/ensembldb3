from collections import defaultdict

import sqlalchemy as sql

from cogent3.core.tree import PhyloNode
from cogent3.util.table import Table
from numpy import empty

from .assembly import location_query
from .database import Database
from .genome import Genome
from .host import get_ensembl_account
from .related_region import RelatedGenes, SyntenicRegions
from .species import Species as _Species
from .util import NoItemError, asserted_one


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2016-, The EnsemblDb3 Project"
__credits__ = ["Gavin Huttley", "Hua Ying", "Jason Merkin"]
__license__ = "BSD"
__version__ = "2021.04.01"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


class Compara(object):
    """comaparison among genomes"""

    def __init__(
        self, species, release, account=None, pool_recycle=None, division=None
    ):
        assert release, "invalid release specified"
        self.release = str(release)
        self.general_release = None
        if account is None:
            account = get_ensembl_account(release=release)
        self._account = account
        self._pool_recycle = pool_recycle
        self._compara_db = None
        sp = sorted(
            [_Species.get_species_name(sp, level="raise") for sp in set(species)]
        )
        self.species = tuple(sp)
        self._genomes = {}
        self._attach_genomes()

        self._species_id_map = None
        self._species_db_map = None
        self._dbid_species_map = None
        self._species_set = None
        self._method_species_link = None
        self.division = division

    def _attach_genomes(self):
        updated_species = []
        for species in self.species:
            attr_name = _Species.get_compara_name(species)
            genome = Genome(
                species=species, release=self.release, account=self._account
            )
            if self.general_release is None:
                self.general_release = genome.general_release

            species = genome.CoreDb.db_name.species
            updated_species.append(species)
            self._genomes[species] = genome
            setattr(self, attr_name, genome)

        self.species = tuple(updated_species)

    def __str__(self):
        my_type = self.__class__.__name__
        return "%s(species=%s; release=%s; connected=%s)" % (
            my_type,
            self.species,
            self.release,
            self.ComparaDb is not None,
        )

    def _connect_db(self):
        if self._compara_db is None:
            # TODO can the connection be all done in init?
            connection = dict(
                account=self._account,
                release=self.release,
                pool_recycle=self._pool_recycle,
            )
            self._compara_db = Database(
                db_type="compara", division=self.division, **connection
            )

    def _get_compara_db(self):
        self._connect_db()
        return self._compara_db

    ComparaDb = property(_get_compara_db)

    @property
    def _dbid_genome_map(self):
        """maps genome_db id to Genome instances"""
        if self._dbid_species_map is not None:
            return self._dbid_species_map

        db_species = dict(
            [(_Species.get_ensembl_db_prefix(n), n) for n in self.species]
        )
        db_prefixes = list(db_species.keys())
        genome_db_table = self.ComparaDb.get_table("genome_db")
        query = sql.select(
            [genome_db_table.c.genome_db_id, genome_db_table.c.name],
            genome_db_table.c.name.in_(db_prefixes),
        )
        records = query.execute().fetchall()
        data = {r[0]: self._genomes[db_species[r[1]]] for r in records}
        self._dbid_species_map = data
        return self._dbid_species_map

    def get_species_tree(self, just_members=True):
        """returns the species tree

        Arguments:
        ----------
          - just_members: limits tips to just members of self
        """
        # grab the Ensembl species tree root ID
        sptr = self.ComparaDb.get_table("species_tree_root")
        condition = sptr.c.label == "Ensembl"
        query = sql.select([sptr.c.root_id], whereclause=condition)
        records = query.execute().fetchall()
        assert len(records) == 1, records
        root_id = records[0]["root_id"]

        # get the tree nodes
        sptn = self.ComparaDb.get_table("species_tree_node")
        condition = sql.and_(sptn.c.root_id == root_id)
        query = sql.select([sptn], whereclause=condition)
        records = query.execute().fetchall()

        # get the genome db -> name map
        gen_db = self.ComparaDb.get_table("genome_db")
        db_ids = [r["genome_db_id"] for r in records]
        query = sql.select(
            [gen_db.c.genome_db_id, gen_db.c.name],
            whereclause=gen_db.c.genome_db_id.in_(db_ids),
        )
        id_name = dict(query.execute().fetchall())
        for id_, name in id_name.items():
            name = _Species.get_species_name(name)
            id_name[id_] = None if name == "None" else name

        nodes = {}
        parents = defaultdict(list)
        for record in records:
            parent_id = record["parent_id"]
            node_id = record["node_id"]
            length = record["distance_to_parent"]
            name = record["node_name"]
            gen_dbid = record["genome_db_id"]
            n = None if gen_dbid is None else id_name[gen_dbid]
            name = name if n is None else n
            node = PhyloNode(length=length, name=name)
            nodes[node_id] = node
            parents[parent_id].append(node)

        root = None
        for parent, value in parents.items():
            if parent not in nodes:
                node = PhyloNode(name="root")
                nodes[parent] = node

            node = nodes[parent]
            for child in parents[parent]:
                child.parent = node

            if len(value) == 1:
                root = node

        # convert tip-names to match genome db names
        if just_members:
            root = root.get_sub_tree(self.species, tipsonly=True)

        return root

    def _get_species_set(self):
        if self._species_set is not None:
            return self._species_set
        # we make sure the species set contains all species
        species_set_table = self.ComparaDb.get_table("species_set")
        query = sql.select(
            [species_set_table],
            species_set_table.c.genome_db_id.in_(list(self._dbid_genome_map.keys())),
        )
        species_sets = {}
        for record in query.execute():
            gen_id = record["genome_db_id"]
            sp_set_id = record["species_set_id"]
            if sp_set_id in species_sets:
                species_sets[sp_set_id].update([gen_id])
            else:
                species_sets[sp_set_id] = set([gen_id])

        expected = set(self._dbid_genome_map.keys())
        species_set_ids = [
            sp_set
            for sp_set, gen_id in list(species_sets.items())
            if expected <= gen_id
        ]

        self._species_set = species_set_ids
        return self._species_set

    species_set = property(_get_species_set)

    def _get_method_link_species_set(self):
        if self._method_species_link is not None:
            return self._method_species_link

        method_link_table = self.ComparaDb.get_table("method_link")
        query = sql.select(
            [method_link_table],
            method_link_table.c["class"].like("%" + "alignment" + "%"),
        )
        methods = query.execute().fetchall()
        method_link_ids = dict([(r["method_link_id"], r) for r in methods])
        method_link_species_table = self.ComparaDb.get_table("method_link_species_set")
        query = sql.select(
            [method_link_species_table],
            sql.and_(
                method_link_species_table.c.species_set_id.in_(self.species_set),
                method_link_species_table.c.method_link_id.in_(
                    list(method_link_ids.keys())
                ),
            ),
        )
        records = query.execute().fetchall()
        # store method_link_id, type, species_set_id,
        # method_link_species_set.name, class
        header = [
            "method_link_species_set_id",
            "method_link_id",
            "species_set_id",
            "align_method",
            "align_clade",
        ]
        rows = []
        for record in records:
            ml_id = record["method_link_id"]
            sp_set_id = record["species_set_id"]
            ml_sp_set_id = record["method_link_species_set_id"]
            clade_name = record["name"]
            aln_name = method_link_ids[ml_id]["type"]
            rows += [[ml_sp_set_id, ml_id, sp_set_id, aln_name, clade_name]]

        if rows == []:
            rows = empty((0, len(header)))

        t = Table(
            header=header,
            data=rows,
            space=2,
            index_name="method_link_species_set_id",
            title="Align Methods/Clades",
        )
        self._method_species_link = t
        return t

    method_species_links = property(_get_method_link_species_set)

    def get_related_genes(
        self, gene_region=None, stableid=None, relationship=None, DEBUG=False
    ):
        """returns a RelatedGenes instance.

        Arguments:
            - gene_region: a Gene instance
            - stableid: ensembl stable_id identifier
            - relationship: the types of related genes sought"""
        assert gene_region is not None or stableid is not None, "No identifier provided"

        # TODO understand why this has become necessary to suppress warnings
        # in SQLAlchemy 0.6
        relationship = None if relationship is None else str(relationship)

        stableid = stableid or gene_region.stableid

        if self.general_release > 75:
            mem_name = "gene_member"
            mem_id = "gene_member_id"
            frag_strand = "dnafrag_strand"
        else:
            mem_name = "member"
            mem_id = "member_id"
            frag_strand = "chr_strand"

        member_table = self.ComparaDb.get_table(mem_name)
        homology_member_table = self.ComparaDb.get_table("homology_member")
        homology_table = self.ComparaDb.get_table("homology")
        # homolog.c.gene_tree_root_id "The root_id of the gene tree from which
        # the homology is derived"
        member_ids = sql.select(
            [
                member_table.c[mem_id],
                member_table.c.taxon_id,
                member_table.c.genome_db_id,
            ],
            member_table.c.stable_id == str(stableid),
        )

        member_records = member_ids.execute()
        if DEBUG:
            print("member_records", member_records)

        if not member_records:
            return None
        member_record = asserted_one(member_records)
        member_id = member_record[mem_id]

        # in case query gene_region is not provided
        if gene_region is None:
            ref_genome = self._dbid_genome_map[member_record["genome_db_id"]]
            gene_region = ref_genome.get_gene_by_stableid(stableid)

        homology_ids = sql.select(
            [homology_member_table.c.homology_id, homology_member_table.c[mem_id]],
            homology_member_table.c[mem_id] == member_id,
        )
        homology_ids = [r["homology_id"] for r in homology_ids.execute()]
        if not homology_ids:
            return None

        if DEBUG:
            print("1 - homology_ids", homology_ids)

        condition = homology_table.c.homology_id.in_(homology_ids)
        if relationship is not None:
            condition = sql.and_(
                condition, homology_table.c.description == relationship
            )
        homology_records = sql.select(
            [
                homology_table.c.homology_id,
                homology_table.c.description,
                homology_table.c.method_link_species_set_id,
                homology_table.c.gene_tree_root_id,
            ],
            condition,
        )

        homology_ids = []
        gene_tree_roots = set()
        for r in homology_records.execute():
            homology_ids.append(
                (r["homology_id"], (r["description"], r["method_link_species_set_id"]))
            )
            gene_tree_roots.update([r["gene_tree_root_id"]])
        homology_ids = dict(homology_ids)

        if DEBUG:
            print("2 - homology_ids", homology_ids)
        if not homology_ids:
            return None

        ortholog_ids = sql.select(
            [homology_member_table.c[mem_id], homology_member_table.c.homology_id],
            homology_member_table.c.homology_id.in_(list(homology_ids.keys())),
        )

        ortholog_ids = dict(
            [(r[mem_id], r["homology_id"]) for r in ortholog_ids.execute()]
        )

        if DEBUG:
            print("ortholog_ids", ortholog_ids)
        if not ortholog_ids:
            return None

        gene_set = sql.select(
            [
                member_table.c.gene_member_id,
                member_table.c.stable_id,
                member_table.c.dnafrag_strand,
                member_table.c.genome_db_id,
                homology_member_table.c.homology_id,
            ],
            sql.and_(
                member_table.c[mem_id].in_(list(ortholog_ids.keys())),
                member_table.c.genome_db_id.in_(list(self._dbid_genome_map.keys())),
                member_table.c.gene_member_id == homology_member_table.c.gene_member_id,
                homology_member_table.c.homology_id.in_(list(homology_ids.keys())),
                member_table.c.stable_id != stableid,
            ),
        )  # exclude query gene
        related, stableids = defaultdict(list), set()
        for record in gene_set.execute():
            homid = record["homology_id"]
            # stableid has been taken by the query gene
            sid = record["stable_id"]
            assert sid not in stableids  # no repeated record for the same gene
            stableids.update([sid])

            genome = self._dbid_genome_map[record["genome_db_id"]]
            gene = genome.get_gene_by_stableid(sid)
            assert gene.location.strand == record[frag_strand]
            # mid is method_link_species_set_id
            reltype, mid = homology_ids[homid]
            related[reltype].append(gene)
        if not related:
            return None

        for reltype in related:
            genes = related[reltype] + [gene_region]
            yield RelatedGenes(
                self, genes, relationship=reltype, gene_tree_root=gene_tree_roots
            )

    def _get_dnafrag_id_for_coord(self, coord):
        """returns the dnafrag_id for the coordnate"""
        dnafrag_table = self.ComparaDb.get_table("dnafrag")
        genome_db_table = self.ComparaDb.get_table("genome_db")

        # column renamed between versions
        prefix = coord.genome.species.lower()
        if int(self.release) > 58:
            prefix = _Species.get_ensembl_db_prefix(prefix)

        query = sql.select(
            [dnafrag_table.c.dnafrag_id, dnafrag_table.c.coord_system_name],
            sql.and_(
                dnafrag_table.c.genome_db_id == genome_db_table.c.genome_db_id,
                genome_db_table.c.name == prefix,
                dnafrag_table.c.name == str(coord.coord_name),
            ),
        )
        try:
            record = asserted_one(query.execute().fetchall())
            dnafrag_id = record["dnafrag_id"]
        except NoItemError:
            raise RuntimeError("No DNA fragment identified")
        return dnafrag_id

    def _get_genomic_align_blocks_for_dna_frag_id(
        self, method_clade_id, dnafrag_id, coord
    ):
        genomic_align_table = self.ComparaDb.get_table("genomic_align")
        query = sql.select(
            [
                genomic_align_table.c.genomic_align_id,
                genomic_align_table.c.genomic_align_block_id,
            ],
            sql.and_(
                genomic_align_table.c.method_link_species_set_id == method_clade_id,
                genomic_align_table.c.dnafrag_id == dnafrag_id,
            ),
        )
        query = location_query(
            genomic_align_table,
            coord.ensembl_start,
            coord.ensembl_end,
            start_col="dnafrag_start",
            end_col="dnafrag_end",
            query=query,
        )

        return query.execute().fetchall()

    def _get_joint_genomic_align_dnafrag(self, genomic_align_block_id):
        genomic_align_table = self.ComparaDb.get_table("genomic_align")
        dnafrag_table = self.ComparaDb.get_table("dnafrag")
        query = sql.select(
            [
                genomic_align_table.c.genomic_align_id,
                genomic_align_table.c.genomic_align_block_id,
                genomic_align_table.c.dnafrag_start,
                genomic_align_table.c.dnafrag_end,
                genomic_align_table.c.dnafrag_strand,
                dnafrag_table,
            ],
            sql.and_(
                genomic_align_table.c.genomic_align_block_id == genomic_align_block_id,
                genomic_align_table.c.dnafrag_id == dnafrag_table.c.dnafrag_id,
                dnafrag_table.c.genome_db_id.in_(list(self._dbid_genome_map.keys())),
            ),
        )
        return query.execute().fetchall()

    def get_syntenic_regions(
        self,
        species=None,
        coord_name=None,
        start=None,
        end=None,
        strand=1,
        ensembl_coord=False,
        region=None,
        align_method=None,
        align_clade=None,
        method_clade_id=None,
    ):
        """returns a SyntenicRegions instance

        Arguments:
            - species: the species name
            - coord_name, start, end, strand: the coordinates for the region
            - ensembl_coord: whether the coordinates are in Ensembl form
            - region: a region instance or a location, in which case the
              coord_name etc .. arguments are ignored
            - align_method, align_clade: the alignment method and clade to use
              Note: the options for this instance can be found by printing
              the method_species_links attribute of this object.
            - method_clade_id: over-rides align_method/align_clade. The entry
              in method_species_links under method_link_species_set_id
        """
        assert (
            align_method and align_clade
        ) or method_clade_id, (
            "Must specify (align_method & align_clade) or method_clade_id"
        )
        if method_clade_id is None:
            for row in self.method_species_links:
                if (
                    align_method.lower() in row["align_method"].lower()
                    and align_clade.lower() in row["align_clade"].lower()
                ):
                    method_clade_id = row["method_link_species_set_id"]

        if method_clade_id is None:
            raise RuntimeError(
                "Invalid align_method[%s] or align_clade "
                "specified[%s]" % (align_method, align_clade)
            )

        if region is None:
            ref_genome = self._genomes[_Species.get_species_name(species)]
            region = ref_genome.make_location(
                coord_name=coord_name,
                start=start,
                end=end,
                strand=strand,
                ensembl_coord=ensembl_coord,
            )
        elif hasattr(region, "location"):
            region = region.location

        # make sure the genome instances match
        ref_genome = self._genomes[region.genome.species]
        if ref_genome is not region.genome:
            # recreate region from our instance
            region = ref_genome.make_location(
                coord_name=region.coord_name,
                start=region.start,
                end=region.end,
                strand=region.strand,
            )

        ref_dnafrag_id = self._get_dnafrag_id_for_coord(region)
        blocks = self._get_genomic_align_blocks_for_dna_frag_id(
            method_clade_id, ref_dnafrag_id, region
        )

        for block in blocks:
            genomic_align_block_id = block["genomic_align_block_id"]
            # we get joint records for these identifiers from
            records = self._get_joint_genomic_align_dnafrag(genomic_align_block_id)
            members = []
            ref_location = None
            for record in records:
                genome = self._dbid_genome_map[record.genome_db_id]
                # we have a case where we getback different coordinate system
                # results for the ref genome. We keep only those that match
                # the coord_name of region

                if genome is region.genome and record.name == region.coord_name:
                    # this is the ref species and we adjust the ref_location
                    # for this block
                    diff_start = record.dnafrag_start - region.ensembl_start
                    shift_start = [0, diff_start][diff_start > 0]
                    diff_end = record.dnafrag_end - region.ensembl_end
                    shift_end = [diff_end, 0][diff_end > 0]
                    try:
                        ref_location = region.resized(shift_start, shift_end)
                    except ValueError:
                        # we've hit some ref genome fragment that matches
                        # but whose coordinates aren't right
                        continue
                elif genome is region.genome:
                    continue
                members += [(genome, record)]
            assert ref_location is not None, "Failed to make the reference" " location"
            yield SyntenicRegions(
                self,
                members,
                ref_location=ref_location,
                method_clade_id=method_clade_id,
            )

    def get_distinct(self, property_type):
        """returns the Ensembl data-bases distinct values for the named
        property_type.

        Arguments:
            - property_type: valid values are relationship"""
        property_type = property_type.lower()
        db = self.ComparaDb
        property_map = {
            "relationship": ("homology", "description"),
            "clade": ("method_link_species_set", "name"),
        }
        if property_type not in property_map:
            raise RuntimeError(f"ERROR: Unknown property type: {property_type}")
        table_name, column = property_map[property_type]
        return list(db.get_distinct(table_name, column))
