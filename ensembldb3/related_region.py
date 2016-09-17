from pprint import pprint
from collections import defaultdict
import warnings

import sqlalchemy as sql
from cogent3 import DNA, LoadTree
from cogent3.core.tree import PhyloNode
from cogent3.core.alignment import SequenceCollection, Alignment, Aligned
from cogent3.parse import cigar

from .util import LazyRecord, asserted_one, NoItemError
from .assembly import location_query
from .species import Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"


class _RelatedRegions(LazyRecord):
    # a basic related region, capable of providing the sequences
    # obtaining SyntenicRegions -- for getting aligned blocks -- is delegated
    # to compara
    type = None

    def __init__(self):
        super(_RelatedRegions, self).__init__()

    def __str__(self):
        # temporary string method, just to demo correct assembly
        # TODO StableID and description
        my_type = self.__class__.__name__

        data = list(map(repr, self.members))
        data.insert(0, '%s(' % my_type)
        data.append(')')
        return "\n\t".join(data)

    def get_seq_collection(self, feature_types=None, where_feature=None):
        """returns a SequenceCollection instance of the unaligned sequences"""
        seqs = []
        for member in self.members:
            if feature_types:
                seq = member.get_annotated_seq(feature_types, where_feature)
            else:
                seq = member.seq
            if seq is None:
                continue
            seqs.append((seq.name, seq))
        return SequenceCollection(data=seqs, moltype=DNA)

    def getseq_lengths(self):
        """returns a vector of lengths"""
        return [len(member) for member in self.members]

    def get_species_set(self):
        """returns the latin names of self.Member species as a set"""
        return set([m.location.species for m in self.members
                    if m.location is not None])


class RelatedGenes(_RelatedRegions):
    type = 'related_genes'

    def __init__(self, compara, members, relationships, gene_tree_root=None):
        super(RelatedGenes, self).__init__()
        self.compara = compara
        self.members = tuple(m for m in members if m.location is not None)
        self.relationships = relationships
        # gene_tree_root is the id from homology
        self._gene_tree_root = gene_tree_root

    def __str__(self):
        my_type = self.__class__.__name__

        display = ['%s:' % my_type,
                   ' relationships=%s' % self.relationships]
        display += ['  %s' % m for m in self.members]
        return '\n'.join(display)

    def __repr__(self):
        return self.__str__()

    def get_max_cds_lengths(self):
        """returns the vector of maximum cds lengths from member transcripts"""
        return [max(member.get_cds_lengths()) for member in self.members]
    
    def get_tree(self, just_members=True):
        """returns the gene tree with tip names as gene stableid's
        
        Arguments:
        ----------
          - just_members: limits tips to just members of self
        """
        gtrn = self.compara.ComparaDb.get_table("gene_tree_node")
        condition = gtrn.c.root_id == list(self._gene_tree_root)[0]
        query = sql.select([gtrn.c.node_id, gtrn.c.parent_id,
                            gtrn.c.distance_to_parent, gtrn.c.seq_member_id],
                           whereclause=condition)
        records = query.execute().fetchall()
        
        # get the gene stable IDs, via join of seq_member with gene_member
        # on seq_member_id
        seqmem_ids = [r['seq_member_id'] for r in records]
        seqmem = self.compara.ComparaDb.get_table("seq_member")
        genmem = self.compara.ComparaDb.get_table("gene_member")
        condition = genmem.c.gene_member_id == seqmem.c.gene_member_id
        joined = seqmem.join(genmem,
                             genmem.c.gene_member_id == seqmem.c.gene_member_id)
        query = sql.select([joined.c.seq_member_seq_member_id,
                            joined.c.gene_member_stable_id]).\
            where(joined.c.seq_member_seq_member_id.in_(seqmem_ids)).select_from(joined)
        
        gene_ids = dict(query.execute().fetchall())
        nodes = {}
        parents = defaultdict(list)
        for record in records:
            parent_id = record['parent_id']
            node_id = record['node_id']
            length = record['distance_to_parent']
            name = gene_ids.get(record['seq_member_id'], None)
            node = PhyloNode(length=length, name=name)
            nodes[node_id] = node
            parents[parent_id].append(node)
        
        root = None
        for parent in parents:
            if parent not in nodes:
                node = PhyloNode(name='root')
                nodes[parent] = node
            
            node = nodes[parent]
            for child in parents[parent]:
                child.parent = node
            if len(parents[parent]) == 1:
                root = node
        
        if just_members:
            stableids = [g.stableid for g in self.members]
            root = root.get_sub_tree(stableids)
        
        return root
        


class SyntenicRegion(LazyRecord):
    """a class that takes the genome, compara instances and is used to build
    Aligned sequences for Ensembl multiple alignments"""

    def __init__(self, parent, genome, identifiers_values, am_ref_member,
                 location=None):
        # create with method_link_species_set_id, at least, in
        # identifiers_values
        super(SyntenicRegion, self).__init__()
        self.parent = parent
        self.compara = parent.compara
        self.genome = genome
        self.am_ref_member = am_ref_member
        self.aln_map = None
        self.aln_loc = None
        self._make_map_func = [self._make_map_from_ref,
                               self._make_ref_map][am_ref_member]

        if location is not None:
            if hasattr(location, 'location'):  # likely to be a feature region
                region = location
            else:
                region = genome.get_region(region=location)
            self._cached['region'] = region

        for identifier, value in list(dict(identifiers_values).items()):
            self._cached[identifier] = value

    def __len__(self):
        return len(self._get_cached_value('region', self._make_map_func))

    def _get_location(self):
        region = self._get_cached_value('region', self._make_map_func)
        location = None if region is None else region.location
        return location

    location = property(_get_location)

    def _get_region(self):
        region = self._get_cached_value('region', self._make_map_func)
        return region

    region = property(_get_region)

    def _get_cigar_record(self):
        genomic_align_table = \
            self.parent.compara.ComparaDb.get_table('genomic_align')
        query = sql.select([genomic_align_table.c.cigar_line],
                           genomic_align_table.c.genomic_align_id ==
                           self._cached['genomic_align_id'])
        record = asserted_one(query.execute())
        self._cached['cigar_line'] = record['cigar_line']
        return record

    def _get_cigar_line(self):
        return self._get_cached_value('cigar_line', self._get_cigar_record)

    cigar_line = property(_get_cigar_line)

    def _make_ref_map(self):
        if self.aln_map and self.aln_loc is not None:
            return

        ref_record = self._cached
        record_start = ref_record['dnafrag_start']
        record_end = ref_record['dnafrag_end']
        record_strand = ref_record['dnafrag_strand']

        block_loc = self.genome.make_location(coord_name=ref_record['name'],
                                             start=record_start,
                                             end=record_end,
                                             strand=record_strand,
                                             ensembl_coord=True)

        ref_location = self.parent.ref_location
        relative_start = ref_location.start - block_loc.start
        relative_end = relative_start + len(ref_location)
        if block_loc.strand != 1:
            relative_start = len(block_loc) - relative_end
            relative_end = relative_start + len(ref_location)

        aln_map, aln_loc = cigar.slice_cigar(self.cigar_line, relative_start,
                                             relative_end, by_align=False)

        self.aln_map = aln_map
        self.aln_loc = aln_loc
        region_loc = ref_location.copy()
        region_loc.strand = block_loc.strand
        region = self.genome.get_region(region=region_loc)
        self._cached['region'] = region

    def _make_map_from_ref(self):
        # this is the 'other' species
        if self.aln_loc and self.aln_map is not None:
            return
        
        record = self._cached
        aln_map, aln_loc = cigar.slice_cigar(self.cigar_line,
                                             self.parent.cigar_start,
                                             self.parent.cigar_end,
                                             by_align=True)
        if not aln_loc:
            self._cached['region'] = self.NULL_VALUE
            return
        
        self.aln_map = aln_map
        self.aln_loc = aln_loc  # probably unnecesary to store??

        # we make a loc for the aligned region
        block_loc = self.genome.make_location(coord_name=record['name'],
                                             start=record['dnafrag_start'],
                                             end=record['dnafrag_end'],
                                             strand=record[
                                                 'dnafrag_strand'],
                                             ensembl_coord=True)
        relative_start = aln_loc[0]
        relative_end = aln_loc[1]
        # new location with correct length ensembl_start
        loc = block_loc.copy()
        loc.end = loc.start + (relative_end - relative_start)

        if block_loc.strand != 1:
            shift = len(block_loc) - relative_end
        else:
            shift = relative_start
        loc = loc.shifted(shift)
        region = self.genome.get_region(region=loc)
        self._cached['region'] = region

    def _make_aligned(self, feature_types=None, where_feature=None):
        if self.aln_loc is None or self.aln_map is None:  # is this required?
            self._make_map_func()
        region = self._cached['region']
        if region is None:
            self._cached['aligned_seq'] = None
            return
        if feature_types:
            seq = region.get_annotated_seq(feature_types, where_feature)
        else:
            seq = region.seq

        # we get the seq objects to allow for copying of their annotations
        gapped_seq = Aligned(self.aln_map, seq)

        self._cached['aligned_seq'] = gapped_seq

    def _get_aligned_seq(self):
        aligned = self._get_cached_value('aligned_seq', self._make_aligned)
        return aligned

    aligned_seq = property(_get_aligned_seq)

    def get_annotated_aligned(self, feature_types, where_feature=None):
        """returns aligned seq annotated for the specified feature types"""
        region = self._get_cached_value('region', self._make_map_func)
        if region is None:
            return None
        self._make_aligned(feature_types=feature_types,
                           where_feature=where_feature)
        return self.aligned_seq


class SyntenicRegions(_RelatedRegions):
    type = 'syntenic_regions'

    def __init__(self, compara, members, ref_location, method_clade_id=None):
        super(SyntenicRegions, self).__init__()
        self.compara = compara
        self._method_clade_id = method_clade_id
        mems = []
        ref_member = None
        self.ref_location = ref_location
        for genome, data in members:
            if genome is ref_location.genome:
                ref_member = SyntenicRegion(self, genome, dict(data),
                                            am_ref_member=True, location=ref_location)
            else:
                mem = SyntenicRegion(self, genome, dict(data),
                                           am_ref_member=False)
                mems.append(mem)
                

        assert ref_member is not None, "Can't match a member to ref_location"
        self.ref_member = ref_member
        
        self.members = tuple([ref_member] + mems)
        self.num_members = len(self.members)
        self.aln_loc = None
        self._do_rc = None

    def __str__(self):
        my_type = self.__class__.__name__

        display = ['%s:' % my_type]
        display += ['  %r' % m.location for m in self.members
                    if m.region is not None]
        return '\n'.join(display)

    def __repr__(self):
        return self.__str__()

    def _populate_ref(self):
        """near (don't actually get the sequence) completes construction of
        ref sequence"""
        self.ref_member._make_map_func()
        self._cached['cigar_start'] = self.ref_member.aln_loc[0]
        self._cached['cigar_end'] = self.ref_member.aln_loc[1]

    def _get_rc_state(self):
        """determines whether the ref_member strand is the same as that from
        the align block, if they diff we will rc the alignment, seqs,
        seq_names"""
        if self._do_rc is not None:
            return self._do_rc
        self._populate_ref()
        inferred = self.ref_member._cached['region'].location.strand
        self._do_rc = self.ref_location.strand != inferred
        return self._do_rc

    _rc = property(fget=_get_rc_state)

    def __len__(self):
        return self.cigar_end - self.cigar_start

    def _get_ref_start(self):
        return self._get_cached_value('cigar_start', self._populate_ref)

    cigar_start = property(_get_ref_start)

    def _get_ref_end(self):
        return self._get_cached_value('cigar_end', self._populate_ref)

    cigar_end = property(_get_ref_end)

    def get_alignment(self, feature_types=None, where_feature=None,
                     omit_redundant=True):
        """Arguments:
            - feature_types: annotations to be applied to the returned
              sequences
            - omit_redundant: exclude redundant gap positions"""
        seqs = []
        annotations = {}

        for member in self.members:
            if feature_types:
                seq = member.get_annotated_aligned(feature_types, where_feature)
            else:
                seq = member.aligned_seq
            if seq is None:
                continue
            name = seq.name

            if self._rc:  # names should reflect change to strand
                loc = member.location.copy()
                loc.strand *= -1
                name = str(loc)

            annotations[name] = seq.data.annotations
            seq.name = seq.data.name = name
            seqs += [(name, seq)]

        if seqs is None:
            return None

        aln = Alignment(data=seqs, moltype=DNA)

        if self._rc:
            aln = aln.rc()

        if omit_redundant:
            aln = aln.omit_gap_pos()

        return aln
    
    def get_species_tree(self):
        """returns the species tree for the members
        
        Arguments:
        ----------
        - just_members: limits tips to species that are just members of self
        """
        tree = self.compara.get_species_tree()
        return tree.get_sub_tree(self.get_species_set())
