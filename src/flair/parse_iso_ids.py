
from collections import namedtuple
import re


#IsoformId      iso_id|iso_name|gene_id
#IsoformId      FL:12-1|ENST123|ENSGXXX
#IsoformId      FL:12-2|ENST123|ENSGXXX

#Fusion isoform     FL:13-1|ENSTXXXX--ENSTYYYY|ENGSXXXX--ENSGYYYY
#Fusion isoform     FL:13-1|FL:13-1|ENGSXXXX--ENSGYYYY


#Fusion isoform     FL:13-1:fusion:X--Y|segment1|X
#Fusion isoform     FL:13-1:fusion:X--Y|segment2|Y


#Fusion isoform     FL:13-1|FL:13-1|X--Y    10
#->
#                   FL:13-1|FL:13-1:X--Y|X    10
#                   FL:13-1|FL:13-1:X--Y|Y    10
#                   FL:14-1|FL:14-1|X    20     20/30 = 0.333
#                   FL:15-1|FL:15-1|Y    30     30/40 = .25




#how to keep the isoform


#iso_id: unique identifier for each transcript
#   contains:
#       tag: origin of iso (currently FL or ANNOT (ANNOT only applies to files that are a conversion/subset of the annotated gtf))
#       sj_id: unique id for set of splice junctions contained in the dataset
#       ends_id: ID for set of ends within the space of the unique splice junction
#       add_id: additional ID (fusion or haplotype) - remove??

#iso_name: additional info about isoform (annotated name)
#   formatting/parsing is less constrained, not required to be in a known format
#   can get reset to iso_id

#when one isoform requires multiple ids:
#   fusion: isoA-locus1, isoA-locus2
#   alleles: isoA-a1, isoA-a2

# fusion ID - don't change the isoform ID


# variant-aware-isoform ID - keep separate?
#   can derive from isoform ID, but it is a separate concept
#   variant identifier: scope could be isoform (could also be gene)
#   var:12


class IsoformInfo(namedtuple('IsoformInfo', ('iso_id','iso_name','gene_id'))):
    __slots__ = ()
    # def __new__(cls, gene_id, sj_id, ends_id=1, tag='FL', iso_name=None, add_id=None):
    #     #sjID, endsID, fusionID = cls._parseID(iso_id)
    #     return super(IsoformInfo, cls).__new__(cls, iso_id, iso_name, gene_id)

    def __new__(cls, iso_id, iso_name, gene_id):
        if not iso_name: iso_name = iso_id
        return super(IsoformInfo, cls).__new__(cls, iso_id, iso_name, gene_id)

    def get_sub_ids(self):
        ##FL:12-1
        # ANNOT:31-
        all_ids = [int(x) for x in re.findall(r'\d+', self.iso_id)]
        sjID, endsID = all_ids[:2]
        addID = all_ids[2] if len(all_ids) > 2 else None
        tag = self.iso_id.split(':')[0]
        return tag, sjID, endsID, addID

    @classmethod
    def get_iso_id(cls, tag, sj_id, ends_id, add_id):
        id = f'{tag}:{sj_id}-{ends_id}'
        if add_id: id += f'-{add_id}'
        return id

    @property
    def full_transcript_name(self):
        return self.iso_id + '|' + self.iso_name

    @property
    def short_iso_id(self):
        tag, sjID, endsID, addID = self.get_sub_ids()
        return self.get_iso_id(tag, sjID, endsID, None)

    @property
    def has_iso_name(self):
        return self.iso_id != self.iso_name

    @property
    def fusion_index(self):
        fusion_id = self.get_sub_ids()[-1]
        if fusion_id and fusion_id.is_digit():
            return int(fusion_id)
        else: return None

    @property
    def is_fusion_locus(self):
        fusion_id = self.get_sub_ids()[-1]
        return fusion_id and fusion_id.is_digit()

    @classmethod
    def parse_from_ids(cls, gene_id, sj_id, ends_id=1, tag='FL', iso_name=None, add_id=None):
        iso_id = cls.get_iso_id(tag, sj_id, ends_id, add_id)
        return IsoformInfo(iso_id, iso_name, gene_id)

    # @classmethod
    # def parse_ids(cls, iso_id, iso_name, gene_id):
    #     # tag,sj_id, ends_id, add_id = get_sub_ids(iso_id)
    #
    #     return IsoformInfo(gene_id, sj_id, )

    @classmethod
    def parse_string(cls, name):
        iso_id, iso_name, gene_id = name.split('|')
        return IsoformInfo(iso_id, iso_name, gene_id)

    def __str__(self):
        return f'{self.iso_id}|{self.iso_name}|{self.gene_id}'
