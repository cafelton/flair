
from collections import namedtuple
import re



class IsoformInfo(namedtuple('IsoformInfo', ('iso_id','iso_name','gene_id'))):
    __slots__ = ()
    # def __new__(cls, gene_id, sj_id, ends_id=1, tag='FL', iso_name=None, add_id=None):
    #     #sjID, endsID, fusionID = cls._parseID(iso_id)
    #     return super(IsoformInfo, cls).__new__(cls, iso_id, iso_name, gene_id)

    def __new__(cls, iso_id, iso_name, gene_id):
        if not iso_name: iso_name = iso_id
        return super(IsoformInfo, cls).__new__(cls, iso_id, iso_name, gene_id)

    def get_sub_ids(self):
        p = re.compile(r'\d+')
        all_ids = [int(x) for x in p.findall(self.iso_id)]
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
