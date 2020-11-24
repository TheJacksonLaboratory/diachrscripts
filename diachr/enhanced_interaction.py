from typing import Tuple, List


class EnhancedInteraction:
    """
    chr_a: Chromosome of the first digest of the interaction (e.g. chr1)
    sta_a: Start coordinate of the first digest of the interaction (e.g. 123456)
    end_a: End coordinate of the first digest of the interaction (e.g. 234567)
    syms_a: Comma separated list of gene symbols associated with TSS on the first digest (e.g. HIST1H1E,HIST1H2BD)
    tsss_a: Comma separated list of TSS coordinates on the first digest (e.g. chr6:26158121:+,chr6:26156331:+)
    chr_b, sta_b, end_b, syms_b, tsss_b: Same as for the first digest but for the second digest
    enrichment_pair_tag: Two letter tag indicating the enrichment status of the two digests (e.g.  AA, AI, II)
    strand_pair_tag: Two symbol tag separated by '/' indicating the strands of TSS on first and second digest (e.g. -/-, +/+, -/+, +/-, -/d, ...)
    interaction_category: Interaction category wi
    """
    def __init__(self, chrA:str, staA:int, endA:int, symsA:str, tsssA:str,chrB:str, staB:int, endB:int, symsB:str, tsssB:str,
                 enr_pair_tag:str, strand_pair_tag:str, intera_cat:str, neg_log_p:float, read_pair_left:int, read_pair_right:int,i_dist:int):
        self._chr_a = chrA
        self._start_a = staA
        self._end_a = endA
        self._syms_a = symsA
        self._symbols_a = symsA.split(",")
        self._tsss_a = tsssA
        self._trans_starts_a = tsssA.split(",")
        self._chr_b = chrB
        self._start_b = staB
        self._end_b = endB
        self._syms_b = symsB
        self._symbols_b = symsB.split(",")
        self._tsss_b = tsssB
        self._trans_starts_b = tsssB.split(",")
        if len(enr_pair_tag) != 2:
            raise ValueError("Malformed enrichment pair tag (length is not 2): '%s'" % enr_pair_tag)
        self._enrichment_pair_tag = enr_pair_tag
        self._strand_pair_tag = strand_pair_tag
        self._interaction_category = intera_cat
        self._neg_log_p_value = neg_log_p
        self._read_pair_left = read_pair_left
        self._read_pair_right = read_pair_right
        self._i_dist = i_dist

    @property
    def chr_a(self) -> str:
        return self._chr_a

    @property
    def sta_a(self) -> int:
        return self._start_a

    @property
    def end_a(self) -> int:
        return self._end_a

    @property
    def chr_b(self) -> str:
        return self._chr_b

    @property
    def sta_b(self) -> int:
        return self._start_b

    @property
    def end_b(self) -> int:
        return self._end_b

    @property
    def syms_a(self) -> str:
        return self._syms_a

    @property
    def syms_b(self) -> str:
        return self._syms_b

    @property
    def tsss_a(self) -> str:
        return self._tsss_a

    @property
    def tsss_b(self) -> str:
        return self._tsss_b

    @property
    def enrichment_pair_tag(self) -> str:
        """
        A string such as AA (enriched/enriched) or AI (enriched, unenriched)
        """
        return self._enrichment_pair_tag

    def digest_a_enriched(self) -> bool:
        return self._enrichment_pair_tag[0] == 'A'

    def digest_b_enriched(self) -> bool:
        return self._enrichment_pair_tag[1] == 'A'

    @property
    def strand_pair_tag(self):
        return self._strand_pair_tag

    @property
    def interaction_category(self) -> str:
        return self._interaction_category

    @property
    def neg_log_p_value(self) -> float:
        return self._neg_log_p_value

    @property
    def read_pair_left(self) -> int:
        """
        these are the simple read pairs
        TODO -- deprecate this and replace by simple_count
        """
        return self._read_pair_left

    @property
    def simple_count(self) -> int:
        return self._read_pair_left


    @property
    def read_pair_right(self) -> int:
        """
        these are the twisted read pairs
        TODO -- deprecate this and replace by simple_count
        """
        return self._read_pair_right

    @property
    def twisted_count(self) -> int:
        return self._read_pair_right

    @property
    def i_dist(self) -> int:
        return self._i_dist

    @property
    def rp_total(self) -> int:
        return self._read_pair_right + self._read_pair_left

    @property
    def coordinate_key_a(self) -> str:
        return "%s\t%d\t%d" % (self._chr_a, self._start_a, self._end_a)

    @property
    def coordinate_key_b(self) -> str:
        return "%s\t%d\t%d" % (self._chr_b, self._start_b, self._end_b)

    def set_enrichment_pair_tag(self, tag):
        self._enrichment_pair_tag = tag

    def set_interaction_category(self, cat):
        """
        Set the interaction cateogory to one of 'DI', TODO details
        """
        self._interaction_category = cat

    def in_cis(self) -> bool:
        """
        return true iff both interactions are on the same chromosome.
        """
        return self._chr_a == self._chr_b

    def get_line(self):
        """
        recreate the output line
        chr1:156041392-156055887;chr1:156164786-156166517	108899	DIAI	LAMTOR2,UBQLN4;	0:10	AI	6.93	d/-1	chr1:156053798:-,chr1:156054782:+,chr1:156054726:+;
        """
        fields = []
        # field 0
        chrA = "%s:%d-%d" % (self._chr_a, self._start_a, self._end_a)
        chrB = "%s:%d-%d" % (self._chr_b, self._start_b, self._end_b)
        chroms = "%s;%s" % (chrA, chrB)
        fields.append(chroms)
        # field 1, interaction distance
        fields.append(str(self._i_dist))
        # field 2, interaction_category
        fields.append(self._interaction_category)
        # field 3, symbols
        syms = "%s;%s" % (self._syms_a, self._syms_b)
        fields.append(syms)
        # field 4, readpair counts
        rpairs = "%d:%d" % (self._read_pair_left, self._read_pair_right)
        fields.append(rpairs)
        # field 5, enrichment_pair_tag
        fields.append(self._enrichment_pair_tag)
        # field 6, neg_log_p_value
        fields.append(str(self._neg_log_p_value))
        # field 7, strand_pair_tag
        fields.append(self._strand_pair_tag)
        # field 8, tsss_a, tsss_b 
        tss = "%s;%s" % (self._tsss_a, self._tsss_b)
        fields.append(tss)
        # field 9, 
        return "\t".join(fields)







