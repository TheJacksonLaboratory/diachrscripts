import gzip
from scipy.stats import binom
import numpy as np
from decimal import *
import warnings
warnings.filterwarnings("error")
from decimal import Decimal, getcontext
from math import log, log10
import sys
import warnings

#######################################################################################################################

class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest

    # Attributes
    chromosome = None
    start = None
    end = None
    active = False

    # Initializer

    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    # Methods
    def set_active(self):
        self.active = True

    def is_active(self):
        return self.active

    def get_chromosome(self):
        return self.chromosome

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

#######################################################################################################################

class Interaction:

    # Class to represent a genomic interaction between two restriction digests

    # Attributes

    digest_distance = None  # Distance between the two interacting digests
    n_simple = 0            # Number of simple read pairs
    n_twisted = 0           # Number of twisted read pairs
    cis = None              # Intrachromosomal interaction
    type = None             # Either simple, twisted or undirected
    digest_1 = None
    digest_2 = None

    # Initializer
    # def __init__(self, digest_1, digest_2, n_simple, n_twisted):
    #     self.digest_1 = digest_1
    #     self.digest_2 = digest_2
    #     if digest_1.get_chromosome() == digest_2.get_chromosome():
    #         self.cis = True
    #     else:
    #         self.cis = False
    #     self.n_simple = n_simple
    #     self.n_twisted = n_twisted

    def __init__(self, diachromatic_interaction_line):

        fields = diachromatic_interaction_line.split("\t")
        self.digest_1 = Digest(fields[0], int(fields[1]), int(fields[2]))
        if fields[3] == 'A':
            self.digest_1.set_active()
        self.digest_2 = Digest(fields[4], int(fields[5]), int(fields[6]))
        if fields[7] == 'A':
            self.digest_2.set_active()

        if ':' in fields[8]:  # line regular Diachromatic file
            fields2 = fields[8].split(":")
            self.n_simple = int(fields2[0])
            self.n_twisted = int(fields2[1])
            i_type = "TBD"  # to be determined using binomial P-value
        else: # line from LRT script
            self.n_simple = fields[8]
            self.n_twisted = fields[9]
            i_type = fields[13].rstrip() # has been determined already using LRT
        self.set_interaction_type(i_type)

    # Methods

    def get_digest_distance(self):
        if self.digest_distance == None:
            self.digest_distance = self.digest_2.get_start() - self.digest_1.get_end() # distance between digest ends
        return self.digest_distance

    def is_cis(self):
        return self.cis

    def get_digest_status_pair_flag(self):
        if self.digest_1.is_active():
            state_1 = 'A'
        else:
            state_1 = 'I'
        if self.digest_2.is_active():
            state_2 = 'A'
        else:
            state_2 = 'I'
        category = state_1 + state_2
        return sorted(category)[0]+sorted(category)[1]

    def get_binomial_p_value(self): # ToDo: check this function for small simple and twisted read counts in R
        if self.n_simple < self.n_twisted:
            return 1 - binom.cdf(self.n_twisted-1, self.n_simple + self.n_twisted, 0.5)
        else:
            return 1 - binom.cdf(self.n_simple-1, self.n_simple + self.n_twisted, 0.5)

    def set_interaction_type(self, type):
        """
        This function determines whether this interaction is significantly directed, i.e. simple or twisted, using a
        binomial distribution with p=0.5 and returns a 'S' for simple and 'T' for twisted.
        Otherwise, the function returns an 'U' for undirected or a 'NA', if the sum of twisted and simple is smaller
        than five because interactions with less read pairs cannot be significant using the binomial distribution.
        If the interaction input file was generated by the LRT script, the type of this interaction is set to the
        given type. Note that interactions from LRT cannot have the type 'NA'.

        :return: 'S' (simple),'T' (twisted), 'U' (undirected), 'NA' (not available)
        """
        if(type != 'S' and type != 'T' and type != 'U' and type != "TBD"):
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'TBD' but was " + type + ".")
        else:
            if type == "TBD":
                if self.n_twisted + self.n_simple < 5:
                    type = 'NA' # an interaction with less than 5 read pairs cannot be significant with a binomial P-value threshold of 0.05
                elif self.get_binomial_p_value() <= 0.05:
                    if self.n_twisted < self.n_simple:
                        type = 'S'
                    else:
                        type = 'T'
                else:
                    type = 'U'
            self.type = type

    def get_interaction_type(self):
        return self.type

    def get_first_digest(self):
        return self.digest_1

    def get_second_digest(self):
        return self.digest_2

    def is_cis_long_range(self, d_dist):
        return self.digest_1.get_chromosome() == self.digest_1.get_chromosome() and 10000 < self.get_digest_distance()



#######################################################################################################################

class TSSInfo:

    # Class that represents additional information about TSS beyond coordinates

    # Attributes

    strand = None
    gene_id = None
    gene_symbol = None
    fpkm = None
    expression_category = None

    # Initializer

    def __init__(self, strand, gene_id, gene_symbol):
        self.strand = strand
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol

    # Methods

    def set_fpkm(self, fpkm):
        self.fpkm = fpkm

    def set_expression_category(self, category):
        self.expression_category = category

#######################################################################################################################

class TSSCoordinate:

    # Class that represents a genomic coordinate that corresponds to a TSS for one or more genes

    # Attributes

    chromosome = None
    position = None
    tss_info_dict = None   # stores information for all TSS at this position

    # Initializer

    def __init__(self, chromosome, position, strand, gene_id, gene_symbol):
        self.chromosome = chromosome
        self.position = position
        self.tss_info_dict = {}
        self.tss_info_dict[gene_id] = TSSInfo(strand, gene_id, gene_symbol)

    # Methods

    def append_TSSInfo(self, strand, gene_id, gene_symbol):
        self.tss_info_dict[gene_id] = TSSInfo(strand, gene_id, gene_symbol)

    def get_num_of_genes(self):
        return len(self.tss_info_dict)

    def has_multiple_genes(self):
        return 1 < self.get_num_of_genes()

    def has_tss_on_both_strands(self):
        current_strand = None
        for info in self.tss_info_dict.values():
            if current_strand == None:
                current_strand = info.strand
            elif current_strand != info.strand:
                return True
        return False

    def print_tss_info(self):
        if self.has_tss_on_both_strands():
            print("[INFO] Found TSS on both strands for " + str(self.get_num_of_genes()) + " genes at " + self.chromosome + ":" + str(self.position))
        else:
            print("[INFO] Found TSS for " + str(self.get_num_of_genes()) + " genes at " + self.chromosome + ":" + str(self.position))
        for info in self.tss_info_dict.values():
            print("\tgene_id:" + info.gene_id + "|" + "gene_symbol:" + info.gene_symbol + "|" + "strand:" + info.strand + "|" + "fpkm:" + str(info.fpkm))

#######################################################################################################################

class TSSCoordinateMap:

    # Class that represents all coordinates that have one or more TSS

    # Attributes

    tss_coord_dict = None
    annotation_format = None
    annotation_file_path = None

    fpkm_n_zero = None
    fpkm_upper_first_q = None
    fpkm_upper_second_q = None
    fpkm_upper_third_q = None

    # Initializer

    def __init__(self, annotation_file_path, format):
        self.annotation_format = format
        self.annotation_file_path = annotation_file_path
        if format == "refGene":
            self.tss_coord_dict = {}
            self.parse_ref_gene_file()

    # Methods

    def parse_ref_gene_file(self):

        print("[INFO] Parsing " + self.annotation_format + " annotation file: " + self.annotation_file_path + " ...")

        # open file
        with gzip.open(self.annotation_file_path, mode='rt') as fp:

            # get first line
            line = fp.readline()

            # iterate file
            while line:

                # parse line
                values = line.split('\t')
                gene_id = values[12] # for refGene we use gene symbols as gene ID because refSeq gene IDs correspond to transcripts
                chromosome = values[2]
                strand = values[3]
                if strand == "+":
                    position = values[4]
                else:
                    position = values[5]
                gene_symbol = values[12]

                # construct key
                chr_pos_key = chromosome + ":" + position

                # check whether this coordinate has been seen already
                if chr_pos_key in self.tss_coord_dict.keys():
                    # append TSS info to existing coordinate
                    self.tss_coord_dict[chr_pos_key].append_TSSInfo(strand, gene_id, gene_symbol)
                else:
                    # create new coordinate
                    self.tss_coord_dict[chr_pos_key] = TSSCoordinate(chromosome, position, strand, gene_id, gene_symbol)

                # get next line
                line = fp.readline()

        # close file
        fp.close()

    def analyze_coordinates_and_print_report(self):

        # declare count variables
        n_coord_total = 0
        n_genes_total = 0
        n_coord_multiple_genes = 0
        n_coord_multiple_genes_different_strands = 0

        # iterate coordinates and count
        for coord in self.tss_coord_dict.values():
            n_coord_total += 1
            n_genes_total += coord.get_num_of_genes()
            if 1 < coord.get_num_of_genes():
                n_coord_multiple_genes += 1
                if coord.has_tss_on_both_strands():
                    n_coord_multiple_genes_different_strands += 1

        # print report
        print("\t[INFO] Found " + str(n_coord_total) + " coordinates that host TSS for " + str(n_genes_total) + " genes.")
        print("\t[INFO] Found " + str(n_coord_multiple_genes) + " coordinates that host TSS for more than one gene.")
        print("\t[INFO] Found " + str(n_coord_multiple_genes_different_strands) + " coordinates that host TSS on different strands.")


    def get_num_of_coords_with_miltiple_genes(self):
        num = 0
        for coord in self.tss_coord_dict.values():
            if 1 < coord.get_num_of_genes():
                num += 1
        return num

    def parse_cuffdiff_genes_fpkm_tracking_file(self, cuffdiff_genes_fpkm_tracking_file):

        print("[INFO] Parsing genes.fpkm_tracking file of cuffdiff: " + cuffdiff_genes_fpkm_tracking_file + " ...")

        # read FPKM values to temporary dictionary

        tmp_fpkm_dict = {}

        with open(cuffdiff_genes_fpkm_tracking_file) as fp:

            line = fp.readline() # skip first line
            line = fp.readline()
            while line:
                line_fields = line.split("\t")
                gene_id = line_fields[3]
                fpkm = float(line_fields[9])
                tmp_fpkm_dict[gene_id] = fpkm
                line = fp.readline()

        # determine number of zero FPKMs and quartiles
        tmp_fpkm_list = tmp_fpkm_dict.values()
        self.fpkm_n_zero = tmp_fpkm_list.count(0)
        tmp_fpkm_list = filter(lambda a: a != 0, tmp_fpkm_list) # remove zero FPKMs
        self.fpkm_upper_first_q = round(np.quantile(tmp_fpkm_list, .25),2)
        self.fpkm_upper_second_q = round(np.quantile(tmp_fpkm_list, .50),2)
        self.fpkm_upper_third_q = round(np.quantile(tmp_fpkm_list, .75),2)

        # add FPKM values to TSSInfo
        genes_without_fpkm = []
        for coord in self.tss_coord_dict.values():
            for tss_info in coord.tss_info_dict.values():
                if tss_info.gene_id in tmp_fpkm_dict:
                    tss_info.set_fpkm(tmp_fpkm_dict[tss_info.gene_id])
                else:
                    genes_without_fpkm.append(tss_info.gene_id)

        if 0 < len(genes_without_fpkm):
            print("\t[WARNING] No FPKM for the following " + str(len(genes_without_fpkm)) + " genes: " + " ".join(genes_without_fpkm))

        # print info about zero FPKM and quartiles
        print ("\t[INFO] There were", self.fpkm_n_zero, "genes with zero FPKM.")
        print ("\t[INFO] Upper FPKM limit of the first quartile:", self.fpkm_upper_first_q)
        print ("\t[INFO] Upper FPKM limit of the second quartile:", self.fpkm_upper_second_q)
        print ("\t[INFO] Upper FPKM limit of the third quartile:", self.fpkm_upper_third_q)

    def set_expression_categories(self, q_limit):

        print("[INFO] Setting expression level categories...")

        if q_limit == 1:
            thresh = self.fpkm_upper_first_q
            print("\t[INFO] Will use upper limit of the first quartile (" + str(self.fpkm_upper_first_q)  + ") as FPKM threshold for inactive genes.")
        elif q_limit == 2:
            thresh = self.fpkm_upper_second_q
            print("\t[INFO] Will use upper limit of the second quartile (" + str(self.fpkm_upper_second_q) + ") as FPKM threshold for inactive genes.")
        elif q_limit == 3:
            thresh = self.fpkm_upper_third_q
            print("\t[INFO] Will use upper limit of the third quartile (" + str(self.fpkm_upper_third_q) + ") as FPKM threshold for inactive genes.")
        else:
            raise Exception("\t[FATAL] Threshold must be an upper limits of the first three quartiles. Should be either 1, 2 or 3 but was " + str(q_limit) + ".")

        none_fpkm = 0
        inactive_category = 0
        active_category = 0

        for coord in self.tss_coord_dict.values():
            for tss_info in coord.tss_info_dict.values():
                if tss_info.fpkm == None:
                    tss_info.set_expression_category(None)
                    none_fpkm += 1
                elif tss_info.fpkm == 0 or tss_info.fpkm < thresh:
                    tss_info.set_expression_category(0)
                    inactive_category +=1
                else:
                    tss_info.set_expression_category(1)
                    active_category += 1

        print("\t[WARNING] There are " + str(none_fpkm) + " without FPKM. Category was set to 'None'.")
        print("\t[INFO] " + str(inactive_category) + " genes were categorized as inactive.")
        print("\t[INFO] " + str(active_category) + " genes were categorized as active.")

    def has_key(self, key):
        return key in self.tss_coord_dict

    def get_coord_category(self, key):

        coord_expression_categories = [] # note that there may be multiple TSS a given coordinate

        if key in self.tss_coord_dict:
            for tss_info in self.tss_coord_dict[key].tss_info_dict.values():
                coord_expression_categories.append(tss_info.expression_category)
        else:
            return -1 # no TSS annotated at this position

        if 1 in coord_expression_categories:
            return 1 # one active coord is enough to make the digest active
        elif 0 in coord_expression_categories:
            return 0
        else:
            return tss_info.gene_id # no FPKM for annotated TSS

    def get_coord_strand(self, coord_key):

        coord_strands = []  # note that there may be multiple TSS a given coordinate

        if coord_key in self.tss_coord_dict:
            for tss_info in self.tss_coord_dict[coord_key].tss_info_dict.values():
                coord_strands.append(tss_info.strand)
        else:
            return -1

        if '-' in coord_strands and '+' in coord_strands:
            return 'd'
        elif '-' in coord_strands:
            return '-'
        elif '+' in coord_strands:
            return '+'
        else:
            raise Exception("[FATAL] Invalid strand symbol. Should be either '-' or '+' but was " + str(coord_strands) + ".")



#######################################################################################################################

class PairKeyDict:

    states = []

    def __init__(self, states):
        self.states = states
        self.pair_dict = {}
        for i in states:
            for j in states:
                key = i + "/" + j
                self.pair_dict[key] = 0

    def print_pair_key_dict(self):
        for i in self.states:
            for j in self.states:
                key = i + "/" + j
                print(key + "\t" + str(self.pair_dict[key]))


#######################################################################################################################

def get_fraction(value, values):
    """
    :param value: A single value.
    :param values: A list of values that contains 'value'.
    :return: A string representing a formatted number with two decimal places.
    """
    if 0 < sum(values):
        return str(Decimal(float(1.0*value/sum(values))).quantize(Decimal(10) ** -2))
    else:
        return str(0.00)

def get_string_formatted_fraction(numerator, denominator):
    if 0 < denominator or 0 > denominator:
        return str(Decimal(float(1.0 * numerator / denominator)).quantize(Decimal(10) ** -2))
    else:
        return str(0.00)

def get_x_inter_binomial_log10_p_value(k, n, x):
    """
    Calculate a P-value using a binomial distribution for observing k out of n interactions with a total number of x
    read pairs that have either zero simple or zero twisted read pairs, given an equal probability of 0.5 for simple
    and twisted read pairs.

    For instance, a single interaction with zero simple and two twisted read pairs corresponds to sequence of two
    independent Bernoulli trials and in each trial the probability for drawing a twisted read pair is 0.5.
    Therefore, the probability to draw two twisted read pairs in a sequence of two trials is 0.5^2,
    and the probability to draw three twisted read pairs in a sequence of three trials is 0.5^3 and so on,
    i.e. p=0.5^x and since we do not distinguish between zero simple and zero twisted interactions, we have 2*0.5^x.

    This function uses a binomial distribution binom(n,k,p) in order to calculate the P-value.

    :param k: Number of interactions with either zero simple or zero twisted read pairs
    :param n: Total number of interactions with a fixed number of x read pairs
    :param x: Total number of read pairs in n interactions
    :return: Negative decadic logarithm of the binomial P-value and error code
    """
    if n==0: # no interaction, no P-value
        return "NA", 1
    #print("-binom.logsf(" + str(k) + ", " + str(n) + ", 2*(0.5 ** " + str(x) + "))/log(10)")
    try:
        return -binom.logsf(k, n, 2*(0.5 ** x))/log(10), 0
    except RuntimeWarning: # underflow
       return -log(sys.float_info.min * sys.float_info.epsilon)/log(10), 2 # return smallest possible float