import os
import math
import gzip

from .enhanced_interaction_parser import EnhancedInteractionParser
from .diachr_util import find_minimum_powered_n



class InteractionClassifier:


    def __init__(self, interaction_file:str, prefix: str, p_threshold:float = 0.001):
        if not os.path.exists(interaction_file):
            raise FileNotFoundError("Could not find enhanced interaction file at ", interaction_file)
        self._enhanced_interaction_file = interaction_file
        self._prefix = prefix
        self._p_threshold = p_threshold
        print("[INFO] " + "Input parameters")
        print("\t[INFO] Analysis for: " + self._prefix)
        print("\t[INFO] Interaction file: " +  self._enhanced_interaction_file)
        print("\t[INFO] --p-value-cutoff: " + str(self._p_threshold))
        # P-value threshold
        self._neg_log_p_val_thresh =  (-1) * math.log(self._p_threshold)
        self._n_interaction_total = 0
        # Number of directed interactions and involved digests
        self._dir_inter_num = 0
        self._dir_inter_involved_dig_num = None
        # Number of undirected interactions and involved digests
        self._undir_inter_num = 0
        self._undir_inter_involved_dig_num = None
        # Number of underpower interactions
        self._indef_inter_num = 0
        # Number of exclusive undirected interactions
        self._undir_exc_inter_num = 0
        # Number of exclusive undirected interactions
        self._undir_inc_inter_num = 0
        # Set containing all digests involved in directed interactions
        self._dir_dig_set = set()
        # Set containing all digests involved in undirected interactions
        self._undir_dig_set = set()
        # Set containing all digests involved in exclusive undirected interactions
        self._undir_exc_dig_set = set()
        # Set containing all digests involved in not exclusive undirected interactions
        self._undir_inc_dig_set = set()
        # Smallest n that required for significance given the P-value threshold
        self._indef_n, _ = find_minimum_powered_n(self._p_threshold, verbose = True)
        self._einteractions = None # initiialize
       

    def identify_directed_interactions(self):
        """
        input the data and identify directed interactions based on the p-value threshold
        """
        print("[INFO] 1st pass: Identify directed interactions (DI) based on P-value threshold ...")

        print("\t[INFO] Iterating enhanced interaction file ...")
        parser = EnhancedInteractionParser(self._enhanced_interaction_file)
        self._einteractions = parser.parse()
        self._n_interaction_total = 0
        for ei in self._einteractions:
            self._n_interaction_total += 1
            if self._n_interaction_total % 1000000 == 0:
                print("\t\t[INFO]", self._n_interaction_total, "interactions processed ...")
            # Add digest of directed interactions to digest set
            if self._neg_log_p_val_thresh < ei.neg_log_p_value:
                self._dir_inter_num += 1
                self._dir_dig_set.add(ei.chr_a + "\t" + str(ei.sta_a) + "\t" + str(ei.end_a))
                self._dir_dig_set.add(ei.chr_b + "\t" + str(ei.sta_b) + "\t" + str(ei.end_b))
        print("\t[INFO] ... done; ingested %d extended interactions." % self._n_interaction_total)

    def identify_undirected_interactions(self):
        """
        Identify inclusive and exclusive undirected interactions (UII and UIE)
        """
        print("\t[INFO] Iterating enhanced interaction file ...")
        if self._einteractions is None or len(self._einteractions) == 0:
            raise ValueError("einteractions not initialized, did you run identify_directed_interactions?")
        n_interaction_total = 0
        for ei in self._einteractions :
            n_interaction_total += 1
            if n_interaction_total % 1000000 == 0:
                print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

            # Skip and count indefinable interactions
            if ei.rp_total < self._indef_n:
                self._indef_inter_num += 1
                continue
            # Print line with directed interaction to file (field 3 will be 'DI')
            if self._neg_log_p_val_thresh < ei.neg_log_p_value:
                ei.set_interaction_category("DI")
                continue
            # All remaining interactions must be undirected
            self._undir_inter_num += 1
            # Get digest coordinates used as keys for digest sets
            d1_coords = ei.coordinate_key_a
            d2_coords = ei.coordinate_key_b
            self._undir_dig_set.add(d1_coords)
            self._undir_dig_set.add(d2_coords)
            # Identify exclusive undirected interactions
            if d1_coords not in self._dir_dig_set and d2_coords not in self._dir_dig_set:
                self._undir_exc_dig_set.add(d1_coords)
                self._undir_exc_dig_set.add(d2_coords)
                ei.set_interaction_category("UIE")
                self._undir_exc_inter_num += 1
            else:
                self._undir_inc_dig_set.add(d1_coords)
                self._undir_inc_dig_set.add(d2_coords)
                ei.set_interaction_category("UII")
                self._undir_inc_inter_num += 1

    def write_to_file(self):
        # Prepare stream for output of filtered interactions annotated with respect to exclusive undirected interactions
        ei_file_output = self._prefix + "_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
        fh = gzip.open(ei_file_output, 'wt')
        for ei in self._einteractions:
            fh.write(ei.get_line() + "\n")
        fh.close()


    def output_statistics(self):
        """ 
        Output statistics about interaction and digest sets
        """
        print("[INFO] Finish: Output statistics about interaction and digest sets ...")
        print("\t[INFO] Interaction statistics")
        n_interaction_total = len(self._einteractions)
        print("\t\t[INFO] Total number of processed interactions: %d." % n_interaction_total)
        print("\t\t[INFO] Smallest n that can yield a significant interaction: %d." % self._indef_n)
        print("\t\t[INFO] Number of underpowered interactions (discarded): %d."  % self._indef_inter_num)
        indef_inter_percentage = "{0:.2f}".format(100 * self._indef_inter_num / n_interaction_total)
        print("\t\t[INFO] Percentage of indefinable interactions: " + str(indef_inter_percentage) + "%")
        print("\t\t[INFO] Number of directed interactions (DI): " + str(self._dir_inter_num))
        dir_inter_percentage = "{0:.2f}".format(100 * self._dir_inter_num / n_interaction_total)
        print("\t\t[INFO] Percentage of directed interactions: " + str(dir_inter_percentage) + "%")
        print("\t\t[INFO] Number of undirected interactions (UI): " + str(self._undir_inter_num))
        undir_inter_percentage = "{0:.2f}".format(100 * self._undir_inter_num / n_interaction_total)
        print("\t\t[INFO] Percentage of undirected interactions: " + str(undir_inter_percentage) + "%")
        print("\t\t\t[INFO] Number of exclusive undirected interactions (UIE): %d." % self._undir_exc_inter_num)
        undir_exc_inter_percentage = "{0:.2f}".format(100 * self._undir_exc_inter_num / n_interaction_total)
        print("\t\t\t[INFO] Percentage of exclusive undirected interactions: " + str(undir_exc_inter_percentage) + "%")
        print("\t\t\t[INFO] Number of inclusive undirected interactions (UII): %d."  % self._undir_inc_inter_num)
        undir_inc_inter_percentage = "{0:.2f}".format(100 * self._undir_inc_inter_num / n_interaction_total)
        print("\t\t\t[INFO] Percentage of inclusive undirected interactions: " + str(undir_inc_inter_percentage) + "%")
        print("\t[INFO] Digest and connectivity statistics")
        union_dir_undir_inter_involved_dig_num = len(self._dir_dig_set.union(self._undir_dig_set))
        print("\t\t[INFO] Total number of digests involved in DI and UI: %d." % union_dir_undir_inter_involved_dig_num)
        union_dir_undir_inter_connectivity = "{0:.2f}".format(1 - (union_dir_undir_inter_involved_dig_num / (2 * self._dir_inter_num + self._undir_inter_num)))
        print("\t\t[INFO] Connectivity factor for DI and UI: " + str(union_dir_undir_inter_connectivity))
        dir_inter_involved_dig_num = len(self._dir_dig_set)
        print("\t\t[INFO] Number of digests that are involved in DI: " + str(dir_inter_involved_dig_num))
        dir_inter_connectivity = "{0:.2f}".format(1 - (dir_inter_involved_dig_num / (2 * self._dir_inter_num)))
        print("\t\t[INFO] Connectivity factor for DI: " + str(dir_inter_connectivity))
        print()
        undir_inter_involved_dig_num = len(self._undir_dig_set)
        print("\t\t[INFO] Number of digests that are involved in UI: " + str(undir_inter_involved_dig_num))
        undir_inter_connectivity = "{0:.2f}".format(1 - (self._undir_inter_involved_dig_num / (2 * self._undir_inter_num)))
        print("\t\t[INFO] Connectivity factor for UI: " + str(undir_inter_connectivity))
        print()
        undir_exc_inter_involved_dig_num = len(self._undir_exc_dig_set)
        print("\t\t\t[INFO] Number of digests that are involved UIE: " + str(undir_exc_inter_involved_dig_num))
        undir_exc_inter_connectivity = "{0:.2f}".format(1 - (undir_exc_inter_involved_dig_num / (2 * self._undir_exc_inter_num)))
        print("\t\t\t[INFO] Connectivity factor for UIE: " + str(undir_exc_inter_connectivity))
        print()
        undir_inc_inter_involved_dig_num = len(self._undir_inc_dig_set)
        print("\t\t\t[INFO] Number of digests that are involved in UII: " + str(undir_inc_inter_involved_dig_num))
        undir_inc_inter_connectivity = "{0:.2f}".format(1 - (undir_inc_inter_involved_dig_num / (2 * self._undir_inc_inter_num)))
        print("\t\t\t[INFO] Connectivity factor for UII: " + str(undir_inc_inter_connectivity))
        print()
        print("\t[INFO] Intersects between digests")
        print()
        dir_dig_undir_exc_dig_intersect_num = len(self._dir_dig_set.intersection(self._undir_exc_dig_set))
        print("\t\t[INFO] Intersect of digests involved in DI and digest involved in UIE contains " + str(dir_dig_undir_exc_dig_intersect_num) + " digests. Sanity check: Must be 0!")
        dir_dig_undir_inc_dig_intersect_num = len(self._dir_dig_set.intersection(self._undir_inc_dig_set))
        print("\t\t[INFO] Intersect of digests involved in DI and digest involved in UII contains " + str(dir_dig_undir_inc_dig_intersect_num) + " digests.")
        dir_dig_undir_inc_dig_intersect_percentage = "{0:.2f}".format(100 * dir_dig_undir_inc_dig_intersect_num / len(self._dir_dig_set))
        print("\t\t[INFO] Percentage of digests involved in DI that are also involved in UII: " + str(dir_dig_undir_inc_dig_intersect_percentage) + "%.")
        print()
        tab_file_stats_output = self._prefix + "_stats_di_uii_and_uie.tsv"
        fh = open(tab_file_stats_output, 'wt')
        fh.write(
    "out_prefix" + "\t" +                               # Prefix for output

    "p_value_threshold" + "\t" +                        # P-value threshold used to define directed interactions

    "n_interaction_total" + "\t" +                      # Total number of interactions

    "indef_n" + "\t" +                                  # Smallest n that can yield a significant interaction given the P-value threshold
    "indef_pv" + "\t" +                                 # P-value that corresponds to indefinable n
    "indef_inter_num" + "\t" +                          # Number of indefinable interactions
    "indef_inter_percentage" + "\t" +                   # Percentage of indefinable interactions

    "dir_inter_num" + "\t" +                            # Number of DI
    "dir_inter_percentage" + "\t" +                     # Percentage of DI

    "undir_inter_num" + "\t" +                          # Number of UI
    "undir_inter_percentage" + "\t" +                   # Percentage of UI

    "undir_exc_inter_num" + "\t" +                      # Number of UIE
    "undir_exc_inter_percentage" + "\t" +               # Percentage of UIE

    "undir_inc_inter_num" + "\t" +                      # Number of UII
    "undir_inc_inter_percentage" + "\t" +               # Percentage of UII

    "union_dir_undir_inter_involved_dig_num" + "\t" +   # Number of digests involved in DI and UI
    "union_dir_undir_inter_connectivity" + "\t" +       # Connectivity factor for DI and UI

    "dir_inter_involved_dig_num" + "\t" +               # Number of digests involved in DI
    "dir_inter_connectivity" + "\t" +                   # Connectivity factor for DI

    "undir_inter_involved_dig_num" + "\t" +             # Number of digests involved in UI
    "undir_inter_connectivity" + "\t" +                 # Connectivity factor for UI

    "undir_exc_inter_involved_dig_num" + "\t" +         # Number of digests involved in UIE
    "undir_exc_inter_connectivity" + "\t" +             # Connectivity factor for UIE

    "undir_inc_inter_involved_dig_num" + "\t" +         # Number of digests involved in UII
    "undir_inc_inter_connectivity" + "\t" +             # Connectivity factor for UII

    "dir_dig_undir_exc_dig_intersect_num" + "\t" +      # Number of digests involved in DI and UIE (sanity check: Must be 0!)
    "dir_dig_undir_inc_dig_intersect_num" + "\t" +      # Number of digests involved in DI and UII
    "dir_dig_undir_inc_dig_intersect_percentage" +      # Percentage of digests involved in DI and UII

    "\n"
    )
        fh.write(

    str(self._prefix) + "\t" +

    str(self._p_threshold) + "\t" +

    str(n_interaction_total) + "\t" +

    str(self._indef_n) + "\t" +
    str(self._indef_pv) + "\t" +
    str(self._indef_inter_num) + "\t" +
    str(self._indef_inter_percentage) + "\t" +

    str(self._dir_inter_num) + "\t" +
    str(self._dir_inter_percentage) + "\t" +

    str(self._undir_inter_num) + "\t" +
    str(undir_inter_percentage) + "\t" +

    str(self._undir_exc_inter_num) + "\t" +
    str(undir_exc_inter_percentage) + "\t" +

    str(self._undir_inc_inter_num) + "\t" +
    str(undir_inc_inter_percentage) + "\t" +

    str(union_dir_undir_inter_involved_dig_num) + "\t" +
    str(union_dir_undir_inter_connectivity) + "\t" +

    str(dir_inter_involved_dig_num) + "\t" +
    str(dir_inter_connectivity) + "\t" +

    str(undir_inter_involved_dig_num) + "\t" +
    str(undir_inter_connectivity) + "\t" +

    str(undir_exc_inter_involved_dig_num) + "\t" +
    str(undir_exc_inter_connectivity) + "\t" +

    str(undir_inc_inter_involved_dig_num) + "\t" +
    str(undir_inc_inter_connectivity) + "\t" +

    str(dir_dig_undir_exc_dig_intersect_num) + "\t" +
    str(dir_dig_undir_inc_dig_intersect_num) + "\t" +
    str(dir_dig_undir_inc_dig_intersect_percentage) +

    "\n"
    )
        fh.close()


        