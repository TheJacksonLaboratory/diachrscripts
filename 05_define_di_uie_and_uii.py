#!/usr/bin/env python

"""
This script takes an enhanced interaction file, supplemented with digest associated gene symbols and TSS as well as
directionality P-values, and segments the interactions based on a P-value threshold into the following subsets:

   1. Indefinable interactions that cannot be significant with the chosen P-value threshold
   2. DI - Directed interactions
   3. UI - Undirected interactions
      i. UIE - Exclusive undirected interactions. Subset of UI. No digest interacts with a digest involved in DI
      ii. UII - Inclusive undirected interactions. Subset of UI. UI without UIE

The script implements a two pass approach:

   1. Iterate interactions in order to identify DI and involved digests.
   2. Iterate interactions a second time and use information from the first pass in order to identify UIE.

DI, UII and UIE interactions are written to a file in enhanced interaction format with the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_exclusive_ui.tsv.gz'

Column 3 of the created file contains DI or UIE. Indefinable interactions are discarded.

Finally, summary statistics about interaction and associated digest sets as well as connectivity within interaction
subsets are reported and written to a tab separated file with the name:

   '<OUT_PREFIX>_stats_exclusive_ui.tsv'

containing one header line and one line with the corresponding values:

   1. Prefix for output

   2. Total number of interactions

   3. Smallest n that can yield a significant interaction given the P-value threshold
   4. Number of indefinable interactions
   5. Percentage of indefinable interactions

   6. Number of DI
   7. Percentage of DI

   8. Number of UI
   9. Percentage of UI

  10. Number of UIE
  11. Percentage of UIE

  12. Number of UII
  13. Percentage of UII

  14. Number of digests involved in DI and UI
  15. Connectivity factor for DI and UI

  16. Number of digests involved in DI
  17. Connectivity factor for DI

  18. Number of digests involved in UI
  19. Connectivity factor for UI

  20. Number of digests involved in UIE
  21. Connectivity factor for UIE

  22. Number of digests involved in UII
  23. Connectivity factor for UII

  24. Number of digests involved in DI and UIE (sanity check: Must be 0!)
  25. Number of digests involved in DI and UII
  26. Percentage of digests involved in DI and UII
"""


import argparse
import math
import gzip
import diachrscripts_toolkit as dclass
from collections import defaultdict


### Parse command line
######################

parser = argparse.ArgumentParser(description='Identify directed and exclusive undirected interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--p-value-threshold', help='P-value threshold for directed interactions.', default=0.001)


parser.add_argument('-m','--usemod', help="Use new module", dest='usemod', action='store_true')

args = parser.parse_args()

out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
p_value_threshold = float(args.p_value_threshold)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --p-value-cutoff: " + str(p_value_threshold))

if args.usemod:
    print("[INFO] Using module")

    exit(1)

# Explore indefinable n
# import numpy as np
# for pc in np.arange(0.00001, 0.01, 0.00001):
#     n_1, pv_n_1 = dclass.find_indefinable_n(pc, False)
#     print(str(pc) + "\t" + str(pv_n_1) + "\t" + str(n_1))
#
# exit(1)

### Prepare variables, data structures and streams for output files
###################################################################

# P-value threshold
neg_log_p_val_thresh =  (-1) * math.log(p_value_threshold)

# Number of directed interactions and involved digests
dir_inter_num = 0
dir_inter_involved_dig_num = None

# Number of undirected interactions and involved digests
undir_inter_num = 0
undir_inter_involved_dig_num = None

# Number of indefinable interactions
indef_inter_num = 0

# Number of exclusive undirected interactions
undir_exc_inter_num = 0

# Number of exclusive undirected interactions
undir_inc_inter_num = 0

# Set containing all digests involved in directed interactions
dir_dig_set = set()

# Set containing all digests involved in undirected interactions
undir_dig_set = set()

# Set containing all digests involved in exclusive undirected interactions
undir_exc_dig_set = set()

# Set containing all digests involved in not exclusive undirected interactions
undir_inc_dig_set = set()

# Smallest n that required for significance given the P-value threshold
indef_n, indef_pv = dclass.find_indefinable_n(p_value_threshold, verbose = True)

# Prepare stream for output of filtered interactions annotated with respect to exclusive undirected interactions
enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_with_di_uii_and_uie.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')

# Prepare stream for output of statistics
tab_file_stats_output = out_prefix + "_stats_di_uii_and_uie.tsv"
tab_stream_stats_output = open(tab_file_stats_output, 'wt')


### 1st pass: Identifiy directed interactions (DI) based on P-value threshold
#############################################################################

print("[INFO] 1st pass: Identify directed interactions (DI) based on P-value threshold ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Add digest of directed interactions to digest set
        if neg_log_p_val_thresh < neg_log_p_value:
            dir_inter_num += 1
            dir_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            dir_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")


### 2nd pass: Identifiy exclusive undirected interactions and associated digests
################################################################################

print("[INFO] 2nd pass: Identifiy inclusive and exclusive undirected interactions (UII and UIE) ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Skip and count indefinable interactions
        if rp_total < indef_n:
            indef_inter_num += 1
            line = fp.readline()
            continue

        # Print line with directed interaction to file (field 3 will be 'DI')
        if neg_log_p_val_thresh < neg_log_p_value:
            enhanced_interaction_stream_output.write(dclass.set_interaction_category_in_enhanced_interaction_line(line, "DI") + "\n")
            line = fp.readline()
            continue

        # All remaining interactions must be undirected
        undir_inter_num += 1

        # Get digest coordinates used as keys for digest sets
        d1_coords = chr_a + "\t" + str(sta_a) + "\t" + str(end_a)
        d2_coords = chr_b + "\t" + str(sta_b) + "\t" + str(end_b)

        # Add digests of undirected interactions
        undir_dig_set.add(d1_coords)
        undir_dig_set.add(d2_coords)

        # Identify exclusive undirected interactions
        if d1_coords not in dir_dig_set and d2_coords not in dir_dig_set:
            undir_exc_dig_set.add(d1_coords)
            undir_exc_dig_set.add(d2_coords)
            # Print line with exclusive undirected interaction to file (field 3 will be 'UIE')
            enhanced_interaction_stream_output.write(dclass.set_interaction_category_in_enhanced_interaction_line(line, "UIE") + "\n")
            undir_exc_inter_num += 1
        else:
            undir_inc_dig_set.add(d1_coords)
            undir_inc_dig_set.add(d2_coords)
            # Print line with inclusive undirected interaction to file (field 3 will be 'UII')
            enhanced_interaction_stream_output.write(dclass.set_interaction_category_in_enhanced_interaction_line(line, "UII") + "\n")
            undir_inc_inter_num += 1

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")

enhanced_interaction_stream_output.close()


### Finish: Output statistics about interaction and digest sets
###############################################################

print("[INFO] Finish: Output statistics about interaction and digest sets ...")

print()

print("\t[INFO] Interaction statistics")

print()

print("\t\t[INFO] Total number of processed interactions: " + str(n_interaction_total))

print()

print("\t\t[INFO] Smallest n that can yield a significant interaction: " + str(indef_n))

print()

print("\t\t[INFO] Number of indefinable interactions (discarded): " + str(indef_inter_num))
indef_inter_percentage = "{0:.2f}".format(100 * indef_inter_num / n_interaction_total)
print("\t\t[INFO] Percentage of indefinable interactions: " + str(indef_inter_percentage) + "%")

print()

print("\t\t[INFO] Number of directed interactions (DI): " + str(dir_inter_num))
dir_inter_percentage = "{0:.2f}".format(100 * dir_inter_num / n_interaction_total)
print("\t\t[INFO] Percentage of directed interactions: " + str(dir_inter_percentage) + "%")

print()

print("\t\t[INFO] Number of undirected interactions (UI): " + str(undir_inter_num))
undir_inter_percentage = "{0:.2f}".format(100 * undir_inter_num / n_interaction_total)
print("\t\t[INFO] Percentage of undirected interactions: " + str(undir_inter_percentage) + "%")

print()

print("\t\t\t[INFO] Number of exclusive undirected interactions (UIE): " + str(undir_exc_inter_num))
undir_exc_inter_percentage = "{0:.2f}".format(100 * undir_exc_inter_num / n_interaction_total)
print("\t\t\t[INFO] Percentage of exclusive undirected interactions: " + str(undir_exc_inter_percentage) + "%")

print()

print("\t\t\t[INFO] Number of inclusive undirected interactions (UII): " + str(undir_inc_inter_num))
undir_inc_inter_percentage = "{0:.2f}".format(100 * undir_inc_inter_num / n_interaction_total)
print("\t\t\t[INFO] Percentage of inclusive undirected interactions: " + str(undir_inc_inter_percentage) + "%")

print()

print("\t[INFO] Digest and connectivity statistics")

print()

union_dir_undir_inter_involved_dig_num = len(dir_dig_set.union(undir_dig_set))
print("\t\t[INFO] Total number of digests involved in DI and UI: " + str(union_dir_undir_inter_involved_dig_num))
union_dir_undir_inter_connectivity = "{0:.2f}".format(1 - (union_dir_undir_inter_involved_dig_num / (2 * dir_inter_num + undir_inter_num)))
print("\t\t[INFO] Connectivity factor for DI and UI: " + str(union_dir_undir_inter_connectivity))

print()

dir_inter_involved_dig_num = len(dir_dig_set)
print("\t\t[INFO] Number of digests that are involved in DI: " + str(dir_inter_involved_dig_num))
dir_inter_connectivity = "{0:.2f}".format(1 - (dir_inter_involved_dig_num / (2 * dir_inter_num)))
print("\t\t[INFO] Connectivity factor for DI: " + str(dir_inter_connectivity))

print()

undir_inter_involved_dig_num = len(undir_dig_set)
print("\t\t[INFO] Number of digests that are involved in UI: " + str(undir_inter_involved_dig_num))
undir_inter_connectivity = "{0:.2f}".format(1 - (undir_inter_involved_dig_num / (2 * undir_inter_num)))
print("\t\t[INFO] Connectivity factor for UI: " + str(undir_inter_connectivity))

print()

undir_exc_inter_involved_dig_num = len(undir_exc_dig_set)
print("\t\t\t[INFO] Number of digests that are involved UIE: " + str(undir_exc_inter_involved_dig_num))
undir_exc_inter_connectivity = "{0:.2f}".format(1 - (undir_exc_inter_involved_dig_num / (2 * undir_exc_inter_num)))
print("\t\t\t[INFO] Connectivity factor for UIE: " + str(undir_exc_inter_connectivity))

print()

undir_inc_inter_involved_dig_num = len(undir_inc_dig_set)
print("\t\t\t[INFO] Number of digests that are involved in UII: " + str(undir_inc_inter_involved_dig_num))
undir_inc_inter_connectivity = "{0:.2f}".format(1 - (undir_inc_inter_involved_dig_num / (2 * undir_inc_inter_num)))
print("\t\t\t[INFO] Connectivity factor for UII: " + str(undir_inc_inter_connectivity))

print()

print("\t[INFO] Intersects between digests")

print()

dir_dig_undir_exc_dig_intersect_num = len(dir_dig_set.intersection(undir_exc_dig_set))
print("\t\t[INFO] Intersect of digests involved in DI and digest involved in UIE contains " + str(dir_dig_undir_exc_dig_intersect_num) + " digests. Sanity check: Must be 0!")
dir_dig_undir_inc_dig_intersect_num = len(dir_dig_set.intersection(undir_inc_dig_set))
print("\t\t[INFO] Intersect of digests involved in DI and digest involved in UII contains " + str(dir_dig_undir_inc_dig_intersect_num) + " digests.")
dir_dig_undir_inc_dig_intersect_percentage = "{0:.2f}".format(100 * dir_dig_undir_inc_dig_intersect_num / len(dir_dig_set))
print("\t\t[INFO] Percentage of digests involved in DI that are also involved in UII: " + str(dir_dig_undir_inc_dig_intersect_percentage) + "%.")

print()

tab_stream_stats_output.write(

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
tab_stream_stats_output.write(

    str(out_prefix) + "\t" +

    str(p_value_threshold) + "\t" +

    str(n_interaction_total) + "\t" +

    str(indef_n) + "\t" +
    str(indef_pv) + "\t" +
    str(indef_inter_num) + "\t" +
    str(indef_inter_percentage) + "\t" +

    str(dir_inter_num) + "\t" +
    str(dir_inter_percentage) + "\t" +

    str(undir_inter_num) + "\t" +
    str(undir_inter_percentage) + "\t" +

    str(undir_exc_inter_num) + "\t" +
    str(undir_exc_inter_percentage) + "\t" +

    str(undir_inc_inter_num) + "\t" +
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

tab_stream_stats_output.close()
