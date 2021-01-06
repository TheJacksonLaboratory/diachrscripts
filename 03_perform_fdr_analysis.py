#!/usr/bin/env python

"""
This takes a path to a Diachromatic interaction file and determines a P-value threshold that corresponds to a
given FDR threshold.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""


import argparse
import math
from scipy.stats import binom
from collections import defaultdict
import numpy as np


from diachr.fdr_analysis import FdrAnalysis
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.binomial_model import BinomialModel


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine a P-value threshold that corresponds to a given FDR threshold.')
parser.add_argument('-p', '--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('-f', '--fdr-threshold', help='Use this switch to estimate a P-value cutoff that corresponds to a given FDR threshold.', default=0.25)
parser.add_argument('--p-val-c-min', help='Smallest P-value cutoff.', default=0.00025)
parser.add_argument('--p-val-c-max', help='Largest P-value cutoff.', default=0.05)
parser.add_argument('--p-val-step-size', help='P-value step size.', default=0.00025)

parser.add_argument('-m','--usemod', help="Use new module", dest='usemod', action='store_true')

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_file = args.interaction_file
fdr_threshold = float(args.fdr_threshold)
p_val_c_min = float(args.p_val_c_min)
p_val_c_max = float(args.p_val_c_max)
p_val_step_size = float(args.p_val_step_size)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --diachromatic-interaction-file: " + interaction_file)
print("\t[INFO] --fdr-threshold: " + str(fdr_threshold))
print("\t[INFO] --p-val-c-min: " + str(p_val_c_min))
print("\t[INFO] --p-val-c-max: " + str(p_val_c_max))
print("\t[INFO] --p-val-step-size: " + str(p_val_step_size))


if args.usemod:
    print("[INFO] Using FdrAnalysis module")
    fdra = FdrAnalysis(diachromatic_interaction_file=interaction_file, fdr_threshold=fdr_threshold, p_val_c_min=p_val_c_min,
                       p_val_c_max=p_val_c_max, step_size=p_val_step_size, prefix=out_prefix)
    fdra.calculate_fdr_threshild()
    exit(10)


### Define auxiliary functions
##############################

# Get object that returns P-values
p_values = BinomialModel()

def get_pvals_permuted_counts():
    """
    This function generates and returns a list of randomized binomial P-values.

    :return: List of randomized P-values
    """

    # Init list of randomized P-values
    pvals_permuted_counts = []

    # Iterate dictionary with numbers of interactions foreach read pair number n
    for n, i_num in N_DICT.items():

        # Generate random simple read pair counts for current n
        simple_count_list = list(binom.rvs(n, p = 0.5, size = i_num))

        for simple_count in simple_count_list:

            # Get twisted count
            twisted_count = n - simple_count

            # Get binomial P-value and add to list
            pv_nln = p_values.get_binomial_nnl_p_value(simple_count, twisted_count)
            pvals_permuted_counts.append(pv_nln)

    return pvals_permuted_counts


### Prepare variables, data structures and streams for output files
###################################################################

# Dictionary that stores the numbers of interactions with n read pairs
N_DICT = {}

# Total number of input interactions
n_interaction = 0

# List of P-values for observed interactions
p_val_o_list = []

# Prepare stream for output of results
txt_file_name_results = out_prefix + "_fdr_analysis_results.txt"
txt_file_stream_results = open(txt_file_name_results, 'wt')
txt_file_stream_results.write("OUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O" + "\n")


### Start execution
###################

# Set random seed
np.random.seed(42)

# Get list of Diachromatic interaction objects
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(interaction_file)
d_inter_list = interaction_set.interaction_list

# Iterate Dichromatic interaction objects
n_interaction = 0
for d_inter in d_inter_list:
    n_interaction += 1
    if n_interaction % 100000 == 0:
        print("\t[INFO] Processed " + str(n_interaction) + " interaction objects ...")

    # Get negative natural logarithm of P-value of a two-sided binomial test
    nnl_p_value = float("{:.2f}".format(p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)))
    p_val_o_list.append(nnl_p_value)

    # Add the sum of simple and twisted read pair counts to dictionary that will be used for randomization
    if d_inter.rp_total in N_DICT:
        N_DICT[d_inter.rp_total] += 1
    else:
        N_DICT[d_inter.rp_total] = 1

print("[INFO] Total number of interactions: {}".format(n_interaction))

print("[INFO] Getting sorted lists of P-values ...")

# Get list of randomized P-values
p_val_r_list = get_pvals_permuted_counts()

# Sort P-value lists
p_val_o_list.sort(key=float, reverse=True)
p_val_r_list.sort(key=float, reverse=True)

print("[INFO] Estimating FDR for increasing P-value thresholds ...")

# Estimate FDR for increasing P-value thresholds from original and randomized P-value lists
S_o = 0 # start index in sorted list
S_r = 0 # start index in sorted list
for pc in np.arange(p_val_c_min, p_val_c_max, p_val_step_size):

    # Transform P-value to negative natural logarithm (nnl)
    pc_nnl = - math.log(pc)

    # Get number of significant interactions for original P-values
    #S_o = sum(pc_nnl < pv_nnl for pv_nnl in p_val_o_list) # slow
    while S_o < len(p_val_o_list) and p_val_o_list[S_o] > pc_nnl:
        S_o += 1

    # Get number of significant interactions for randomized P-values
    #S_r = sum(pc_nnl < pv_nnl for pv_nnl in p_val_r_list) # slow
    while S_r < len(p_val_r_list) and p_val_r_list[S_r] > pc_nnl:
        S_r += 1

    # Estimate FDR
    fdr = S_r / S_o

    # Write results for this threshold to file
    txt_file_stream_results.write(out_prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(S_r) + "\t" + str(S_o) + "\n")

    # Print results for this threshold to screen
    print("\t" + out_prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(S_r) + "\t" + str(S_o))

    # Keep track of the largest P-value that satisfies the FDR threshold
    if fdr < fdr_threshold:
        fdr_last = fdr
        S_o_last = S_o
        S_r_last = S_r
        pc_last = pc

txt_file_stream_results.close()

# Print results for the largest P-value that satisfies the FDR threshold to the screen
print()
print("\tOUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O")
print("\t" + out_prefix + "\t" + str(fdr_last) + "\t" + str(pc_last) + "\t" + str(S_r_last) + "\t" + str(S_o_last))
