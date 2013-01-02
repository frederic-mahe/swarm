#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2013/01/01"
__version__ = "$Revision: 3.0"

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser

#**********************************************************************#
#                                                                      #
#                            Functions                                 #
#                                                                      #
#**********************************************************************#

def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Find clusters of nearly identical sequences in giant
    barcoding projects."""
    
    parser = OptionParser(usage="usage: %prog --input_file filename --threshold integer",
        description = desc,
        version = "%prog version 1.0")

    parser.add_option("-i", "--input_file",
        metavar = "<FILENAME>",
        action = "store",
        dest = "input_file",
        help = "set <FILENAME> as input file. Ungapped lowercase fasta file (acgt only).")

    parser.add_option("-t", "--threshold",
        metavar = "<THRESHOLD>",
        action = "store",
        dest = "threshold",
        type = "int",
        default = 7,
        help = "set <THRESHOLD> for the swarm building.")

    (options, args) = parser.parse_args()
    
    return options.input_file, options.threshold


def needleman_wunsch(seqA, seqB):
    """
    Global pairwise alignment algorithm with a linear gap
    penalty. Code adapted from Wikipedia's Needleman-Wunsch
    pseudo-code
    (https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
    """
    # Substitution matrix (S), mismatch penalty (p) and gap opening
    # penalty (d) used in Swarm (transformed +5/-4/+12/-4 model).
    d = 7
    S = {'a': {'a': 0, 'g': 3, 'c': 3, 't': 3},
         'g': {'a': 3, 'g': 0, 'c': 3, 't': 3},
         'c': {'a': 3, 'g': 3, 'c': 0, 't': 3},
         't': {'a': 3, 'g': 3, 'c': 3, 't': 0}}

    # Initialize array F with zeroes (and set outer lines and columns
    # with cumulative penalty value)
    A = seqA
    B = seqB
    I = range(len(A) + 1)
    J = range(len(B) + 1)
    F = [[0 for j in J] for i in I]
    for i in I:
        F[i][0] = d * i
    for j in J:
        F[0][j] = d * j
	
    ## Scoring
    #
    # Replace S mapping with a simple if statement (faster)? If the
    # matrix is symetrical (delete == insert penalty), there is no
    # need to compute both? In think I am wrong.
    for i in I[1:]:
	for j in J[1:]:
            match = F[i-1][j-1] + S[A[i-1]][B[j-1]]
            delete = F[i-1][j] + d
            insert = F[i][j-1] + d
            # Use max() if F is a similarity matrix 
            F[i][j] = min(match, insert, delete)

    score = F[-1][-1]
    return score


#**********************************************************************#
#                                                                      #
#                              Body                                    #
#                                                                      #
#**********************************************************************#

if __name__ == '__main__':
    
    # Parse command line options.
    input_file, threshold = option_parse()

    # Substitution matrix (S), mismatch penalty (p) and gap opening
    # penalty (d) used in Swarm (transformed +5/-4/+12/-4 model).
    d = 7
    p = 3

    # Scoring system and length differences. Swarm model is based on a
    # mismatch penalty of 3 (p) and a gap opening penalty of 7 (d)
    # if 0 < threshold < 7: max_length_difference = 0
    # elif 7 <= threshold < 10: max_length_difference = 1
    # elif 10 <= threshold < 13: max_length_difference = 2
    # elif 13 <= threshold < 16: max_length_difference = 3
    # elif 16 <= threshold < 19: max_length_difference = 4
    # elif 19 <= threshold < 21: max_length_difference = 5
    # else: max_length_difference = 100
    # That can be expressed like that:
    if threshold < d:
        max_length_difference = 0
    else:
        max_length_difference = ((threshold - d) / p) + 1

    # Scoring system and max number of mismatches
    # 0: identical sequences
    # 3: 1 mismatch
    # 6: 2 mismatches
    # 7: 1 gap of length 1
    # 9: 3 mismatches
    # 10: 1 gap of length 2, or 1 gap of length 1 + 1 mismatch
    # 12: 4 mismatches
    # 13: 1 gap of length 3, or 1 gap of length 1 + 2 mismatches, or 1 gap of length 2 + 1 mismatch
    # 14: 2 gaps of length 1 each
    # 15: 5 mismatches
    # 16: 1 gap of length 4, or 1 gap of length 1 + 3 mismatches, or 1 gap of length 2 + 2 mismatches, or 1 gap of length 3 + 1 mismatch 
    max_number_of_mismatches = threshold / p

    # Build a list of sequences and count nucleotide occurences
    input_format = "fasta"
    nucleotides = ("a","c","g","t")
    records_list = list()
    status = list()
    with open(input_file, "rU") as input_file:
        records = SeqIO.parse(input_file, input_format)
        records_list = [(record.id.split("_")[0],
                         record.id.split("_")[1],
                         len(record.seq),
                         str(record.seq),
                         [record.seq.lower().count(n) for n in nucleotides])
                        for record in records]
        status = [True] * len(records_list)
        
    # Start swarming
    while True:
        try:
            # Search the next master seed
            i = status.index(True)
        except ValueError:
            # All sequences have been treated
            break

        # Master seed
        swarm = [records_list[i][0]]
        status[i] = False

        # List remaining non-swarmed sequences (boolean test)
        comparisons = [j for j, state in enumerate(status) if state]
        
        # Deal with the last remaining item
        if not comparisons:
            print(" ".join(swarm), file=sys.stdout)
            break

        # Candidates list is initialized for each major seed
        candidates = [(j, needleman_wunsch(records_list[i][3], records_list[j][3])) for j in comparisons]

        # Parse candidates and select sons
        firstseeds = [j for j, d in candidates if d <= threshold]
        swarm.extend([records_list[j][0] for j in firstseeds])
        for j in firstseeds: status[j] = False

        # Loop other the first subseeds (if they exists)
        all_subseeds = list()
        if firstseeds:
            all_subseeds.append(firstseeds)
            # Loop other the list of firstseeds lists
            for k, subseeds in enumerate(all_subseeds):
                nextseeds = list()
                frontier = (k + 2) * threshold
                # I could replace that for loop with a list
                # comprehension, but the risk is to perform un-needed
                # amplicon comparisons, as the status list is not
                # updated regularly.
                for l in subseeds:
                    # Candidates already have been filtered for "left
                    # hand" reads. Do not treat already assigned
                    # sequences. Do not compare sequences we know to
                    # be too distant from the seed. Do not compare
                    # sequences with a length difference greater than
                    # the threshold. Do not compare sequences if their
                    # nucleotide profiles are too divergent.
                    hits = [j for j, d in candidates
                            if status[j]
                            and d <= frontier
                            and abs(cmp(records_list[l][2], records_list[j][2])) <= max_length_difference
                            and sum([abs(cmp(couple[0], couple[1])) for couple in zip(records_list[l][4], records_list[j][4])]) <= 2 * max_number_of_mismatches - abs(cmp(records_list[l][2], records_list[j][2]))
                            and needleman_wunsch(records_list[l][3], records_list[j][3]) <= threshold]
                    if hits:
                        nextseeds.extend(hits)
                        swarm.extend([records_list[j][0] for j in hits])
                        for j in hits: status[j] = False
                # Stop condition
                if nextseeds:
                    all_subseeds.append(nextseeds)
                else:
                    # No new subseeds
                    print(" ".join(swarm), file=sys.stdout)
                    break
        else:
            # Deal with isolated sequences
            print(" ".join(swarm), file=sys.stdout)
        
sys.exit(0)

# Profiling
# ---------
#
# python -m cProfile swarm.py -i file.fas > tmp
#
# Expected results
# ----------------
#
# time python swarm.py -i AF091148.fas | awk '{print NF}' | sort -nr | uniq -c
#       1 878
#       1 423
#       1 15
#       1 5
#       1 4
#       4 3
#       5 2
#      56 1
