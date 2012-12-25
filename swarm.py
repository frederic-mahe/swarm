#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@sb-roscoff.fr>"
__date__ = "2012/05/14"
__version__ = "$Revision: 1.0"

import sys
from editdist import distance
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
        help = "set <FILENAME> as input file. Ungapped fasta file.")

    parser.add_option("-t", "--threshold",
        metavar = "<THRESHOLD>",
        action = "store",
        dest = "threshold",
        type = "int",
        default = 1,
        help = "set <THRESHOLD> for the clique building.")

    (options, args) = parser.parse_args()
    
    return options.input_file, options.threshold

#**********************************************************************#
#                                                                      #
#                              Body                                    #
#                                                                      #
#**********************************************************************#

if __name__ == '__main__':
    
    # Parse command line options.
    input_file, threshold = option_parse()

    # Build a list of sequences and initialize lists
    input_format = "fasta"
    records_list = list()
    status = list()
    distances = list()
    with open(input_file, "rU") as input_file:
        records = SeqIO.parse(input_file, input_format)
        records_list = tuple([(record.id.split("_")[0], record.id.split("_")[1], len(record.seq), str(record.seq)) for record in records])
        status = [True for i in xrange(len(records_list))]

    # Start cliquing
    while True:
        try:
            # Search the next master seed
            i = status.index(True)
        except ValueError:
            # All sequences have been treated
            break

        # Master seed
        father_and_sons = dict()
        individual_results = dict()
        clique = [records_list[i][0]]
        status[i] = False

        # List remaining non-cliqued sequences
        comparisons = tuple([j for j, state in enumerate(status) if state is True])
        
        # Deal with the last remaining item
        if not comparisons:
            print(" ".join(clique), file=sys.stdout)
            break

        # Candidates list is initialized for each major seed
        candidates = tuple([(j, distance(records_list[i][3], records_list[j][3])) for j in comparisons])
        individual_results[records_list[i][0]] = dict([[j, d] for j, d in candidates if d > threshold])

        # Parse candidates and select sons
        firstseeds = [j for j, d in candidates if d <= threshold]
        clique.extend([records_list[j][0] for j in firstseeds])
        for j in firstseeds:
            status[j] = False
            father_and_sons[records_list[j][0]] = records_list[i][0]

        # Loop other the first subseeds (if they exists)
        all_subseeds = list()
        if firstseeds:
            all_subseeds.append(firstseeds)
            # Loop other the list of firstseeds lists
            for k, subseeds in enumerate(all_subseeds):
                nextseeds = list()
                frontier = (2 + k) * threshold
                # I could replace that for loop with a list
                # comprehension, but the risk is to perform un-needed
                # distance computations, as the status list is not
                # updated regularly .
                for l in subseeds:
                    # Candidates already have been filtered for "left
                    # hand" reads. Do not treat already assigned
                    # sequences. Do not compare sequences we know to
                    # be too distant. Do not compare sequences with a
                    # length difference greater than the threshold.
                    # Keep if its below the threshold.
                    #
                    # Update the candidate list with data from
                    # sons. We have to climb back the tree of
                    # relations up to the sommit. We do it only if the
                    # father is not the root (k != 0 is interpreted as
                    # k is True)
                    if k > 0:
                        candidates_list = list()
                        son = records_list[l][0]
                        # Loop until you reach the root
                        while True:
                            try:
                                father = father_and_sons[son]
                                candidates_list.append(individual_results[father])
                                son = father
                            except KeyError:
                                # The record has no father
                                break
                        # Create an updated candidate list by merging
                        # the individual results (if the list is not
                        # empty)
                        if candidates_list:
                            # The last list is the root results dict
                            updated_candidates = candidates_list.pop()
                            # Reverse so we can go from the oldest to
                            # the newest results
                            candidates_list.reverse()
                            # Convert dict to tuples
                            for d in candidates_list:
                                updated_candidates.update(d)
                        # Convert dict to tuples (sorted)
                        updated_candidates = updated_candidates.items()
                    else:
                        updated_candidates = candidates
                    # Compute distances
                    distances = [(j, distance(records_list[l][3], records_list[j][3])) for j, d in updated_candidates
                            if status[j]
                            and d <= frontier
                            and abs((records_list[l][2] - records_list[j][2])) <= threshold]
                    # Extract new sons
                    hits = [j for j, d in distances if d <= threshold]
                    nextseeds.extend(hits)
                    clique.extend([records_list[j][0] for j in hits])
                    for j in hits:
                        status[j] = False
                        father_and_sons[records_list[j][0]] = records_list[l][0]
                    # Store results for future sons
                    individual_results[records_list[l][0]] = dict([[j, d] for j, d in distances if d > threshold])
                # Stop condition
                if nextseeds:
                    all_subseeds.append(nextseeds)
                else:
                    # No new subseeds
                    print(" ".join(clique), file=sys.stdout)
                    # print(father_and_sons, file=sys.stdout)
                    break
        else:
            # Deal with isolated sequences (no father)
            print(" ".join(clique), file=sys.stdout)
        
sys.exit(0)

# Notes
#
# Profiling (python -m cProfile DNA_clique.py -i file.fas > tmp)
# shows that more than 80% of the computing time is spent on
# editdist.distance. If I want to accelerate the algorithm, I only
# have two options:
#
# - use a faster distance function (80 µs/comparison),
# - reduce the number of sequence comparisons.
#
# I don't know how to do that yet. The improvement possibilities in
# the python part are now close to zero.
# 
# - should I use arrays: NO
# - I use tuples and tuples of # tuples when possible (less memory, faster handling)
# - Comparing length may be a waste of time: NO
#
# Future improvements
# -------------------
#
# time python DNA_clique.py -i AF091148.fas > tmp && awk '{print NF}' tmp | sort -nr | uniq -c
#       1 878
#       1 423
#       1 15
#       1 5
#       1 4
#       4 3
#       5 2
#      56 1
