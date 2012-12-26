#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@sb-roscoff.fr>"
__date__ = "2012/05/14"
__version__ = "$Revision: 2.0"

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
        help = "set <FILENAME> as input file. Ungapped lowercase fasta file (acgt only).")

    parser.add_option("-t", "--threshold",
        metavar = "<THRESHOLD>",
        action = "store",
        dest = "threshold",
        type = "int",
        default = 1,
        help = "set <THRESHOLD> for the swarm building.")

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

    # Build a list of sequences and count nucleotide occurences
    input_format = "fasta"
    records_list = list()
    status = list()
    distances = list()
    with open(input_file, "rU") as input_file:
        records = SeqIO.parse(input_file, input_format)
        records_list = [(record.id.split("_")[0], record.id.split("_")[1], len(record.seq), str(record.seq), [record.seq.count(nuc) for nuc in ("a","c","g","t")]) for record in records]
        status = [True for i in xrange(len(records_list))]

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

        # List remaining non-swarmed sequences
        comparisons = [j for j, state in enumerate(status) if state is True]
        
        # Deal with the last remaining item
        if not comparisons:
            print(" ".join(swarm), file=sys.stdout)
            break

        # Candidates list is initialized for each major seed
        candidates = [(j, distance(records_list[i][3], records_list[j][3])) for j in comparisons]

        # Parse candidates and select sons
        firstseeds = [j for j, d in candidates if d <= threshold]
        swarm.extend([records_list[j][0] for j in firstseeds])
        for j in firstseeds:
            status[j] = False

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
                    # length difference greater than the threshold. Do
                    # not compare sequences if their nucleotide
                    # profile Keep if its below the threshold.
                    hits = [j for j, d in candidates
                            if status[j] is True
                            and d <= frontier
                            and abs((records_list[l][2] - records_list[j][2])) <= threshold
                            and sum([abs(cmp(couple[0], couple[1])) for couple in zip(records_list[l][4], records_list[j][4])]) <= 2 * threshold
                            and distance(records_list[l][3], records_list[j][3]) <= threshold]
                    nextseeds.extend(hits)
                    swarm.extend([records_list[j][0] for j in hits])
                    for j in hits:
                        status[j] = False
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
