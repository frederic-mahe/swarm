#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2012/12/28"
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
        default = 1,
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
    # Penalty matrix and gap penalty (d) used in Swarm (transformed
    # +5/-4/+12/-4 model).
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
    for i in I[1:]:
	for j in J[1:]:
            match = F[i-1][j-1] + S[A[i-1]][B[j-1]]
            delete = F[i-1][j] + d
            insert = F[i][j-1] + d
            # Use max() for a similarity matrix 
            F[i][j] = min(match, insert, delete)

    ## Traceback
    #
    # The algorithm progresses backward. Appending to a string or to
    # the begining of a list is slow. Use a list instead of a string,
    # reverse the list and convert it back to a string at the end.
    AlignmentA = list()
    AlignmentB = list()
    i = len(A)
    j = len(B)

    while i > 0 and j > 0:
        Score = F[i][j]
        ScoreDiag = F[i - 1][j - 1]
        ScoreUp = F[i][j - 1]
        ScoreLeft = F[i - 1][j]
        if Score == ScoreDiag + S[A[i-1]][B[j-1]]:
            AlignmentA.append(A[i-1])
            AlignmentB.append(B[j-1])
            i -= 1
            j -= 1
        elif Score == ScoreLeft + d:
            AlignmentA.append(A[i-1])
            AlignmentB.append("-")
            i -= 1
        elif Score == ScoreUp + d:
            AlignmentA.append("-")
            AlignmentB.append(B[j-1])
            j -= 1
        else:
            print("Something went really bad!", A, B, sep="\n", file=sys.stderr)
            sys.exit(-1)
    # Deal with overhanging 5' parts.
    while i:
        AlignmentA.append(A[i-1])
        AlignmentB.append("-")
        i -= 1
    while j:
        AlignmentA.append("-")
        AlignmentB.append(B[j-1])
        j -= 1
    # Pythonic way?
    # if i:
    #     AlignmentA.append(A[i-1:len(i)])
    #     AlignmentB.append("-" * i)
    # if j:
    #     AlignmentA.append("-" * j)
    #     AlignmentB.append(B[j-1:len(j)])
        

    ## Similarity or Dissimilarity (reverse and convert back the lists
    # to strings)
    AlignmentA = "".join(AlignmentA[::-1])
    AlignmentB = "".join(AlignmentB[::-1])
    lenA = len(AlignmentA)
    lenB = len(AlignmentB)

    if lenA >= lenB:
        sim1, sim2 = AlignmentA, AlignmentB
        len0 = lenA
    else:
        sim1, sim2 = AlignmentB, AlignmentA
        len0 = lenB

    # matches = len([True for k in range(len0) if sim1[k] is sim2[k]])
    mismatches = len([False for k in range(len0) if sim1[k] is not sim2[k]])
    # score = abs(F[-1][-1])

    ## Visualize
    # visual = list()
    # for k in range(len0):
    #     if sim1[k] is sim2[k]:
    #         visual.append("|")
    #     else:
    #         visual.append(" ")
    # print(AlignmentA, "".join(visual), AlignmentB, sep="\n", file=sys.stdout)

    # similarity = 100.0 * matches / len0
    # dissimilarity = 100.0 * mismatches / len0
    
    # return score
    return mismatches


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

        # List remaining non-swarmed sequences
        comparisons = [j for j, state in enumerate(status) if state is True]
        
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
                frontier = (2 + k) * threshold
                # I could replace that for loop with a list
                # comprehension, but the risk is to perform un-needed
                # amplicon comparisons, as the status list is not
                # updated regularly .
                for l in subseeds:
                    # Candidates already have been filtered for "left
                    # hand" reads. Do not treat already assigned
                    # sequences. Do not compare sequences we know to
                    # be too distant from the seed. Do not compare
                    # sequences with a length difference greater than
                    # the threshold. Do not compare sequences if their
                    # nucleotide profiles are too divergent.
                    hits = [j for j, d in candidates
                            if status[j] is True
                            and d <= frontier
                            and abs(cmp(records_list[l][2], records_list[j][2])) <= threshold
                            and sum([abs(cmp(couple[0], couple[1])) for couple in zip(records_list[l][4], records_list[j][4])]) <= 2 * threshold - abs(cmp(records_list[l][2], records_list[j][2]))
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
