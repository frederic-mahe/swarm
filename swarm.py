#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects. That version works with score values (not the number of
    mismatches). Should I go back to counting mismatches? It saves the
    backtracking, but complexifies the threshold computation.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2013/02/24"
__version__ = "$Revision: 4.0"

import sys
from Bio import SeqIO
from optparse import OptionParser

#******************************************************************************#
#                                                                              #
#                                  Functions                                   #
#                                                                              #
#******************************************************************************#


def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Find clusters of nearly identical amplicons in giant
    amplicon-based environmental or clinical projects."""

    parser = OptionParser(usage="usage: %prog --input_file filename --threshold integer",
                          description=desc,
                          version="%prog version 1.0")

    parser.add_option("-i", "--input_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="input_file",
                      help="set <FILENAME> as input file. Ungapped fasta file (acgt only, case insensitive).")

    parser.add_option("-t", "--threshold",
                      metavar="<THRESHOLD>",
                      action="store",
                      dest="threshold",
                      type="int",
                      default=1,
                      help="set <THRESHOLD> for the swarm building.")

    (options, args) = parser.parse_args()
    return options.input_file, options.threshold


def parse_input_file(input_file):
    """
    Build a list of amplicons and count nucleotide occurences
    """
    input_format = "fasta"
    nucleotides = ("a", "c", "g", "t")
    with open(input_file, "rU") as input_file:
        amplicons = [(record.id.split("_")[0],
                      record.id.split("_")[1],
                      len(record.seq),
                      str(record.seq),
                      [record.seq.lower().count(n) for n in nucleotides])
                     for record in SeqIO.parse(input_file, input_format)]
        status = [True] * len(amplicons)
    return amplicons, status


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

    # Compute score for each cell
    for i in I[1:]:
        for j in J[1:]:
            match = F[i-1][j-1] + S[A[i-1]][B[j-1]]
            delete = F[i-1][j] + d
            insert = F[i][j-1] + d
            # Use max() for a similarity matrix
            F[i][j] = min(match, insert, delete)

    # Trackback
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
    # Deal with overhanging 5' parts (replace a while loop)
    if i:
        AlignmentA.extend(A[i-1:])
        AlignmentB.extend("-" * i)
    if j:
        AlignmentA.extend("-" * j)
        AlignmentB.extend(B[j-1:])

    # Mismatches
    AlignmentA = "".join(AlignmentA[::-1])
    AlignmentB = "".join(AlignmentB[::-1])
    lenA = len(AlignmentA)
    lenB = len(AlignmentB)

    if lenA != lenB:
        print("Something went really bad!", lenA, lenB, file=sys.stderr)

    mismatches = len([False for k in xrange(lenA)
                      if AlignmentA[k] is not AlignmentB[k]])
    return mismatches


def find_kseeds(candidates, status, frontier, amplicons, threshold):
    """
    For a given subseed, find all k-seeds with threshold or less
    differences with the subseed. Candidates already have been
    filtered for "left hand" reads. Do not treat already assigned
    amplicons. Do not compare amplicons we know to be too distant from
    the seed. Do not compare amplicons with a length difference
    greater than the threshold. Do not compare amplicons if their
    nucleotide profiles are too divergent.
    """
    hits = [j for j, d in candidates
            if status[j] and
            d <= frontier and
            sum([abs(cmp(couple[0], couple[1])) for couple in zip(amplicons[l][4], amplicons[j][4])]) <= 2 * threshold - abs(cmp(amplicons[l][2], amplicons[j][2])) and
            needleman_wunsch(amplicons[l][3], amplicons[j][3]) <= threshold]
    return hits


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    # Parse command line options.
    input_file, threshold = option_parse()

    # Build a list of amplicons and count nucleotide occurences
    amplicons, status = parse_input_file(input_file)

    # Start swarming
    while True:
        try:
            # Search the next master seed (first amplicon with status
            # set to "True")
            i = status.index(True)
            swarm = [amplicons[i][0]]
            status[i] = False
        except ValueError:
            # All amplicons have been treated, stop the process
            break

        # List remaining non-swarmed amplicons (status set to True)
        comparisons = [j for j, status[j] in enumerate(status) if status[j]]

        # No amplicon remains, stop the process
        if not comparisons:
            print(" ".join(swarm), file=sys.stdout)
            break

        # Compute and store pairwise distances between the master seed
        # and all remaining amplicons ("candidates")
        candidates = [(j, needleman_wunsch(amplicons[i][3], amplicons[j][3]))
                      for j in comparisons]

        # Parse candidates and select the first generation of subseeds
        firstseeds = [j for j, d in candidates if d <= threshold]
        swarm.extend([amplicons[j][0] for j in firstseeds])
        for j in firstseeds:
            status[j] = False

        # Loop other the first subseeds (if they exist)
        all_subseeds = list()
        if firstseeds:
            all_subseeds.append(firstseeds)
            # Loop other the list of firstseeds lists
            for k, subseeds in enumerate(all_subseeds):
                nextseeds = list()
                frontier = (k + 2) * threshold
                for l in subseeds:
                    hits = find_kseeds(candidates, status, frontier,
                                       amplicons, threshold)
                    if hits:
                        nextseeds.extend(hits)
                        swarm.extend([amplicons[j][0] for j in hits])
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
            # Deal with isolated amplicons
            print(" ".join(swarm), file=sys.stdout)

sys.exit(0)

# Profiling
# ---------
#
# python -m cProfile swarm.py -i ../examples/AF091148.fas > tmp
#
# Expected results
# ----------------
#
# time python swarm.py -i ../examples/AF091148.fas | awk '{print NF}' | sort -nr | uniq -c
#  1 878
#  1 423
#  1 15
#  1 5
#  1 4
#  4 3
#  5 2
# 56 1
