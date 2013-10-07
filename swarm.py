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
import itertools
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


def compute_kmer_vectors(amplicons, k):
    """
    Produce a kmer-vector representation for all amplicons
    """
    # Create a list of all possible kmers
    all_kmers = ["".join(kmer)
                 for kmer in itertools.product('acgt', repeat=k)]
    # Represent each amplicon as a 4^k-vector of booleans
    vectors = list()
    for amplicon in amplicons:
        sequence = amplicon[3]
        amplicon_kmers = set([sequence[i:i+k]
                              for i in range(len(sequence) - k)])
        vector = [kmer in amplicon_kmers for kmer in all_kmers]
        vectors.append(vector)
    return vectors


def compare_vectors(seed_vector, candidate_vectors):
    """
    XOR and POPCNT the seed vector against all candidate vectors
    """
    popcnts = []
    for j, candidate_vector in candidate_vectors:
        xoring = [status1 ^ status2
                  for status1, status2 in zip(seed_vector, candidate_vector)]
        popcnt = xoring.count(True)
        popcnts.append((j, popcnt))
    return popcnts


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


def find_kseeds(candidates, amplicons, threshold, i):
    """
    For a given subseed, find all k-seeds with threshold or less
    differences with the subseed. Candidates already have been
    filtered for "left hand" reads. Do not treat already assigned
    amplicons. Do not compare amplicons with a length difference
    greater than the threshold. Do not compare amplicons if their
    nucleotide profiles are too divergent.
    """
    # I do not use the nucleotide composition comparisons
    # sum([abs(cmp(couple[0], couple[1])) for couple in zip(amplicons[i][4], amplicons[j][4])]) <= 2 * threshold - abs(cmp(amplicons[i][2], amplicons[j][2]))
    hits = [j for j in candidates
            if needleman_wunsch(amplicons[i][3], amplicons[j][3]) <= threshold]
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

    ## Create a list of all possible kmers
    k = 5
    vectors = compute_kmer_vectors(amplicons, k)

    # Torbjørn formula is mindiff = (popcnt + 2*k - 1)/(2*k). I don't
    # understand why. I double-checked the maximum acceptable popcnt
    # is: d * k *2

    # Start swarming
    for i, amplicon in enumerate(amplicons):

        if not status[i]:  # Skip amplicons already swarmed
            continue

        swarm = [amplicons[i][0]]
        status[i] = False

        candidate_vectors = [(j, vectors[j]) for j, status[j] in enumerate(status) if status[j]]
        if not candidate_vectors:  # No more amplicons, stop swarm
            print(" ".join(swarm), file=sys.stdout)
            break

        # Search sub-seed candidates
        popcnts = compare_vectors(vectors[i], candidate_vectors)
        candidates_nw = [j for j, popcnt in popcnts if popcnt <= 2 * threshold * k]
        if not candidates_nw:  # subseed has no sons
            print(" ".join(swarm), file=sys.stdout)
            continue
        hits = find_kseeds(candidates_nw, amplicons, threshold, i)
        if not hits:  # Singleton, go to the next master seed
            print(" ".join(swarm), file=sys.stdout)
            continue
        swarm.extend([amplicons[j][0] for j in hits])
        for hit in hits:
            status[hit] = False

        # Work on subseeds
        all_subseeds = [hits]
        for subseeds in all_subseeds:
            nextseeds = list()
            for l in subseeds:
                status[l] = False
                # To break superswarms, I just had to allow only
                # candidates on the right side of l (i.e. at l+1)
                candidate_vectors = [(j, vectors[j]) for j, status[j] in enumerate(status) if status[j]]
                if not candidate_vectors:  # No more amplicons, stop swarm
                    print(" ".join(swarm), file=sys.stdout)
                    break
                popcnts = compare_vectors(vectors[l], candidate_vectors)
                candidates_nw = [j for j, popcnt in popcnts if popcnt <= 2 * threshold * k]
                if not candidates_nw:  # subseed has no sons
                    continue
                hits = find_kseeds(candidates_nw, amplicons, threshold, l)
                if hits:
                    nextseeds.extend(hits)
                    swarm.extend([amplicons[hit][0] for hit in hits])
                    for hit in hits:
                        status[hit] = False
            if not nextseeds:  # No new subseeds
                print(" ".join(swarm), file=sys.stdout)
                break
            all_subseeds.append(nextseeds)

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
