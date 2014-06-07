#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2013/10/07"
__version__ = "$Revision: 6.0"

import sys
import itertools
from Bio import SeqIO
from bitarray import bitarray, bitdiff
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
    Build a list of amplicons
    """
    input_format = "fasta"
    with open(input_file, "rU") as input_file:
        amplicons = [(record.id.split("_")[0],
                      record.id.split("_")[1],
                      len(record.seq),
                      str(record.seq))
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
        vector = bitarray([kmer in amplicon_kmers for kmer in all_kmers])
        vectors.append(vector)
    return vectors


def compare_vectors(seed_vector, candidate_vectors, max_kmer_diff):
    """
    XOR and POPCNT the seed vector against all candidate vectors
    (bitdiff does the same as (a ^ b).count(), but is more memory
    efficient)
    """
    candidates = [j for j, candidate_vector in candidate_vectors
                  if bitdiff(seed_vector, candidate_vector) <= max_kmer_diff]
    return candidates


def pairwise_alignment(seed, candidate):
    """
    Fast pairwise alignment (search only for a single mutation)
    """
    seed_length = len(seed)
    candidate_length = len(candidate)
    # Eliminate sequence too long or too short
    diff = abs(seed_length - candidate_length)
    if diff > 1:
        mismatches = diff
        return mismatches
    # Identify the longest sequence (in practice, turn all
    # insertion cases into deletions)
    if seed_length >= candidate_length:
        query, subject = seed, candidate
    else:
        query, subject = candidate, seed
    # Compare the sequences character by character
    length = len(query)
    stop_forward = length
    for i in xrange(0, length - 1, 1):
        if query[i] != subject[i]:
            stop_forward = i
            break
    stop_reverse = stop_forward
    for i in xrange(1, length - stop_forward, 1):
        if query[-i] != subject[-i]:
            stop_reverse = length - i
            break
    # Do we detect a single mutation or insertion-deletion?
    if stop_forward == stop_reverse:
        mismatches = 1
    else:
        mismatches = 2
    return mismatches


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
    max_kmer_diff = threshold * k * 2

    # Start swarming
    for i, amplicon in enumerate(amplicons):

        if not status[i]:  # Skip amplicons already swarmed
            continue

        swarm = [amplicons[i][0]]
        status[i] = False

        # List remaining amplicons
        candidate_vectors = [(j, vectors[j]) for j, status[j] in enumerate(status) if status[j]]
        if not candidate_vectors:  # No more amplicons, stop swarm
            print(" ".join(swarm), file=sys.stdout)
            break

        # Search sub-seed candidates
        candidates_nw = compare_vectors(vectors[i], candidate_vectors, max_kmer_diff)
        if not candidates_nw:  # subseed has no sons
            print(" ".join(swarm), file=sys.stdout)
            continue
        hits = [j for j in candidates_nw
                if pairwise_alignment(amplicons[i][3], amplicons[j][3]) <= threshold]
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
                candidates_nw = compare_vectors(vectors[l], candidate_vectors, max_kmer_diff)
                if not candidates_nw:  # subseed has no sons
                    continue
                hits = [j for j in candidates_nw
                        if pairwise_alignment(amplicons[l][3], amplicons[j][3]) <= threshold]
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

# The next step is to introduce 3-mer and 4-mer, to evaluate the number
# of 5-mer comparisons avoided. I also tried using 6-mers and it yielded
# a 10 s gain (34 s to run). The future is maybe to compute 6-mer
# vectors or 7-mers vectors when needed to keep avoiding the cost of a
# Needleman-Wunsch alignment.

# Finally, I need a faster Needleman-Wunsch implementation. I can do my
# own C version with python bindings, or wait for Torbjørn's version
# prepared by Umer.
