#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/06/07"
__version__ = "$Revision: 6.0"

import sys
from Bio import SeqIO
from optparse import OptionParser

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Find clusters of nearly identical amplicons in giant
    amplicon-based environmental or clinical projects."""

    parser = OptionParser(usage="usage: %prog --input_file filename",
                          description=desc,
                          version="%prog version 1.0")

    parser.add_option("-i", "--input_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="input_file",
                      help="set <FILENAME> as input fasta file.")

    (options, args) = parser.parse_args()
    return options.input_file


def parse_input_file(input_file):
    """
    Build a list of amplicons
    """
    input_format = "fasta"
    amplicons = dict()
    order = list()
    with open(input_file, "rU") as input_file:
        for record in SeqIO.parse(input_file, input_format):
            amplicons[str(record.seq)] = [record.id,
                                          len(record.seq),
                                          True]
            order.append(str(record.seq))
    return amplicons, order


def produce_microvariants(seq):
    """
    For a given sequence, produce without repetition all possible
    micro-variants with one difference (mutation, insertion,
    deletion).

    Some micro-variants are identical to the mother sequence: remove
    them.
    """
    nucleotides = ["a", "c", "g", "t"]
    seq = list(seq)
    length = len(seq)
    microvariants = list()
    # Insertions
    for i in xrange(0, length + 1, 1):
        for nuc in nucleotides:
            microvariants.append("".join(seq[0:i] + [nuc] + seq[i:]))
    # Mutations
    for i in xrange(0, length, 1):
        for nuc in nucleotides:
            tmp = seq[:]
            tmp[i] = nuc
            microvariants.append("".join(tmp))
    # Deletions
    for i in xrange(0, length, 1):
        tmp = seq[:]
        del tmp[i]
        microvariants.append("".join(tmp))
    # Cleaning
    microvariants = list(set(microvariants))
    seq = "".join(seq)
    del microvariants[microvariants.index(seq)]
    return microvariants


def main():
    """
    """
    # Parse command line options.
    input_file = option_parse()

    # Build a list of amplicons and count nucleotide occurences
    amplicons, order = parse_input_file(input_file)

    # Start swarming
    for seed in order:

        amplicon = amplicons[seed]
        if not amplicon[2]:  # Skip amplicons already swarmed
            continue

        swarm = [amplicon[0]]
        amplicons[seed][2] = False

        # Create micro-variants
        microvariants = produce_microvariants(seed)

        # Which of these microvariants are in our dataset?
        hits = [microvariant for microvariant in microvariants
                if microvariant in amplicons
                and amplicons[microvariant][2]]  # WARNING! for the
                                                 # post-processing,
                                                 # these hits would
                                                 # need to be sorted
                                                 # by decreasing
                                                 # abundance!

        # Add them to the swarm (if any)
        if not hits:  # Singleton, go to the next master seed
            print(" ".join(swarm), file=sys.stdout)
            continue
        swarm.extend([amplicons[j][0] for j in hits])
        for hit in hits:
            amplicons[hit][2] = False

        # Work on subseeds
        all_subseeds = [hits]
        for subseeds in all_subseeds:
            nextseeds = list()
            for subseed in subseeds:
                amplicons[subseed][2] = False
                # Search for k-seeds
                microvariants = produce_microvariants(subseed)
                hits = [microvariant for microvariant in microvariants
                        if microvariant in amplicons
                        and amplicons[microvariant][2]]   # WARNING!
                                                          # for the
                                                          # post-processing,
                                                          # these hits
                                                          # would need
                                                          # to be
                                                          # sorted by
                                                          # decreasing
                                                          # abundance!
                if not hits:  # subseed has no sons
                    continue
                nextseeds.extend(hits)
                swarm.extend([amplicons[hit][0] for hit in hits])
                for hit in hits:
                    amplicons[hit][2] = False
            if not nextseeds:  # No new subseeds
                print(" ".join(swarm), file=sys.stdout)
                break
            all_subseeds.append(nextseeds)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

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



# The remaining questions are:
# - how fast can we produce the microvariants (we have to do it for (n-1) of the n sequences)?
# - what are the most efficient mapping implementations to get items?

# The number of mapping queries will be (n - 1) * 8 *
# average__sequence_length (example: 33 billion lookups for TARA V9). Hence, we have a O(n) complexity!
