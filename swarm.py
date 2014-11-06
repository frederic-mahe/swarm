#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/11/06"
__version__ = "$Revision: 7.0"

import sys
from Bio import SeqIO
from operator import itemgetter
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

    parser.add_option("-b", "--breaker",
                      action="store_true",
                      dest="output_pairwise",
                      default=False,
                      help="Output pairwise relations.")

    (options, args) = parser.parse_args()
    return options.input_file, options.output_pairwise


def parse_input_file(input_file):
    """
    Build a list of amplicons
    """
    input_format = "fasta"
    amplicons = dict()
    order = list()
    with open(input_file, "rU") as input_file:
        for record in SeqIO.parse(input_file, input_format):
            # Store 0) amplicon_id, 1) amplicon abundance, 2) amplicon status
            amplicons[str(record.seq)] = [record.id,
                                          int(record.id.split("_")[1]),
                                          True]
            order.append(str(record.seq))
    return amplicons, order


def produce_microvariants(seq):
    """
    For a given sequence, produce without repetition all possible
    micro-variants with one difference (mutation, insertion,
    deletion).
    """
    nucleotides = ("a", "c", "g", "t")
    seq = list(seq)
    length = len(seq)
    microvariants = list()
    # Insertions
    for i in xrange(0, length, 1):
        tmp = seq[:]
        tmp.insert(i, "")  # insert once (costly)
        for nuc in nucleotides:  # change four times (cheap)
            if tmp[i+1] is not nuc:
                tmp[i] = nuc
                microvariants.append("".join(tmp))
    # Insertions at the last position
    for nuc in nucleotides:
        tmp = seq[:] + [nuc]
        microvariants.append("".join(tmp))
    # Mutations and deletions
    for i in xrange(0, length, 1):
        tmp = seq[:]
        initial = tmp[i]
        for nuc in nucleotides:
            if initial is not nuc:  # Avoid useless mutations
                tmp[i] = nuc
                microvariants.append("".join(tmp))
        tmp[i] = initial  # Restore the initial sequence
        # Deletions
        try:
            if tmp[i] is not tmp[i+1]:  # delete at the end of homopolymers
                del tmp[i]
                microvariants.append("".join(tmp))
        except IndexError:  # Deletion at the last position
            del tmp[i]
            microvariants.append("".join(tmp))
    return microvariants


def main():
    """
    Load and parse input fasta file and clusterize it.
    """
    # Parse command line options.
    input_file, output_pairwise = option_parse()

    # Build a list of amplicons and count nucleotide occurences
    amplicons, order = parse_input_file(input_file)

    # Start swarming
    for seed in order:

        amplicon = amplicons[seed]
        if not amplicon[2]:  # Skip amplicons already swarmed
            continue

        # Seed id, abundance and status
        swarm = [amplicon[0]]
        amplicons[seed][2] = False

        # Create micro-variants
        microvariants = produce_microvariants(seed)

        # Which of these microvariants are in our dataset? 
        hits = [(microvariant, amplicons[microvariant][1])
                for microvariant in microvariants
                if microvariant in amplicons
                and amplicons[microvariant][2]]  # No check abundance here

        # Isolated seed? close the swarm
        if not hits:
            print(" ".join(swarm), file=sys.stdout)
            continue

        # Sort by decreasing abundance and remove abundance values
        hits.sort(key=itemgetter(1, 0), reverse=True)
        hits = map(itemgetter(0), hits)

        # Add them to the swarm and update their status
        swarm.extend([amplicons[hit][0] for hit in hits])
        for hit in hits:
            amplicons[hit][2] = False

        if output_pairwise:  # Swarm breaker option activated
            for hit in hits:
                print("@@", amplicon[0], amplicons[hit][0], "1",
                      sep="\t", file=sys.stderr)

        # Work on subseeds
        all_subseeds = [hits]
        for subseeds in all_subseeds:
            nextseeds = list()
            for subseed in subseeds:
                subseed_abundance = amplicons[subseed][1]
                amplicons[subseed][2] = False
                # Search for k-seeds (discard hits with higher abundance value)
                microvariants = produce_microvariants(subseed)
                hits = [(microvariant, amplicons[microvariant][1])
                        for microvariant in microvariants
                        if microvariant in amplicons
                        and amplicons[microvariant][2]
                        and amplicons[microvariant][1] <= subseed_abundance]

                if not hits:  # subseed has no son
                    continue

                # Sort by decreasing abundance and remove abundance values
                hits.sort(key=itemgetter(1, 0), reverse=True)
                hits = map(itemgetter(0), hits)
                nextseeds.extend(hits)

                # Add them to the swarm and update their status
                swarm.extend([amplicons[hit][0] for hit in hits])
                for hit in hits:
                    amplicons[hit][2] = False

                if output_pairwise:  # Swarm breaker option activated
                    for hit in hits:
                        print("@@", amplicon[0], amplicons[hit][0], "1",
                              sep="\t", file=sys.stderr)

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
