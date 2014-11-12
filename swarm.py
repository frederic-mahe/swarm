#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Find clusters of nearly identical sequences in giant barcoding
    projects.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/11/10"
__version__ = "$Revision: 8.0"

import os
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
            seq = str(record.seq).lower()  # Convert all sequences to lowercase
            # Store 0) amplicon_id, 1) amplicon abundance, 2) amplicon
            # status, 3) swarm mass, 4) swarm seed id
            amplicons[seq] = [record.id,
                              int(record.id.split("_")[1]),
                              True,
                              0,
                              ""]
            order.append(seq)
    return amplicons, order


def output_swarms(input_file, order, amplicons, swarms, d):
    """
    Write swarms to a file
    """
    # Create new file name
    extension = "_" + str(d) + ".swarms_2_prototype"
    output_file = os.path.splitext(os.path.abspath(input_file))[0] + extension
    with open(output_file, "w") as output_file:
        for seed in order:
            seed_id = amplicons[seed][0]
            if seed_id in swarms:
                print(" ".join(swarms[seed_id]), file=output_file)
    return


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

    # Build a dict of amplicons and a list to preserve input order
    # (assuming decreasing abundance)
    amplicons, order = parse_input_file(input_file)

    # Store each swarm created in a dictionary: swarm[seed] = list(amplicons)
    swarms = dict()

    # Phase 1: d = 1 clustering -----------------------------------------------

    # Start swarming
    for seed in order:

        seed_id, seed_abundance, seed_status = amplicons[seed][0:3]
        if not seed_status:  # Skip amplicons already swarmed
            continue

        # Seed id and status
        swarm = [seed_id]
        amplicons[seed][2] = False  # could be seed_status = False
        amplicons[seed][4] = seed_id  # point to itself

        # Create micro-variants
        microvariants = produce_microvariants(seed)

        # Which of these microvariants are in our dataset?
        hits = [(microvariant, amplicons[microvariant][1])
                for microvariant in microvariants
                if microvariant in amplicons
                and amplicons[microvariant][2]]  # No need to check
                                                 # abundance here

        # Isolated seed? close the swarm
        if not hits:
            swarms[swarm[0]] = swarm
            # Add mass information
            amplicons[seed][3] = seed_abundance
            continue

        # Sort by decreasing abundance
        hits.sort(key=itemgetter(1, 0), reverse=True)

        # Add them to the swarm and update their status
        swarm.extend([amplicons[hit[0]][0] for hit in hits])
        for hit in hits:
            hit_seq = hit[0]
            amplicons[hit_seq][2] = False
            amplicons[hit_seq][4] = seed_id  # point to swarm seed
            if output_pairwise:  # Swarm breaker option activated
                print("@@", seed_id, amplicons[hit_seq][0], "1",
                      sep="\t", file=sys.stderr)

        # Work on subseeds (also save a list of hits along the way)
        all_hits = [(seed, seed_abundance)]
        all_subseeds = [hits]
        all_hits.extend(hits)
        for subseeds in all_subseeds:
            nextseeds = list()
            for subseed in subseeds:
                subseed_seq, subseed_abundance = subseed[0], subseed[1]
                # Update subseed status
                amplicons[subseed_seq][2] = False
                # Produce all microvariants of subseed
                microvariants = produce_microvariants(subseed_seq)

                # Which of these microvariants are in our dataset?
                # (discard hits with higher abundance value)
                hits = [(microvariant, amplicons[microvariant][1])
                        for microvariant in microvariants
                        if microvariant in amplicons
                        and amplicons[microvariant][2]
                        and amplicons[microvariant][1] <= subseed_abundance]

                if not hits:  # subseed has no son
                    continue

                # Sort by decreasing abundance
                hits.sort(key=itemgetter(1, 0), reverse=True)
                nextseeds.extend(hits)  # HITS ARE NOT GLOBALLY SORTED
                                        # BY DECREASING
                                        # ABUNDANCE. POSSIBLE SOURCE
                                        # OF CLUSTERING VARIATION?

                # Add hits to the swarm and update their status
                swarm.extend([amplicons[hit[0]][0] for hit in hits])
                all_hits.extend(hits)
                for hit in hits:
                    hit_seq = hit[0]
                    amplicons[hit_seq][2] = False
                    amplicons[hit_seq][4] = seed_id  # point to swarm seed
                    if output_pairwise:  # Swarm breaker option activated
                        print("@@", amplicons[subseed_seq][0],
                              amplicons[hit_seq][0], "1",
                              sep="\t", file=sys.stderr)

            if not nextseeds:  # No new subseeds, end of the all_subseeds list
                swarms[swarm[0]] = swarm  # Memorize the swarm
                mass = sum([hit[1] for hit in all_hits])
                for hit in all_hits:  # Set swarm mass value for each amplicon
                    amplicons[hit[0]][3] = mass
                break
            all_subseeds.append(nextseeds)

    # Output swarms (d = 1)
    output_swarms(input_file, order, amplicons, swarms, 1)

    # Phase 2: d = 2 clustering -----------------------------------------------

    # The process seems to detect the same candidates many times. I
    # need to mark the sequences already seen to get a more precise
    # evaluation of the dynamic of recruitment. After the first phase,
    # all amplicons are marked as False. I can re-use that data
    # structure to my advantage.

    # Number of OTUs
    n_OTUs = len(swarms)

    for seed in order:
        captured_OTUs = 0
        seed_id, seed_abundance, seed_status, swarm_mass, main_seed = amplicons[seed]

        # Do not deal with rare amplicons
        if seed_abundance == 2:
            break

        # Create micro-variants
        microvariants = produce_microvariants(seed)

        # Which of these microvariants are not in our dataset?
        fails = [microvariant
                 for microvariant in microvariants
                 if microvariant not in amplicons]

        for fail in fails:
            microvariants_2d = produce_microvariants(fail)

            # Which of these microvariants are not in our dataset?
            # OTUs of size 2 or less are dust.
            hits = [(microvariant, amplicons[microvariant])
                    for microvariant in microvariants_2d
                    if microvariant in amplicons
                    and amplicons[microvariant][3] <= 2
                    and amplicons[microvariant][2] is False]

            # Add hits to the swarm to which the seed belongs
            for hit in hits:
                amplicons[hit[0]][2] = True
                hit_id, hit_pointer = hit[1][0], hit[1][4]
                if hit_pointer in swarms:
                    swarms[main_seed] += swarms[hit_pointer]
                    del swarms[hit_pointer]
                    captured_OTUs += 1

        # Basic stats
        n_OTUs -= captured_OTUs
        print(main_seed, seed_id, seed_abundance,
              len(seed), len(microvariants),
              len(fails), captured_OTUs, n_OTUs,
              file=sys.stderr)

    # Output swarms (d = 2)
    output_swarms(input_file, order, amplicons, swarms, 2)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)

## TODO
# the algorithm could be simplified instead of adding the mass values
# to the amplicons, I added it to the swarms dict. It would cost a
# look up in the amplicons dict first to get the swarm seed, then a
# look up in the swarms dict to get the mass of the swarm.
# 
# During the second step, the mass values should be updated?? Not sure
# it is usefull.
