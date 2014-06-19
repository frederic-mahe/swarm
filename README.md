# swarm #

A robust and fast clustering method for amplicon-based studies. This
is an experimental python version of
[swarm](https://github.com/torognes/swarm "swarm public
repository"). The objective of that version is to identify the fastest
way to perform the clustering when the number of differences *d*
between amplicons is set to 1. That *d* parameters is set to 1 by
default in the official swarm version, as it yields high-resolution
clustering results, and as datasets grow in size, using *d* = 1 is
arguably the best decision.

With current clustering algorithms, the computation complexity is
**quadratic**. Multiplying by 10 the size of the dataset increases by a
factor 100 the computation time.

In swarm, using a fixed *d* = 1 value allows a radical change in the
algorithm design (while of course remaining exact). That new
algorithm, implemented in the python script presented here, has a
fantastic property: it has a **linear** computation
complexity. Increasing the dataset 10 times only increases the
computation time by a factor of 10. A major change in scalability.

To give some perspective, the monothreaded python script is already as
fast as the C++ multithreaded implementation of swarm on mid-size
amplicon datasets (appr. 1 million unique amplicons). The python
script is more than 10 times faster on large amplicon datasets (32
million unique amplicons). On an extremely large dataset (154 million
unique amplicons!), the script takes only 6 days to run, where the C++
version would take several months on a 16-cores computer.

Early tests with a C version of the key-function of the new algorithm
show a 10 times speed-up. We can confidently estimate that a smart
C/C++ re-implementation of the algorithm will bring a very significant
speed-up.

In conclusion, with that new swarm algorithm, fast, exact and accurate
partitioning of extremely large datasets, intractable with current
clustering methods suddenly becomes an easy computation task.

## Warning ##

Tested with python 2.7.3

## Quick start ##

To get basic usage and help, use the following command:

```
python swarm.py -h
```
