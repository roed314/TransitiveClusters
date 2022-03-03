#!/usr/bin/env python

# This script is used to collect the hashes computed by ComputeHashes.m into clusters.

import os
from collections import defaultdict

opj = os.path.join

clusters = defaultdict(list)
for nTt in os.listdir(opj("DATA", "trun")):
    with open(opj("DATA", "trun", nTt)) as F:
        ordhsh = F.read().strip()
    clusters[ordhsh].append(nTt)
for ordhsh, cluster in clusters.items():
    if len(cluster) == 1:
        fname = opj("DATA", "hash_unique", ordhsh)
    else:
        fname = opj("DATA", "active", ordhsh)
    with open(fname), "w") as F:
        _ = F.write("\n".join(cluster) + "\n")
