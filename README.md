# TransitiveClusters
Code for dividing up transitive groups (which are classified up to isomorphism as permutation groups) into isomorphism classes as abstrac groups.

Overall approach
================

Given a set of transitive permutation groups, we begin by computing a hash function for each, invariant under isomorphism of abetract groups.
This hash incorporates the following information:

- If G has order contained in the database of small groups (roughly, order at most 2000),
we use the small group id as the hash.
- If G is abelian, we use the abelian invariants for the hash.
- Otherwise, for each maximal subgroup (up to conjugacy), we compute the small group id or abelian invariants
(if applicable), or the orders and sizes of conjugacy classes (if not).

The combination of order and hash value divide up the input groups into clusters,
with groups in different clusters definitely not isomorphic, and groups within
the same cluster possibly isomorphic.  We then apply several strategies to refine these clusters
by finding isomorphisms between groups in them or splitting them apart into smaller clusters.

0. Of course, if a cluster has size 1 it is finalized.
1. Using a randomized strategy, we search for subgroups of groups in the cluster with trivial core and specified index
(determined by the set of degrees within the cluster).  The corresponding permutation action on cosets
gives another representation of the original group as a transitive group, which may be non-isomorphic
as a permutation group thus providing an isomorphism with another group in the cluster.
This strategy leaves the set of clusters unchanged, but finds isomorphisms within a cluster.
If all groups within a cluster are found to be isomorphic to each other, we can finalize the cluster.
It is implemented in `RandomizedMerge.m`.
2. We compute additional invariants for groups in the cluster, attempting to distinguish them.
In particular, we compute hashes for their Sylow subgroups, derived subgroups, minimal normal subgroups
and the corresponding normalizers/quotients, together with the degrees of the irreducible characters.
This process does not find isomorphisms within a cluster, but can break clusters apart into smaller pieces
(possibly finalizing them if all groups within the smaller piece are known to be isomorphic).
It is implemented in `RefineHashes.m`.
3. If neither of these strategies are successful, we attempt pairwise isomorphism tests between groups
within the cluster using Magma's IsIsomorphic function.  This process is generally slower
than the previous two approaches, but can both find isomorphisms within a cluster and break clusters
into smaller pieces.  It is implemented in `FullIsoTest.m`.
4. Another fallback is to find all core-free subgroups of a given order rather than searching randomly.
This is usually the slowest approach, but the additional information of the number of "siblings" of each type
is useful in some applications.  It is implemented in `FullSibSplit.m`.
5. If these are too slow, `IterativeSubgroupSetup.m` can be used to create scripts that allow for easily parallelizing the call to `Subgroups` (for example, by finding subgroups of index 4 inside subgroups of index 2 rather than subgroups of index 8).

Usage
=====

The overall input is a file `DATA/hash.todo`, containing a list of transitive labels (like 44T1780) to consider,
one on each line.  In the original use case, we took the set of all transitive groups contained
within the LMFDB (degree up to 47 skipping those of degree 32 with order between 512 and 4*10^10)
that were outside the scope of the small group database.  Note that it was difficult to compute the hash for
some of the larger groups, but we were able to finalize them early because they were the only group with
a given order (like 47!) and thus finding a hash value was unnecessary.

Clusters are stored in separate files:
- the file name indicates identifying information about the cluster: the order, hash value, refined hash value (optional), or first label within the cluster (used by FullIsoTest.m and FullSibSplit.m since a refined hash may not be sufficient).
- Each row of the file consists of groups known to be isomorphic to each other, listed as space separated transitive labels (like 44T1780)
- Different rows of the cluster may or may not be isomorphic.  Once a cluster has been reduced to a single row it is "finalized" and moved out of the active folder and into a final folder.

Before beginning, you should create the following folders.  The following folders are where finalized clusters will be stored:
```
DATA/hash_unique/ -- clusters finalized because there is only one group with a given order and hash
DATA/merge_finished/ -- clusters that were finalizeded by merging (finding isomorphisms within the cluster until all groups belonged to a single isomorphism class)
DATA/refined_unique/ -- clusters that were finalized when computing refined hash values
DATA/isotest_finished/ -- clusters that were finished using Magma's IsIsomorphic function
DATA/sibs_finished/ -- clusters that were finished using Magma's Subgroups function
```
The following folders are used to hold information about sibling counts
```
DATA/sibs_with_count/ -- additional information computed giving sibling information; not in the standard cluster format
DATA/last.selfsib/ -- used for storing generators of core-free subgroups found by `IterativeSubgroupSetup.m` that yield the same transitive label (and thus do not help with finding isomorphisms within a cluster, but are useful if trying to assemble sibling data
DATA/last.osib/ -- used for storing generators of core-free subgroups found by `IterativeSubgroupSetup.m` that yield a different transitive label, thus producing an isomorphism within the cluster.
```
The following folders are used by some of the scripts as a place to hold intermediate results:
```
DATA/trun/ -- files with basic hash values are stored here
DATA/active/ -- active clusters are stored here
DATA/inactive/ -- clusters are moved here when they are split into pieces by RefineHashes.m (contents can be deleted if there are no problems)
DATA/tmp/ -- used when modifying files in the active directory; should be empty after any run
DATA/merge_check/ -- used for recording runs of RandomizedMerge.m
DATA/last/ -- used for holding the scripts created by `IterativeSubgroupSetup.m`
DATA/lastdone/ -- used for monitoring whether the scripts created by `IterativeSubgroupSetup.m` have been run to completion
```
In addition, you should create folders to store timings:
```
DATA/trun.timings/ -- timings for computing the initial hash values
DATA/merge.timings/ -- timings for merging (only recorded on successful runs)
DATA/refining.timings/ -- timings for computing refined hash values
DATA/isotest.timings/ -- timings for running IsIsomorphic on pairs of groups within a cluster
DATA/sibs.timings/ -- timings for running Subgroups on groups in a cluster
DATA/last.timings/ -- timings for running Subgroups iteratively
```

1. After creating the input file `DATA/hash.todo`, to compute the initial hash values, run something like
```
cat DATA/hash.todo | parallel -j128 magma nTt:="{1}" ComputeHashes.m
```
2. Once complete, assemble the hash values into an `active` folder as follows:
```
./AssembleHashes.py
```
3. Do a pass with randomized isomorphism testing to attempt to merge some clusters:
```
ls DATA/active | parallel -j128 --timeout 600 magma hsh:="{1}" timeout:=600 RandomizedMerge.m
```
4. For the clusters that remain, compute refined hashes to try to split them up:
```
ls DATA/active | parallel -j128 --timeout 14400 magma hsh:="{1}" RefineHashes.m
```
5. Do another pass with randomized isomorphism testing:
```
ls DATA/active | parallel -j128 --timeout 1800 magma hsh:="{1}" timeout:=1800 RandomizedMerge.m
```
6. Try using IsIsomorphic (possibly modifying the script to skip large clusters that remain):
```
ls DATA/active | parallel -j128 --timeout 3600 --memfree 100G --joblog DATA/isotest.log magma hsh:="{1}" FullIsoTest.m
```
7. Try using Subgroups:
```
ls DATA/active | parallel -j128 --timeout 3600 --memfree 100G --joblog DATA/sibs.log magma hsh:="{1}" FullSibSplit.m
```
8. Try some more randomized isomorphism testing:
```
ls DATA/active | parallel -j128 --timeout 3600 magma hsh:="{1}" timeout:=3600 RandomizedMerge.m
```
9. If clusters remain, use `IterativeSubgroupSetup.m`, which hasn't been automated and still needs a bit of attention.

The commands above were run on a server with 128 physical cores and 2TB of RAM; they should be modified as appropriate to take your resources into account.

Assumptions
===========

Some of the code may assume that all transitive groups with a given order and hash are listed in the initial `DATA/hash.todo`, since this was the case in our application.