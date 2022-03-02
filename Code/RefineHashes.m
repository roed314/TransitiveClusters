/*****************************************************************************************************
This file is used for attempting to split up a cluster of groups having the same basic hash
using more refined hash functions.

USAGE:

ls DATA/active | parallel -j128 --timeout 14400 magma hsh:="{1}" RefineHashes.m

INPUT:

You should provide a variable `hsh` at the command line: the order and the hash value, separated by a period.

OUTPUT:

It will write new files with an additional hash value appended to the filename, either in the
DATA/refined_unique/ directory (for groups that have a distinct hash from any other row in the input), or
DATA/active (for clusters that still need further division).  Afterward, it will move the input file
to the DATA/inactive directory.

In addition, it will write timings to the DATA/refining.timings/ directory.
*****************************************************************************************************/

AttachSpec("hashspec");

SetColumns(0);
// Check to see that this is a simple hash, rather than a refined hash
pieces := Split(hsh, ".");
if #pieces ne 2 then
    print hsh, "is not simple";
    exit;
end if;
nTt_lookup := AssociativeArray();
nTts := [];
for cluster in Split(Read("DATA/active/" * hsh), "\n") do
    first := Split(cluster, " ")[1];
    nTt_lookup[first] := cluster;
    Append(~nTts, first);
end for;
t0 := Cputime();
collator := DistinguishingHashes(nTts);
PrintFile("DATA/refining.timings/" * hsh, Sprint(Cputime() - t0));
for hashes -> gps in collator do
    newhsh := hsh * "." * Sprint(CollapseIntList(hashes));
    folder := (#gps eq 1) select "refined_unique/" else "active/";
    PrintFile("DATA/" * folder * newhsh, Join([nTt_lookup[gp] : gp in gps], "\n"));
end for;
System("mv DATA/active/" * hsh * " DATA/inactive/");
exit;
