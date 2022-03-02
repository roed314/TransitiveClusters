/*****************************************************************************************************
This file uses Magma's Subgroups function to find all core-free subgroups of an appropriate index,
and the TransitiveGroupIdentification function to determine the transitive labels of groups within the
isomorphism class.  It is often slower than FullIsoTest.m and can take a lot of memory for large groups,
but it does provide additional information (the counts of the number of siblings).

USAGE:

ls DATA/active | parallel -j24 --timeout 3600 --memfree 100G --joblog DATA/sibs.log magma hsh:="{1}" FullSibSplit.m

INPUT:

You should provide a variable `hsh` at the command line: the order and hash value, separated by a period.

For simplicity, the code currently assumes that the input clusters all consist of transitive groups of a single degree
(except that in some cases every row has a group of degree 40 even though some also have groups of degree 20).
If this assumption is not satisfied it will currently raise an assertion error, but it should be fairly
easy to remove this hypothesis.

OUTPUT:

Every time it finds a complete isomorphism class from among the rows of the input file,
it will write the result to a file in the DATA/sibs_finished directory.  When all input rows have been
assigned, it will delete the input file.  Timings and memory usage are written to the DATA/isotest.timings/
folder, and information on the sibling counts are written to the DATA/sibs_with_count/ folder.
*****************************************************************************************************/

AttachSpec("hashspec");

SetColumns(0);
cluster_lookup := AssociativeArray();
first_lookup := AssociativeArray();
G_lookup := AssociativeArray();
nTts := [];
ns := {};
activefile := "DATA/active/" * hsh;
for cluster in Split(Read(activefile), "\n") do
    first := Split(cluster, " ")[1];
    cluster_lookup[first] := cluster;
    for nTt in Split(cluster, " ") do
        pair := [StringToInteger(c) : c in Split(nTt, "T")];
        first_lookup[pair] := first;
    end for;
    n, t := Explode([StringToInteger(c) : c in Split(first, "T")]);
    G_lookup[first] := TransitiveGroup(n, t);
    Append(~nTts, first);
    // The active clusters can all be handled by considering a single index: 36, 40 or 44.  Some of the 40 clusters have degree 20 siblings mixed in.
    if n eq 20 then
        for nTt in Split(cluster, " ") do
            if nTt[1..3] eq "40T" then
                n := 40;
                break;
            end if;
        end for;
    end if;
    Include(~ns, n);
end for;

assert #ns eq 1;

while #nTts gt 0 do
    t0 := Cputime();
    label := nTts[1];
    G := G_lookup[label];
    S := Subgroups(G : IndexEqual:=n);
    cnt1 := #S;
    S := [H`subgroup : H in S | #Core(G, H`subgroup) eq 1];
    cnt2 := #S;
    ts := [TransitiveGroupIdentification(Image(CosetAction(G, H))) : H in S];
    sibs := AssociativeArray();
    for t in ts do
        if not IsDefined(sibs, [n,t]) then
            first := first_lookup[[n,t]];
            Exclude(~nTts, first);
            cluster := cluster_lookup[first];
            for nTt in Split(cluster, " ") do
                sibs[[StringToInteger(c) : c in Split(nTt, "T")]] := 0;
            end for;
        end if;
        sibs[[n,t]] +:= 1;
    end for;
    // Since this process may get killed, we want to write output now
    print "Writing progress";
    mem := GetMaximumMemoryUsage();
    if mem lt 2^30 then
        mem := Sprintf("%oMB", RealField(4)!(mem / 2^20));
    else
        mem := Sprintf("%oGB", RealField(4)!(mem / 2^30));
    end if;
    PrintFile("DATA/sibs.timings/" * hsh, Sprintf("Subs(%o) -> %o -> %o -> %o in %os using %o", label, cnt1, cnt2, #sibs, Cputime() - t0, mem));
    sibs := [<k, v> : k -> v in sibs];
    Sort(~sibs);
    withcount := Join([Sprintf("%oT%o:%o", x[1][1], x[1][2], x[2]) : x in sibs], " ");
    nocount := Join([Sprintf("%oT%o", x[1][1], x[1][2]) : x in sibs], " ");
    PrintFile("DATA/sibs_finished/" * hsh * "." * label, nocount);
    PrintFile("DATA/sibs_with_count/" * hsh * "." * label, withcount);
    if #nTts gt 1 then
        tmp := "DATA/tmp/" * hsh;
        PrintFile(tmp, Join([cluster_lookup[nTt] : nTt in nTts], "\n"));
        System("mv " * tmp * " " * activefile); // mv is atomic
    else
        if #nTts eq 1 then
            print "Writing last";
            PrintFile("DATA/sibs.timings/" * hsh, "One cluster left");
            PrintFile("DATA/sibs_finished/" * hsh * "." * nTts[1], cluster_lookup[nTts[1]]);
        end if;
        System("rm " * activefile);
        break;
    end if;
end while;
exit;
