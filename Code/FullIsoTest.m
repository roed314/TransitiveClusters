/*****************************************************************************************************
This file uses Magma's IsIsomorphic function to determine how a cluster of groups is actually
divided into isomorphism classes.  It is usually faster than FullSibSplit.m, but not always
(and it does not provide sibling counts).

USAGE:

ls DATA/active | parallel -j64 --timeout 3600 --memfree 100G --joblog DATA/isotest.log magma hsh:="{1}" FullIsoTest.m

INPUT:

You should provide a variable `hsh` at the command line: the order and hash value, separated by a period.

OUTPUT:

Every time it finds a complete isomorphism class from among the rows of the input file,
it will write the result to a file in the DATA/isotest_finished directory.  When all input rows have been
assigned, it will delete the input file.  Timings are written to the DATA/isotest.timings/ folder.
*****************************************************************************************************/

AttachSpec("hashspec");

SetColumns(0);
nTt_lookup := AssociativeArray();
G_lookup := AssociativeArray();
nTts := [];
activefile := "DATA/active/" * hsh;
for cluster in Split(Read(activefile), "\n") do
    first := Split(cluster, " ")[1];
    nTt_lookup[first] := cluster;
    G_lookup[first] := StringToGroup(first);
    Append(~nTts, first);
end for;

//if #nTts gt 6 then
//    // Skip bigger clusters
//    print "Exiting since", #nTts, "clusters";
//    exit;
//end if;

while #nTts gt 0 do
    t0 := Cputime();
    label := nTts[1];
    G := G_lookup[label];
    cluster := [nTt_lookup[label]];
    nTts := nTts[2..#nTts];
    i := 1;
    while i le #nTts do
        now := nTts[i];
        H := G_lookup[now];
        print "Checking isomorphism", label, now;
        if IsIsomorphic(G, H) then
            Append(~cluster, nTt_lookup[now]);
            Remove(~nTts, i);
        else
            i +:= 1;
        end if;
    end while;
    // Since this process may get killed, we want to write output now
    print "Writing progress";
    PrintFile("DATA/isotest.timings/" * hsh, Sprintf("%o size %o(%o), %o", label, #cluster, &+[#Split(x, " ") : x in cluster], Cputime() - t0));
    // Can't use label since it might be too big in general
    if #label gt 30 then
        label := Sprint(Hash(label[1..#label div 2])) * Sprint(Hash(label[(#label div 2) + 1..#label]));
    end if;
    PrintFile("DATA/isotest_finished/" * hsh * "." * label, Join(cluster, " "));
    if #nTts gt 0 then
        tmp := "DATA/tmp/" * hsh;
        PrintFile(tmp, Join([nTt_lookup[nTt] : nTt in nTts], "\n"));
        System("mv " * tmp * " " * activefile); // mv is atomic
    else
        System("rm " * activefile);
    end if;
end while;
exit;
