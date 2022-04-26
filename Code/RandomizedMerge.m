/*****************************************************************************************************
This file applies a randomized strategy to find isomorphisms between groups in a cluster.
It is based on the strategy described in RandomizedTMerge.m, adapted to work for arbitrary groups.

Given a set of groups together with a bound B, we search for subgroups H_i < G_i with trivial core and index at most B.  If one is found with index m, we construct the corresponding permutation representation P_i in S_m and compute the CycleHash h, storing the triple <H_i, P_i, G_i> under the key [m, h].

When a collision is found between <H_i, G_i> and <H_j, G_j> for G_i != G_j, we can then test whether P_i is conjugate to P_j inside S_m.  If so, we've found an isomorphism.

In the current implementation, we are only interested in whether groups are isomorphic; it would be easy to actually recover the isomorphism.

USAGE:

ls DATA/active | parallel -j64 --timeout 900 magma -b hsh:="{1}" timeout:=900 RandomizedMerge.m

INPUT:

You should provide a variable `hsh` at the command line, corresponding to a file in the DATA/active/ directory.
This file should consist of lines of space separated group descriptions (as accepted by StringToGroup);
all groups on the same line should be isomorphic, and groups on different lines may or may not be isomorphic.

Since the strategry is randomized, it's common to want to rerun multiple times, yet it's also helpful to
have some indication of progress.  To this end, you must also specify `timeout` on the command line.
If the active file has already been run with this value of `timeout` it will exit immediately;
otherwise it will write the value to a corresponding file in the DATA/merge_check directory.

OUTPUT:

If the merging is completely successful and all groups in the input file are found to be isomorphic,
a file will be created in the DATA/merge_finished directory with a single line containing all group labels,
and the input file will be DELETED from DATA/active.
If not, the program will keep running until terminated (GNU parallel with a timeout is a good way to do this).
Either way, partial progress will be recorded by UPDATING the input file in DATA/active, merging
lines as appropriate.
Since the mixture of overwriting the input and automatic termination has the potential for data corruption,
this is achieved by writing to a DATA/tmp/ directory, then moving the result into place.

Timing information is output to DATA/merge.timings/
*****************************************************************************************************/

AttachSpec("hashspec");
SetColumns(0);

function RandomCorefreeSubgroup(G, H, B : max_tries:=20)
    // Try to find a subgroup of G properly containing H with index at most B in G and trivial core.
    // H should have trivial core, and G should be nontrivial.
    // May fail, in which case H is returned
    H0 := H;
    tries := 1;
    while tries le max_tries do
        g := Random(G);
        if g in H then
            continue; // H has trivial core so has index at least 2, thus this should happen rarely
        end if;
        K := sub<G|H,g>;
        if #Core(G, K) eq 1 then
            if Index(G, K) le B then
                return K; // success
            end if;
            H := K;
            tries := 0;
        else
            // try to take powers of g, since otherwise it's hard to get lower order elements
            D := Divisors(Order(g));
            for m in D[2..#D-1] do
                h := g^m;
                if h in G then
                    continue;
                end if;
                K := sub<G|H,h>;
                if #Core(G, K) eq 1 then
                    if Index(G, K) le B then
                        return K; // success
                    end if;
                    H := K;
                    tries := 0;
                    break;
                end if;
            end for;
        end if;
        tries +:= 1;
    end while;
    return H0; // failure
end function;

print hsh;
activefile := "DATA/active/" * hsh;
file_exists, ifile := OpenTest("DATA/active/" * hsh, "r");
if not file_exists then
    print "File for", hsh, "does not exist!";
    exit;
end if;
cfname := "DATA/merge_check/" * hsh;
file_exists, cfile := OpenTest(cfname, "r");
if file_exists then
    if timeout in Split(Read(cfile)) then
        exit;
    end if;
end if;
PrintFile(cfname, timeout);
N := StringToInteger(Split(hsh, ".")[1]);

descs := [];
groups := [];
for s in Split(Read(ifile), "\n") do
    Append(~descs, s);
    labels := Split(s, " ");
    // We pick a random description in case the particular representation matters (can improve on another run of the script)
    label := Random(labels);
    Append(~groups, StringToGroup(label));
end for;

active := [1..#groups]; // As we find isomorphisms, we'll choose only one i from each pair of isomorphic groups
clusters := [[i] : i in active];
whichcluster := AssociativeArray();
for i in active do
    whichcluster[i] := i;
end for;

function WriteProgress(active, clusters, i0, j0, t0)
    PrintFile("DATA/merge.timings/" * hsh, Sprintf("%o=%o in %o", descs[i0], descs[j0], Cputime() - t0));
    lines := Join([Join([descs[i] : i in C], " ") : C in clusters | #C gt 0], "\n");
    if #active eq 1 then
        PrintFile("DATA/merge_finished/" * hsh, lines);
        System("rm " * activefile);
        exit;
    else
        tmp := "DATA/tmp/" * hsh;
        PrintFile(tmp, lines);
        System("mv " * tmp * " " * activefile); // mv is atomic
    end if;
    return Cputime();
end function;

t0 := Cputime();
progress_ctr := 0;
known_corefree := AssociativeArray(); // holds known corefree subgroups as starting points
best_deg := N;
for i in active do
    G := groups[i];
    known_corefree[i] := AssociativeArray();
    known_corefree[i][N] := [sub<G|>];
    // We start by generating some corefree subgroups in an attempt to find a good bound for when to compute cycle hashes (since these will be expensive in very large degree)
    tries := 0;
    B := N;
    while tries lt 100 and (tries lt 10 or B gt 1024) do
        d := (tries mod 2 eq 0) select B else Random(Keys(known_corefree[i]));
        H := Random(known_corefree[i][d]);
        K := RandomCorefreeSubgroup(G, H, d); // any improvement on H is okay here.
        if #H ne #K then
            d := Index(G, K);
            if IsDefined(known_corefree[i], d) then
                Append(~known_corefree[i][d], K);
            else
                known_corefree[i][d] := [K];
            end if;
            B := Min(B, d);
        end if;
        tries +:= 1;
    end while;
    best_deg := Min(best_deg, B);
end for;
reps := AssociativeArray();
for i in active do
    if IsDefined(known_corefree[i], best_deg) then
        G := groups[i];
        for K in known_corefree[i][best_deg] do
            P := Image(CosetAction(G, K));
            h := CycleHash(P);
            found := false;
            if IsDefined(reps, <best_deg, h>) then
                for old in reps[<best_deg, h>] do
                    j, Q := Explode(old);
                    if best_deg lt 48 or IsConjugate(Sym(best_deg), P, Q) then
                        found := true;
                        if i ne j then
                            // A nontrivial isomorphism
                            i0, j0 := Explode([Min(i,j), Max(i,j)]);
                            clusters[i0] := Sort(clusters[i0] cat clusters[j0]);
                            clusters[j0] := [];
                            Exclude(~active, j0);
                            t0 := WriteProgress(active, clusters, i0, j0, t0);
                        end if;
                        break;
                    end if;
                end for;
                if not found then
                    Append(~reps[<best_deg, h>], <i, P>);
                end if;
            else
                reps[<best_deg, h>] := [<i, P>];
            end if;
        end for;
    end if;
end for;
PrintFile("DATA/merge.timings/" * hsh, Sprintf("Setup done in %o, finished %o", Cputime() - t0, #groups - #active));
while true do
    i := Random(active);
    G := groups[i];
    // Want to weight lower degrees more heavily, so we do the selection manually
    keys := Sort([x : x in Keys(known_corefree[i])]);
    r := Random(Integers()!(#keys * (#keys + 1) / 2 - 1));
    for i in [1..#keys] do
        r -:= i;
        if r lt 0 then
            d := keys[#keys+1 - i];
            break;
        end if;
    end for;
    H := Random(known_corefree[i][d]);
    K := RandomCorefreeSubgroup(G, H, best_deg);
    if #H ne #K then
        d := Index(G, K);
        if IsDefined(known_corefree[i], d) then
            Append(~known_corefree[i][d], K);
        // Might need to remove items here to prevent memory usage growth, but not going to worry for now.
        else
            known_corefree[i][d] := [K];
        end if;
        P := Image(CosetAction(G, K));
        h := CycleHash(P);
        found := false;
        if IsDefined(reps, <d, h>) then
            for old in reps[<d, h>] do
                j, Q := Explode(old);
                if d lt 48 or IsConjugate(Sym(d), P, Q) then
                    found := true;
                    if i ne j then
                        // A nontrivial isomorphism
                        i0, j0 := Explode([Min(i,j), Max(i,j)]);
                        clusters[i0] := Sort(clusters[i0] cat clusters[j0]);
                        clusters[j0] := [];
                        Exclude(~active, j0);
                        t0 := WriteProgress(active, clusters, i0, j0, t0);
                    end if;
                    break;
                end if;
            end for;
            if not found then
                Append(~reps[<d, h>], <i, P>);
            end if;
        else
            reps[<best_deg, h>] := [<i, P>];
        end if;
    end if;
end while;
