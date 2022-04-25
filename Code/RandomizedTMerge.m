/*****************************************************************************************************
This file applies a randomized strategy to find isomorphisms between transitive groups in a cluster.
As an example, 40T194505 and 40T195727 have the same hash (in fact, they end up abstractly isomorphic)
If we compute a random index-40 core-free subgroup H of 40T194505, we can use TransitiveGroupIdentification
on the induced permutation representation.  Some of the time, we'll get 195727 out, thus exhibiting an
isomorphism between the two permutation groups.  This randomized strategy is often much faster than
falling back on IsIsomorphic.

We will iteratively build up an index-40 core-free subgroup as follows.  At each stage, H is the
subgroup under construction, N is its normalizer and G the ambient group.  We start with H = {1}.
1. If H != N, choose a random g in N and set K = <H, g>.
   a) If the index of K is not divisible by 40, replace g with an appropriate power and set K = <H, g>
   b) While [N_G(K) : K] is prime, set K = N_G(K) (since this is the only way to enlarge K within its normalizer)
   c) If K has trivial core, set H = K and proceed; otherwise choose another random g.
2. If H = N or step 1 has failed too many times, instead just choose random g in G until <H, g> is core-free
   with index divisible by 40.

The current version can probably be improved; it failed to find isomorphisms in some cases that should
have been possible:
40T271235 = 40T272233
40T271277 = 40T272110
44T1780 = 44T1783
40T308155 = 40T308167
40T308164 = 40T308166

USAGE:

ls DATA/active | parallel -j96 --timeout 900 magma hsh:="{1}" timeout:=900 RandomizedTMerge.m

INPUT:

You should provide a variable `hsh` at the command line, corresponding to a file in the DATA/active/ directory.
This file should consist of lines of space separated transitive groups (such as 44T1780);
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

SMALL_TRIES := 40;
function SmallJump(G, H, N, M)
    // G -- ambient group
    // H -- current subgroup
    // N -- the normalizer of H
    // M -- the desired order of a subgroup of G
    // output -- a corefree subgroup K, H < K < G, with the order of K dividing M
    //           if no such subgroup found, returns H
    //        -- the normalizer of K
    for ctr in [1..SMALL_TRIES] do
        g := Random(N);
        K := sub<G|H,g>;
        if GCD(#K, M) ne #H then
            toofar := #K div GCD(M, #K);
            if toofar gt 1 then
                g := g^toofar;
                K := sub<G|H,g>;
            end if;
            // Expand K repeatedly when forced:
            // if the index of K in its normalizer is prime
            N1 := Normalizer(G, K);
            C := Core(G, K);
            while IsPrime(Index(N1, K)) do
                if #K eq M and #C eq 1 then
                    return K, N1;
                end if;
                K := N1;
                N1 := Normalizer(G, K);
                C := Core(G, K);
            end while;
            if #C eq 1 then
                //print "Small jump successful, #K", #K, "#N", #N1;
                return K, N1;
            end if;
        end if;
    end for;
    //print "Small jump failed";
    return H, N;
end function;

BIG_TRIES := 400;
function BigJump(G, H, N, M)
    // G -- ambient group
    // H -- current subgroup
    // N -- the normalizer of H
    // M -- the desired order of a subgroup of G
    // output -- a corefree subgroup K, H < K < G, with the order of K dividing M
    //           if no such subgroup found, returns H
    //        -- the normalizer of K
    for ctr in [1..BIG_TRIES] do
        g := Random(G);
        K := sub<G|H,g>;
        if #K gt #H and IsDivisibleBy(M, #K) and #Core(G, K) eq 1 then
            N := Normalizer(G, K);
            //print "Big jump successful, #K", #K, "#N", #N;
            return K, N;
        end if;
    end for;
    //print "Big jump failed";
    return H, N;
end function;

function RandomCorelessSubgroup(G, m : max_tries:=40)
    // If m was chosen incorrectly, there may be no coreless subgroups of index m
    // Given that this is just part of a loop below where we try different G and m,
    // we just give up after a certain number of tries
    assert IsDivisibleBy(#G, m);
    M := #G div m;
    tries := 0;
    while true do
        // start at the trivial group
        H := sub<G|>;
        N := G;
        while true do
            K, N := SmallJump(G, H, N, M);
            if #K eq #H then
                // small jump failed
                K, N := BigJump(G, H, N, M);
                if #K eq #H then
                    //big jump also failed, so start anew
                    //print "Restarting";
                    tries +:= 1;
                    if tries ge max_tries then
                        //printf "Giving up after %o tries\n", tries;
                        return -1, tries;
                    end if;
                    break;
                end if;
            end if;
            if #K eq M then
                //printf "Done after %o tries\n", tries;
                return K, tries;
            end if;
            H := K;
        end while;
    end while;
end function;

SetColumns(0);

print hsh;
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
lookup := AssociativeArray();
groups := [];
degrees := [];
i := 1;
for s in Split(Read(ifile)) do
    labels := Split(s, " ");
    nTts := [[StringToInteger(m) : m in Split(label, "T")] : label in labels];
    Append(~groups, nTts);
    Append(~degrees, {nTt[1] : nTt in nTts});
    for nTt in nTts do
        lookup[nTt] := i;
    end for;
    i +:= 1;
end for;
// We actually need the degrees of the OTHER groups
degrees := [&join[degrees[j] : j in [1..#degrees] | j ne i] : i in [1..#degrees]];

active := [1..#groups]; // As we find isomorphisms, we'll choose only one i from each pair of isomorphic groups
t0 := Cputime();
progress_ctr := 0;
max_tries := 40;
failures := 0;
total_restarts := 0;
while true do
    progress_ctr +:= 1;
    i := Random(active);
    n, t := Explode(groups[i][1]);
    G := TransitiveGroup(n, t);
    m := Random(degrees[i]);
    H, restarts := RandomCorelessSubgroup(G, m : max_tries:=max_tries);
    total_restarts +:= restarts;
    if H cmpeq -1 then
        failures +:= 1;
        continue;
    end if;
    s := TransitiveGroupIdentification(Image(CosetAction(G, H)));
    j := lookup[[m, s]];
    if i ne j then
        i, j := Explode([Min(i,j), Max(i,j)]);
        PrintFile("DATA/merge.timings/" * hsh, Sprintf("%oT%o=%oT%o %o %o %o %o", groups[i][1][1], groups[i][1][2], groups[j][1][1], groups[j][1][2], progress_ctr, failures, total_restarts, Cputime() - t0));
        groups[i] := Sort(groups[i] cat groups[j]);
        groups[j] := [];
        for k -> v in lookup do
            if v eq j then
                lookup[k] := i;
            end if;
        end for;
        Exclude(~active, j);
        lines := Join([Join([Sprintf("%oT%o", tgp[1], tgp[2]) : tgp in groups[k]], " ") : k in active], "\n");
        activefile := "DATA/active/" * hsh;
        if #active eq 1 then
            PrintFile("DATA/merge_finished/" * hsh, lines);
            System("rm " * activefile);
            print "Finished!";
            break;
        else
            tmp := "DATA/tmp/" * hsh;
            PrintFile(tmp, lines);
            System("mv " * tmp * " " * activefile); // mv is atomic
        end if;
        printf "Isomorphism found after %o tries and %o seconds, %o clusters remaining\n", progress_ctr, Cputime() - t0, #active;
        progress_ctr := 0;
        t0 := Cputime();
        failures +:= 0;
        total_restarts := 0;
        max_tries := 40;
    elif progress_ctr mod 10 eq 0 then
        printf "Not yet successful; on loop %o\n", progress_ctr;
        max_tries +:= 1;
    end if;
end while;

exit;
