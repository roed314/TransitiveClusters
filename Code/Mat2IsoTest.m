// There were several 2x2 matrix cases where IsIsomophic and RandomizedMerge.m both failed.  This script takes advantage of the structure of matrix groups to handle some of these cases.

SetColumns(0);
AttachSpec("hashspec");

N := StringToInteger(Split(hsh, ".")[1]);
activefile := "DATA/active/" * hsh;
file_exists, ifile := OpenTest("DATA/active/" * hsh, "r");
if not file_exists then
    print "File for", hsh, "does not exist!";
    exit;
end if;
descs := Split(Read(ifile), "\n");
groups := [StringToGroup(Split(desc, " ")[1]) : desc in descs];
if &and[IsAbelian(G) : G in groups] then
    // all abelian, so we can just use abelian invariants
    // The hash should be enough to distinguish in this case, but we double check
    invs := {AbelianInvariants(G) : G in groups};
    assert #invs eq 1;
    PrintFile("DATA/mat2_finished/" * hsh, Join(descs, " "));
    System("rm " * activefile);
    exit;
end if;
ps := {};
if &and[Type(G) eq GrpMat : G in groups] then
    ps := {@ #CoefficientRing(G) : G in groups @};
end if;
if &and[Type(G) eq GrpMat and Degree(G) eq 2 and #ps eq 1 and IsPrime(ps[1]) : G in groups] then
    p := ps[1];
    k := GF(p);
    allborel := true;
    As := [];
    for G in groups do
        Agens := [];
        m1101 := G![1,1,0,1];
        found := false;
        for g in Generators(G) do
            if g eq m1101 then
                found := true;
            else
                if Eltseq(g)[2] ne 0 or Eltseq(g)[3] ne 0 then
                    found := false;
                    break;
                end if;
                Append(~Agens, g);
            end if;
        end for;
        if not found then
            allborel := false;
            break;
        end if;
        Append(~As, sub<G|Agens>);
    end for;
    if allborel then
        // Let T be the diagonal subgroup of G; G = <m1101> \rtimes T, and K the centralizer of m1101 in T.
        // I claim that it suffices to check that the Ks are direct summands of the Ts and are all isomorphic.
        // Since p is prime, Aut(<m1101>) is cyclic of order p-1, and G is a semidirect product of <m1101> and T.
        // This semidirect product structure is determined by the conjugation homomorphism T -> Aut(<m1101>).
        // If K is a semidirect summand, then T is isomorphic to the direct product of K and the cyclic image T
        // in Aut(<m1101>).  If the Ks are isomorphic, the cyclic complements automatically are by counting.
        Ks := [Centralizer(A, m1101) : A in As];
        invs := {AbelianInvariants(K) : K in Ks};
        // If K has index p-1 in A then the sequence is automatically a direct summand since T has exponent dividing p-1.
        if #invs eq 1 and (Index(As[1], Ks[1]) eq p-1 or &and[#Complements(As[i], Ks[i]) gt 0 : i in [1..#As]]) then
            // Success!
            PrintFile("DATA/mat2_finished/" * hsh, Join(descs, " "));
            System("rm " * activefile);
            exit;
        end if;
    end if;
end if;
P := PrimeDivisors(N);
if #P gt 1 and &and[IsNormal(G, SylowSubgroup(G, p)) : G in groups, p in P] then
    // direct product of Sylow subgroups, so we try using IsIsomorphic on the Sylow subgroups
    active := [1..#groups];
    while #active gt 0 do
        t0 := Cputime();
        cluster := [descs[active[1]]];
        label := Split(descs[active[1]], " ")[1];
        G := groups[active[1]];
        SylG := [SylowSubgroup(G, p) : p in P];
        active := active[2..#active];
        m := 1;
        while m le #active do
            j := active[m];
            H := groups[j];
            SylH := [SylowSubgroup(H, p) : p in P];
            if &and[IsIsomorphic(SylG[k], SylH[k]) : k in [1..#P]] then
                Append(~cluster, descs[j]);
                Remove(~active, m);
            else
                m +:= 1;
            end if;
        end while;
        // This process may get killed, so we write output
        print "Writing progress";
        PrintFile("DATA/sylow.timings/" * hsh, Sprintf("%o size %o(%o), %o", label, #cluster, &+[#Split(x, " ") : x in cluster], Cputime() - t0));
        // Can't use label since it might be too big in general
        if #label gt 30 then
            label := Sprint(Hash(label[1..#label div 2])) * Sprint(Hash(label[(#label div 2) + 1..#label]));
        end if;
        PrintFile("DATA/mat2_finished/" * hsh * "." * label, Join(cluster, " "));
        if #active gt 0 then
            tmp := "DATA/tmp/" * hsh;
            PrintFile(tmp, Join([descs[j] : j in active], "\n"));
            System("mv " * tmp * " " * activefile); // mv is atomic
        else
            System("rm " * activefile);
        end if;
    end while;
    exit;
end if;
if &and[Type(G) eq GrpMat : G in groups] then
    t0 := Cputime();
    Zs := [G meet Center(GL(Degree(G), CoefficientRing(G))) : G in groups];
    zids := {IdentifyGroup(Z) : Z in Zs};
    Cs := [Complements(groups[i], Zs[i]) : i in [1..#groups]];
    if #zids eq 1 and &and[#C gt 0 : C in Cs] then
        // Each G is the direct product of Z and G/Z, and all the Zs are isomorphic
        active := [1..#groups];
        Qs := [groups[i] / Zs[i] : i in [1..#groups]];
        while #active gt 0 do
            cluster := [descs[active[1]]];
            label := Split(descs[active[1]], " ")[1];
            G := Qs[active[1]];
            active := active[2..#active];
            m := 1;
            while m le #active do
                j := active[m];
                H := Qs[j];
                if IsIsomorphic(G, H) then
                    Append(~cluster, descs[j]);
                    Remove(~active, m);
                else
                    m +:= 1;
                end if;
            end while;
            // This process may get killed, so we write output
            print "Writing progress";
            PrintFile("DATA/proj.timings/" * hsh, Sprintf("%o size %o(%o), %o", label, #cluster, &+[#Split(x, " ") : x in cluster], Cputime() - t0));
            // Can't use label since it might be too big in general
            if #label gt 30 then
                label := Sprint(Hash(label[1..#label div 2])) * Sprint(Hash(label[(#label div 2) + 1..#label]));
            end if;
            PrintFile("DATA/mat2_finished/" * hsh * "." * label, Join(cluster, " "));
            if #active gt 0 then
                tmp := "DATA/tmp/" * hsh;
                PrintFile(tmp, Join([descs[j] : j in active], "\n"));
                System("mv " * tmp * " " * activefile); // mv is atomic
            else
                System("rm " * activefile);
            end if;
        end while;
        exit;
    end if;
end if;

// None of the above solutions worked, so we print something to describe which case we're in.
typs := [];
for row in descs do
    rowout := [];
    for desc in Split(row, " ") do
        if "Mat" in desc then
            Append(~rowout, Split(desc, "Mat")[1] * "Mat");
        else
            Append(~rowout, Split(desc, "Perm")[1] * "Perm");
        end if;
    end for;
    Append(~typs, rowout);
end for;
if &and[#r eq 1 : r in typs] then
    typs := [r[1] : r in typs];
end if;
print "Failure:", hsh, typs;
for i in [1..#groups] do
    G := groups[i];
    if Type(G) eq GrpMat and Degree(G) eq 2 then
        print hsh, i;
        R := CoefficientRing(G);
        N := #R;
        if IsAbelian(G) then
            print "abelian";
        elif IsPrime(N) then
            k := GF(N);
            Agens := [];
            found := false;
            m1101 := G![1,1,0,1];
            for g in Generators(G) do
                if g eq m1101 then
                    found := true;
                else
                    Append(~Agens, g);
                end if;
            end for;
            A := sub<G|Agens>;
            if found and IsAbelian(A) then
                // Borel
                C := AbelianGroup([N-1]);
                ahom := hom<A -> C | [<a, Log(k!Eltseq(m1101^a)[2])*C.1> : a in Generators(A)]>;
                iord := #Image(ahom);
                if iord eq N-1 then
                    print "ok";
                else
                    K := Kernel(ahom);
                    comps := Complements(A, K);
                    //lift := (((N-1) div iord) * C.1) @@ ahom;
                    if #comps gt 0 then
                       print "split Borel";
                    else
                        print "nonsplit Borel";
                    end if;
                end if;
            else
                print "not Borel";
            end if;
        else
            Z := G meet Center(GL(2, Integers(N)));
            Append(~Zids, IdentifyGroup(Z));
            C := Complements(G, Z);
            if #C gt 0 then
                print "projective quotient splits (thus direct product)";
            else
                print "nonsplit";
                /*p, e := Explode(Factorization(N)[1]);
                desc := descs[i];
                amb, gens := Split(desc, "Mat");
                red := StringToGroup(Sprintf("2,%oMat%o", p^(e-1), gens));
                rhom := hom<G -> red | [<g, red!Eltseq(g)> : g in Generators(G)]>;*/
            end if;
        end if;
    end if;
end for;
