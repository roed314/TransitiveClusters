// There were several 2x2 matrix cases where IsIsomophic and RandomizedMerge.m both failed.  This script takes advantage of the structure of matrix groups to handle these cases.

SetColumns(0);
AttachSpec("hashspec");

activefile := "DATA/active/" * hsh;
file_exists, ifile := OpenTest("DATA/active/" * hsh, "r");
if not file_exists then
    print "File for", hsh, "does not exist!";
    exit;
end if;
descs := [Split(x, " ")[1] : x in Split(Read(ifile), "\n")];
groups := [StringToGroup(desc) : desc in descs];
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
                if #Image(ahom) eq N-1 then
                    print "ok";
                else
                    lift := (((N-1) div #Image(ahom)) * C.1) @@ ahom;
                    if Order(lift) eq Order(N-1)
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
            C := Complements(G, Z);
            if #C gt 0 then
                print "split";
            else
                print "nonsplit";
            end if;
        end if;
    end if;
end for;
