/*****************************************************************************************************
This file splits up the computation of Subgroups (which was too slow and memory intensive in some cases)
into multiple files.  It is aimed at five clusters that did not succumb to the other methods, each of which
had order of the form 2^k*p for p=5 or p=11, and where the desired subgroups had index 40 or 44.

In this case, it suffices to compute the index 4 or 8 subgroups of the 2-Sylow subgroup, and we can
parallelize this task by first finding the index 2 subgroups (of which there were either 511 or 4095 in practice)
and then finding index 2 or 4 subgroups of those and scanning for those with trivial core.

This file sets up the computation by finding these initial maximal subgroups, and for each one writing
a magma file that executes the rest of the job.

USAGE:

magma nTt:=44T1780 IterativeSubgroupSetup.m
ls DATA/last/ | parallel -j 96 magma {1}

INPUT:

Give a transitive group label at the command line to this file

OUTPUT:

It will write a bunch of magma files to the folder DATA/last/ which can then be run to produce output in
DATA/last.selfsib/ and DATA/last.osib/.  Output in the first case is for generators of core-free subgroups of the
original group that yield the same transitive label, while in the second case those that give a different
label.  The output in the second case give isomorphisms with other transitive groups.  Note that these scripts
do not delete files from DATA/active/ or create finalized clusters.

Timing information is printed to files in DATA/last.timings/
*****************************************************************************************************/

SetColumns(0);

n, t := Explode([StringToInteger(c) : c in Split(nTt, "T")]);
G := TransitiveGroup(n, t);
P := SylowSubgroup(G, 2);
t0 := Cputime();
Hs := [H`subgroup : H in Subgroups(P : IndexEqual := 2)];
msg := Sprintf("Computed %o index 2 subgroups of %o in %o", #Hs, nTt, Cputime() - t0);
print msg;
PrintFile("DATA/last.timings", msg);
if n eq 40 then
    m := 4;
else
    m := 2;
end if;
for i in [1..#Hs] do
    slab := nTt * "." * Sprint(i);
    F := "DATA/last/" * slab * ".m";
    PrintFile(F, "SetColumns(0);");
    PrintFile(F, Sprintf("checkfile := \"DATA/lastdone/%o\";", slab));
    PrintFile(F, "file_exists, ifile := OpenTest(checkfile, \"r\");");
    PrintFile(F, "if file_exists then");
    PrintFile(F, "    print \"Already done\";");
    PrintFile(F, "    exit;");
    PrintFile(F, "end if;");
    PrintFile(F, Sprintf("G := TransitiveGroup(%o, %o);", n, t));
    PrintFile(F, Sprintf("H := %m;", Hs[i]));
    PrintFile(F, "t0 := Cputime();");
    PrintFile(F, Sprintf("S := [K`subgroup : K in Subgroups(H : IndexEqual:=%o)];", m));
    PrintFile(F, "PrintFile(\"DATA/last.timings\", Sprintf(\"Found %o subgroups within " * slab * " in %o\", #S, Cputime() - t0));");
    PrintFile(F, "t0 := Cputime();");
    PrintFile(F, "c := 0;");
    PrintFile(F, "for K in S do");
    PrintFile(F, "    if #Core(G, K) eq 1 then");
    PrintFile(F, "        c +:= 1;");
    PrintFile(F, "        t := TransitiveGroupIdentification(Image(CosetAction(G, K)));");
    PrintFile(F, Sprintf("        if t eq %o then", t));
    PrintFile(F, "            PrintFile(\"DATA/last.selfsib/"*nTt*"\", Sprintf(\"%o\", Generators(K)));");
    PrintFile(F, "        else");
    PrintFile(F, "            PrintFile(\"DATA/last.osib/"*nTt*"\", Sprintf(\"%o\", Generators(K)));");
    PrintFile(F, "        end if;");
    PrintFile(F, "    end if;");
    PrintFile(F, "end for;");
    PrintFile(F, "PrintFile(\"DATA/last.timings\", Sprintf(\"Finished checking %o subgroups within " * slab * " in %o\", #S, Cputime() - t0));");
    PrintFile(F, "PrintFile(checkfile, Sprint(c));");
    PrintFile(F, "exit;");
end for;
exit;
