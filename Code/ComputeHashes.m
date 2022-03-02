/*****************************************************************************************************
This file is used for computing a hash value for a transitive group using the basic hash function
defined in Hash.m.

USAGE:

cat DATA/hash.todo | parallel -j128 magma nTt:="{1}" ComputeHashes.m

INPUT:

You should provide a variable `nTt` at the command line, such as 44T1780.

OUTPUT:

It will write the order and hash value to a file named for the transitive group in the DATA/trun/ directory,
and the time spent to a file in the DATA/trun.timings/ directory.
*****************************************************************************************************/

AttachSpec("hashspec");

SetColumns(0);
nTt := Split(nTt, " ")[1];
n, t := Explode([StringToInteger(c) : c in Split(nTt, "T")]);
G := TransitiveGroup(n, t);
t0 := Cputime();
hsh := hash(G);
PrintFile("DATA/trun.timings/" * nTt, Sprint(Cputime() - t0));
PrintFile("DATA/trun/" * nTt, Sprintf("%o.%o", Order(G), hash(G)));

exit;
