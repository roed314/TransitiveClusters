AttachSpec("hashspec");

// cat DATA/hashclusters/hash.todo | parallel -j128 magma nTt:="{1}" AddTHashes2.m

SetColumns(0);
// Split off the sibling data
nTt := Split(nTt, " ")[1];
n, t := Explode([StringToInteger(c) : c in Split(nTt, "T")]);
G := TransitiveGroup(n, t);
t0 := Cputime();
hsh := hash(G);
PrintFile("DATA/hashclusters/trun.times/" * nTt, Sprint(Cputime() - t0));
PrintFile("DATA/hashclusters/trun/" * nTt, Sprintf("%o.%o", Order(G), hash(G)));

exit;
