/*****************************************************************************************************
This file is used for groups where the hash value is too hard to compute.  Instead, it computes the easyhash for line of input
and writes the possibilities that have the same easyhash as the first line.

USAGE:

ls Problems | parallel -j64 magma infile:="Problems/{1}" outfile:="Solutions/{1}" ComputeEasyHashes.m

INPUT:

You should provide an infile and outfile as variables at the command line

OUTPUT:

It will write a line consiting of the group description for each input description that has the same easyhash as the first line.
*****************************************************************************************************/

SetColumns(0);
AttachSpec("hashspec");

inputs := Split(Read(infile), "\n");
gps := [StringToGroup(desc) : desc in inputs];
hshs := [EasyHash(G) : G in gps];
matches := [i : i in [1..#inputs] | hshs[i] eq hshs[1]];
for i in matches do
    PrintFile(outfile, inputs[i]);
end for;
exit;
