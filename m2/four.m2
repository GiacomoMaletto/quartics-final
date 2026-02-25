needsPackage "Topcom";

A = matrix{
    {0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0},
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4}};

at = topcomAllTriangulations(A);

atht = new MutableHashTable;

for i from 0 to #at-1 do (
    atht#(at_i//sort/sort) = true;
);

symmetries = {
    {0, 5, 9, 12, 14, 1, 6, 10, 13, 2, 7, 11, 3, 8, 4},
    {4, 3, 2, 1, 0, 8, 7, 6, 5, 11, 10, 9, 13, 12, 14},
    {4, 8, 11, 13, 14, 3, 7, 10, 12, 2, 6, 9, 1, 5, 0},
    {14, 12, 9, 5, 0, 13, 10, 6, 1, 11, 7, 2, 8, 3, 4},
    {14, 13, 11, 8, 4, 12, 10, 7, 3, 9, 6, 2, 5, 1, 0}
};

for i from 0 to #at-1 do (
    tr = at_i;
    if atht#tr then (
        symms = apply(symmetries, s -> sort(tr / (t -> sort (s_t)))) - set {tr};
        for s in symms do (
            atht#s = false;
        );
    );
);

den = method();
den (ZZ) := n -> 1;
den (QQ) := q -> if q == 0 then 1 else denominator q;

out = openOut "four.csv";

for tr in keys(atht) do (
    if atht#tr then (
        htQQ = topcomRegularTriangulationWeights(A, tr);
        ht = htQQ * lcm (htQQ / den) / (n -> lift(n, ZZ));
        out << toString tr << toString ht << endl << flush;
    );
);

out << close