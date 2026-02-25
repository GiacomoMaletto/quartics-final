needsPackage "Topcom";

A = matrix {{0,1,2,3,0,1,2,0,1,0},{0,0,0,0,1,1,1,2,2,3}};

at = topcomAllTriangulations(A);

atht = new MutableHashTable;

for i from 0 to #at-1 do (
    atht#(at_i//sort/sort) = true;
);

symmetries = {
    {9,7,4,0,8,5,1,6,2,3},
    {3,6,8,9,2,5,7,1,4,0},
    {0,4,7,9,1,5,8,2,6,3},
    {3,2,1,0,6,5,4,8,7,9},
    {9,8,6,3,7,5,2,4,1,0}
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

out = openOut "three.csv";

for tr in keys(atht) do (
    if atht#tr then (
        htQQ = topcomRegularTriangulationWeights(A, tr);
        ht = htQQ * lcm (htQQ / den) / (n -> lift(n, ZZ));
        out << toString tr << toString ht << endl << flush;
    );
);

out << close