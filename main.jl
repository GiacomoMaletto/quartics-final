using Test
using Graphs

include("dyck.jl")
include("my_graphs.jl")

@inline function setindexsymmetric!(S, i, j, v)
    S[i, j] = v
    S[j, i] = v
end

@inline function remove_consecutive(v)
    return [v[1]; v[2:end][v[2:end].!=v[1:end-1]]]
end

@inline function reversebool(b::Bool, A)
    return b ? reverse(A) : A
end

@inline function product1d(x...)
    rect = Iterators.product(x...)
    return reshape(collect(rect), length(rect))
end

# Combinatorial cellular complex
# d: vector indicating the degree of each curve (0 if it's just a curve)
# nV: number of vertices
# E: vector of tuples (s, d, i) where s is the source, d the destination, i the type of each edge
# F: vector of vectors of tuples (s, e) where e is the edge and s the sign (1==true: opposite orientation)
struct CellComplex
    d::Vector{Int}
    nV::Int
    E::Vector{@NamedTuple{s::Int, d::Int, i::Int}}
    F::Vector{Vector{@NamedTuple{s::Bool, e::Int}}}
end

Base.:(==)(cc1::CellComplex, cc2::CellComplex) = cc1.d == cc2.d && cc1.nV == cc2.nV && cc1.E == cc2.E && cc1.F == cc2.F
Base.hash(cc::CellComplex) = hash((cc.d, cc.nV, cc.E, cc.F))

function source(e)
    return e.s
end

function source(s, e)
    return s ? e.d : e.s
end

function destination(e)
    return e.d
end

function destination(s, e)
    return s ? e.s : e.d
end

function iscellcomplex(cc::CellComplex)::Bool
    for (e, (_, _, i)) in enumerate(cc.E)
        if !(i in eachindex(cc.d) || i == 0)
            return false
        end
        if !(sum([count(se -> se.e == e, F) for F in cc.F]) == 2)
            return false
        end
    end

    for f in cc.F
        for ((s1, e1), (s2, e2)) in zip(f, circshift(f, -1))
            if destination(s1, cc.E[e1]) != source(s2, cc.E[e2])
                return false
            end
        end
    end

    return true
end

# Combinatorial arrangement
# (n, W, T) data of a curve D of degree d intersecting transversely the cellular complex cc
# n[e]: the number of intersections of C with the edge e
# W[f]: the Dyck word of the face f
# T: vector of tuples (g, T) where T is a rooted tree through a Dyck word, with ground g
struct NWT
    cc::CellComplex
    d::Int
    n::Vector{Int}
    W::Vector{Vector{Bool}}
    T::Vector{@NamedTuple{g::Int, T::Vector{Bool}}}
end

Base.:(==)(nwt1::NWT, nwt2::NWT) = nwt1.cc == nwt2.cc && nwt1.d == nwt2.d && nwt1.n == nwt2.n && nwt1.W == nwt2.W && nwt1.T == nwt2.T
Base.hash(nwt::NWT) = hash((nwt.cc, nwt.d, nwt.n, nwt.W, nwt.T))

function f_ground_sequence(Ws::Vector{Vector{Bool}})::Vector{Vector{Int}}
    i = 1
    fgs = []
    for W in Ws
        push!(fgs, dyck_sequence(W, i))
        i += div(length(W), 2) + 1
    end
    return fgs
end

@inline function f_ground_sequence(nwt::NWT)::Vector{Vector{Int}}
    return f_ground_sequence(nwt.W)
end

# output: figs[(f, i)] returns the ground sequence at the i-th edge e of f,
# with the same orientation as e
function fi_ground_sequence(nwt::NWT)::Dict{Tuple{Int,Int},Vector{Int}}
    fgs = f_ground_sequence(nwt)

    figs = Dict()
    for (f, F) in enumerate(nwt.cc.F)
        i = 1
        for (j, (s, e)) in enumerate(F)
            figs[(f, j)] = reversebool(s, fgs[f][i:(i+nwt.n[e])])
            i += nwt.n[e]
        end
    end

    return figs
end

# compute the groundmatrix of a nwt
function groundmatrix(nwt)::Vector{Matrix{Bool}}
    cc = nwt.cc

    fgs = f_ground_sequence(nwt)
    ng = length(unique(vcat(fgs...)))
    gmat = [zeros(Bool, ng, ng) for _ in 1:(length(cc.d)+1)]

    for seq in fgs
        for (g1, g2) in zip(seq[1:end-1], seq[2:end])
            setindexsymmetric!(gmat[end], g1, g2, 1)
        end
    end

    efi = [Tuple{Int,Int}[] for _ in eachindex(cc.E)]
    for (f, F) in enumerate(cc.F)
        for (i, (s, e)) in enumerate(F)
            push!(efi[e], (f, i))
        end
    end
    @assert(all(length.(efi) .== 2))

    figs = fi_ground_sequence(nwt)

    for e in eachindex(nwt.n)
        ((f1, i1), (f2, i2)) = efi[e]
        for (g1, g2) in zip(figs[(f1, i1)], figs[(f2, i2)])
            setindexsymmetric!(gmat[cc.E[e].i], g1, g2, 1)
        end
    end

    return gmat
end

function floatingforest(nwt)
    ng = length(unique(vcat(f_ground_sequence(nwt)...)))
    graph = SimpleGraph(ng)

    vss = Vector{Int}[]
    for (g, T) in nwt.T
        rt, root = dyck_to_rooted_tree(T)
        push!(vss, [g, nv(graph) + root])
        graph = my_union(graph, rt)
    end

    graph, _ = my_merge_vertices_list(graph, vss)
    return graph
end

function isNWT(nwt::NWT)::Bool
    cc = nwt.cc
    if !iscellcomplex(cc)
        return false
    end

    if length(nwt.n) != length(cc.E)
        return false
    end

    for (i, d) in enumerate(cc.d)
        if d != 0 && nwt.d != 0
            es = findall(E -> E.i == i, cc.E)
            n_t = sum(nwt.n[es])
            if n_t > d * nwt.d || (d * nwt.d - n_t) % 2 != 0
                return false
            end
        end
    end

    for (f, F) in enumerate(cc.F)
        es = getproperty.(F, :e)
        n_f = sum(nwt.n[es])
        if n_f % 2 != 0 || n_f != length(nwt.W[f])
            return false
        end
    end

    if !issorted(getproperty.(nwt.T, :g))
        return false
    end
    gs = unique(vcat(f_ground_sequence(nwt)...))
    for (g, T) in nwt.T
        if !(g in gs)
            return false
        end
        if T != forest_to_dyck(dyck_to_forest(T))
            return false
        end
    end

    return true
end

# input:
# - in_nwt: the input nwt, where in_cc = in_nwt.cc
# - out_cc: cc of the output; must have the same vertices, a subset of the edges, and the merged faces
# - Imap: if an edge has type i in out_cc, it has type Imap in in_cc
# - Emap: edge e in out_cc.E is edge Emap[e] in in_cc.E
# - Fmap: the i-th border of the face f in out_cc.F is the in_i-th border of the face in_f in in_cc.F, where Fmap[f][i]=(in_f,in_i)
# output:
# - out_nwt: obtained by removing some edges from in_nwt
function forget_edges(in_nwt::NWT, out_cc::CellComplex, Imap::Vector{Int}, Emap::Vector{Int}, Fmap::Vector{Vector{Tuple{Int,Int}}})::NWT
    in_cc = in_nwt.cc

    (in_cc.nV == out_cc.nV) || error("Not the same number of vertices.")

    (length(out_cc.E) == length(Emap) == length(unique(Emap))) || error("Not a subset of the edges.")
    for (e, E) in enumerate(out_cc.E)
        in_E = in_cc.E[Emap[e]]
        (E.s == in_E.s && E.d == in_E.d && Imap[E.i] == in_E.i) || error("Edges do not match.")
    end

    for (f, F) in enumerate(out_cc.F)
        (length(Fmap[f]) == length(F)) || error("Not enough borders.")
    end

    in_figs = fi_ground_sequence(in_nwt)

    gmat = groundmatrix(in_nwt)

    removed_i = setdiff(eachindex(in_cc.d), Imap)
    removed_matrix = reduce((m1, m2) -> m1 .|| m2, gmat[removed_i])
    removed_components = connected_components(SimpleGraph(removed_matrix))
    removed_map = Dict{Int,Int}()
    for (i, c) in enumerate(removed_components)
        for g in c
            removed_map[g] = i
        end
    end

    merged_fgs = Vector{Int}[]
    for f in eachindex(out_cc.F)
        push!(merged_fgs, [])
        for ((in_f, in_i), (out_s, out_e)) in zip(Fmap[f], out_cc.F[f])
            (in_cc.F[in_f][in_i].e == Emap[out_e]) || error("Fmap format not correct.")
            append!(merged_fgs[end], [removed_map[g] for g in reversebool(out_s, in_figs[(in_f, in_i)])])
        end
    end
    merged_fgs = remove_consecutive.(merged_fgs)
    out_D = sequence_dyck.(merged_fgs)

    merged_gs = vcat(merged_fgs...)
    out_gs = vcat(f_ground_sequence(out_D)...)
    gmap = Dict(zip(out_gs, merged_gs))

    unused_ground = [g for g in 1:size(gmat[1], 1) if !(removed_map[g] in merged_gs)]

    flforest = floatingforest(in_nwt)
    for ug in unused_ground
        for neigh in findall(gmat[end][ug, :])
            add_edge!(flforest, ug, neigh)
        end
    end
    flforest, vmap = my_merge_vertices_list(flforest, removed_components)

    out_T = []
    for (ng, mg) in gmap
        D = rooted_tree_to_dyck(flforest, vmap[removed_components[mg][1]])
        if !isempty(D)
            push!(out_T, (g=ng, T=D))
        end
    end

    return NWT(out_cc, in_nwt.d, in_nwt.n[Emap], out_D, out_T)
end

# input:
# - in_nwt: the input nwt, where in_cc = in_nwt.cc
# - out_cc: cc of the output; must have the a subset of the vertices, the merged edges, and the same faces
# - Vmap: vertex v in out_cc is vertex Vmap[v] in in_cc
# - Emap: edge e in out_cc is made by merging edges Emap[e] in in_cc
# output:
# - out_nwt: obtained by removing some vertices from in_nwt
function forget_vertices(in_nwt::NWT, out_cc::CellComplex, Vmap::Vector{Int}, Emap::Vector{Vector{Int}})
    in_cc = in_nwt.cc

    (length(out_cc.nV) == length(Vmap) == length(unique(Vmap))) || error("Not a subset of the vertices.")
    for v in setdiff(1:in_cc.nV, Vmap)
        n = count(==(v), getproperty.(in_cc.E, :s)) + count(==(v), getproperty.(in_cc.E, :d))
        (n in [0, 2]) || error("Vertex not occurring 0 or 2 times.")
    end

    (in_cc.d == out_cc.d) || error("Not the same number of types.")

    out_E = [sum(in_nwt.n[Emap[e]]) for e in eachindex(out_cc.E)]

    return NWT(out_cc, in_nwt.d, out_E, in_nwt.W, in_nwt.T)
end

# the ground faces are indexed as in the indexing used in nwt.T
function groundcc(in_nwt::NWT)::CellComplex
    in_cc = in_nwt.cc

    outd = [in_cc.d; in_nwt.d]

    outnV = in_cc.nV
    outE = @NamedTuple{s::Int, d::Int, i::Int}[]
    ine_oute = Vector{Int}[]
    for (ine, inE) in enumerate(in_cc.E)
        ine_outV = [inE.s; (outnV+1):(outnV+in_nwt.n[ine]); inE.d]
        outnV += in_nwt.n[ine]

        push!(ine_oute, [])
        for (v1, v2) in zip(ine_outV[1:end-1], ine_outV[2:end])
            push!(outE, (v1, v2, inE.i))
            push!(ine_oute[ine], length(outE))
        end
    end

    outF = Vector{@NamedTuple{s::Bool, e::Int}}[]
    for (inf, inF) in enumerate(in_cc.F)
        B = Vector{@NamedTuple{s::Bool, e::Int}}[[]]
        for (ins, ine) in inF
            push!(B[end], (ins, ine_oute[ine][1]))
            for oute in ine_oute[ine][2:end]
                push!(B, [(ins, oute)])
            end
        end
        pushfirst!(B[1], pop!(B)[1])

        W = in_nwt.W[inf]
        @assert(length(B) == length(W))
        mabel = dyck_mabel(W, length(outE) + 1)
        for m in unique(mabel)
            (i, j) = findall(==(m), mabel)
            push!(outE, (outE[B[i][end].e].d, outE[B[j][end].e].d, length(outd)))
        end

        start = collect(eachindex(W))
        while !isempty(start)
            push!(outF, [])
            level = 0
            walker = start[1]
            while walker <= length(W)
                if level == 0
                    setdiff!(start, walker)
                    append!(outF[end], B[walker])
                    push!(outF[end], (!W[walker], mabel[walker]))
                    if !W[walker]
                        break
                    end
                end
                level += W[walker] ? 1 : (-1)
                walker += 1
            end
        end
    end

    return CellComplex(outd, outnV, outE, outF)
end

# Let cc be the cellular complex on which nwt is based;
# Suppose cc is the projective plane and iline is the index of a line C_iline on cc;
# then curve_type(nwt, iline) is the topological type of the curve nwt
function curve_type(nwt::NWT, iline::Int)::Vector{Bool}

end

# Generate all NWT with specified cc, d, n and no trees
function allW(cc::CellComplex, d::Int, n::Vector{Int})
    fn = [sum(n[getproperty(F, :e)]) for F in cc.F]
    (fn .% 2 .== 0) || error("Does not add up to an even number.")

    for W in product1d(all_dyck_words.(fn))
        push!(result, NWT(cc, d, n, W, []))
    end
    return result
end

# given a CellComplex cc and a closed fpath fp,
# returns the possible new arrangements
function add_line(incc::CellComplex, fp::Vector{Int})::Vector{NWT}
    result = []
    # assume fp = [f1, ..., fk, f1]
    # then edgepaths is a vector whose elements are tuples (e12, e23, ..., ek1)
    # where ei{i+1} is an edge between fi and f{i+1}
    inFe = [getproperty.(F, :e) for F in eachindex(incc.F)]

    for epath in product1d([intersect(inFe[f1], inFe[f2]) for (f1, f2) in zip(fp[1:end-1], fp[2:end])]...)
        n = [0 for _ in incc.E]
        for e in epath
            n[e] += 1
        end
        for nwt in allW(incc, 1, n)
            if curve_type(nwt) == []
                push!(result, nwt)
            end
        end
    end

    return result

    #     # newV are the labels of the vertices of the new arrangement
    #     # and can be either the original vertices 1, ..., nv
    #     # or new vertices (e, i) where e is the edge over which the new vertex lie
    #     # and i is an increasing index;
    #     # newEV keeps track of the indices of the new vertices lying over the old edges,
    #     # so that newEV[e] = [i1, ..., ik] and newV[ij] = (e, j)
    #     newV = Vector{Union{Int,Tuple{Int,Int}}}(collect(1:incc.nV))
    #     newEV = [Int[] for _ in eachindex(cc.E)]
    #     for (e, k) in counter(epath)
    #         push!(newEV[e], ((length(newV)+1):(length(newV)+k))...)
    #         push!(newV, [(e, i) for i in 1:k]...)
    #     end

    #     # in this way, we get
    #     # path2     =   [f1,  f2,  f3, ..., f{k-1}, fk]
    #     # edgepath2 = [ek1, e12, e23, ...,   e{k-1}k, ek1]
    #     path2 = fp[1:end-1]
    #     epath2 = [epath[end], epath...]

    #     vmatchings = [Vector{Tuple{Int,Int}}[] for _ in eachindex(incc.F)]
    #     for (f, F) in enumerate(incc.F)
    #         p = [sort((epath2[i], epath2[i+1])) for i in findall(==(f), path2)]
    #         vs = [reversebool(s, newEV[e]) for (s, e) in F]
    #         vmatchings[f] = generate_matchings(inFe[f], p, vs)
    #     end
    # end

end

### TESTS ###

three_lines = CellComplex(
    [1, 1, 1],
    3,
    [(2, 3, 1), (3, 1, 2), (1, 2, 3), (2, 3, 1), (3, 1, 2), (1, 2, 3)],
    [[(0, 1), (0, 2), (0, 3)], [(0, 1), (0, 5), (0, 6)], [(0, 4), (0, 2), (0, 6)], [(0, 4), (0, 5), (0, 3)]])

two_lines_three_points = CellComplex(
    [1, 1],
    3,
    [(2, 3, 1), (1, 2, 2), (2, 3, 1), (1, 2, 2)],
    [[(0, 1), (1, 3), (1, 4), (0, 2)], [(0, 1), (1, 3), (1, 2), (0, 4)]])

two_lines = CellComplex(
    [1, 1],
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(0, 1), (0, 2)], [(0, 1), (1, 2)]]
)

one_line = CellComplex(
    [1],
    1,
    [(1, 1, 1)],
    [[(0, 1), (0, 1)]]
)

zero_lines = CellComplex(
    [],
    0,
    [],
    [[]]
)

fine = NWT(
    three_lines,
    4,
    [4, 1, 3, 0, 1, 1],
    [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 0], [1, 0], [1, 1, 0, 0]],
    [(1, [1, 0, 1, 0]), (10, [1, 0])])

nonviro = NWT(
    three_lines,
    4,
    [3, 4, 1, 1, 0, 3],
    [[1, 0, 1, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 1, 0, 0], [1, 0]],
    [])

marge = NWT(
    three_lines,
    0,
    [4, 7, 3, 2, 1, 3],
    [[1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 0, 0], [1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0], [1, 0, 1, 1, 0, 0]],
    [(1, [1, 0]), (5, [1, 0]), (17, [1, 0]), (18, [1, 0])]
)

@testset "CellComplex" begin
    @test iscellcomplex(three_lines)
    @test iscellcomplex(two_lines)
    @test iscellcomplex(one_line)
    @test iscellcomplex(zero_lines)
end

@testset "isNWT" begin
    @test isNWT(fine)
    @test isNWT(nonviro)
    @test isNWT(marge)
end

@testset "remove_edges" begin
    marge_e = NWT(
        two_lines_three_points,
        0,
        [4, 3, 2, 3],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_ev = NWT(
        two_lines,
        0,
        [6, 6],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_eve = NWT(
        one_line,
        0,
        [6],
        [[1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    fine_e = NWT(
        two_lines_three_points,
        4,
        [4, 3, 0, 1],
        [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 1, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 0])]
    )

    @test forget_edges(marge, two_lines_three_points, [1, 3], [1, 3, 4, 6], [[(1, 1), (3, 1), (3, 3), (1, 3)], [(2, 1), (4, 1), (4, 3), (2, 3)]]) == marge_e
    @test forget_vertices(marge_e, two_lines, [2], [[1, 3], [4, 2]]) == marge_ev
    @test forget_edges(marge_ev, one_line, [2], [2], [[(2, 2), (1, 2)]]) == marge_eve

    @test forget_edges(fine, two_lines_three_points, [1, 3], [1, 3, 4, 6], [[(1, 1), (3, 1), (3, 3), (1, 3)], [(2, 1), (4, 1), (4, 3), (2, 3)]]) == fine_e
end

@testset "nonline" begin
    nonlinecc = CellComplex([1, 2], 2, [(1, 2, 1), (2, 1, 1), (2, 1, 2), (2, 1, 2)], [[(0, 2), (1, 3), (0, 2), (1, 4)], [(0, 1), (0, 3)], [(0, 1), (0, 4)]])

    @test iscellcomplex(nonlinecc)

    nonlinenwt = NWT(nonlinecc, 0, [3, 2, 3, 3], [[1, 1, 0, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]], [])

    @test isNWT(nonlinenwt)

    nonlinenwt2 = NWT(nonlinecc, 0, [0, 2, 0, 0], [[1, 0, 1, 0], [], []], [])

    @test isNWT(nonlinenwt2)
end