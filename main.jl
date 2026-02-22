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
    k::Int
    nV::Int
    E::Vector{@NamedTuple{s::Int, d::Int, i::Int}}
    F::Vector{Vector{@NamedTuple{s::Bool, e::Int}}}
end

Base.:(==)(cc1::CellComplex, cc2::CellComplex) = cc1.k == cc2.k && cc1.nV == cc2.nV && cc1.E == cc2.E && cc1.F == cc2.F
Base.hash(cc::CellComplex) = hash((cc.k, cc.nV, cc.E, cc.F))

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
        if !(i in 1:cc.k)
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

# a projective cell complex is a complex whose realization is the real projective plane, and whose first curve is a line
function isprojective(cc::CellComplex)::Bool
    return true
end

# Combinatorial arrangement
# (n, W, T) data of a curve D of degree d intersecting transversely the cellular complex cc
# n[e]: the number of intersections of C with the edge e
# W[f]: the Dyck word of the face f
# T: vector of tuples (g, T) where T is a rooted tree through a Dyck word, with ground g
struct NWT
    cc::CellComplex
    n::Vector{Int}
    W::Vector{Vector{Bool}}
    T::Vector{@NamedTuple{g::Int, T::Vector{Bool}}}
end

Base.:(==)(nwt1::NWT, nwt2::NWT) = nwt1.cc == nwt2.cc && nwt1.n == nwt2.n && nwt1.W == nwt2.W && nwt1.T == nwt2.T
Base.hash(nwt::NWT) = hash((nwt.cc, nwt.n, nwt.W, nwt.T))

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
    gmat = [zeros(Bool, ng, ng) for _ in 1:(cc.k+1)]

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

function curvematrix(nwt)
    grmatrix = groundmatrix(nwt)
    ng = size(grmatrix[1], 1)
    flforest = floatingforest(nwt)
    cmatrix = adjacency_matrix(flforest)
    cmatrix[1:ng, 1:ng] += grmatrix[end]
    return cmatrix
end

function isNWT(nwt::NWT)::Bool
    cc = nwt.cc
    if !iscellcomplex(cc)
        return false
    end

    if length(nwt.n) != length(cc.E)
        return false
    end

    # for (i, d) in enumerate(cc.d)
    #     if d != 0 && nwt.d != 0
    #         es = findall(E -> E.i == i, cc.E)
    #         n_t = sum(nwt.n[es])
    #         if n_t > d * nwt.d || (d * nwt.d - n_t) % 2 != 0
    #             return false
    #         end
    #     end
    # end

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

    removed_i = setdiff(1:in_cc.k, Imap)
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

    return NWT(out_cc, in_nwt.n[Emap], out_D, out_T)
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

    (in_cc.k == out_cc.k) || error("Not the same type.")

    out_E = [sum(in_nwt.n[Emap[e]]) for e in eachindex(out_cc.E)]

    return NWT(out_cc, out_E, in_nwt.W, in_nwt.T)
end

# the ground faces are indexed as in the indexing used in nwt.T
function groundcc(in_nwt::NWT)::CellComplex
    in_cc = in_nwt.cc

    outk = in_cc.k + 1

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
            push!(outE, (outE[B[i][end].e].d, outE[B[j][end].e].d, outk))
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

    return CellComplex(outk, outnV, outE, outF)
end

# given a nwt on a projective CellComplex, return the topological type of the curve (b, T)
# where b==true iff there is a pseudoline, and T is the rooted tree representing the regions of the curve,
# where the root is the unique region whose interior of the closure is nonorientable
function curve_type(nwt::NWT)::Tuple{Bool,Vector{Bool}}
    gm = groundmatrix(nwt)
    n = size(gm[1], 1)

    cm = curvematrix(nwt)
    N = size(cm, 1)

    lm = falses(N, N)
    lm[1:n, 1:n] += gm[1]

    vss = connected_components(SimpleGraph(max.(gm[2:end-1]...)))
    lc_g, lc_vmap = my_merge_vertices_list(SimpleGraph(cm), vss)
    l_g, l_vmap = my_merge_vertices_list(SimpleGraph(lm), vss)
    @assert(lc_vmap == l_vmap)

    c_g, c_vmap = my_merge_vertices_list(lc_g, connected_components(l_g))

    l_ext = [v for v in 1:nv(l_g) if has_edge(l_g, v, v)]
    c_ext = [v for v in 1:nv(c_g) if has_edge(c_g, v, v)]
    if isempty(l_ext)
        (length(c_ext) == 1) || error("There should be a unique self-loop in the odd case")
        rem_edge!(c_g, c_ext[1], c_ext[1])
        return (true, rooted_tree_to_dyck(c_g, c_ext[1]))
    else
        (length(c_ext) == 0) || error("There should be no self-loops in the even case")
        return (false, rooted_tree_to_dyck(c_g, only(unique(c_vmap[l_ext]))))
    end
end

# # Generate all NWT with specified cc, n and no trees
# function allW(cc::CellComplex, n::Vector{Int})::Vector{NWT}
#     result = NWT[]
#     fn = [sum(n[getproperty.(F, :e)]) for F in cc.F]
#     (all(fn .% 2 .== 0)) || error("Does not add up to an even number.")

#     for W in product1d((all_dyck_words.(fn))...)
#         push!(result, NWT(cc, n, [W...], []))
#     end
#     return result
# end

function add_line(incc::CellComplex, fp::Vector{Int})::Vector{NWT}
    result = []
    inFe = [getproperty.(F, :e) for F in incc.F]

    for edgepath in product1d([intersect(inFe[f1], inFe[f2]) for (f1, f2) in zip(fp[1:end-1], fp[2:end])]...)
        edgepath2 = [edgepath[end], edgepath...]

        n = [0 for _ in incc.E]
        for e in edgepath
            n[e] += 1
        end

        fW = [[] for _ in eachindex(incc.F)]
        for (f, F) in enumerate(incc.F)
            p = sort([sort((edgepath2[i], edgepath2[i+1])) for i in findall(==(f), fp[1:end-1])])
            es = vcat([fill(e, n[e]) for (s, e) in F]...)

            for W in all_dyck_words(length(p))
                if p == sort([sort((es[i], es[j])) for (i, j) in dyck_pairs(W)])
                    push!(fW[f], W)
                end
            end
        end

        for Ws in product1d(fW...)
            nwt = NWT(incc, n, [Ws...], [])
            if curve_type(nwt) == (true, [])
                push!(result, nwt)
            end
        end
    end

    return result
end

### TESTS ###

three_lines = CellComplex(
    3,
    3,
    [(2, 3, 1), (3, 1, 2), (1, 2, 3), (2, 3, 1), (3, 1, 2), (1, 2, 3)],
    [[(0, 1), (0, 2), (0, 3)], [(0, 1), (0, 5), (0, 6)], [(0, 4), (0, 2), (0, 6)], [(0, 4), (0, 5), (0, 3)]])

two_lines_three_points = CellComplex(
    2,
    3,
    [(2, 3, 1), (1, 2, 2), (2, 3, 1), (1, 2, 2)],
    [[(0, 1), (1, 3), (1, 4), (0, 2)], [(0, 1), (1, 3), (1, 2), (0, 4)]])

two_lines = CellComplex(
    2,
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(0, 1), (0, 2)], [(0, 1), (1, 2)]]
)

one_line = CellComplex(
    1,
    1,
    [(1, 1, 1)],
    [[(0, 1), (0, 1)]]
)

zero_lines = CellComplex(
    0,
    0,
    [],
    [[]]
)

fine = NWT(
    three_lines,
    [4, 1, 3, 0, 1, 1],
    [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 0], [1, 0], [1, 1, 0, 0]],
    [(1, [1, 0, 1, 0]), (10, [1, 0])])

nonviro = NWT(
    three_lines,
    [3, 4, 1, 1, 0, 3],
    [[1, 0, 1, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 1, 0, 0], [1, 0]],
    [])

marge = NWT(
    three_lines,
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
        [4, 3, 2, 3],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_ev = NWT(
        two_lines,
        [6, 6],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_eve = NWT(
        one_line,
        [6],
        [[1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    fine_e = NWT(
        two_lines_three_points,
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
    nonlinecc = CellComplex(2, 2, [(1, 2, 1), (2, 1, 1), (2, 1, 2), (2, 1, 2)], [[(0, 2), (1, 3), (0, 2), (1, 4)], [(0, 1), (0, 3)], [(0, 1), (0, 4)]])

    @test iscellcomplex(nonlinecc)

    nonlinenwt = NWT(nonlinecc, [3, 2, 3, 3], [[1, 1, 0, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]], [])

    @test isNWT(nonlinenwt)

    nonlinenwt2 = NWT(nonlinecc, [0, 2, 0, 0], [[1, 0, 1, 0], [], []], [])

    @test isNWT(nonlinenwt2)
end