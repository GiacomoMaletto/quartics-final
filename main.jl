using Test
using Graphs

@inline function setindexsymmetric!(S, i, j, v)
    S[i, j] = v
    S[j, i] = v
end

# Combinatorial cellular complex
# T: vector indicating the degree of each curve (0 if it's just a curve)
# nV: number of vertices
# E: vector of tuples (s, d, t) where s is the source, d the destination, t the type of each edge
# F: vector of vectors of tuples (e, s) where e is the edge and s the sign (1==true: opposite orientation)
struct CellComplex
    T::Vector{Int}
    nV::Int
    E::Vector{@NamedTuple{s::Int, d::Int, t::Int}}
    F::Vector{Vector{@NamedTuple{e::Int, s::Bool}}}
end

Base.:(==)(cc1::CellComplex, cc2::CellComplex) = cc1.T == cc2.T && cc1.nV == cc2.nV && cc1.E == cc2.E && cc1.F == cc2.F
Base.hash(cc::CellComplex) = hash([cc.T, cc.nV, cc.E, cc.F])

# (n, D, F) data of a curve C of degree d intersecting the cellular complex cc
# n[e]: the number of intersections of C with the edge e
# D[f]: the Dyck word of the face f
# F: vector of tuples (g, F) where F is a forest expressed through a Dyck word, with ground g
struct NDF
    cc::CellComplex
    d::Int
    n::Vector{Int}
    D::Vector{Vector{Bool}}
    F::Vector{@NamedTuple{g::Int, F::Vector{Bool}}}
end

Base.:(==)(ndf1::NDF, ndf2::NDF) = ndf1.cc == ndf2.cc && ndf1.d == ndf2.d && ndf1.n == ndf2.n && ndf1.D == ndf2.D && ndf1.F == ndf2.F
Base.hash(ndf::NDF) = hash([ndf.cc, ndf.d, ndf.n, ndf.D, ndf.F])

function source(e)
    return e.s
end

function source(e, s)
    return s ? e.d : e.s
end

function destination(e)
    return e.d
end

function destination(e, s)
    return s ? e.s : e.d
end

function iscellcomplex(cc::CellComplex)::Bool
    for (e, (_, _, t)) in enumerate(cc.E)
        if !(t in eachindex(cc.T) || t == 0)
            return false
        end
        if !(count(f -> e in getproperty.(f, :e), cc.F) in [1, 2])
            return false
        end
    end

    for f in cc.F
        for ((e1, s1), (e2, s2)) in zip(f, circshift(f, -1))
            if destination(cc.E[e1], s1) != source(cc.E[e2], s2)
                return false
            end
        end
    end

    return true
end

# function istwosided(cc::CellComplex)::Bool
#     if !iscellcomplex(cc)
#         return false
#     end

#     for e in eachindex(cc.E)
#         if !(count(f -> e in getproperty.(f, :e), cc.F) == 2)
#             return false
#         end
#     end

#     return true
# end

function edge_faces(cc::CellComplex)
    return [findall(f -> e in getproperty.(f, :e), cc.F) for e in eachindex(cc.E)]
end

# ChatGPT-generated
function dyck_sequence(D::Vector{Bool}, start_label)
    out = Int[]
    push!(out, start_label)            # start at root labelled 1
    stack = [start_label]              # current path (top = current node)
    next_label = start_label + 1

    for b in D
        if b == 1
            push!(stack, next_label)   # create and move to new child
            push!(out, next_label)
            next_label += 1
        elseif b == 0
            pop!(stack)                # move up (must be valid Dyck word)
            push!(out, stack[end])     # record the parent's label
        else
            throw(ArgumentError("D must contain only 0 or 1"))
        end
    end
    return out
end

function sequence_dyck(seq::AbstractVector{<:Integer})
    m = length(seq)
    m >= 1 || throw(ArgumentError("borseq must be nonempty"))
    # borseq[1] == 1 || throw(ArgumentError("first label must be 1"))

    # quick type-correct output
    # dyck = Vector{Bool}(undef, max(0, m - 1))
    dyck = Bool[]

    # Step 1: infer bits by comparing consecutive borseq
    for i in 1:m-1
        a, b = seq[i], seq[i+1]
        if b == a
            throw(ArgumentError("invalid transition: two consecutive equal regions")) # continue
        end
        push!(dyck, b > a)
    end

    # @assert(dyck_sequence(dyck, seq[1]) == seq)

    return dyck
end

function f_ground_sequence(Ds::Vector{Vector{Bool}})::Vector{Vector{Int}}
    i = 1
    fgs = []
    for D in Ds
        push!(fgs, dyck_sequence(D, i))
        i += div(length(fgs[end]), 2) + 1
    end
    return fgs
end

function f_ground_sequence(ndf::NDF)::Vector{Vector{Int}}
    return f_ground_sequence(ndf.D)
end

function reversebool(A, b::Bool)
    return b ? reverse(A) : A
end

function fe_ground_sequence(ndf::NDF)::Dict{Tuple{Int,Int},Vector{Int}}
    fgs = f_ground_sequence(ndf)

    fegs = Dict()
    for (f, F) in enumerate(ndf.cc.F)
        i = 1
        for (e, s) in F
            fegs[(f, e)] = reversebool(fgs[f][i:(i+ndf.n[e])], s)
            i += ndf.n[e]
        end
    end

    return fegs
end

function fi_ground_sequence(ndf::NDF)::Dict{Tuple{Int,Int},Vector{Int}}
    fgs = f_ground_sequence(ndf)

    figs = Dict()
    for (f, F) in enumerate(ndf.cc.F)
        i = 1
        for (j, (e, s)) in enumerate(F)
            figs[(f, j)] = reversebool(fgs[f][i:(i+ndf.n[e])], s)
            i += ndf.n[e]
        end
    end

    return figs
end

# ChatGPT-generated
function dyck_mabel(word, start_label)
    labels = Int[]
    stack = Int[]
    next_label = start_label

    for x in word
        if x == 1 || x === true
            push!(stack, next_label)
            push!(labels, next_label)
            next_label += 1
        else
            isempty(stack) && error("Invalid Dyck word: too many closes")
            id = pop!(stack)
            push!(labels, id)
        end
    end

    isempty(stack) || error("Invalid Dyck word: unclosed opens")
    return labels
end

# ChatGPT-generated
function mabel_dyck(labels)::Vector{Bool}
    seen = Set{Int}()
    dyck = Vector{Int}(undef, length(labels))

    for (i, x) in pairs(labels)
        if x ∈ seen
            dyck[i] = 0   # close
        else
            push!(seen, x)
            dyck[i] = 1   # open
        end
    end

    return dyck
end

# ChatGPT-generated
function remabel(v)
    mapping = Dict{Int,Int}()
    next_id = 0
    out = similar(v)

    for (i, x) in pairs(v)
        if !haskey(mapping, x)
            next_id += 1
            mapping[x] = next_id
        end
        out[i] = mapping[x]
    end

    return out
end

function f_mabel(Ds::Vector{Vector{Bool}})::Vector{Vector{Int}}
    i = 1
    fgs = []
    for D in Ds
        push!(fgs, dyck_mabel(D, i))
        i += div(length(fgs[end]), 2) + 1 # i += div(length(fgs[end]), 2)
    end
    return fgs
end

function f_mabel(ndf::NDF)::Vector{Vector{Int}}
    return f_mabel(ndf.D)
end

# function f_mabel(ndf::NDF)
#     i = 1
#     fgs = []
#     for D in ndf.D
#         push!(fgs, dyck_mabel(D, i))
#         i += div(length(fgs[end]), 2) + 1 # i += div(length(fgs[end]), 2)
#     end
#     return fgs
# end

function fe_mabel(ndf::NDF)::Dict{Tuple{Int,Int},Vector{Int}}
    fm = f_mabel(ndf)

    fem = Dict()
    for (f, F) in enumerate(ndf.cc.F)
        i = 1
        for (e, s) in F
            fem[(f, e)] = reversebool(fm[f][i:(i+ndf.n[e]-1)], s)
            i += ndf.n[e]
        end
    end

    return fem
end

# the root is always the first vertex
# ChatGPT-generated, some alterations by me
function dyck_to_rooted_tree(D::Vector{Bool})
    g = SimpleGraph(1)
    stack = Int[1]

    for b in D
        if b
            v = nv(g) + 1
            add_vertex!(g)
            if isempty(stack)
                error("Invalid Dyck word: too many closes")
            else
                add_edge!(g, stack[end], v)
            end
            push!(stack, v)
        else
            pop!(stack)
        end
    end

    (stack == [1]) || error("Invalid Dyck word: unclosed opens")
    return g, 1
end

# ChatGPT-generated, some alterations by me
function rooted_tree_to_dyck(g::SimpleGraph, root::Int)
    D = Bool[]

    for w in neighbors(g, root)
        _dfs_dyck!(D, g, w, root)
    end

    return forest_to_dyck(dyck_to_forest(D))
end
function _dfs_dyck!(D::Vector{Bool}, g::SimpleGraph, v::Int, parent::Int)
    push!(D, true)   # open parenthesis

    for w in neighbors(g, v)
        if w != parent
            _dfs_dyck!(D, g, w, v)
        end
    end

    push!(D, false)  # close parenthesis
end

function union2(g::SimpleGraph, h::SimpleGraph)
    return SimpleGraph(blockdiag(adjacency_matrix(g), adjacency_matrix(h)))
end

# input: a graph g and a list of vertices vs
# output: the graph gprime obtained by merging every vertex in vs to minimum(vs) and otherwise keeping the same order,
# and a list vmap corresponding to the morphism from g to gprime
function merge_vertices1(g, vs)
    vs = unique(vs)
    if length(vs) <= 1
        return g, collect(vertices(g))
    end
    gprime = merge_vertices(g, vs)
    vmap = Int[]
    vs = sort(vs)
    i = -1
    for v in 1:nv(g)
        if v in vs
            push!(vmap, vs[1])
            i += 1
        else
            push!(vmap, v - max(0, i))
        end
    end
    return gprime, vmap
end

# input: a graph g and a list of lists of vertices vs
# output: the graph obtained by successively merging the vertices in vs[1], vs[2], ...
# and a list vmap corresponding to the morphism from the original graph to the new one
function merge_vertices_list(g, vss)
    vmap_final = collect(1:nv(g))
    for i in eachindex(vss)
        g, vmap = merge_vertices1(g, vss[i])
        vss = map(vs -> unique([vmap[v] for v in vs]), vss)
        vmap_final = [vmap[x] for x in vmap_final]
    end

    return g, vmap_final
end

# compute the groundmatrix of a ndf over a two-sided cellular complex
function groundmatrix(ndf)::Vector{Matrix{Bool}}
    cc = ndf.cc

    fgs = f_ground_sequence(ndf)
    ng = length(unique(vcat(fgs...)))
    gmat = [zeros(Bool, ng, ng) for _ in 1:(length(cc.T)+1)]

    for seq in fgs
        for (g1, g2) in zip(seq[1:end-1], seq[2:end])
            setindexsymmetric!(gmat[end], g1, g2, 1)
        end
    end

    efis = [Tuple{Int,Int,Bool}[] for _ in eachindex(cc.E)]
    for (f, F) in enumerate(cc.F)
        for (i, (e, s)) in enumerate(F)
            push!(efis[e], (f, i, s))
        end
    end
    @assert(all(length.(efis) .== 2))

    figs = fi_ground_sequence(ndf)

    for e in eachindex(ndf.n)
        ((f1, i1, s1), (f2, i2, s2)) = efis[e]
        seq1 = reversebool(figs[(f1, i1)], s1)
        seq2 = reversebool(figs[(f2, i2)], s2)
        for (g1, g2) in zip(seq1, seq2)
            setindexsymmetric!(gmat[cc.E[e].t], g1, g2, 1)
        end
    end

    return gmat

    # all(==(2), length.(edge_faces(cc))) || error("Not two-sided!")
    # fgs = f_ground_sequence(ndf)
    # ng = length(unique(vcat(fgs...)))
    # gmat = zeros(Int, ng, ng)
    # nt = length(cc.T) + 1
    # for seq in fgs
    #     for (g1, g2) in zip(seq[1:end-1], seq[2:end])
    #         setindexsymmetric!(gmat, g1, g2, nt)
    #     end
    # end

    # fes = fe_ground_sequence(ndf)

    # ef = [Tuple{Int,Bool}[] for _ in eachindex(cc.E)]
    # for (f, F) in enumerate(cc.F)
    #     for (e, s) in F
    #         push!(ef[e], (f, s))
    #     end
    # end
    # @assert(all(length.(ef) .== 2))

    # for e in eachindex(ndf.n)
    #     ((f1, s1), (f2, s2)) = ef[e]
    #     s1 = reversebool(fes[(f1, e)], s1)
    #     s2 = reversebool(fes[(f2, e)], s2)
    #     for (g1, g2) in zip(s1, s2)
    #         setindexsymmetric!(gmat, g1, g2, cc.E[e].t)
    #     end
    # end

    # return gmat
end

function floatingforest(ndf)
    ng = length(unique(vcat(f_ground_sequence(ndf)...)))
    graph = SimpleGraph(ng)

    vss = []
    for (g, F) in ndf.F
        rt, root = dyck_to_rooted_tree(F)
        push!(vss, [g, nv(graph) + root])
        graph = union2(graph, rt)
    end

    graph, _ = merge_vertices_list(graph, vss)
    return graph
end

# compute the regionmatrix of a ndf over a two-sided cellular complex
function regionmatrix(ndf)
    all(==(2), length.edge_faces(cc)) || error("Not two-sided!")
    grmatrix = groundmatrix(ndf)
    ng = size(grmatrix, 1)
    flforest = floatingforest(ndf)
    nt = length(ndf.cc.T) + 1
    rmatrix = nt * adjacency_matrix(flforest)
    rmatrix[1:ng, 1:ng] += grmatrix
    return rmatrix
end

# ChatGPT-generated, some alterations by me
function dyck_to_forest(D::Vector{Bool})
    stack = Vector{Vector{Any}}()
    push!(stack, Any[])  # this will be the forest

    for b in D
        if b  # open parenthesis
            push!(stack, Any[])
        else  # close parenthesis
            node = pop!(stack)
            push!(stack[end], node)
            sort!(stack[end])
        end
    end

    return only(stack)
end

# ChatGPT-generated
function forest_to_dyck(forest::Vector{Any})
    D = Bool[]
    for node in forest
        append!(D, node_to_dyck(node))
    end
    return D
end
function node_to_dyck(node::Vector{Any})
    D = Bool[1]
    for child in node
        append!(D, node_to_dyck(child))
    end
    push!(D, 0)
    return D
end

function isNDF(ndf::NDF)::Bool
    cc = ndf.cc
    if !iscellcomplex(cc)
        return false
    end

    if length(ndf.n) != length(cc.E)
        return false
    end

    for (t, d) in enumerate(cc.T)
        if d != 0 && ndf.d != 0
            es = findall(E -> E.t == t, cc.E)
            n_t = sum(ndf.n[es])
            if n_t > d * ndf.d || (d * ndf.d - n_t) % 2 != 0
                return false
            end
        end
    end

    for (f, F) in enumerate(cc.F)
        es = getproperty.(F, :e)
        n_f = sum(ndf.n[es])
        if n_f % 2 != 0 || n_f != length(ndf.D[f])
            return false
        end
    end

    if !issorted(getproperty.(ndf.F, :g))
        return false
    end
    gs = unique(vcat(f_ground_sequence(ndf)...))
    for (g, F) in ndf.F
        if !(g in gs)
            return false
        end
        if F != forest_to_dyck(dyck_to_forest(F))
            return false
        end
    end

    return true
end

function remove_consecutive(v)
    return [v[1]; v[2:end][v[2:end].!=v[1:end-1]]]
end

# input:
# - in_ndf: the input ndf, where in_cc = in_ndf.cc must be two-sided
# - out_cc: cc of the output, and must have the same vertices, a subset of the edges, and the merged faces
# - Tmap: if an edge has type t in out_cc, it has type Tmap in in_cc
# - Emap: edge e in out_cc.E is edge Emap[e] in in_cc.E
# - Fmap: the i-th border of the face f in out_cc.F comes from the face Fmap[f][i] of in_cc.F
# output:
# - out_ndf: obtained by removing some edges from in_ndf
function remove_edges(in_ndf::NDF, out_cc::CellComplex, Tmap::Vector{Int}, Emap::Vector{Int}, Fmap::Vector{Vector{Int}})::NDF
    in_cc = in_ndf.cc
    fs = edge_faces(in_cc)
    all(==(2), length.(fs)) || error("Not two-sided.")

    (in_cc.nV == out_cc.nV) || error("Not the same number of vertices.")

    (length(out_cc.E) == length(Emap) == length(unique(Emap))) || error("Not a subset of the edges.")
    for (e, E) in enumerate(out_cc.E)
        in_E = in_cc.E[Emap[e]]
        (E.s == in_E.s && E.d == in_E.d && Tmap[E.t] == in_E.t) || error("Edges do not match.")
    end

    for (f, F) in enumerate(out_cc.F)
        (length(Fmap[f]) == length(F)) || error("Not enough borders.")
    end

    in_fegs = fe_ground_sequence(in_ndf)

    gmat = groundmatrix(in_ndf)

    removed_t = setdiff(eachindex(in_cc.T), Tmap)
    removed_matrix = reduce((m1, m2) -> m1 .|| m2, gmat[removed_t])
    removed_components = connected_components(SimpleGraph(removed_matrix))
    removed_map = Dict{Int,Int}()
    for (i, c) in enumerate(removed_components)
        for g in c
            removed_map[g] = i
        end
    end

    # unified_grounds = Dict([(g, g) for g in vcat(f_ground_sequence(in_ndf)...)])
    # for e in findall(E -> !(E.t in Tmap), in_cc.E)
    #     for (g1, g2) in zip(in_fegs[(fs[e][1], e)], in_fegs[(fs[e][2], e)])
    #         unified_grounds[g2] = g1
    #     end
    # end

    # println(f_ground_sequence(in_ndf))
    # println(unified_grounds)

    merged_fgs = Vector{Int}[]
    for (f, F) in enumerate(out_cc.F)
        push!(merged_fgs, [])
        for (in_f, (e, s)) in zip(Fmap[f], F)
            append!(merged_fgs[end], [removed_map[g] for g in reversebool(in_fegs[(in_f, Emap[e])], s)])
        end
    end
    merged_fgs = remove_consecutive.(merged_fgs)
    out_D = sequence_dyck.(merged_fgs)

    merged_gs = vcat(merged_fgs...)
    out_gs = vcat(f_ground_sequence(out_D)...)
    gmap = Dict(zip(out_gs, merged_gs))

    unused_ground = [g for g in 1:size(gmat[1], 1) if !(removed_map[g] in merged_gs)]

    # println(removed_map)
    # println(size(gmat[1], 1))
    # println(removed_components)
    # println(merged_gs)
    # println(unused_ground)

    # setdiff(values(unified_grounds), merged_gs)
    # gmat = groundmatrix(in_ndf)
    flforest = floatingforest(in_ndf)
    for ug in unused_ground
        for neigh in findall(gmat[end][ug, :])
            add_edge!(flforest, ug, neigh)
        end
    end
    flforest, vmap = merge_vertices_list(flforest, removed_components)

    out_F = []
    for (ng, mg) in gmap
        D = rooted_tree_to_dyck(flforest, vmap[removed_components[mg][1]])
        if !isempty(D)
            push!(out_F, (g=ng, F=D))
        end
    end

    return NDF(out_cc, in_ndf.d, in_ndf.n[Emap], out_D, out_F)
end

# removes vertices, provided they are used by either 0 or 2 edges
function remove_vertices(in_ndf::NDF, out_cc::CellComplex, Vmap::Vector{Int}, Emap::Vector{Vector{Int}})
    in_cc = in_ndf.cc

    (length(out_cc.nV) == length(Vmap) == length(unique(Vmap))) || error("Not a subset of the vertices.")
    for v in setdiff(1:in_cc.nV, Vmap)
        n = count(==(v), getproperty.(in_cc.E, :s)) + count(==(v), getproperty.(in_cc.E, :d))
        (n in [0, 2]) || error("Vertex not occurring 0 or 2 times.")
    end

    (in_cc.T == out_cc.T) || error("Not the same number of types.")

    # (length(Emap) == length(out_cc.E)) || error("Not enough elements in Emap.")
    # for in_es in Emap
    #     for ((e1, s1), (e2, s2)) in zip(in_es[1:end-1], in_es[2:end])
    #         if destination(in_cc.E[e1], s1) != source(in_cc.E[e2], s2) || in_cc.E[e1].t != in_cc.E[e2].t
    #             error("Edges in Emap not agreeing.")
    #         end
    #         if destination(in_cc.E[e1], s1) in Vmap
    #             error("Merging two edges whose common vertex hasn't been removed.")
    #         end
    #     end
    # end

    out_E = [sum(in_ndf.n[Emap[e]]) for e in eachindex(out_cc.E)]

    return NDF(out_cc, in_ndf.d, out_E, in_ndf.D, in_ndf.F)
end

three_lines = CellComplex(
    [1, 1, 1],
    3,
    [(2, 3, 1), (3, 1, 2), (1, 2, 3), (2, 3, 1), (3, 1, 2), (1, 2, 3)],
    [[(1, 0), (2, 0), (3, 0)], [(1, 0), (5, 0), (6, 0)], [(4, 0), (2, 0), (6, 0)], [(4, 0), (5, 0), (3, 0)]])

two_lines_three_points = CellComplex(
    [1, 1],
    3,
    [(2, 3, 1), (1, 2, 2), (2, 3, 1), (1, 2, 2)],
    [[(1, 0), (3, 1), (4, 1), (2, 0)], [(1, 0), (3, 1), (2, 1), (4, 0)]])

two_lines = CellComplex(
    [1, 1],
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(1, 0), (2, 0)], [(1, 0), (2, 1)]]
)

one_line = CellComplex(
    [1],
    1,
    [(1, 1, 1)],
    [[(1, 0), (1, 0)]]
)

zero_lines = CellComplex(
    [],
    0,
    [],
    [[]]
)

fine = NDF(
    three_lines,
    4,
    [4, 1, 3, 0, 1, 1],
    [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 0], [1, 0], [1, 1, 0, 0]],
    [(1, [1, 0, 1, 0]), (10, [1, 0])])

nonviro = NDF(
    three_lines,
    4,
    [3, 4, 1, 1, 0, 3],
    [[1, 0, 1, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 1, 0, 0], [1, 0]],
    [])

marge = NDF(
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

    # @test forget_curve(three_lines, 2) == CellComplex([1, 1], )
end

@testset "isNDF" begin
    @test isNDF(fine)
    @test isNDF(nonviro)
    @test isNDF(marge)
end

@testset "forest" begin
    @test dyck_to_forest(Bool[1, 0, 1, 1, 0, 0]) == [[], [[]]]
    @test dyck_to_forest(Bool[1, 1, 0, 0, 1, 0]) == [[], [[]]]
    @test forest_to_dyck(Any[Any[], Any[Any[]]]) == Bool[1, 0, 1, 1, 0, 0]
    @test forest_to_dyck([]) == Bool[]
end

@testset "merge" begin
    g = SimpleGraph(Edge.([
        (1, 3), (1, 10),
        (2, 3),
        (3, 5), (3, 6), (3, 9),
        (4, 8), (4, 10),
        (5, 10),
        (6, 10),
        (7, 10),
        (8, 9)
    ]))

    gprime, vmap = merge_vertices_list(g, [[9, 4], [5, 3, 7]])

    @test gprime == SimpleGraph(Edge.([
        (1, 3), (1, 7),
        (2, 3),
        (3, 4), (3, 5), (3, 7),
        (4, 6), (4, 7),
        (5, 7)
    ]))

    @test vmap == [1, 2, 3, 4, 3, 5, 3, 6, 4, 7]
end

@testset "remove_edges" begin
    marge_e = NDF(
        two_lines_three_points,
        0,
        [4, 3, 2, 3],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_ev = NDF(
        two_lines,
        0,
        [6, 6],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_eve = NDF(
        one_line,
        0,
        [6],
        [[1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    fine_e = NDF(
        two_lines_three_points,
        4,
        [4, 3, 0, 1],
        [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 1, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 0])]
    )

    @test remove_edges(marge, two_lines_three_points, [1, 3], [1, 3, 4, 6], [[1, 3, 3, 1], [2, 4, 4, 2]]) == marge_e
    @test remove_vertices(marge_e, two_lines, [2], [[1, 3], [4, 2]]) == marge_ev
    @test remove_edges(marge_ev, one_line, [2], [2], [[2, 1]]) == marge_eve

    @test remove_edges(fine, two_lines_three_points, [1, 3], [1, 3, 4, 6], [[1, 3, 3, 1], [2, 4, 4, 2]]) == fine_e
end