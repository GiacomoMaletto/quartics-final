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

function twosides(cc::CellComplex)
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
function groundmatrix(ndf)
    cc = ndf.cc
    all(==(2), length.(twosides(cc))) || error("Not two-sided!")
    fgs = f_ground_sequence(ndf)
    ng = length(unique(vcat(fgs...)))
    gmat = zeros(Int, ng, ng)
    nt = length(cc.T) + 1
    for seq in fgs
        for (g1, g2) in zip(seq[1:end-1], seq[2:end])
            setindexsymmetric!(gmat, g1, g2, nt)
        end
    end

    fes = fe_ground_sequence(ndf)

    ef = [Tuple{Int,Bool}[] for _ in eachindex(cc.E)]
    for (f, F) in enumerate(cc.F)
        for (e, s) in F
            push!(ef[e], (f, s))
        end
    end
    @assert(all(length.(ef) .== 2))

    for e in eachindex(ndf.n)
        ((f1, s1), (f2, s2)) = ef[e]
        s1 = reversebool(fes[(f1, e)], s1)
        s2 = reversebool(fes[(f2, e)], s2)
        for (g1, g2) in zip(s1, s2)
            setindexsymmetric!(gmat, g1, g2, cc.E[e].t)
        end
    end

    return gmat
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
    all(==(2), length.twosides(cc)) || error("Not two-sided!")
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

# ndf = marge
# fmap = [[(1, 1, 0), (3, 4, 1), (3, 6, 1), (1, 3, 0)], [(2, 1, 0), (4, 4, 1), (4, 3, 1), (2, 6, 0)]]
# emap = ([1, 3], [[1, 3, 3, 1], [2, 4, 4, 2]])

# nT = [1, 3]
# nE = [1, 3, 4, 6]
# nB = [[1, 3, 3, 1], [2, 4, 4, 2]]
function remove_edges(ndf::NDF, ncc::CellComplex, nT::Vector{Int}, nE::Vector{Int}, nB::Vector{Vector{Int}})::NDF
    fs = twosides(ndf.cc)
    all(==(2), length.(fs)) || error("Not two-sided!")

    fem = fe_mabel(ndf)
    fegs = fe_ground_sequence(ndf)

    unified_mabel = Dict([(l, l) for l in vcat(f_mabel(ndf)...)])
    for e in findall(E -> !(E.t in nT), ndf.cc.E)
        for (l1, l2) in zip(fem[(fs[e][1], e)], fem[(fs[e][2], e)])
            unified_mabel[l2] = l1
        end
    end

    a = []
    b = []
    for nf in eachindex(ncc.F)
        push!(a, [])
        push!(b, [])
        for (f, (ne, s)) in zip(nB[nf], ncc.F[nf])
            e = nE[ne]
            append!(a[end], [unified_mabel[l] for l in reversebool(fem[(f, e)], s)])
            append!(b[end], reversebool(fegs[(f, e)], s))
        end
    end
    nD = mabel_dyck.(remabel.(a))

    println(b)
    println(f_ground_sequence(nD))

    # used_mabel = unique(sort(vcat(a...)))
    # unused_mabel = setdiff(values(unified_mabel), used_mabel)
    # used_ground = used_mabel .+ 1
    # unused_ground = unused_mabel .+ 1
    unused_old_mabel = setdiff(values(unified_mabel), vcat(a...))
    unused_old_ground = unused_old_mabel .+ 1

    # println(a)
    # println(f_mabel(nD))

    nt = length(ndf.cc.T) + 1
    gmat = groundmatrix(ndf)
    flforest = floatingforest(ndf)
    for ug in unused_ground
        for un in findall(==(nt), gmat[ug, :])
            add_edge!(flforest, ug, un)
        end
    end
    nF = []
    for g in used_ground
        D = rooted_tree_to_dyck(flforest, g)
        if !isempty(D)
            push!(nF, (g=g, F=D))
        end
    end


    return NDF(ncc, ndf.d, ndf.n[nE], nD, nF)
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

three_two_lines_emap = ([1, 3], [[1, 3, 3, 1], [2, 4, 4, 2]])

two_lines = CellComplex(
    [1, 1],
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(1, 0), (2, 0)], [(1, 0), (2, 1)]]
)

two_lines_vmap = [2]

one_line = CellComplex(
    [1],
    1,
    [(1, 1, 1)],
    [[(1, 0)]]
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
    [(1, [1, 0]), (5, [1, 0]), (17, [1, 0])]
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