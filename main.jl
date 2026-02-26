using Test
using Graphs
using LinearAlgebra: I
using JSON3, CodecZlib
using Combinatorics: multiset_combinations

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

@inline function distinct_pairs(xs)
    return [(xs[i], xs[j]) for i in eachindex(xs) for j in (i+1):length(xs)]
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

function write_cc_to_jsonl!(io, cc)
    rec = Dict(
        "k" => cc.k,
        "nV" => cc.nV,
        "E" => cc.E,
        "F" => cc.F)
    JSON3.write(io, rec)
    write(io, '\n')
end

function read_cc_from_jsonl_line(line)
    rec = JSON3.read(line, CellComplex)
    return CellComplex(rec.k, rec.nV, rec.E, rec.F)
end

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

Base.show(io::IO, nwt::NWT) = print(io, (nwt.n, nwt.W, nwt.T))

struct NWT_from_jsonl
    n::Vector{UInt8}
    W::Vector{Vector{Bool}}
    T::Vector{@NamedTuple{g::UInt8, T::Vector{Bool}}}
    itriang::Int
    sign::Vector{Bool}
end

function write_nwt_to_jsonl!(io, nwt, itriang, sign)
    rec = Dict(
        "n" => UInt8.(nwt.n),
        "W" => nwt.W,
        "T" => @NamedTuple{g::UInt8, T::Vector{Bool}}.(nwt.T),
        "itriang" => itriang,
        "sign" => sign)
    JSON3.write(io, rec)
    write(io, '\n')
end

function read_nwt_from_jsonl_line(line, cc)
    rec = JSON3.read(line, NWT_from_jsonl)
    return NWT(cc, rec.n, rec.W, rec.T), rec.itriang, rec.sign
end

function emptyNWT(cc::CellComplex)::NWT
    return NWT(cc, [0 for _ in eachindex(cc.E)], [[] for _ in eachindex(cc.F)], [])
end

function f_ground_sequence(Ws::Vector{Vector{Bool}})::Vector{Vector{Int}}
    i = 0
    fgs = []
    for W in Ws
        push!(fgs, dyck_to_sequence(W) .+ i)
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
function groundmatrix(nwt)::Vector{BitMatrix}
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

    @assert(all(length.(efi) .== 1 .|| length.(efi) .== 2))

    figs = fi_ground_sequence(nwt)

    for e in eachindex(nwt.n)
        if length(efi[e]) == 2
            ((f1, i1), (f2, i2)) = efi[e]
            for (g1, g2) in zip(figs[(f1, i1)], figs[(f2, i2)])
                setindexsymmetric!(gmat[cc.E[e].i], g1, g2, 1)
            end
        end
    end

    return gmat
end

function floatingforest(nwt)::SimpleGraph
    ng = length(unique(vcat(f_ground_sequence(nwt)...)))
    graph = SimpleGraph(ng)

    vss = Vector{Int}[]
    for (g, T) in nwt.T
        rt, root = dyck_to_rtgraph(T)
        push!(vss, [g, nv(graph) + root])
        graph = my_union(graph, rt)
    end

    graph, _ = my_merge_vertices_list(graph, vss)
    return graph
end

function totalmatrix(nwt::NWT)::Vector{BitMatrix}
    gm = groundmatrix(nwt)
    n = size(gm[1], 1)

    fm = adjacency_matrix(floatingforest(nwt))
    N = size(fm, 1)

    tm = [falses(N, N) for _ in eachindex(gm)]
    for i in eachindex(gm)
        tm[i][1:n, 1:n] .= gm[i]
    end
    tm[end] += fm

    return tm
end

function isNWT(nwt::NWT)::Bool
    cc = nwt.cc
    if !iscellcomplex(cc)
        return false
    end

    if length(nwt.n) != length(cc.E)
        return false
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
# - Vmap: the i-th vertex of out_cc is the Vmap[i]-th vertex of in_cc
# - Imap: if an edge has type i in out_cc, it has type Imap[i] in in_cc
# - Emap: edge e in out_cc.E is edge Emap[e] in in_cc.E
# - Fmap: the i-th border of the face f in out_cc.F is the in_i-th border of the face in_f in in_cc.F, where Fmap[f][i]=(in_f,in_i)
# output:
# - out_nwt: obtained by removing some edges from in_nwt
function forget_edges(in_nwt::NWT, out_cc::CellComplex, Emap::Vector{Int}, Fmap::Vector{Vector{Tuple{Bool,Int,Int}}})::NWT
    in_cc = in_nwt.cc

    # (in_cc.nV == out_cc.nV == length(Vmap)) || error("Not the same number of vertices.")

    # (length(out_cc.E) == length(Emap) == length(unique(Emap))) || error("Not a subset of the edges.")
    # for (e, E) in enumerate(out_cc.E)
    #     in_E = in_cc.E[Emap[e]]
    #     (Vmap[E.s] == in_E.s && Vmap[E.d] == in_E.d && Imap[E.i] == in_E.i) || error("Edges do not match.")
    # end

    for (f, F) in enumerate(out_cc.F)
        (length(Fmap[f]) == length(F)) || error("Not enough borders.")
    end

    in_figs = fi_ground_sequence(in_nwt)

    gmat = groundmatrix(in_nwt)

    removed_i = setdiff(1:in_cc.k, getproperty.(in_cc.E[Emap], :i))
    removed_matrix = reduce((m1, m2) -> m1 .|| m2, gmat[removed_i]; init=falses(size(gmat[1])))
    removed_components = connected_components(SimpleGraph(removed_matrix))
    removed_map = Dict{Int,Int}()
    for (i, c) in enumerate(removed_components)
        for g in c
            removed_map[g] = i
        end
    end

    merged_fgs = Vector{Int}[]
    for f in eachindex(out_cc.F)
        # println("f: ", f)
        push!(merged_fgs, [])
        for ((in_s, in_f, in_i), (out_s, out_e)) in zip(Fmap[f], out_cc.F[f])
            # println("in_f: ", in_f)
            # println("in_i: ", in_i)
            # println("out_s: ", out_s)
            # println("out_e: ", out_e)
            (in_cc.F[in_f][in_i].e == Emap[out_e]) || error("Fmap format not correct.")
            # println("removed_map: ", removed_map)
            # println("in_figs: ", in_figs)
            append!(merged_fgs[end], [removed_map[g] for g in reversebool(xor(in_s, out_s), in_figs[(in_f, in_i)])])
            # println("merged_fgs: ", merged_fgs)
        end
    end
    merged_fgs = remove_consecutive.(merged_fgs)
    # println("merged_fgs: ", merged_fgs)
    out_D = sequence_to_dyck.(merged_fgs)

    merged_gs = vcat(merged_fgs...)
    # println("merged_gs: ", merged_gs)
    # println("out_D: ", out_D)
    out_gs = vcat(f_ground_sequence(out_D)...)
    # println("out_gs: ", out_gs)
    gmap = Dict(zip(out_gs, merged_gs))
    # println("gmap: ", gmap)

    unused_ground = [g for g in 1:size(gmat[1], 1) if !(removed_map[g] in merged_gs)]

    flforest = floatingforest(in_nwt)
    for ug in unused_ground
        for neigh in findall(gmat[end][ug, :])
            add_edge!(flforest, ug, neigh)
        end
    end
    flforest, vmap = my_merge_vertices_list(flforest, removed_components)
    # println("removed_components: ", removed_components)
    # println("vmap", vmap)
    out_T = []
    for (ng, mg) in gmap
        D = rtgraph_to_dyck((flforest, vmap[removed_components[mg][1]]))
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
        W = in_nwt.W[inf]
        B = Vector{@NamedTuple{s::Bool, e::Int}}[[]]
        for (ins, ine) in inF
            outes = reversebool(ins, ine_oute[ine])
            push!(B[end], (ins, outes[1]))
            for oute in outes[2:end]
                push!(B, [(ins, oute)])
            end
        end
        if isempty(W)
            push!(outF, only(B))
            continue
        end
        B[1] = [pop!(B); B[1]]
        @assert(length(B) == length(W))

        mabel = dyck_to_mabel(W) .+ length(outE)
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
function curve_type(nwt::NWT)::Tuple{Tuple{Bool,Vector{Bool}},Vector{Int}}
    tm = totalmatrix(nwt)

    vss = length(tm) > 2 ? connected_components(SimpleGraph(max.(tm[2:end-1]...))) : Vector{Int}[]
    lc_g, lc_vmap = my_merge_vertices_list(SimpleGraph(tm[end]), vss)
    l_g, l_vmap = my_merge_vertices_list(SimpleGraph(tm[1]), vss)
    @assert(lc_vmap == l_vmap)

    c_g, c_vmap = my_merge_vertices_list(lc_g, connected_components(l_g))

    l_ext = [v for v in 1:nv(l_g) if has_edge(l_g, v, v)]
    c_ext = [v for v in 1:nv(c_g) if has_edge(c_g, v, v)]
    if isempty(l_ext)
        (length(c_ext) == 1) || error("There should be a unique self-loop in the odd case")
        rem_edge!(c_g, c_ext[1], c_ext[1])
        w, wmap = rt_to_dyck(sort_rt(rtgraph_to_rt(c_g, c_ext[1])))
        return (true, w), wmap[c_vmap[lc_vmap[axes(tm[1], 1)]]]
    else
        (length(c_ext) == 0) || error("There should be no self-loops in the even case")
        w, wmap = rt_to_dyck(sort_rt(rtgraph_to_rt(c_g, only(unique(c_vmap[l_ext])))))
        return (false, w), wmap[c_vmap[lc_vmap[axes(tm[1], 1)]]]
    end
end

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
                if p == sort([sort((es[i], es[j])) for (i, j) in dyck_to_pairs(W)])
                    push!(fW[f], W)
                end
            end
        end

        for Ws in product1d(fW...)
            nwt = NWT(incc, n, [Ws...], [])
            if curve_type(nwt)[1] == (true, [])
                push!(result, nwt)
            end
        end
    end

    return result
end

@inline function liftedindex(N, D, n, d)
    acc = [1; accumulate(*, D .+ 1)[1:end-1]]
    return n + N * sum(d .* acc)
end

function liftedgraph(tm::Vector{BitMatrix}, D)
    k = length(tm)
    N = size(tm[1], 1)

    lv = product1d(1:N, [0:d for d in D]...)
    le = zeros(Bool, length(lv), length(lv))

    for i in 1:k
        for r1 in 1:N, r2 in r1:N
            if tm[i][r1, r2]
                start_t = tuple([0 for d in D]...)
                delta_t = tuple(I[1:k, i]...)
                end_t = tuple([d for d in D]...) .- delta_t

                for t in product1d((start_t[i]:end_t[i] for i in 1:k)...)
                    le[liftedindex(N, D, r1, t), liftedindex(N, D, r2, t .+ delta_t)] = 1
                    le[liftedindex(N, D, r2, t), liftedindex(N, D, r1, t .+ delta_t)] = 1
                end
            end
        end
    end

    return lv, le
end

function is_bezout_basic(nwt::NWT, D::Vector{Int})
    if D[end] == 0
        return true
    end
    for (i, d) in enumerate(D[1:end-1])
        if d != 0
            es = findall(E -> E.i == i, nwt.cc.E)
            n_t = sum(nwt.n[es])
            if n_t > d * D[end] || (d * D[end] - n_t) % 2 != 0
                return false
            end
        end
    end

    return true
end

# function linedistance_matmul(A::Matrix{Vector{NTuple{k,Int}}}, B::Matrix{Vector{NTuple{k,Int}}}, D::Vector{Int})::Matrix{Vector{NTuple{k,Int}}} where k
#     m, n = size(A)
#     n2, p = size(B)
#     @assert n == n2

#     C = [Vector{NTuple{k,Int}}() for _ in 1:m, _ in 1:p]
#     for i in 1:m, j in 1:p
#         C[i, j] = [x for x in unique!(vcat([[[x .+ y for x in A[i, k], y in B[k, j]]...] for k in 1:n]...))
#                    if all(x .<= D)]
#     end
#     return C
# end

@inline function linedistance_matmul(A, B, D, k, m, n, p)
    C = [NTuple{k,Int}[] for _ in 1:m, _ in 1:p]
    for i in 1:m, j in 1:p
        buf = NTuple{k,Int}[]
        for l in 1:n
            for x in A[i, l], y in B[l, j]
                s = x .+ y
                if all(s .<= D)
                    push!(buf, s)
                end
            end
        end
        C[i, j] = unique!(buf)
    end
    return C
end

@inline function minimal_elements(v)
    filter(t1 -> !any(t2 -> t2 !== t1 && all(t2 .<= t1), v), v)
end

@inline function feasible(s, D)
    for i in eachindex(s)
        s[i] > D[i] && return false
        (D[i] - s[i]) % 2 != 0 && return false
    end
    return true
end

@inline function distinct_pairs_sums(v)
    [v[i] .+ v[j] for i in eachindex(v) for j in i+1:length(v)]
end

function linedistance_matrix(nwt::NWT, D::Vector{Int})
    tm = totalmatrix(nwt)
    k = length(tm)
    N = size(tm[1], 1)

    A = [Vector{NTuple{k,Int}}() for _ in 1:N, _ in 1:N]
    for i in 1:k
        for (j, v) in enumerate(tm[i])
            if v
                push!(A[j], tuple(I[1:k, i]...))
            end
        end
    end

    B = A

    powers = [Set{NTuple{k,Int}}(A[i, j]) for i in 1:N, j in 1:N]
    for i in 1:N
        push!(powers[i, i], ntuple(_ -> 0, k))
    end

    for _ in 1:sum(D)
        B = linedistance_matmul(A, B, D, k, N, N, N)
        for i in 1:N, j in 1:N
            union!(powers[i, j], B[i, j])
        end
    end

    return [minimal_elements(filter(s -> feasible(s, D), distinct_pairs_sums(collect(p)))) for p in powers]
end

function is_bezout_order_zero(nwt::NWT, D::Vector{Int})
    if !is_bezout_basic(nwt, [D; 0])
        return false
    end

    return all([!isempty(v) for v in linedistance_matrix(nwt, D)])
end

function is_bezout_order_one(nwt::NWT, D::Vector{Int})
    if !is_bezout_order_zero(nwt, D)
        return false
    end

    tm = groundmatrix(nwt)
    N = size(tm[1], 1)
    dps = distinct_pairs(axes(tm[1], 1))

    lv, le = liftedgraph(tm, D)
    lg = SimpleDiGraph(le)
    ds = [dijkstra_shortest_paths(lg, liftedindex(N, D, f, [0 for d in D]), allpaths=true) for f in 1:N]
    gcc = groundcc(nwt)

    while !isempty(dps)
        (f1, f2) = first(dps)

        start_f = liftedindex(N, D, f1, [0 for d in D])
        for middle_f in findall(i -> lv[i][1] == f2 && ds[f1].dists[i] < typemax(Int), eachindex(lv))
            for end_f in [liftedindex(N, D, f1, ds) for ds in product1d([collect(d:(-2):m) for (d, m) in zip(D, lv[middle_f][2:end])]...)]
                end2_t = lv[end_f][2:end] .- lv[middle_f][2:end]
                end2_f = liftedindex(N, D, f1, end2_t)
                if ds[f2].dists[end2_f] < typemax(Int)
                    start2_f = liftedindex(N, D, f2, [0 for d in D])
                    paths1 = all_shortest_paths(start_f, middle_f, ds[f1].predecessors)
                    paths2 = all_shortest_paths(start2_f, end2_f, ds[f2].predecessors)
                    paths = [first.(lv[[path1[1:end-1]; path2]]) for (path1, path2) in Iterators.product(paths1, paths2)]

                    for path in paths
                        for refined_nwt in add_line(gcc, path)
                            if is_bezout_order_zero(refined_nwt, [D; 1])
                                setdiff!(dps, [(lv[f1][1], lv[f2][1]) for (f1, f2) in distinct_pairs(path)])
                                @goto boo_ok
                            end
                        end
                    end
                end
            end
        end

        return false

        @label boo_ok
    end

    return true
end

function nonneg_sum(N, n)
    if n == 1
        return [[N]]
    end
    result = Vector{Vector{Int}}()
    for k in 0:N
        for tail in nonneg_sum(N - k, n - 1)
            push!(result, [k; tail])
        end
    end
    return result
end

function all_n(cc::CellComplex, D::Vector{Int})::Vector{Vector{Int}}
    result = []

    ie = [Int[] for _ in 1:cc.k]
    for (e, E) in enumerate(cc.E)
        push!(ie[E.i], e)
    end

    for n_i in product1d([(d*D[end]):(-2):0 for d in D[1:end-1]]...)
        for n_ie in product1d([nonneg_sum(n_i[i], length(ie[i])) for i in 1:cc.k]...)
            n_e = [0 for _ in eachindex(cc.E)]
            for i in 1:cc.k
                for (e, n) in zip(ie[i], n_ie[i])
                    n_e[e] = n
                end
            end
            for (f, F) in enumerate(cc.F)
                n_f = sum([n_e[e] for (s, e) in F])
                if n_f % 2 != 0
                    @goto next_n_ie
                end
            end
            push!(result, n_e)

            @label next_n_ie
        end
    end

    return result
end

function _d6(i::Int, j::Int)::Tuple{Bool,Vector{Bool}}
    if j == 0
        (false, vcat([[1, 0] for _ in 1:i]...))
    else
        (false, [vcat([[1, 0] for _ in 1:i]...); 1; vcat([[1, 0] for _ in 1:j]...); 0])
    end
end

curve_type_list = Vector{Tuple{Bool,Vector{Bool}}}[
    [(true, [])],
    [(false, []), (false, [1, 0])],
    [(true, []), (true, [1, 0])],
    [
        (false, []), (false, [1, 0]), (false, [1, 0, 1, 0]), (false, [1, 0, 1, 0, 1, 0]),
        (false, [1, 0, 1, 0, 1, 0, 1, 0]), (false, [1, 1, 0, 0])
    ],
    [
        (true, []), (true, [1, 0]), (true, [1, 0, 1, 0]), (true, [1, 1, 0, 0]), (true, [1, 0, 1, 0, 1, 0]),
        (true, [1, 0, 1, 0, 1, 0, 1, 0]), (true, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]), (true, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0])
    ],
    [
        [_d6(o, 0) for o in 0:9]..., # <n>
        [_d6(i, o - 1 - i) for o in 2:9 for i in 0:(o-2)]..., # <n ⊔ 1<m>>
        _d6(10, 0), _d6(8, 1), _d6(5, 4), _d6(4, 5), _d6(1, 8), _d6(0, 9),
        _d6(9, 1), _d6(5, 5), _d6(1, 9),
        (false, [1, 1, 1, 0, 0, 0]),
    ]
]

# Constraint satisfaction problem
function all_sections(vvs, c)
    n = length(vvs)
    results = Vector{Vector{Any}}()
    domains = [copy(vvs[i]) for i in 1:n]  # mutable working domains

    function backtrack(i, current, domains)
        if i > n
            push!(results, copy(current))
            return
        end
        for v in domains[i]
            # Check consistency with already-chosen elements
            if all(c(current[j], v) for j in 1:i-1)
                # Forward checking: prune future domains
                new_domains = [copy(domains[l]) for l in 1:n]
                for l in i+1:n
                    filter!(w -> c(v, w), new_domains[l])
                end
                # Only recurse if no domain is empty
                if all(!isempty(new_domains[l]) for l in i+1:n)
                    push!(current, v)
                    backtrack(i + 1, current, new_domains)
                    pop!(current)
                end
            end
        end
    end

    backtrack(1, [], domains)
    return results
end

function all_bezout_zero(cc::CellComplex, D::Vector{Int})::Vector{NWT}
    result = []
    i = 0
    for n in all_n(cc, D)
        n_f = [sum([n[e] for (s, e) in F]) for F in cc.F]
        for W in product1d([all_dyck_words(div(n_f[f], 2)) for f in eachindex(cc.F)]...)
            nwt0 = NWT(cc, n, [W...], [])

            # if nwt0.W != [[1, 1, 1, 0, 1, 0, 0, 0], [1, 1, 1, 1, 0, 0, 0, 0]]
            #     continue
            # end

            @assert(is_bezout_basic(nwt0, D))

            ct, ct_vmap = curve_type(nwt0)

            if !(ct in curve_type_list[D[end]])
                continue
            end

            if !(is_bezout_order_zero(nwt0, D))
                continue
            end

            ctregions = [findall(==(v), ct_vmap) for v in unique!(sort(ct_vmap))]
            ld = [minimum(last.(v)) for v in linedistance_matrix(nwt0, D)]

            # println("ld: ", ld)

            for big in curve_type_list[D[end]]
                bigg, _ = dyck_to_rtgraph(big[2])
                compls = rtree_nonisomorphic_embeddings(big[2], ct[2])
                for compl in compls
                    cut = difference(bigg, my_induced_subgraph(bigg, compl))
                    attachments = []
                    depths = []
                    diameters = []
                    for v in compl
                        depth = rtgraph_depth((cut, v))
                        if depth > 0
                            diameter = rtgraph_diameter((cut, v))
                            push!(attachments, v)
                            push!(depths, depth)
                            push!(diameters, diameter)
                        end
                    end
                    # println("attachments: ", attachments)
                    # println("depths: ", depths)
                    # println("diameters: ", diameters)
                    # println("cut: ", collect(edges(cut)))
                    feasable = []
                    for (i, v) in enumerate(attachments)
                        push!(feasable, [(i, r) for r in ctregions[v] if maximum(ld[r, :]) <= D[end] - depth[i] && ld[r, r] <= D[end] - diameters[i]])
                    end
                    for section in all_sections(feasable, (((i1, r1), (i2, r2)) -> ld[r1, r2] <= D[end] - depths[i1] - depths[i2]))
                        T = [(r, rtgraph_to_dyck((cut, attachments[i]))) for (i, r) in section]
                        push!(result, NWT(cc, n, [W...], sort(T, by=first)))
                        i += 1
                        if i % 10000 == 0
                            println(i)
                        end
                    end
                end
            end
        end
    end

    return result
end

# function viro(points::Vector{Tuple{Int,Int}}, triangles::Vector{Vector{Int}}, signs::Vector{Int})::NWT

# end

function symmetry_yzx(nwt)
    return forget_edges(nwt, three_lines, [2, 3, 1, 5, 6, 4],
        Vector{Tuple{Bool,Int,Int}}[[(0, 1, 2), (0, 1, 3), (0, 1, 1)], [(0, 3, 2), (0, 3, 3), (0, 3, 1)], [(0, 4, 2), (0, 4, 3), (0, 4, 1)], [(0, 2, 2), (0, 2, 3), (0, 2, 1)]])
end

function symmetry_yxz(nwt)
    return forget_edges(nwt, three_lines, [2, 1, 3, 5, 4, 6],
        Vector{Tuple{Bool,Int,Int}}[[(1, 1, 2), (1, 1, 1), (1, 1, 3)], [(1, 3, 2), (1, 3, 1), (1, 3, 3)], [(1, 2, 2), (1, 2, 1), (1, 2, 3)], [(1, 4, 2), (1, 4, 1), (1, 4, 3)]])
end

function symmetry_xYZ(nwt)
    return forget_edges(nwt, three_lines, [1, 5, 6, 4, 2, 3],
        Vector{Tuple{Bool,Int,Int}}[[(0, 2, 1), (0, 2, 2), (0, 2, 3)], [(0, 1, 1), (0, 1, 2), (0, 1, 3)], [(0, 4, 1), (0, 4, 2), (0, 4, 3)], [(0, 3, 1), (0, 3, 2), (0, 3, 3)]])
end

function symmetry_XyZ(nwt)
    return forget_edges(nwt, three_lines, [4, 2, 6, 1, 5, 3],
        Vector{Tuple{Bool,Int,Int}}[[(0, 3, 1), (0, 3, 2), (0, 3, 3)], [(0, 4, 1), (0, 4, 2), (0, 4, 3)], [(0, 1, 1), (0, 1, 2), (0, 1, 3)], [(0, 2, 1), (0, 2, 2), (0, 2, 3)]])
end

function symmetries(nwt::NWT)::Vector{NWT}
    syms = NWT[nwt]
    syms = [syms; symmetry_yzx.(syms); symmetry_yzx.(symmetry_yzx.(syms))]
    syms = [syms; symmetry_yxz.(syms)]
    syms = [syms; symmetry_xYZ.(syms)]
    syms = [syms; symmetry_XyZ.(syms)]
    return syms
end

function strictsymmetries(nwt::NWT)::Vector{NWT}
    syms = NWT[nwt]
    syms = [syms; symmetry_xYZ.(syms)]
    syms = [syms; symmetry_XyZ.(syms)]
    return syms
end

function remove_symmetries(nwts::Vector{NWT}, symm_function)::Vector{NWT}
    dict = Dict{NWT,Bool}()
    for nwt in nwts
        dict[nwt] = true
    end

    for (nwt, v) in dict
        if v
            for symm in setdiff(symm_function(nwt), [nwt])
                dict[symm] = false
            end
        end
    end

    return [nwt for (nwt, v) in dict if v]
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

    @test forget_edges(marge, two_lines_three_points, [1, 3, 4, 6], Vector{Tuple{Bool,Int,Int}}[[(0, 1, 1), (0, 3, 1), (0, 3, 3), (0, 1, 3)], [(0, 2, 1), (0, 4, 1), (0, 4, 3), (0, 2, 3)]]) == marge_e
    @test forget_vertices(marge_e, two_lines, [2], [[1, 3], [4, 2]]) == marge_ev
    @test forget_edges(marge_ev, one_line, [2], Vector{Tuple{Bool,Int,Int}}[[(0, 2, 2), (0, 1, 2)]]) == marge_eve

    @test forget_edges(fine, two_lines_three_points, [1, 3, 4, 6], Vector{Tuple{Bool,Int,Int}}[[(0, 1, 1), (0, 3, 1), (0, 3, 3), (0, 1, 3)], [(0, 2, 1), (0, 4, 1), (0, 4, 3), (0, 2, 3)]]) == fine_e
end

@testset "nonline" begin
    nonlinecc = CellComplex(2, 2, [(1, 2, 1), (2, 1, 1), (2, 1, 2), (2, 1, 2)], [[(0, 2), (1, 3), (0, 2), (1, 4)], [(0, 1), (0, 3)], [(0, 1), (0, 4)]])

    @test iscellcomplex(nonlinecc)

    nonlinenwt = NWT(nonlinecc, [3, 2, 3, 3], [[1, 1, 0, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]], [])

    @test isNWT(nonlinenwt)

    nonlinenwt2 = NWT(nonlinecc, [0, 2, 0, 0], [[1, 0, 1, 0], [], []], [])

    @test isNWT(nonlinenwt2)
end