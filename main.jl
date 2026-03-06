using Test
using Graphs
using LinearAlgebra: I
using JSON3, CodecZlib
using DataStructures: counter
using Combinatorics: multiset_combinations

include("dyck.jl")
include("my_graphs.jl")

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
# k: number of curves
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
Base.hash(nwt::NWT) = hash(([nwt.cc.k, nwt.cc.nV, nwt.cc.E, nwt.cc.F], nwt.n, nwt.W, nwt.T))

Base.show(io::IO, nwt::NWT) = print(io, (nwt.n, nwt.W, nwt.T))

struct NWT_from_jsonl
    n::Vector{UInt8}
    W::Vector{Vector{Bool}}
    T::Vector{@NamedTuple{g::UInt8, T::Vector{Bool}}}
    itriang::Int32
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
            gmat[end][g1, g2] = 1
            gmat[end][g2, g1] = 1
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
                gmat[cc.E[e].i][g1, g2] = 1
                gmat[cc.E[e].i][g2, g1] = 1
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

    for (f, F) in enumerate(out_cc.F)
        (length(Fmap[f]) == length(F)) || error("Not enough borders.")
    end

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

    in_figs = fi_ground_sequence(in_nwt)

    merged_fgs = Vector{Int}[]
    for f in eachindex(out_cc.F)
        push!(merged_fgs, [])
        for ((in_s, in_f, in_i), (out_s, out_e)) in zip(Fmap[f], out_cc.F[f])
            (in_cc.F[in_f][in_i].e == Emap[out_e]) || error("Fmap format not correct.")
            append!(merged_fgs[end], [removed_map[g] for g in reversebool(xor(in_s, out_s), in_figs[(in_f, in_i)])])
        end
    end
    merged_fgs = remove_consecutive.(merged_fgs)
    out_D = sequence_to_dyck.(merged_fgs)

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
        D = rtgraph_to_dyck((flforest, vmap[removed_components[mg][1]]))
        if !isempty(D)
            push!(out_T, (g=ng, T=D))
        end
    end

    return NWT(out_cc, in_nwt.n[Emap], out_D, sort(out_T, by=first))
end

# input:
# - in_nwt: the input nwt, where in_cc = in_nwt.cc
# - out_cc: cc of the output; must have the a subset of the vertices, the merged edges, and the same faces
# - Emap: edge e in out_cc is made by merging edges Emap[e] in in_cc
# output:
# - out_nwt: obtained by removing some vertices from in_nwt
function forget_vertices(in_nwt::NWT, out_cc::CellComplex, Emap::Vector{Vector{Int}})
    out_n = [sum(in_nwt.n[Emap[e]]) for e in eachindex(out_cc.E)]

    return NWT(out_cc, out_n, in_nwt.W, in_nwt.T)
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
        w, wmap = rt_to_dyck(sort_rt(rtgraph_to_rt((c_g, c_ext[1]))))
        return (true, w), [wmap[v] for v in c_vmap[lc_vmap[axes(tm[1], 1)]]]
    else
        (length(c_ext) == 0) || error("There should be no self-loops in the even case")
        w, wmap = rt_to_dyck(sort_rt(rtgraph_to_rt((c_g, only(unique(c_vmap[l_ext]))))))
        return (false, w), [wmap[v] for v in c_vmap[lc_vmap[axes(tm[1], 1)]]]
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
    [
        (true, [])],
    [
        (false, []),
        (false, [1, 0])],
    [
        (true, []),
        (true, [1, 0])],
    [
        (false, []),
        (false, [1, 0]),
        (false, [1, 0, 1, 0]),
        (false, [1, 0, 1, 0, 1, 0]),
        (false, [1, 0, 1, 0, 1, 0, 1, 0]),
        (false, [1, 1, 0, 0])
    ],
    [
        (true, []),
        (true, [1, 0]),
        (true, [1, 1, 0, 0]),
        (true, [1, 0, 1, 0, 1, 0]),
        (true, [1, 0, 1, 0, 1, 0, 1, 0]),
        (true, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]),
        (true, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]),
        (true, [1, 0, 1, 0])
    ],
    [
        [_d6(o, 0) for o in 0:10]..., # <n>
        [_d6(i, o - 1 - i) for o in 2:9 for i in 0:(o-2)]..., # <n ⊔ 1<m>>
        _d6(10, 0),
        _d6(8, 1),
        _d6(5, 4),
        _d6(4, 5),
        _d6(1, 8),
        _d6(0, 9),
        _d6(9, 1),
        _d6(5, 5),
        _d6(1, 9),
        (false, [1, 1, 1, 0, 0, 0])]]

function all_nonnested_floating(nwt::NWT, D::Vector{Int})
    result = NWT[]

    ld = [minimum(last.(v)) for v in linedistance_matrix(nwt, D)]
    N = size(ld, 1)
    _, ct_vmap = curve_type(nwt)
    d = D[end]
    maxevenovals = div(d^2 - 3d + 4, 2) - ((d % 2 == 0) ? 0 : 1)
    maxmissingovals = maxevenovals - (length(unique(ct_vmap)) - 1)

    floating = Int[]
    for r in 1:N
        if maximum(ld[r, :]) <= d - 2
            if ld[r, r] <= d - 4
                append!(floating, fill(r, maxmissingovals))
            else
                push!(floating, r)
            end
        end
    end

    floatingVector = zeros(Int, N)
    for i in 0:maxmissingovals
        for c in multiset_combinations(floating, i)
            for j1 in 1:i
                for j2 in j1+1:i
                    if ld[c[j1], c[j2]] > d - 4
                        @goto skip
                    end
                end
            end

            fill!(floatingVector, 0)
            for v in c
                floatingVector[v] += 1
            end

            T = Tuple{Int,Vector{Bool}}[]
            for (i, t) in enumerate(floatingVector)
                if t > 0
                    push!(T, (i, vcat([[1, 0] for _ in 1:t]...)))
                end
            end

            push!(result, NWT(nwt.cc, nwt.n, nwt.W, T))

            @label skip
        end
    end

    return result
end

function all_floatless_bezout_zero(cc::CellComplex, D::Vector{Int})::Vector{NWT}
    result = []

    for n in all_n(cc, D)
        n_f = [sum([n[e] for (s, e) in F]) for F in cc.F]
        for W in product1d([all_dyck_words(div(n_f[f], 2)) for f in eachindex(cc.F)]...)
            nwt0 = NWT(cc, n, [W...], [])

            @assert(is_bezout_basic(nwt0, D))

            ct, _ = curve_type(nwt0)

            if !(ct in curve_type_list[D[end]])
                continue
            end

            if !(is_bezout_order_zero(nwt0, D))
                continue
            end

            push!(result, nwt0)
        end
    end

    return result
end

function viro_patchworking!(
    deg::Int,
    itriang::Int,
    vertices::Vector{Tuple{Int,Int}},
    triang::Vector{Tuple{Int,Int,Int}},
    dict_nwt::Dict{NWT,Tuple{Int,Vector{Bool}}})

    function viro_edge_type(s, d)
        S = vertices[s]
        D = vertices[d]
        if S[1] == D[1] == 0
            return 1
        end
        if S[2] == D[2] == 0
            return 2
        end
        if S[1] + S[2] == D[1] + D[2] == deg
            return 3
        end
        return 4
    end

    edges = @NamedTuple{s::Int, d::Int, i::Int}[]
    faces = Vector{@NamedTuple{s::Bool, e::Int}}[]
    for tr in triang
        new_edges = [(tr[1], tr[2]), (tr[2], tr[3]), (tr[3], tr[1])]
        new_signs = [ne[1] > ne[2] for ne in new_edges]
        new_indices = []
        for (ns, ne) in zip(new_signs, new_edges)
            (s, d) = reversebool(ns, ne)
            index = findfirst(e -> (e.s, e.d) == (s, d), edges)
            if isnothing(index)
                push!(edges, (s, d, viro_edge_type(s, d)))
                index = length(edges)
            end
            push!(new_indices, index)
        end
        push!(faces, collect(zip(new_signs, new_indices)))
    end

    cc = CellComplex(4, length(vertices), edges, faces)

    x_edges = sort([e for (e, E) in enumerate(edges) if E.i == 1], by=(e -> edges[e].s), rev=true)
    y_edges = sort([e for (e, E) in enumerate(edges) if E.i == 2], by=(e -> edges[e].s))
    z_edges = sort([e for (e, E) in enumerate(edges) if E.i == 3], by=(e -> edges[e].s))
    Emap = [x_edges; y_edges; z_edges]
    Fmap = [vcat([[(false, f, i) for (f, F) in enumerate(faces) for (i, (s, e)) in enumerate(F) if e == be]
                  for be in Emap]...)]

    one_F = [[(s=(edges[e].i == 1), e=out_e) for (out_e, e) in enumerate(Emap)]]
    one_cc = CellComplex(3, length(vertices), cc.E[Emap], one_F)

    fixedvertices = [findfirst(v -> v .% 2 == s, vertices) for s in [(0, 0), (1, 0), (0, 1), (1, 1)]]
    filter!(!isnothing, fixedvertices)
    fixedvertices = fixedvertices[1:min(3, length(fixedvertices))]
    unfixedvertices = setdiff(eachindex(vertices), fixedvertices)

    for unfixedsign in Iterators.product([false:true for _ in eachindex(unfixedvertices)]...)
        sign_xyz = [false for _ in vertices]
        for (i, s) in zip(unfixedvertices, unfixedsign)
            sign_xyz[i] = s
        end
        sign_xYZ = [xor(s, x % 2 != 0) for ((x, y), s) in zip(vertices, sign_xyz)]
        sign_XyZ = [xor(s, y % 2 != 0) for ((x, y), s) in zip(vertices, sign_xyz)]
        sign_XYz = [xor(s, x % 2 != 0) for ((x, y), s) in zip(vertices, sign_XyZ)]

        pr_n = [
            sum([(sign_xyz[E.s] == sign_xyz[E.d]) ? 0 : 1 for E in edges[x_edges]]),
            sum([(sign_xyz[E.s] == sign_xyz[E.d]) ? 0 : 1 for E in edges[y_edges]]),
            sum([(sign_xyz[E.s] == sign_xyz[E.d]) ? 0 : 1 for E in edges[z_edges]]),
            sum([(sign_XYz[E.s] == sign_XYz[E.d]) ? 0 : 1 for E in edges[x_edges]]),
            sum([(sign_xYZ[E.s] == sign_xYZ[E.d]) ? 0 : 1 for E in edges[y_edges]]),
            sum([(sign_XyZ[E.s] == sign_XyZ[E.d]) ? 0 : 1 for E in edges[z_edges]])]

        pr_W = Vector{Bool}[]
        pr_T = @NamedTuple{g::Int, T::Vector{Bool}}[]
        nground = 0
        for sign in [sign_xyz, sign_xYZ, sign_XyZ, sign_XYz]
            n = Int[(sign[E.s] == sign[E.d]) ? 0 : 1 for E in edges]
            n_f = [sum([n[e] for (s, e) in F]) for F in cc.F]
            Ws = product1d([all_dyck_words(div(n_f[f], 2)) for f in eachindex(cc.F)]...)
            nwt = NWT(cc, n, [only(Ws)...], [])

            one_nwt = forget_edges(nwt, one_cc, Emap, Fmap)
            push!(pr_W, one_nwt.W[1])
            append!(pr_T, [(g + nground, T) for (g, T) in one_nwt.T])
            nground += div(length(pr_W[end]), 2) + 1
        end

        pr_nwt = NWT(lines_xyz, pr_n, pr_W, pr_T)
        if !haskey(dict_nwt, pr_nwt)
            for symm in symmetries(pr_nwt)
                dict_nwt[symm] = (itriang, sign_xyz)
            end
        end
    end
end

function symmetry_yzx(nwt)
    return forget_edges(nwt, lines_xyz, [2, 3, 1, 5, 6, 4],
        Vector{Tuple{Bool,Int,Int}}[[(0, 1, 2), (0, 1, 3), (0, 1, 1)], [(0, 3, 2), (0, 3, 3), (0, 3, 1)], [(0, 4, 2), (0, 4, 3), (0, 4, 1)], [(0, 2, 2), (0, 2, 3), (0, 2, 1)]])
end

function symmetry_yxz(nwt)
    return forget_edges(nwt, lines_xyz, [2, 1, 3, 5, 4, 6],
        Vector{Tuple{Bool,Int,Int}}[[(1, 1, 2), (1, 1, 1), (1, 1, 3)], [(1, 3, 2), (1, 3, 1), (1, 3, 3)], [(1, 2, 2), (1, 2, 1), (1, 2, 3)], [(1, 4, 2), (1, 4, 1), (1, 4, 3)]])
end

function symmetry_xYZ(nwt)
    return forget_edges(nwt, lines_xyz, [1, 5, 6, 4, 2, 3],
        Vector{Tuple{Bool,Int,Int}}[[(0, 2, 1), (0, 2, 2), (0, 2, 3)], [(0, 1, 1), (0, 1, 2), (0, 1, 3)], [(0, 4, 1), (0, 4, 2), (0, 4, 3)], [(0, 3, 1), (0, 3, 2), (0, 3, 3)]])
end

function symmetry_XyZ(nwt)
    return forget_edges(nwt, lines_xyz, [4, 2, 6, 1, 5, 3],
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

function symmetries(nwts::Vector{NWT})::Vector{NWT}
    return unique!(vcat(symmetries.(nwts)...))
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

function translate_up(nwt::NWT)::NWT
    @assert(nwt.cc == lines_xyz)

    nwt_xz = forget_edges(
        nwt,
        lines_xz_3points,
        [1, 3, 4, 6],
        [
            [(false, 1, 1), (false, 3, 1), (false, 3, 3), (false, 1, 3)],
            [(false, 2, 1), (false, 4, 1), (false, 4, 3), (false, 2, 3)]])

    return NWT(
        lines_xyz,
        [
            nwt_xz.n[1] + nwt_xz.n[3],
            nwt_xz.n[4],
            nwt_xz.n[2],
            0,
            nwt_xz.n[2],
            nwt_xz.n[4]],
        [
            nwt_xz.W[1],
            nwt_xz.W[2],
            vcat(fill(true, nwt_xz.n[4]), fill(false, nwt_xz.n[4])),
            vcat(fill(true, nwt_xz.n[2]), fill(false, nwt_xz.n[2])),
        ],
        nwt_xz.T)
end

### TESTS ###

lines_xyz = CellComplex(
    3,
    3,
    [(2, 3, 1), (3, 1, 2), (1, 2, 3), (2, 3, 1), (3, 1, 2), (1, 2, 3)],
    [[(0, 1), (0, 2), (0, 3)], [(0, 1), (0, 5), (0, 6)], [(0, 4), (0, 2), (0, 6)], [(0, 4), (0, 5), (0, 3)]])

lines_xz_3points = CellComplex(
    2,
    3,
    [(2, 3, 1), (1, 2, 2), (2, 3, 1), (1, 2, 2)],
    [[(0, 1), (1, 3), (1, 4), (0, 2)], [(0, 1), (1, 3), (1, 2), (0, 4)]])

lines_xz = CellComplex(
    2,
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(0, 1), (0, 2)], [(0, 1), (1, 2)]]
)

lines_yz_3points = CellComplex(
    2,
    3,
    [(3, 1, 1), (1, 2, 2), (3, 1, 1), (1, 2, 2)],
    [[(1, 1), (0, 3), (0, 4), (1, 2)], [(1, 1), (0, 3), (0, 2), (1, 4)]])

lines_yz = CellComplex(
    2,
    1,
    [(1, 1, 1), (1, 1, 2)],
    [[(0, 1), (0, 2)], [(0, 1), (1, 2)]]
)

line_z_point_010 = CellComplex(
    1,
    1,
    [(1, 1, 1)],
    [[(0, 1), (0, 1)]]
)

fine = NWT(
    lines_xyz,
    [4, 1, 3, 0, 1, 1],
    [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 0], [1, 0], [1, 1, 0, 0]],
    [(1, [1, 0, 1, 0]), (10, [1, 0])])

nonviro = NWT(
    lines_xyz,
    [3, 4, 1, 1, 0, 3],
    [[1, 0, 1, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 1, 0, 0], [1, 0]],
    [])

marge = NWT(
    lines_xyz,
    [4, 7, 3, 2, 1, 3],
    [[1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 0, 0], [1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0], [1, 0, 1, 1, 0, 0]],
    [(1, [1, 0]), (5, [1, 0]), (17, [1, 0]), (18, [1, 0])]
)

@testset "CellComplex" begin
    @test iscellcomplex(lines_xyz)
    @test iscellcomplex(lines_xz)
    @test iscellcomplex(line_z_point_010)
end

@testset "isNWT" begin
    @test isNWT(fine)
    @test isNWT(nonviro)
    @test isNWT(marge)
end

@testset "remove_edges" begin
    marge_e = NWT(
        lines_xz_3points,
        [4, 3, 2, 3],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_ev = NWT(
        lines_xz,
        [6, 6],
        [[1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    marge_eve = NWT(
        line_z_point_010,
        [6],
        [[1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0]],
        [(1, [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0])]
    )

    fine_e = NWT(
        lines_xz_3points,
        [4, 3, 0, 1],
        [[1, 1, 0, 0, 1, 0, 1, 0], [1, 0, 1, 1, 0, 1, 0, 0]],
        [(1, [1, 0, 1, 0, 1, 0])]
    )

    @test forget_edges(marge, lines_xz_3points, [1, 3, 4, 6], Vector{Tuple{Bool,Int,Int}}[[(0, 1, 1), (0, 3, 1), (0, 3, 3), (0, 1, 3)], [(0, 2, 1), (0, 4, 1), (0, 4, 3), (0, 2, 3)]]) == marge_e
    @test forget_vertices(marge_e, lines_xz, [[1, 3], [4, 2]]) == marge_ev
    @test forget_edges(marge_ev, line_z_point_010, [2], Vector{Tuple{Bool,Int,Int}}[[(0, 2, 2), (0, 1, 2)]]) == marge_eve

    @test forget_edges(fine, lines_xz_3points, [1, 3, 4, 6], Vector{Tuple{Bool,Int,Int}}[[(0, 1, 1), (0, 3, 1), (0, 3, 3), (0, 1, 3)], [(0, 2, 1), (0, 4, 1), (0, 4, 3), (0, 2, 3)]]) == fine_e
end

@testset "nonline" begin
    nonlinecc = CellComplex(2, 2, [(1, 2, 1), (2, 1, 1), (2, 1, 2), (2, 1, 2)], [[(0, 2), (1, 3), (0, 2), (1, 4)], [(0, 1), (0, 3)], [(0, 1), (0, 4)]])

    @test iscellcomplex(nonlinecc)

    nonlinenwt = NWT(nonlinecc, [3, 2, 3, 3], [[1, 1, 0, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]], [])

    @test isNWT(nonlinenwt)

    nonlinenwt2 = NWT(nonlinecc, [0, 2, 0, 0], [[1, 0, 1, 0], [], []], [])

    @test isNWT(nonlinenwt2)
end