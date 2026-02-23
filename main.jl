using Test
using Graphs
using LinearAlgebra: I

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

# function curvematrix(nwt)
#     grmatrix = groundmatrix(nwt)
#     ng = size(grmatrix[1], 1)
#     flforest = floatingforest(nwt)
#     cmatrix = adjacency_matrix(flforest)
#     cmatrix[1:ng, 1:ng] += grmatrix[end]
#     return cmatrix
# end

function totalmatrix(nwt::NWT)
    gm = groundmatrix(nwt)
    n = size(gm[1], 1)

    flforest = floatingforest(nwt)
    N = nv(flforest)

    tm = [falses(N, N) for _ in eachindex(gm)]
    for i in eachindex(gm)
        tm[i][1:n, 1:n] .= gm[i]
    end
    tm[end] += adjacency_matrix(flforest)

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
            # println("outF: ", outF)
            continue
        end
        B[1] = [pop!(B); B[1]]
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
    # gm = groundmatrix(nwt)
    # n = size(gm[1], 1)

    # cm = curvematrix(nwt)
    # N = size(cm, 1)

    # lm = falses(N, N)
    # lm[1:n, 1:n] += gm[1]
    tm = totalmatrix(nwt)

    vss = connected_components(SimpleGraph(max.(tm[2:end-1]...)))
    lc_g, lc_vmap = my_merge_vertices_list(SimpleGraph(tm[end]), vss)
    l_g, l_vmap = my_merge_vertices_list(SimpleGraph(tm[1]), vss)
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


# pseudoline_multable = [
#     0 0 0 0 0 0 0 0 0;
#     0 0 7 6 1 8 0 0 0;
#     0 7 0 5 2 0 8 0 0;
#     0 6 5 0 3 0 0 8 0;
#     0 1 2 3 4 5 6 7 8;
#     0 8 0 0 5 0 0 0 0;
#     0 0 8 0 6 0 0 0 0;
#     0 0 0 8 7 0 0 0 0;
#     0 0 0 0 8 0 0 0 0
# ]

# pseudoline_weight = [0, 1, 1, 1, 0, 2, 2, 2, 3]

# function pseudoline_mul!(u, v, w)
#     n = size(v, 1)
#     for i in 1:n, j in 1:n
#         u[i, j] = maximum(pseudoline_multable[v[i, k]+1, w[k, j]+1] for k in 1:n)
#     end
#     return u
# end

# function pseudoline_min_path(u)
#     mp = 100
#     n = length(u)
#     for i in 1:n, j in i+1:n
#         if u[i] > 0 && u[j] > 0 && abs(u[i] - u[j]) == 4
#             mp = min(mp, (i - 1) - pseudoline_weight[u[i]+1] + (j - 1) - pseudoline_weight[u[j]+1])
#         end
#     end
#     return mp
# end

# @inline function linedistance_add(a::Vector{NTuple{k,Int}}, b::Vector{NTuple{k,Int}}) where k
#     return unique!(vcat(a, b))
# end

# @inline function linedistance_mul(a::Vector{NTuple{k,Int}}, b::Vector{NTuple{k,Int}}) where k
#     out = Vector{NTuple{k,Int}}()
#     for x in a, y in b
#         push!(out, ntuple(i -> x[i] + y[i], K))
#     end
#     return unique!(out)
# end

# @inline function linedistance_filter(D, a::Vector{NTuple{k,Int}})
#     return [x for x in a if x <= D]
# end

# function linedistance_matmul(A::Matrix{Vector{NTuple{k,Int}}}, B::Matrix{Vector{NTuple{k,Int}}}, D::Vector{Int})::Matrix{Vector{NTuple{k,Int}}} where k
#     D = tuple(D...)
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

# function linedistance_matrix(nwt::NWT, D::Vector{Int})
#     tm = totalmatrix(nwt)
#     k = length(tm)
#     N = size(tm[1], 1)

#     A = [Vector{NTuple{k,Int}}() for _ in 1:N, _ in 1:N]
#     for i in 1:k
#         for (j, v) in enumerate(tm[i])
#             if v
#                 push!(A[j], tuple(I[1:k, i]...))
#             end
#         end
#     end

#     B = deepcopy(A)
#     powers = [Vector{NTuple{k,Int}}([tuple([0 for _ in 1:k]...); A[i, j]]) for i in 1:N, j in 1:N]
#     for _ in 1:sum(D)
#         B = linedistance_matmul(A, B, D)
#         powers = unique!.(vcat.(powers, B))
#     end

#     for v in powers
#         for (a, b) in distinct_pairs(v)
#             c = a .+ b
#             if all(c .<= D .&& (D .- c) .% 2 .== 0)
#                 @goto ld_next
#             end
#         end
#         return false
#         @label ld_next
#     end

#     return true
# end

# Generate all NWT with specified cc, n and no trees
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

# function face_graph(cc::CellComplex, D::Vector{Int})::Tuple{Vector{NTuple{cc.k + 1,Int}},Matrix{Int}}
#     @assert(length(D) == cc.k)
#     faceV = product1d(eachindex(cc.F), [0:d for d in D]...)
#     faceE = zeros(Int, length(faceV), length(faceV))

#     efs = [[] for _ in cc.E]
#     for (f, F) in enumerate(cc.F)
#         for (s, e) in F
#             push!(efs[e], f)
#         end
#     end
#     @assert(all(length.(efs) .== 2))

#     for (e, (s, d, i)) in enumerate(cc.E)
#         start_t = tuple([0 for d in D]...)
#         delta_t = tuple(I[1:cc.k, i]...)
#         end_t = tuple([d for d in D]...) .- delta_t

#         (f1, f2) = efs[e]

#         for t in product1d((start_t[i]:end_t[i] for i in 1:cc.k)...)
#             if1 = findfirst(==(tuple(f1, t...)), faceV)
#             if2 = findfirst(==(tuple(f2, (t .+ delta_t)...)), faceV)
#             faceE[if1, if2] = 1

#             if1 = findfirst(==(tuple(f2, t...)), faceV)
#             if2 = findfirst(==(tuple(f1, (t .+ delta_t)...)), faceV)
#             faceE[if1, if2] = 1
#         end
#     end

#     return faceV, faceE
# end

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

function is_bezout_order_zero(nwt::NWT, D::Vector{Int})
    tm = totalmatrix(nwt)
    N = size(tm[1], 1)
    dps = distinct_pairs(axes(tm, 1))

    lv, le = liftedgraph(tm, D)
    lg = SimpleDiGraph(le)
    fws = floyd_warshall_shortest_paths(lg)
    ep = enumerate_paths(fws)

    while !isempty(dps)
        (f1, f2) = first(dps)
        start_f = liftedindex(N, D, f1, [0 for d in D])
        ends_f = [liftedindex(N, D, f1, ds) for ds in product1d([collect(d:(-2):0) for d in D]...)]

        for middle_f in findall(i -> lv[i][1] == f2 && fws.dists[start_f, i] < typemax(Int), eachindex(lv))
            for end_f in ends_f
                if fws.dists[middle_f, end_f] < typemax(Int)
                    path1 = ep[start_f][middle_f]
                    path2 = ep[middle_f][end_f]
                    path = [path1[1:end-1]; path2]

                    setdiff!(dps, [(lv[f1][1], lv[f2][1]) for (f1, f2) in distinct_pairs(path)])

                    @goto boz_ok
                end
            end
        end

        return false

        @label boz_ok
    end

    return true
end

# function regiongraph(cc::CellComplex, T::Vector{@NamedTuple{g::Int, T::Vector{Bool}}}, D::Vector{Int})::Tuple{Vector{NTuple{cc.k + 1,Int}},Matrix{Int}}
#     @assert(length(D) == cc.k)
#     faceV = product1d(eachindex(cc.F), [0:d for d in D]...)
#     faceE = zeros(Int, length(faceV), length(faceV))

#     efs = [[] for _ in cc.E]
#     for (f, F) in enumerate(cc.F)
#         for (s, e) in F
#             push!(efs[e], f)
#         end
#     end
#     @assert(all(length.(efs) .== 2))

#     for (e, (s, d, i)) in enumerate(cc.E)
#         start_t = tuple([0 for d in D]...)
#         delta_t = tuple(I[1:cc.k, i]...)
#         end_t = tuple([d for d in D]...) .- delta_t

#         (f1, f2) = efs[e]

#         for t in product1d((start_t[i]:end_t[i] for i in 1:cc.k)...)
#             if1 = findfirst(==(tuple(f1, t...)), faceV)
#             if2 = findfirst(==(tuple(f2, (t .+ delta_t)...)), faceV)
#             faceE[if1, if2] = 1

#             if1 = findfirst(==(tuple(f2, t...)), faceV)
#             if2 = findfirst(==(tuple(f1, (t .+ delta_t)...)), faceV)
#             faceE[if1, if2] = 1
#         end
#     end

#     return faceV, faceE
# end

# function is_bezout_order_zero(cc::CellComplex, D::Vector{Int})
#     dps = distinct_pairs(eachindex(cc.F))

#     faceV, faceE = face_graph(cc, D)
#     faceG = SimpleDiGraph(faceE)
#     fws = floyd_warshall_shortest_paths(faceG)
#     ep = enumerate_paths(fws)

#     while !isempty(dps)
#         (f1, f2) = first(dps)
#         start_f = findfirst(==(tuple(f1, [0 for d in D]...)), faceV)
#         ends_f = [findfirst(==(tuple(f1, ds...)), faceV) for ds in product1d([collect(d:(-2):0) for d in D]...)]

#         for middle_f in findall(i -> faceV[i][1] == f2 && fws.dists[start_f, i] < typemax(Int), eachindex(faceV))
#             for end_f in ends_f
#                 if fws.dists[middle_f, end_f] < typemax(Int)
#                     path1 = ep[start_f][middle_f]
#                     path2 = ep[middle_f][end_f]
#                     path = [path1[1:end-1]; path2]

#                     setdiff!(dps, [(faceV[f1][1], faceV[f2][1]) for (f1, f2) in distinct_pairs(path)])

#                     @goto boz_ok
#                 end
#             end
#         end

#         return false

#         @label boz_ok
#     end

#     return true
# end

function emptyNWT(cc::CellComplex)::NWT
    return NWT(cc, [0 for _ in eachindex(cc.E)], [[] for _ in eachindex(cc.F)], [])
end

function is_bezout_order_n(cc::CellComplex, D::Vector{Int}, n::Int)
    if n == 0
        return is_bezout_order_zero(emptyNWT(cc), D)
    end

    if !is_bezout_order_n(cc, D, n - 1)
        return false
    end

    tm = totalmatrix(nwt)
    N = size(tm[1], 1)
    dps = distinct_pairs(axes(tm, 1))

    lv, le = liftedgraph(tm, D)
    lg = SimpleDiGraph(le)
    ds = [dijkstra_shortest_paths(lg, liftedindex(N, D, f, [0 for d in D]), allpaths=true) for f in eachindex(cc.F)]

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
                        for refined_nwt in add_line(cc, path)
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

# function is_bezout_order_one(cc::CellComplex, D::Vector{Int})
#     if !is_bezout_order_zero(cc, D)
#         return false
#     end
#     dps = distinct_pairs(eachindex(cc.F))

#     faceV, faceE = face_graph(cc, D)
#     faceG = SimpleDiGraph(faceE)
#     ds = [dijkstra_shortest_paths(faceG, findfirst(==(tuple(f, [0 for d in D]...)), faceV), allpaths=true) for f in eachindex(cc.F)]

#     while !isempty(dps)
#         (f1, f2) = first(dps)

#         start_f = findfirst(==((f1, [0 for d in D]...)), faceV)
#         for middle_f in findall(i -> faceV[i][1] == f2 && ds[f1].dists[i] < typemax(Int), eachindex(faceV))
#             for end_f in [findfirst(==(tuple(f1, ds...)), faceV) for ds in product1d([collect(d:(-2):m) for (d, m) in zip(D, faceV[middle_f][2:end])]...)]
#                 end2_t = faceV[end_f][2:end] .- faceV[middle_f][2:end]
#                 end2_f = findfirst(==(tuple(f1, end2_t...)), faceV)
#                 if ds[f2].dists[end2_f] < typemax(Int)
#                     start2_f = findfirst(==((f2, [0 for d in D]...)), faceV)
#                     paths1 = all_shortest_paths(start_f, middle_f, ds[f1].predecessors)
#                     paths2 = all_shortest_paths(start2_f, end2_f, ds[f2].predecessors)
#                     paths = [first.(faceV[[path1[1:end-1]; path2]]) for (path1, path2) in Iterators.product(paths1, paths2)]

#                     for path in paths
#                         for refined_nwt in add_line(cc, path)
#                             if is_bezout_order_zero(groundcc(refined_nwt), [D...; 1])
#                                 setdiff!(dps, [(faceV[f1][1], faceV[f2][1]) for (f1, f2) in distinct_pairs(path)])
#                                 @goto boo_ok
#                             end
#                         end
#                     end
#                 end
#             end
#         end

#         return false

#         @label boo_ok
#     end

#     return true
# end

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