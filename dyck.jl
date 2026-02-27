using Graphs
using Test

# ChatGPT-generated
function dyck_to_sequence(W::Vector{Bool})::Vector{Int}
    out = Int[]
    label = 1
    push!(out, label)            # start at root labelled 1
    stack = [label]              # current path (top = current node)

    for b in W
        if b == 1
            label += 1
            push!(stack, label)   # create and move to new child
            push!(out, label)
        elseif b == 0
            pop!(stack)                # move up (must be valid Dyck word)
            push!(out, stack[end])     # record the parent's label
        else
            throw(ArgumentError("W must contain only 0 or 1"))
        end
    end
    return out
end

function sequence_to_dyck(seq::Vector{Int})::Vector{Bool}
    seq = [findfirst(==(n), unique(seq)) for n in seq]
    m = length(seq)
    m >= 1 || throw(ArgumentError("borseq must be nonempty"))
    # borseq[1] == 1 || throw(ArgumentError("first label must be 1"))

    # quick type-correct output
    # dyck = Vector{Bool}(undef, max(0, m - 1))
    W = Bool[]

    # Step 1: infer bits by comparing consecutive borseq
    for i in 1:m-1
        a, b = seq[i], seq[i+1]
        if b == a
            throw(ArgumentError("invalid transition: two consecutive equal regions")) # continue
        end
        push!(W, b > a)
    end

    # @assert(dyck_to_sequence(dyck, seq[1]) == seq)

    return W
end

# the root is always the first vertex
# ChatGPT-generated, some alterations by me
function dyck_to_rtgraph(w::Vector{Bool})::Tuple{SimpleGraph,Int}
    g = SimpleGraph(1)
    stack = Int[1]

    for b in w
        if b
            v = nv(g) + 1
            add_vertex!(g)
            (!isempty(stack)) || error("Invalid Dyck word: too many closes")
            add_edge!(g, stack[end], v)
            push!(stack, v)
        else
            pop!(stack)
        end
    end

    (stack == [1]) || error("Invalid Dyck word: unclosed opens")
    return g, 1
end

# ChatGPT-generated, some alterations by me
function rtgraph_to_dyck((g, r)::Tuple{SimpleGraph,Int})::Vector{Bool}
    D = Bool[]

    for w in neighbors(g, r)
        _dfs_dyck!(D, g, w, r)
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


# ChatGPT-generated
function dyck_to_mabel(W::Vector{Bool})::Vector{Int}
    labels = Int[]
    stack = Int[]
    next_label = 1

    for x in W
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
function mabel_to_dyck(mabel::Vector{Int})::Vector{Bool}
    seen = Set{Int}()
    W = Vector{Int}(undef, length(mabel))

    for (i, x) in pairs(mabel)
        if x ∈ seen
            W[i] = 0   # close
        else
            push!(seen, x)
            W[i] = 1   # open
        end
    end

    return W
end

function dyck_to_pairs(W)::Vector{Tuple{Int,Int}}
    stack = []  # To store the positions of unmatched open parentheses
    pairs = []  # To store the pairs of matching parentheses

    for (i, bit) in enumerate(W)
        if bit == 1
            # Push the position of the open parenthesis onto the stack
            push!(stack, i)
        elseif bit == 0
            # Pop the last unmatched open parenthesis and form a pair
            if !isempty(stack)
                open_pos = pop!(stack)
                push!(pairs, (open_pos, i))
            else
                throw(ArgumentError("Invalid Dyck word: Unmatched closing parenthesis."))
            end
        end
    end
    if !isempty(stack)
        throw(ArgumentError("Invalid Dyck word: Unmatched opening parenthesis."))
    end

    return pairs
end

function all_dyck_words(n::Integer)::Vector{Vector{Bool}}
    L = 2n
    out = Vector{Vector{Bool}}()
    buf = Vector{Bool}(undef, L)            # reuse and copy when a word is complete

    function back!(pos::Int, opens::Int, closes::Int)
        if pos > L
            push!(out, copy(buf))
            return
        end
        # you may place an open if we still need more opens
        if opens < n
            buf[pos] = 1
            back!(pos + 1, opens + 1, closes)
        end
        # you may place a close only if it won't break the prefix condition
        if closes < opens
            buf[pos] = 0
            back!(pos + 1, opens, closes + 1)
        end
    end

    back!(1, 0, 0)
    return out
end

function rtgraph_to_rt((g, r)::Tuple{SimpleGraph,Int}, p::Int=0)
    ngbs = setdiff(neighbors(g, r), p)
    return (r, Any[rtgraph_to_rt((g, n), r) for n in ngbs])
end

function unlabel_rt((v, rt))
    return Any[unlabel_rt(br) for br in rt]
end

function sort_rt((v, rt))
    return (v, sort(Any[sort_rt(b) for b in rt], by=unlabel_rt))
end

function rt_branch_to_dyck((v, rt))::Tuple{Vector{Bool},Vector{Int}}
    w = Bool[1]
    vmapinv = Int[v]
    for (w_br, vmapinv_br) in rt_branch_to_dyck.(rt)
        append!(w, w_br)
        append!(vmapinv, vmapinv_br)
    end
    return [w; 0], vmapinv
end

function rt_to_dyck((v, rt))::Tuple{Vector{Bool},Dict{Int,Int}}
    w = Bool[]
    vmapinv = Int[v]
    for (w_br, vmapinv_br) in rt_branch_to_dyck.(rt)
        append!(w, w_br)
        append!(vmapinv, vmapinv_br)
    end
    # vmap = similar(vmapinv)
    # println("vmapinv: ", vmapinv)
    # for (k, v) in enumerate(vmapinv)
    #     vmap[v] = k
    # end
    vmap = Dict(value => key for (key, value) in enumerate(vmapinv))
    return w, vmap
end

function rtree_nonisomorphic_embeddings(bigw::Vector{Bool}, smallw::Vector{Bool})
    bigg, _ = dyck_to_rtgraph(bigw)
    smallg, _ = dyck_to_rtgraph(smallw)

    subg_verts = [first.(ps) for ps in Graphs.Experimental.all_subgraphisomorph(
        bigg, smallg, vertex_relation=((u, v) -> (u == v == 1) || (u != 1 && v != 1)))]

    subg_colors = [fill(2, nv(bigg)) for _ in subg_verts]
    for (i, vs) in enumerate(subg_verts)
        for v in vs
            subg_colors[i][v] = 3
        end
        subg_colors[i][1] = 1
    end
    subgs = collect(eachindex(subg_verts))
    i = 1
    while i <= length(subgs)
        for j in length(subgs):(-1):(i+1)
            if Graphs.Experimental.has_isomorph(bigg, bigg,
                vertex_relation=((u, v) -> subg_colors[subgs[i]][u] == subg_colors[subgs[j]][v]))
                deleteat!(subgs, j)
            end
        end
        i += 1
    end
    nosymms = subg_verts[subgs]
    return nosymms
    # println("nosymms: ", nosymms)
    # diffs = SimpleGraph[]
    # for ns in nosymms
    #     push!(diffs, difference(bigg, my_induced_subgraph(bigg, ns)))
    # end
    # return diffs
end

function rtgraph_depth((g, r)::Tuple{SimpleGraph,Int}, p::Int=0)
    ngbs = setdiff(neighbors(g, r), p)
    if isempty(ngbs)
        return 0
    else
        return 1 + maximum([rtgraph_depth((g, n), r) for n in ngbs])
    end
end

function rtgraph_diameter((g, r)::Tuple{SimpleGraph,Int}, p::Int=0)
    ngbs = setdiff(neighbors(g, r), p)
    if isempty(ngbs)
        return 0
    elseif length(ngbs) == 1
        rtgraph_depth((g, r))
    else
        depths = [rtgraph_depth((g, n), r) for n in ngbs]
        return 2 + maximum([depths[i] + depths[j] for i in eachindex(ngbs) for j in (i+1):length(ngbs)])
    end
end

# function rtgraph_split((g, r)::Tuple{SimpleGraph,Int})::Vector{Tuple{SimpleGraph,Int}}

# end

# claude
# vectors of size N of nonnegative integers that add up to n
function compositions(n::Int, N::Int)::Vector{Vector{Int}}
    results = []
    current = zeros(Int, N)

    function backtrack(remaining, i)
        if i == N
            current[i] = remaining
            push!(results, copy(current))
            return
        end
        for k in 0:remaining
            current[i] = k
            backtrack(remaining - k, i + 1)
        end
    end

    backtrack(n, 1)
    return results
end

### TESTS ###

@testset "dyck" begin
    dyck = Bool[1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0]
    @test dyck_to_sequence(dyck) .+ 10 == [11, 12, 13, 14, 13, 12, 15, 12, 16, 17, 18, 17, 16, 12, 11, 19, 20, 19, 11, 21, 22, 23, 22, 24, 25, 26, 25, 27, 28, 27, 25, 24, 29, 24, 22, 21, 11]
    @test sequence_to_dyck(dyck_to_sequence(dyck)) == dyck
    @test Set(dyck_to_pairs(dyck)) == Set([
        (1, 14), (2, 5), (3, 4), (6, 7), (8, 13), (9, 12), (10, 11), (15, 18), (16, 17),
        (19, 36), (20, 35), (21, 22), (23, 34), (24, 31), (25, 26), (27, 30), (28, 29), (32, 33)])
end

@testset "forest" begin
    @test dyck_to_forest(Bool[1, 0, 1, 1, 0, 0]) == [[], [[]]]
    @test dyck_to_forest(Bool[1, 1, 0, 0, 1, 0]) == [[], [[]]]
    @test forest_to_dyck(Any[Any[], Any[Any[]]]) == Bool[1, 0, 1, 1, 0, 0]
    @test forest_to_dyck([]) == Bool[]
end