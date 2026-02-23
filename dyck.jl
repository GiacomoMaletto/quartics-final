using Test

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

function dyck_pairs(D)::Vector{Tuple{Int,Int}}
    stack = []  # To store the positions of unmatched open parentheses
    pairs = []  # To store the pairs of matching parentheses

    for (i, bit) in enumerate(D)
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

### TESTS ###

@testset "dyck" begin
    dyck = Bool[1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0]
    @test dyck_sequence(dyck, 11) == [11, 12, 13, 14, 13, 12, 15, 12, 16, 17, 18, 17, 16, 12, 11, 19, 20, 19, 11, 21, 22, 23, 22, 24, 25, 26, 25, 27, 28, 27, 25, 24, 29, 24, 22, 21, 11]
    @test sequence_dyck(dyck_sequence(dyck, 11)) == dyck
    @test Set(dyck_pairs(dyck)) == Set([
        (1, 14), (2, 5), (3, 4), (6, 7), (8, 13), (9, 12), (10, 11), (15, 18), (16, 17),
        (19, 36), (20, 35), (21, 22), (23, 34), (24, 31), (25, 26), (27, 30), (28, 29), (32, 33)])
end

@testset "forest" begin
    @test dyck_to_forest(Bool[1, 0, 1, 1, 0, 0]) == [[], [[]]]
    @test dyck_to_forest(Bool[1, 1, 0, 0, 1, 0]) == [[], [[]]]
    @test forest_to_dyck(Any[Any[], Any[Any[]]]) == Bool[1, 0, 1, 1, 0, 0]
    @test forest_to_dyck([]) == Bool[]
end