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