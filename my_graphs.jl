using Graphs

function my_union(g::SimpleGraph, h::SimpleGraph)
    return SimpleGraph(blockdiag(adjacency_matrix(g), adjacency_matrix(h)))
end

# input: a graph g and a list of vertices vs
# output: the graph gprime obtained by merging every vertex in vs to minimum(vs) and otherwise keeping the same order,
# and a list vmap corresponding to the morphism from g to gprime
function my_merge_vertices(g::SimpleGraph, vs::Vector{Int})
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
function my_merge_vertices_list(g::SimpleGraph, vss::Vector{Vector{Int}})
    vmap_final = collect(1:nv(g))
    for i in eachindex(vss)
        g, vmap = my_merge_vertices(g, vss[i])
        vss = map(vs -> unique([vmap[v] for v in vs]), vss)
        vmap_final = [vmap[x] for x in vmap_final]
    end

    return g, vmap_final
end