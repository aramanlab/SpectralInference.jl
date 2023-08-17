module NewickTreeExt

using SpectralInference, NewickTree
using StatsBase: sample, mean

using NewickTree: Node
using NewickTree: prewalk, height, getpath
using NewickTree: setdistance!, setsupport!, support, getpath

"""
    getleafnames(t::NewickTree.Node)

Get names of all the leafs with `t` as ancester
"""
SpectralInference.getleafnames(t::Node) = name.(getleaves(t))

"""
    getleafids(t::NewickTree.Node)

Get IDs of all the leafs with `t` as ancester
"""
SpectralInference.getleafids(t::Node) = id.(getleaves(t))


# how to delete a newick tree node
function Base.delete!(n::Node)
    p = parent(n)
    cs = children(n)
    for c in cs
        push!(p, c)
    end
    delete!(p, n)
end



"""    
    network_distance(n::Node, m::Node)

shortest traversal distance between two nodes

sibling nodes would have network distance = 2
"""
function SpectralInference.network_distance(n::Node, m::Node)
    p1, p2 = getpath(n, m)
    (length(p1) - 1) + (length(p2) - 1)
end

"""    
    network_distances(t::Node)

shortest traversal distance between all leafs with `t` as ancester
sibling nodes would have network distance = 2
"""
function SpectralInference.network_distances(tree::Node)
    leaves = getleaves(tree)
    dists = zeros(length(leaves), length(leaves))
    for j in axes(dists, 2), i in (j+1):lastindex(dists, 1)
        dists[i, j] = network_distance(leaves[i], leaves[j])
        dists[j, i] = dists[i, j]
    end
    return dists
end


"""
    patristic_distance(n::Node, m::Node)

shortest branch length path between two nodes.

sibling nodes `i`, and `j` of parent `p` would have patristic `distance(i, p) + distance(j, p)`
"""
function SpectralInference.patristic_distance(n::Node, m::Node)
    NewickTree.getdistance(n, m)
end

"""    
    patristic_distances(t::Node)

shortest branch length path between all leafs with `t` as ancester

sibling nodes `i`, and `j` of parent `p` would have patristic `distance(i, p) + distance(j, p)`
"""
function SpectralInference.patristic_distances(tree::Node)
    leaves = getleaves(tree)
    dists = zeros(length(leaves), length(leaves))
    for j in axes(dists, 2), i in (j+1):lastindex(dists, 1)
        dists[i, j] = patristic_distance(leaves[i], leaves[j])
        dists[j, i] = dists[i, j]
    end
    return dists
end

""" 
    fscore_precision_recall(reftree, predictedtree) 
fscore, precision, and recall of branches between the two trees
"""
function SpectralInference.fscore_precision_recall(reftree::Node, predtree::Node; β=1.0)
    truesplits = keys(_tally_tree_bifurcations(reftree))
    predsplits = keys(_tally_tree_bifurcations(predtree))
    compsplits = [k for k in truesplits] .== permutedims([k for k in predsplits])
    precision = mean(sum(compsplits; dims=1))
    recall = mean(sum(compsplits; dims=2))
    fscore = (1 + β) * (precision * recall) / (β * precision + recall)
    return fscore, precision, recall
end


"""
    tally_tree_bifurcations(rootnode::AbstractTree, [cntr::Accumulator])

returns tally of leaf groups based on splitting tree at each internal node

assumes nodenames are strings parsible to integers, and returns Accumulator with key
as a string: `1110000` where 1 = belonging to smaller group & 0 = belonging to the larger.
values are the number of times this split is observed. cntr is modified in place so using
the same counter in multiple calls will keep a tally across multiple trees
"""
function _tally_tree_bifurcations(tree::Node, cntr=Dict{BitVector,Int}(); countroot=true)
    leafnames = sort(getleafnames(tree))
    nleaves = length(leafnames)
    for node in prewalk(tree)
        # only internal nodes
        isleaf(node) && continue
        !countroot && isroot(node) && continue
        # get leaves
        key = zeros(Bool, nleaves)
        subsetleaves = indexin(getleafnames(node), leafnames)
        key[subsetleaves] .= true
        # convert to bitstring true := smaller cluster
        keystr = sum(key) <= nleaves / 2 ? key : .!key
        cntr[keystr] = 1 + get(cntr, keystr, 0)
    end
    cntr
end



"""
    cuttree(distfun::Function, tree::NewickTree.Node, θ) 

returns all vector of Nodes where distance to root is
greater than theta `d > θ`.

distance function must take to NewickTree.Node objects and 
compute a scaler distance between them.
"""
function SpectralInference.cuttree(distfun::Function, tree::Node, θ)
    θ ≤ distfun(tree, tree) && return [tree]
    ns = typeof(tree)[]
    # define tree traversal
    function walk!(n, t)
        distn = distfun(n, t)
        distp = distfun(parent(n), t)
        # if edge crosses θ then add current node
        if distp ≤ θ < distn
            push!(ns, n)
            return
            # if node is leaf an within θ add as singlet cluster
        elseif isleaf(n) && distn ≤ θ
            push!(ns, n)
            return
            # otherwise continue taversing tree
        else
            !isleaf(n) && for c in n.children
                walk!(c, t)
            end
            return
        end
    end
    # actually taverse the tree
    for c in tree.children
        walk!(c, tree)
    end
    return ns
end

"""
    mapnodes(fun::Function, tree::NewickTree.Node, args...; filterfun::Function=x->true, kwargs...)

maps function `fun()` across internal nodes of tree.

args and kwargs are passed to `fun()`
"""
function SpectralInference.mapnodes(fun::Function, tree::Node, args...; filterfun::Function=x -> true, kwargs...)
    results = Vector()
    for (i, node) in enumerate(prewalk(tree))
        # only internal nodes
        filterfun(node) || continue
        # run function
        push!(results, fun(node, args...; kwargs...))
    end
    return results
end

"""
    mapinternalnodes(fun::Function, tree::NewickTree.Node, args...; kwargs...)

maps function `fun()` across internal nodes of tree.

args and kwargs are passed to `fun()`
"""
function SpectralInference.mapinternalnodes(fun::Function, tree::Node, args...; kwargs...)
    results = Vector()
    for (i, node) in enumerate(prewalk(tree))
        # only internal nodes
        isleaf(node) && continue
        # run function
        push!(results, fun(node, args...; kwargs...))
    end
    return results
end

"""
    maplocalnodes(fun::Function, tree::NewickTree.Node, args...; kwargs...)

maps function `fun()` across internal nodes of tree conditioned on having
    one direct child that is a leaf.

args and kwargs are passed to `fun()`
"""
function SpectralInference.maplocalnodes(fun::Function, tree::Node, args...; kwargs...)
    results = Vector()
    for (i, node) in enumerate(prewalk(tree))
        # only internal nodes
        isleaf(node) && continue
        # only nodes that have a leaf as child
        all(.!isleaf.(children(node))) && continue
        # run function
        push!(results, fun(node, args...; kwargs...))
    end
    return results
end


"""
    collectiveLCA(nodes)

finds last common ancester of a collection of Nodes
"""
function SpectralInference.collectiveLCA(nodes::AbstractArray{<:NewickTree.Node})
    lca = map(b -> NewickTree.getlca(nodes[1], b), nodes[2:end])
    idx = argmin(NewickTree.height.(lca))
    lca[idx]
end

"""
    as_polytomy(fun::Function, tree::NewickTree.Node)
    as_polytomy!(fun::Function, tree::NewickTree.Node)

removes internal nodes from tree based on `fun()` 
which must return true if node is to be removed

by default removes zero length branches 
(i.e. nodes where distance between child and parent == 0)
"""
function SpectralInference.as_polytomy(fun::Function, tree::Node)
    tree_new = deepcopy(tree)
    as_polytomy!(fun, tree_new)
    tree_new
end

function SpectralInference.as_polytomy!(fun::Function, tree::Node)
    for n in filter(fun, prewalk(tree))
        !isroot(n) && !isleaf(n) && delete!(n)
    end
end

"""
    ladderize!(tree, rev=false)

sorts children of each node by number of leaves descending from the child in ascending order.
"""
function SpectralInference.ladderize!(t; rev=false)
    function walk!(n)
        if isleaf(n)
            return 1
        else
            numleaves = [walk!(c) for c in children(n)]
            n.children .= n.children[sortperm(numleaves, rev=rev)]
            return sum(numleaves)
        end
    end
    walk!(t)
end


"""
    spectral_lineage_encoding(tree::Node, orderedleafnames=getleafnames(tree); filterfun=x->true)

returns vector of named tuples with the id `nodeid` of the node and `sle` a vector of booleans
ordered by `orderedleafnames` where true indicates the leaf descends from the node and false indicates that it does not.
"""
function SpectralInference.spectral_lineage_encoding(tree::Node, orderedleafids=getleafids(tree); filterfun=x -> true)
    map(filter(filterfun, prewalk(tree))) do node
        tmp = falses(length(orderedleafids))
        tmp[indexin(getleafids(node), orderedleafids)] .= true
        nodeid = id(node)
        (; nodeid, sle=tmp)
    end
end

"""
    pairedMI_across_treedepth(metacolumns, metacolumns_ids, tree)
    pairedMI_across_treedepth(metacolumns, metacolumns_ids, compare::Function=(==), tree::Node; ncuts=100, bootstrap=false, mask=nothing)

iterates over each metacolumn and calculates MI between the paired elements of the metacolumn and the paired elements of tree clusters.

Args:
* metacolumns: column iterator, can be fed to `map(metacolumns)`
* metacolumn_ids: ids for each element in the metacolumn. should match the leafnames of the tree, but not necessarily the order.
* tree: NewickTree tree
* compare: function used to calculate similarity of two elements in metacolumn. 
Should be written for each element as it will be broadcast across all pairs.

Returns:
* (; MI, treedepths)
* MI: Vector{Vector{Float64}} MI for each metacolumn and each tree depth
* treedepths Vector{<:Number} tree depth (away from root) for each cut of the tree
"""
function SpectralInference.pairedMI_across_treedepth(metacolumns, metacolumn_ids, tree::Node;
    comparefun::Function=(==), treecut_distancefun::Function=network_distance, ncuts=100, bootstrap=false, mask=nothing)

    clusts, treedepths = clusters_per_cutlevel(treecut_distancefun, tree, ncuts)
    clust_names = getleafnames(tree)
    MI = map(metacolumns) do mcol
        pmcol = comparefun.(mcol, permutedims(mcol))
        _collectMI_across_treedepth(clusts, clust_names, metacolumn_ids, pmcol; bootstrap, mask)
    end
    (; MI, treedepths)
end



"""
    clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)

Returns:
* clusts: vector of cluster-memberships. each value indicates the cluster membership of the leaf at that cut. 
    leaves are ordered in prewalk order within each membership vector
* treedepths: distance from root for each of the `ncuts`
"""
function SpectralInference.clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)
    minmax = extrema(maplocalnodes(distfun, tree, tree))
    treedepths = range(zero(first(minmax)), minmax[2], length=ncuts)
    clusts_nodes = [cuttree(distfun, tree, cut) for cut in treedepths]
    clustmappings = map(c -> getleafnames.(c), clusts_nodes)
    clusts = [Int.(vcat([zeros(length(c)) .+ j for (j, c) in enumerate(clustmapping)]...)) for clustmapping in clustmappings]
    (; clusts, treedepths)
end

function _collectMI_across_treedepth(clusts, clust_names, IDS, ptax; bootstrap=false, mask=nothing)
    uppertriangle = triu(trues(size(ptax)), 1)
    uppertriangle = isnothing(mask) ? uppertriangle : uppertriangle[mask, mask]
    clustorder = indexin(IDS, clust_names)
    # ptax = if isnothing(mask) ptax[uppertriangle] else ptax[mask, mask][uppertriangle] end
    map(clusts) do cids
        ptax, mask, clustorder, uppertriangle
        pcids = cids[clustorder] .== permutedims(cids[clustorder])
        pcids = if isnothing(mask)
            pcids
        else
            pcids[mask, mask]
        end
        wptax = if isnothing(mask)
            ptax
        else
            ptax[mask, mask]
        end
        pcids = pcids[uppertriangle]
        wptax = wptax[uppertriangle]
        if bootstrap
            vals_idx = sample(axes(pcids, 1), length(pcids), replace=true)
            pcids = pcids[vals_idx]
            wptax = wptax[vals_idx]
        end
        empiricalMI(wptax, pcids)
    end
end

using PrecompileTools
@setup_workload begin
    @compile_workload begin
        N = 10
        for T in [Float64, Float32, Int64, Int32]
            m = rand(T, N, N)
            usv = svd(m)
            dij = spectraldistances(usv.U, usv.S, getintervals(usv.S))
            hc = UPGMA_tree(dij)
            tree = readnw(newickstring(hc, string.(1:N)))
            leafnames = getleafnames(tree)
            network_distances(tree)
            patristic_distances(tree)
            fscore_precision_recall(tree, tree) == (1.0, 1.0, 1.0)
            cuttree(network_distance, tree, 0.0)
            collectiveLCA(getleaves(tree))
            as_polytomy(n -> NewickTree.support(n) < 0.5, tree)
            mi, treedepths = pairedMI_across_treedepth((; a=rand(Bool, N), b=rand([1, 0], N)), leafnames, tree; ncuts=5)
            mi, treedepths = pairedMI_across_treedepth((; a=rand(UInt16, N), b=rand(Float64, N)), leafnames, tree; ncuts=5, comparefun=(x, y) -> abs(x - y))
            mi, treedepths = pairedMI_across_treedepth((; a=rand(UInt16, N), b=rand(Float64, N)), leafnames, tree; treecut_distancefun=patristic_distance, ncuts=5, comparefun=(x, y) -> abs(x - y))
            leafids = getleafids(tree)
            spectral_lineage_encoding(tree)
            spectral_lineage_encoding(tree, leafids)
            spectral_lineage_encoding(tree; filterfun=!isleaf)
            spectral_lineage_encoding(tree, leafids; filterfun=!isleaf)
        end
    end
end

end # module