# function definitions are in extensions for various backends of tree implementations
""" returns names of all leaves desended from a node in prewalk order """
function getleafnames end
""" 
    network_distance(leaf_i, leaf_j) 
returns network distance between two leaves 
"""
function network_distance end
""" 
    network_distances(tree::Node) 
returns network distances between all pairs of leaves (leaves are in same order as `getleafnames`)
"""
function network_distances end
""" 
    patristic_distance(leaf_i, leaf_j) 
returns patristic distance between two leaves
"""
function patristic_distance end
"""
    patristic_distances(tree::Node)
returns patristic distances between all pairs of leaves (leaves are in same order as `getleafnames`)
"""
function patristic_distances end
""" 
    fscore_precision_recall(reftree, predictedtree) 
fscore, precision, and recall of branches between the two trees
"""
function fscore_precision_recall end

""" 
    cuttree(tree, θ)

return vector of nodes whose distance from the root are < θ and whose children's distance to the root are > θ
"""
function cuttree end
"""
    mapinternalnodes(f::Function, tree)
maps across all nodes that have children in prewalk order and applies function `f(node)`
"""
function mapinternalnodes end
"""
    maplocalnodes(f::Function, tree)
maps across all nodes that have a leaf as a child in prewalk order and applies function `f(node)`
"""
function maplocalnodes end
"""
    collectiveLCA(treenodes::AbstractVector)
returns last common ancestor of the vector of treenodes
"""
function collectiveLCA end
"""
    as_polytomy!(f::Function, tree)
in-place removal of nodes. function `f` should return `true` if the node is to be removed.
Children of node are attached to the parent of the removed node
"""
function as_polytomy! end
"""
    as_polytomy(f::Function, tree)
Makes a copy of the tree, then removes nodes. function `f` should return `true` if the node is to be removed.
Children of node are attached to the parent of the removed node
"""
function as_polytomy end
"""
    pairedMI_across_treedepth(metacolumns, metacolumns_ids, tree)
    pairedMI_across_treedepth(metacolumns, metacolumns_ids, compare::Function=(==), tree::Node; ncuts=100, bootstrap=false, mask=nothing)

iterates over each metacolumn and calculates MI between the paired elements of the metacolumn and the paired elements of tree clusters.

Args:
* metacolumns: column iterator, can be fed to `map(metacolumns)`
* metacolumn_ids: ids for each element in the metacolumn. should match the leafnames of the tree, but not necessarily the order.
* tree: NewickTree tree
* compare: function used to calculate similarity of two elements in metacolumn. Should be written for each element as it will be broadcast across all pairs.

Returns:
* (; MI, treedepths)
* MI: Vector{Vector{Float64}} MI for each metacolumn and each tree depth
* treedepths Vector{<:Number} tree depth (away from root) for each cut of the tree
"""
function pairedMI_across_treedepth end
"""
    clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)

Returns:
* clusts: vector of cluster-memberships. each value indicates the cluster membership of the leaf at that cut. 
    leaves are ordered in prewalk order within each membership vector
* treedepths: distance from root for each of the `ncuts`
"""
function clusters_per_cutlevel end

"""
    ladderize!(tree, rev=false)
sorts children of each node by number of leaves descending from the child in ascending order.
"""
function ladderize! end
"""
    spectral_lineage_encoding(tree::Node, orderedleafnames=getleafnames(tree); filterfun=x->true)

returns vector of named tuples with the id `nodeid` of the node and `sle` a vector of booleans
ordered by `orderedleafnames` where true indicates the leaf descends from the node and false indicates that it does not.
"""
function spectral_lineage_encoding end