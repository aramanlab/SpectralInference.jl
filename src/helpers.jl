
"""
    zscore by columns by default. set dims=2 for rows
"""
function zscore(X; dims=1)
    ztransform = StatsBase.fit(StatsBase.ZScoreTransform, X; dims=dims); # column zscore
    return StatsBase.transform(ztransform, X)
end

"""
    scaledcumsum(c; dims=1)

cumsum divided by maximum cumulative value
"""
function scaledcumsum(c; dims=1)
    cc=cumsum(c, dims=dims)
    cc./maximum(cc, dims=dims) 
end


"""
    minspaceneeded(n, p; bits=64) = Base.format_bytes(binomial(n,2) * p * bits)

how much memory is needed to store spectral residual trace
"""
minspaceneeded(n,p; bits=64) = Base.format_bytes(binomial(n,2) * p * bits)

"""
    spimtx_spaceneeded(n, p; bits=64) = Base.format_bytes(binomial(n,2) * p * bits)

how much memory is needed to store SPI distance matrix
"""
spimtx_spaceneeded(n; bits=64) = Base.format_bytes(n^2 * bits)


"""
    pairwise(func::Function, m::M) where M<:AbstractMatrix

returns upper offdiagonals of `res[k] = func(i, j)` where `(k, (i,j))` 
are calculated from `enumerate(combinations(1:size(m,2), 2))`
"""
function pairwise(func::Function, m::M) where M<:AbstractMatrix
    result = zeros(binomial(size(m,2),2))
    for (k, (i,j)) in enumerate(combinations(1:size(m,2), 2))
        result[k] = func(m[:,i], m[:,j])
    end
    return result
end


""" 
    nwstr(hc::Hclust[, tiplabels::Vector[<:String]])

convert Hclust to newick tree string
Args:
* hc, `Hclust` object from Clustering package
* tiplabels, `Vector{<:String}` names in same order as distance matrix
"""
nwstr(hc::Hclust; labelinternalnodes=true) = nwstr(hc, string.(1:length(hc.order)); labelinternalnodes)
nwstr(hc::Hclust, tiplabels::Vector{<:Symbol}; labelinternalnodes=true) = nwstr(hc, string.(tiplabels); labelinternalnodes)
function nwstr(hc::Hclust, tiplabels::Vector{<:String}; labelinternalnodes=true)
    r = length(hc.heights)
    _nwstr(view(hc.merges, :, :), view(hc.heights, :), r, r, view(tiplabels, :); labelinternalnodes) * ";"
end
function _nwstr(merges::A, heights::B, i::C, p::C, tiplabels::D; labelinternalnodes=true)::String where {
        A<:AbstractArray{<:Integer}, B<:AbstractVector{<:AbstractFloat},
        C<:Integer, D<:AbstractVector{<:AbstractString}
    }
    j::Int64 = merges[i,1] # left subtree pointer
    k::Int64 = merges[i,2] # right subtree pointer
    a::String = if j < 0 # if tip format tip
            tiplabels[abs(j)] * ':' * @sprintf("%e", heights[i])
        else # recurse and format internal node
            _nwstr(view(merges, :, :), view(heights, :), j, i, view(tiplabels, :); labelinternalnodes)
        end
    b::String = if k < 0 # if tip format tip
            tiplabels[abs(k)] * ':' * @sprintf("%e", heights[i])
        else # recurse and format internal node
            _nwstr(view(merges, :, :), view(heights, :), k, i, view(tiplabels, :); labelinternalnodes)
        end
    nid = labelinternalnodes ? "node" * string(length(heights) + i + 1) : ""
    dist = @sprintf("%e", heights[p] - heights[i])
    _newick_merge_strings(a,b,nid,dist)
end
function _newick_merge_strings(a::S, b::S, n::S, d::S) where S<:String
    '(' *  a *  ',' * b * ')' * n * ':' * d
end