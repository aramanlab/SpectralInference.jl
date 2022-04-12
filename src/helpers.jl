
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