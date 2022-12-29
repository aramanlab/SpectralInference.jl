"""
    scaledcumsum(c; dims=1)

cumsum divided by maximum cumulative value
"""
function scaledcumsum(c; dims=1)
    cc=cumsum(c, dims=dims)
    cc./maximum(cc, dims=dims) 
end


"""
    explainedvariance(s::AbstractVector{<:Number})
"""
function explainedvariance(s::AbstractVector{<:Number})
    s.^2 / sum(s.^2)
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
are calculated from `enumerate(combinations(axes(m,2), 2))`
"""
function pairwise(func::Function, m::M) where M<:AbstractMatrix
    result = zeros(binomial(size(m,2),2))
    for (k, (i,j)) in enumerate(combinations(axes(m,2), 2))
        result[k] = func(m[:,i], m[:,j])
    end
    return result
end


"""
        numpairs2N(x)::Integer
    solve choose(n,k)=x for n
    for numbers around a trillion use a BigInt for x
"""
function numpairs2N(x)
    # solve choose(n,k)=x for n
    n = (1 + sqrt(1 + 8x)) / 2
    # n not whole means that x was not a binomial number
end

"""
    reshape_pairs_to_distance_matrix(pairs::Vector; defaultval=zeros)
take the columnwise vectorized lower diagonal of distance matrix and remake a symetric distance matrix.
"""
function reshape_pairs_to_distance_matrix(pairsvec; defaultval=zeros)
    n = numpairs2N(length(pairsvec))
    n % 1 == 0. || throw(ArgumentError("Your matrix has the wrong number of pairs to be a pairwise combination"))
    distmtx = defaultval(Int(n), Int(n))
    for (k, (i,j)) in enumerate(combinations(1:Int(n), 2))
        distmtx[i,j] = pairsvec[k]
        distmtx[j,i] = pairsvec[k]
    end
    return distmtx
end

"""
    k2ij(k, n)
which pair `(i,j)` produces the `k`th element of combinations(vec), where vec is length `n`
"""
function k2ij(k,n) # column major ordering of lower triangle
    rvLinear = (n*(n-1))÷2-k;
    i = floor(Int, (sqrt(1+8*rvLinear)-1)÷2);
    j = rvLinear - i*(i+1)÷2;
    (n-j, n-(i+1))
end

"""
    ij2k(i,j,n)
with pair `(i,j)` give index `k` to the pairs produced by combinations(vec), where vec is length `n`
"""
function ij2k(i, j, n) # for symetric matrix
    i, j = i < j ? (i, j) : (j, i)
    k = ((n*(n-1))÷2) - ((n-i)*((n-i)-1))÷2 + j - n
end