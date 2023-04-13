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
are calculated from `enumerate(((i, j) for j in axes(m, 2) for i in (j+1):lastindex(m, 2)))`
"""
function pairwise(func::Function, m::M) where M<:AbstractMatrix
    result = zeros(binomial(size(m,2),2))
    for (k, (i,j)) in enumerate(((i, j) for j in axes(m, 2) for i in (j+1):lastindex(m, 2)))
        result[k] = func(m[:,i], m[:,j])
    end
    return result
end


numpairs2N(x) = (1 + sqrt(1 + 8x)) / 2
checkoffdiagonal(d) = numpairs2N(length(d)) % 1 == 0

"""
    squareform(d::AbstractVector, fillvalue=zero(eltype(d)))
    squareform(d::AbstractVector, fillvalue=zero(eltype(d)))

If `d` is a vector, `squareform` checks if it of `n` choose 2 length for integer `n`, then fills the values of a symetric square matrix with the values of `d`.

If `d` is a matrix, `squareform` checks if it is square then fills the values of vector with the lower offdiagonal of matrix `d` in column order form. 

`fillvalue` is the initial value of the produced vector or matrix. Only really apparant in a produced matrix where it will be the values on the diagonal.
"""
function squareform(d::AbstractVector, fillvalue=zero(eltype(d))) 
    checkoffdiagonal(d) || throw(ArgumentError("Vector wrong length, to be square matrix offdiagonals"))
    n = numpairs2N(length(d))
    Dij = fill(fillvalue, Int(n), Int(n))
    for (k,(i,j)) in enumerate(((i, j) for j in axes(Dij, 2) for i in (j+1):lastindex(Dij, 1)))
        Dij[i, j] = d[k]
        Dij[j, i] = d[k]
    end
    Dij
end
function squareform(d::AbstractMatrix, fillvalue=zero(eltype(d)))
    println(size(d))
    (size(d, 1) == size(d, 2)) || throw(ArgumentError("size of d: $(size(d)), is not square"))
    n = binomial(size(d,1), 2)
    Dk = fill(fillvalue, n)
    for (k,(i,j)) in enumerate(((i, j) for j in axes(d, 2) for i in (j+1):lastindex(d, 1)))
        Dk[k] = d[i, j]
    end
    Dk
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