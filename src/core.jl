
struct SPIResult{T}
    usv::SVD{T}
    spimtx::Matrix{T}
    spitree::Hclust
end

struct LSVs end
struct RSVs end

"""
    spiresult(A::Matrix{<:Number})

convenience function for optaining SVD, SPImtx, and SPItree
"""
function spiresult(A::Matrix{<:Number})
    usv = svd(A)
    spimtx = calc_spr_mtx(usv)
    spitree = hclust(spimtx; linkage=:average, branchorder=:optimal)
    return SPIResult(usv,spimtx,spitree)
end


"""
    getintervals(S::Vector{<:Number}; alpha=1.5, q=.75)

finds spectral partitions. Computes log difference between each subsequent singular
value and by default selects the differences that are larger than `1.5 * Q3(differences)`

i.e. finds breaks in the spectrum that explain smaller scales of variance

Args:
* S = singular values of a SVD factorization
* alpha = scalar multiple of `q`
* q = which quantile of log differences to use; by default Q3 

Returns:
* Vector{UnitRange} indices into S corresponding to the spectral partitions

"""
function getintervals(S::Vector{<:Number}; alpha=1.5, q=.75)
    potentialbreaks = abs.(diff(log.(S.+1)))
    θ = alpha * quantile(potentialbreaks, q)
    breaks = findall(potentialbreaks .> θ) .+ 1
    starts, ends = vcat(1, breaks), vcat(breaks.-1, length(S))
    intervals = map((s,e)->s:e, starts, ends)
    return intervals
end


"""
    calc_spi_mtx(A::Matrix; [Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.5, q=.75])
    calc_spi_mtx(usv::SVD; [Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.5, q=.75])
    calc_spi_mtx(usv::SVD[, SPI.LSVs(); Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.5, q=.75])
    calc_spi_mtx(usv::SVD[, SPI.RSVs(); Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.5, q=.75])


computes the cumulative spectral residual distance for spectral phylogenetic inference

```(∑_{p ∈ P} ||UₚΣₚ||₂)²```

where ``P`` are the spectral partitions found with `getintervals`. 

Args:
* A,usv = Matrix or SVD factorization (Matrix is just passed to `svd()` before calculation)
* SPI.Left() computes SPI matrix for LSVs; SPI.Right() computes SPI Matrix for RSVs
* alpha, q are passed to `getintervals()` see its documentation

Returns:
* distance matrix
"""
calc_spi_mtx(A::Matrix{<:Number}; alpha=1.5, q=.75) = calc_spi_mtx(svd(A); alpha, q)
calc_spi_mtx(A::Matrix{<:Number}, ::LSVs; alpha=1.5, q=.75) = calc_spi_mtx(svd(A), LSVs(); alpha, q)
calc_spi_mtx(A::Matrix{<:Number}, ::RSVs; alpha=1.5, q=.75) = calc_spi_mtx(svd(A), RSVs(); alpha, q)

calc_spi_mtx(usv::SVD; alpha=1.5, q=.75) = calc_spi_mtx(usv, LSVs(); alpha, q)
calc_spi_mtx(usv::SVD, ::LSVs; alpha=1.5, q=.75) = calc_spi_mtx(usv.U, usv.S; alpha, q)
calc_spi_mtx(usv::SVD, ::RSVs; alpha=1.5, q=.75) = calc_spi_mtx(Matrix(usv.V), usv.S; alpha, q)
function calc_spi_mtx(vecs::Matrix{<:T}, vals::Vector{<:T}; alpha=1.5, q=.75) where T<:Number
    sprmtx = zeros(size(vecs,1), size(vecs,1))
    for grp in getintervals(vals, alpha=alpha, q=q)
        sprmtx += Distances.pairwise(WeightedEuclidean(vals[grp]), vecs'[grp,:])
    end
    return sprmtx.^2
end


"""
    calc_spi_trace(usv::SVD; alpha=1.5, q=.75)
    calc_spi_trace(usv::SVD[, taxaidxs]; alpha=1.5, q=.75)

calculates spectral residual within each partition of spectrum and each pair of taxa

if `taxaidxs` are provided the `U` matrix is subset and/or reordered based on those indices.

returns matrix where rows are spectral partitions and columns are taxa:taxa pairs
ordered as the upper triangle in rowwise order. 

"""
function calc_spi_trace(usv::SVD; alpha=1.5, q=.75)
    groups = getintervals(usv.S; alpha=alpha, q=q)
    Nsmps = size(usv.U,1)
    r = zeros(binomial(Nsmps, 2), length(groups))
    for (k, (i,j)) in enumerate(combinations(1:Nsmps, 2))
        for (g, grp) in enumerate(groups)
            r[g, k] = Distances.pairwise(WeightedEuclidean(usv.S[grp]), usv.U'[grp, i], usv.U'[grp, j])
        end
    end
    return r
end

function calc_spi_trace(usv::SVD, taxaidxs; alpha=1.5, q=.75)
    groups = getintervals(usv.S; alpha=alpha, q=q)
    Usubset = @view usv.U[taxaidxs,:]
    Nsmps = size(Usubset,1)
    r = zeros(binomial(Nsmps, 2), length(groups))
    for (k, (i,j)) in enumerate(combinations(1:Nsmps, 2))
        for (g, grp) in enumerate(groups)
            r[g, k] = Distances.pairwise(WeightedEuclidean(usv.S[grp]), Usubset'[grp, i], Usubset'[grp, j])
        end
    end
    return r
end

"""
    calc_spi_tree(A[, ids; labelinternalnodes=true])
helper function that immediately returns the newick tree string inferred by SPI
"""
function calc_spi_tree(A; labelinternalnodes=true)
    dij = calc_spi_mtx(A)
    hc = hclust(dij, linkage=:average, branchorder=:optimal)
    nwstr(hc; labelinternalnodes)
end
function calc_spi_tree(A, ids; labelinternalnodes=true)
    dij = calc_spi_mtx(A)
    hc = hclust(dij, linkage=:average, branchorder=:optimal)
    nwstr(hc, ids; labelinternalnodes)
end

"""
    projectinLSV(data::AbstractArray{T}, usv::SVD{T}, [window])

returns estimated left singular vectors (aka: LSV or Û) for new data based on already calculated SVD factorization
"""
projectinLSV(data::AbstractArray{T}, usv::SVD{T}) where T<:Number  = data * usv.V * inv(diagm(usv.S))
projectinLSV(data::AbstractArray{T}, usv::SVD{T}, window) where T<:Number = data * usv.V[:, window] * inv(diagm(usv.S[window]))

"""
    projectinRSV(data::AbstractArray{T}, usv::SVD{T}, [window])

returns estimated transposed right singular vectors (RSV or V̂ᵗ) for new data based on already calculated SVD factorization
"""
projectinRSV(data::AbstractArray{T}, usv::SVD{T}) where T<:Number = inv(diagm(usv.S)) * usv.U' * data
projectinRSV(data::AbstractArray{T}, usv::SVD{T}, window) where T<:Number = inv(diagm(usv.S[window])) * usv.U'[window,:] * data

"""
    projectout(usv::SVD, [window])

recreates original matrix i.e. calculates ``UΣV'`` or if window is included 
creates a spectrally filtered version of the original matrix off of the provided components in `window`.
"""
projectout(usv::SVD) = usv.U * diagm(usv.S) * usv.Vt
projectout(usv::SVD, window) = usv.U[:, window] * diagm(usv.S[window]) * usv.Vt[window, :]


