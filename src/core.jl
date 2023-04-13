
struct SPIResult{T}
    usv::SVD{T}
    spimtx::AbstractMatrix{T}
    spitree::Hclust
end

struct LSVs end
struct RSVs end

"""
    spiresult(A::AbstractMatrix{<:Number})

convenience function for optaining SVD, SPImtx, and SPItree
"""
function spiresult(A::AbstractMatrix{<:Number})
    usv = svd(A)
    spimtx = calc_spr_mtx(usv)
    spitree = hclust(spimtx; linkage=:average, branchorder=:optimal)
    return SPIResult(usv,spimtx,spitree)
end


"""
    getintervals(S::AbstractVector{<:Number}; alpha=1.0, q=0.5)

finds spectral partitions. Computes log difference between each subsequent singular
value and by default selects the differences that are larger than `1.0 * Q2(differences)`

i.e. finds breaks in the spectrum that explain smaller scales of variance

Args:
* S = singular values of a SVD factorization
* alpha = scalar multiple of `q`
* q = which quantile of log differences to use; by default Q2 

Returns:
* AbstractVector{UnitRange} indices into S corresponding to the spectral partitions

"""
function getintervals(S::AbstractVector{<:Number}; alpha=1.0, q=0.5)
    potentialbreaks = abs.(diff(log.(S.+1)))
    θ = alpha * quantile(potentialbreaks, q)
    breaks = findall(potentialbreaks .> θ) .+ 1
    starts, ends = vcat(1, breaks), vcat(breaks.-1, length(S))
    intervals = map((s,e)->s:e, starts, ends)
    return intervals
end


"""
    getintervalsIQR(S::AbstractVector{<:Number}; alpha=1.5, ql=.25, qh=.75)

finds spectral partitions. Computes log difference between each subsequent singular
value and by default selects the differences that are larger than `ALPHA * Q3(differences)`

i.e. finds breaks in the spectrum that explain smaller scales of variance

Args:
* S = singular values of a SVD factorization
* alpha = scalar multiple of `q`
* q = which quantile of log differences to use; by default Q3 

Returns:
* AbstractVector{UnitRange} indices into S corresponding to the spectral partitions

"""
function getintervalsIQR(S::AbstractVector{<:Number}; alpha=1.5, ql=.25, qh=.75)
    potentialbreaks = abs.(diff(log.(S.+1)))
    Q1, Q3 = quantile(potentialbreaks, [ql, qh])
    θ = Q3 + alpha * (Q3 - Q1)
    breaks = findall(potentialbreaks .> θ) .+ 1
    starts, ends = vcat(1, breaks), vcat(breaks.-1, length(S))
    intervals = map((s,e)->s:e, starts, ends)
    return intervals
end


"""
    calc_spi_mtx(A::AbstractMatrix; [Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.0, q=0.5])
    calc_spi_mtx(usv::SVD; [Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.0, q=0.5])
    calc_spi_mtx(usv::SVD[, SPI.LSVs(); Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.0, q=0.5])
    calc_spi_mtx(usv::SVD[, SPI.RSVs(); Nsmps=size(A,1), Nfeats=size(A,2), alpha=1.0, q=0.5])


computes the cumulative spectral residual distance for spectral phylogenetic inference

```(∑_{p ∈ P} ||UₚΣₚ||₂)²```

where ``P`` are the spectral partitions found with `getintervals`. 

Args:
* A,usv = AbstractMatrix or SVD factorization (AbstractMatrix is just passed to `svd()` before calculation)
* SPI.Left() computes SPI matrix for LSVs; SPI.Right() computes SPI matrix for RSVs
* alpha, q are passed to `getintervals()` see its documentation

Returns:
* distance matrix
"""
calc_spi_mtx(A::AbstractMatrix{<:Number}; alpha=ALPHA, q=QUANT) = calc_spi_mtx(svd(A); alpha, q)
calc_spi_mtx(A::AbstractMatrix{<:Number}, ::LSVs; alpha=ALPHA, q=QUANT) = calc_spi_mtx(svd(A), LSVs(); alpha, q)
calc_spi_mtx(A::AbstractMatrix{<:Number}, ::RSVs; alpha=ALPHA, q=QUANT) = calc_spi_mtx(svd(A), RSVs(); alpha, q)

calc_spi_mtx(usv::SVD; alpha=ALPHA, q=QUANT) = calc_spi_mtx(usv, LSVs(); alpha, q)
calc_spi_mtx(usv::SVD, ::LSVs; alpha=ALPHA, q=QUANT) = calc_spi_mtx(usv.U, usv.S; alpha, q)
calc_spi_mtx(usv::SVD, ::RSVs; alpha=ALPHA, q=QUANT) = calc_spi_mtx(AbstractMatrix(usv.V), usv.S; alpha, q)
function calc_spi_mtx(vecs::AbstractMatrix{<:T}, vals::AbstractVector{<:T}; alpha=ALPHA, q=QUANT) where T<:Number
    intervals = getintervals(vals; alpha, q)
    calc_spi_mtx(vecs, vals, intervals)
end
function calc_spi_mtx(vecs::AbstractMatrix{<:T}, vals::AbstractVector{<:T}, intervals::AbstractVector) where T<:Number
    sprmtx = zeros(size(vecs,1), size(vecs,1))
    for grp in intervals
        sprmtx += Distances.pairwise(WeightedEuclidean(vals[grp]), vecs'[grp,:]; dims=2)
    end
    return sprmtx.^2
end


"""
    calc_spi_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)
    calc_spi_trace(vecs, vals, groups)

calculates spectral residual within each partition of spectrum and each pair of taxa

returns matrix where rows are spectral partitions and columns are taxa:taxa pairs
ordered as the upper triangle in rowwise order, or lower triangle in colwise order.

Args:
* method: `calc_spi_trace(vecs, vals, groups)`
    * vecs: either usv.U or usv.V matrix
    * vals: usv.S singular values vector
    * groups: usually calculated with `getintervals(usv.S; alpha=alpha, q=q)`
* method: `calc_spi_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)`     
    * usv: SVD object
    * onrows: switch to calculate SPI on rows (U matrix) or columns (V matrix).
    * groups: if nothing groups are calculated with `getintervals(usv.S; alpha=alpha, q=q)`, 
        otherwise they assume a vector of index ranges `[1:1, 2:3, ...]` to group `usv.S` with. 
    * alpha: passed to `getintervals`
    * q: passed to `getintervals`
"""
function calc_spi_trace(usv::SVD; onrows=true, groups=nothing, alpha=ALPHA, q=QUANT)
    groups = isnothing(groups) ? getintervals(usv.S; alpha=alpha, q=q) : groups
    vecs = onrows ? usv.U : usv.V
    r = zeros(length(groups), binomial(size(vecs, 1), 2))
    calc_spi_trace!(r, vecs, usv.S, groups)
    return r
end

function calc_spi_trace(vecs, vals, groups)
    Nsmps = size(vecs,1)
    r = zeros(length(groups), binomial(Nsmps, 2))
    calc_spi_trace!(r, vecs, vals, groups)
    return r
end

function calc_spi_trace!(r, vecs, vals, groups)
    Nsmps = size(vecs,1)
    for (k,(i,j)) in enumerate(((i, j) for j in 1:Nsmps for i in (j+1):Nsmps))
        for (g, grp) in enumerate(groups)
            r[g, k] = WeightedEuclidean(vals[grp])(vecs'[grp, i], vecs'[grp, j])
        end
    end
end


"""
    calc_spi_tree(A[, ids; labelinternalnodes=true])
helper function that immediately returns the newick tree string inferred by SPI
"""
function calc_spi_tree(A; labelinternalnodes=false)
    dij = calc_spi_mtx(A)
    hc = UPGMA_tree(dij)
    nwstr(hc; labelinternalnodes)
end
function calc_spi_tree(A, ids; labelinternalnodes=false)
    dij = calc_spi_mtx(A)
    hc = UPGMA_tree(dij)
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

"""
    UPGMA_tree(Dij::AbstractMatrix{<:Number})

shorthand for `Clustering.hclust(Dij, linkage=:average, branchorder=:optimal)`
"""
function UPGMA_tree(Dij::AbstractMatrix{<:Number}) 
    Clustering.hclust(Dij, linkage=:average, branchorder=:optimal)
end


"""
    calc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, window)
    calc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, vals::AbstractVector{<:Number}, window)

Calculates pairwise spectral (pearson) correlations for a set of observations. 

Args:
* `vecs`, set of left singular vectors or principal components with observations on rows and components on columns
* `vals`, vector of singular values
* `window`, set of indices of `vecs` columns to compute correlations across

Returns:
* correlation matrix where each pixel is the correlation between a pair of observations
"""
function calc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, window)
    cor(vecs[:, window]')
end
function calc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, vals::AbstractVector{<:Number}, window)
    contributions = vecs[:, window] * diagm(vals[window])
    cor(contributions')
end

""" 
    nwstr(hc::Hclust[, tiplabels::AbstractVector[<:String]])

convert Hclust to newick tree string
Args:
* hc, `Hclust` object from Clustering package
* tiplabels, `AbstractVector{<:String}` names in same order as distance matrix
"""
nwstr(hc::Hclust; labelinternalnodes=false) = nwstr(hc, string.(1:length(hc.order)); labelinternalnodes)
nwstr(hc::Hclust, tiplabels::AbstractVector{<:Symbol}; labelinternalnodes=false) = nwstr(hc, string.(tiplabels); labelinternalnodes)
function nwstr(hc::Hclust, tiplabels::AbstractVector{<:String}; labelinternalnodes=false)
    r = length(hc.heights)
    _nwstr(view(hc.merges, :, :), view(hc.heights, :), r, r, view(tiplabels, :); labelinternalnodes) * ";"
end
function _nwstr(merges::A, heights::B, i::C, p::C, tiplabels::D; labelinternalnodes=false)::String where {
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