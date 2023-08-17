
"""
    getintervals(S::AbstractVector{<:Number}; alpha=1.0, q=0.5)

finds spectral partitions. Computes log difference between each subsequent singular
value and by default selects the differences that are larger than `1.0 * Q2(differences)`

i.e. finds breaks in the spectrum that explain smaller scales of variance

Args:
* S: singular values of a SVD factorization
* alpha: scalar multiple of `q`
* q: which quantile of log differences to use; by default Q2 

Returns:
* AbstractVector{UnitRange} indices into S corresponding to the spectral partitions

"""
function getintervals(S::AbstractVector{<:Number}; alpha=1.0, q=0.5)
    potentialbreaks = abs.(diff(log.(S .+ 1)))
    θ = alpha * quantile(potentialbreaks, q)
    breaks = findall(potentialbreaks .> θ) .+ 1
    starts, ends = vcat(1, breaks), vcat(breaks .- 1, length(S))
    intervals = map((s, e) -> s:e, starts, ends)
    return intervals
end


"""
    getintervals_IQR(S::AbstractVector{<:Number}; alpha=1.5, ql=.25, qh=.75)

finds spectral partitions. Computes log difference between each subsequent singular
value and by default selects the differences that are larger than `1.5 * IQR(differences)`

i.e. finds breaks in the spectrum that explain smaller scales of variance

Args:
* S: singular values of a SVD factorization
* alpha: scalar multiple of `q`
* q: which quantile of log differences to use; by default Q3 

Returns:
* AbstractVector{UnitRange} indices into S corresponding to the spectral partitions

"""
function getintervalsIQR(S::AbstractVector{<:Number}; alpha=1.5, ql=0.25, qh=0.75)
    potentialbreaks = abs.(diff(log.(S .+ 1)))
    Q1, Q3 = quantile(potentialbreaks, [ql, qh])
    θ = Q3 + alpha * (Q3 - Q1)
    breaks = findall(potentialbreaks .> θ) .+ 1
    starts, ends = vcat(1, breaks), vcat(breaks .- 1, length(S))
    intervals = map((s, e) -> s:e, starts, ends)
    return intervals
end


"""
    spectraldistances(A::AbstractMatrix; [onrows=true, alpha=1.0, q=0.5])
    spectraldistances(usv::SVD; [onrows=true, alpha=1.0, q=0.5])
    spectraldistances(vecs::AbstactMatrix, vals::AbstractMatrix, intervals::AbstractVector{<:UnitRange})


computes the cumulative spectral residual distance for spectral phylogenetic inference

```(∑_{p ∈ P} ||UₚΣₚ||₂)²```

where ``P`` are the spectral partitions found with `getintervals`. 

Args:
* A, usv: AbstractMatrix or SVD factorization (AbstractMatrix is passed to `svd()` before calculation)
* onrows: if true will compute spectral distances on the left singular vectors (U matrix), if false will calculate on the right singular vectors or (V matrix)
* alpha, q: are passed to `getintervals()` see its documentation

Returns:
* distance matrix
"""
function spectraldistances(A::AbstractMatrix{<:Number}; onrows=true::Bool, alpha=ALPHA, q=QUANT)
    spectraldistances(svd(A); onrows, alpha, q)
end
function spectraldistances(usv::SVD; onrows=true::Bool, alpha=ALPHA, q=QUANT)
    if onrows
        spectraldistances(usv.U, usv.S; alpha, q)
    else
        spectraldistances(usv.V, usv.S; alpha, q)
    end
end
function spectraldistances(vecs::AbstractMatrix{<:T}, vals::AbstractVector{<:T}; alpha=ALPHA, q=QUANT) where {T<:Number}
    intervals = getintervals(vals; alpha, q)
    spectraldistances(vecs, vals, intervals)
end
function spectraldistances(vecs::AbstractMatrix{<:T}, vals::AbstractVector{<:T}, intervals::AbstractVector) where {T<:Number}
    spimtx = zeros(size(vecs, 1), size(vecs, 1))
    for grp in intervals
        spimtx += Distances.pairwise(WeightedEuclidean(vals[grp]), vecs'[grp, :]; dims=2)
    end
    return spimtx .^ 2
end


"""
    spectraldistances_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)
    spectraldistances_trace(vecs, vals, groups)

calculates spectral residual within each partition of spectrum and each pair of taxa

returns matrix where rows are spectral partitions and columns are taxa:taxa pairs
ordered as the upper triangle in rowwise order, or lower triangle in colwise order.

Args:
* method: `spectraldistances_trace(vecs, vals, groups)`
    * vecs: either usv.U or usv.V matrix
    * vals: usv.S singular values vector
    * groups: usually calculated with `getintervals(usv.S; alpha=alpha, q=q)`
* method: `spectraldistances_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)`     
    * usv: SVD object
    * onrows: true/false switch to calculate spectral distance on rows (U matrix) or columns (V matrix).
    * groups: if nothing groups are calculated with `getintervals(usv.S; alpha=alpha, q=q)`, 
        otherwise they assume a vector of index ranges `[1:1, 2:3, ...]` to group `usv.S` with. 
    * alpha: passed to `getintervals`
    * q: passed to `getintervals`
"""
function spectraldistances_trace(usv::SVD; onrows=true, groups=nothing, alpha=ALPHA, q=QUANT)
    groups = isnothing(groups) ? getintervals(usv.S; alpha=alpha, q=q) : groups
    vecs = onrows ? usv.U : usv.V
    r = zeros(length(groups), binomial(size(vecs, 1), 2))
    spectraldistances_trace!(r, vecs, usv.S, groups)
    return r
end

function spectraldistances_trace(vecs, vals, groups)
    Nsmps = size(vecs, 1)
    r = zeros(length(groups), binomial(Nsmps, 2))
    spectraldistances_trace!(r, vecs, vals, groups)
    return r
end

function spectraldistances_trace!(r, vecs, vals, groups)
    Nsmps = size(vecs, 1)
    for (k, (i, j)) in enumerate(((i, j) for j in 1:Nsmps for i in (j+1):Nsmps))
        for (g, grp) in enumerate(groups)
            r[g, k] = WeightedEuclidean(vals[grp])(vecs'[grp, i], vecs'[grp, j])
        end
    end
end


"""
    projectinLSV(data::AbstractArray{T}, usv::SVD{T}, [window])

returns estimated left singular vectors (aka: LSV or Û) for new data based on already calculated SVD factorization
"""
projectinLSV(data::AbstractArray, usv::SVD) = data * usv.V * inv(diagm(usv.S))
projectinLSV(data::AbstractArray, usv::SVD, window) = data * usv.V[:, window] * inv(diagm(usv.S[window]))

"""
    projectinRSV(data::AbstractArray, usv::SVD, [window])

returns estimated transposed right singular vectors (RSV or V̂ᵗ) for new data based on already calculated SVD factorization
"""
projectinRSV(data::AbstractArray, usv::SVD) = inv(diagm(usv.S)) * usv.U' * data
projectinRSV(data::AbstractArray, usv::SVD, window) = inv(diagm(usv.S[window])) * usv.U'[window, :] * data

"""
    projectout(usv::SVD, [window])

recreates original matrix i.e. calculates ``UΣV'`` or if window is included 
creates a spectrally filtered version of the original matrix off of the provided components in `window`.

i.e., `usv.U[:, window] * diagm(usv.S[window]) * usv.Vt[window, :]`
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
* vecs: set of left or right singular vectors with observations/features on rows and spectral components on columns
* vals: vector of singular values
* window: set of indices of `vecs` columns to compute correlations across

Returns:
* correlation matrix where each pixel is the correlation between a pair of observations
"""
function spectralcorrelations(vecs::AbstractMatrix{<:Number}, window)
    cor(vecs[:, window]')
end
function spectralcorrelations(vecs::AbstractMatrix{<:Number}, vals::AbstractVector{<:Number}, window)
    contributions = vecs[:, window] * diagm(vals[window])
    cor(contributions')
end

""" 
    newickstring(hc::Hclust[, tiplabels::AbstractVector[<:String]])

convert Hclust to newick tree string

Args:
* hc: `Hclust` object from Clustering package
* tiplabels: `AbstractVector{<:String}` names in same order as distance matrix

Returns:
* [newick tree](https://en.wikipedia.org/wiki/Newick_format) formated string
"""
newickstring(hc::Hclust; labelinternalnodes=false) = newickstring(hc, string.(1:length(hc.order)); labelinternalnodes)
newickstring(hc::Hclust, tiplabels::AbstractVector{<:Symbol}; labelinternalnodes=false) = newickstring(hc, string.(tiplabels); labelinternalnodes)
function newickstring(hc::Hclust, tiplabels::AbstractVector{<:AbstractString}; labelinternalnodes=false)
    r = length(hc.heights)
    _newickstring(view(hc.merges, :, :), view(hc.heights, :), r, r, view(tiplabels, :); labelinternalnodes) * ";"
end
function _newickstring(merges::A, heights::B, i::C, p::C, tiplabels::D; labelinternalnodes=false)::String where {
    A<:AbstractArray{<:Integer},B<:AbstractVector{<:AbstractFloat},
    C<:Integer,D<:AbstractVector{<:AbstractString}
}
    j::Int64 = merges[i, 1] # left subtree pointer
    k::Int64 = merges[i, 2] # right subtree pointer
    a::String = if j < 0 # if tip format tip
        tiplabels[abs(j)] * ':' * @sprintf("%e", heights[i])
    else # recurse and format internal node
        _newickstring(view(merges, :, :), view(heights, :), j, i, view(tiplabels, :); labelinternalnodes)
    end
    b::String = if k < 0 # if tip format tip
        tiplabels[abs(k)] * ':' * @sprintf("%e", heights[i])
    else # recurse and format internal node
        _newickstring(view(merges, :, :), view(heights, :), k, i, view(tiplabels, :); labelinternalnodes)
    end
    nid = labelinternalnodes ? "node" * string(length(heights) + i + 1) : ""
    dist = @sprintf("%e", heights[p] - heights[i])
    _newick_merge_strings(a, b, nid, dist)
end
function _newick_merge_strings(a::S, b::S, n::S, d::S) where {S<:String}
    '(' * a * ',' * b * ')' * n * ':' * d
end