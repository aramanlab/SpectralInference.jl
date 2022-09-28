### MI functions ###

"""
    empiricalMI(a::AbstractVector{<:Float}, b::AbstractVector{<:Float}[, nbins=50, normalize=false])

computes empirical MI from identity of ``H(a) + H(b) - H(a,b)``. where
``H := -sum(p(x)*log(p(x))) + log(Δ)``
the ``+ log(Δ)`` corresponds to the log binwidth and unbiases the entropy estimate from binwidth choice.
estimates are roughly stable from ``32`` (``32^2 ≈ 1000`` total bins) to size of sample. going from a small undersestimate to a small overestimate across that range.
We recommend choosing the `sqrt(mean(1000, samplesize))` for `nbins` argument, or taking a few estimates across that range and averaging.

Args:
* a, vecter of length N
* b, AbstractVector of length N
* nbins, number of bins per side, use 1000 < nbins^2 < length(a) for best results
* base, base unit of MI (defaults to nats with base=ℯ)
* normalize, bool, whether to normalize with mi / mean(ha, hb)

Returns:
* MI
"""
function empiricalMI(a::AbstractVector{T}, b::AbstractVector{T}; nbins=50, base=ℯ, normalize=false) where T<:AbstractFloat

    # num samples marginal then total
    N = length(a)
    breaks_a = range(minimum(a), maximum(a), length=nbins)
    breaks_b = range(minimum(b), maximum(b), length=nbins)

    ab_counts = fit(Histogram, (a, b), (breaks_a, breaks_b)).weights
    a_counts = sum(ab_counts; dims=1)
    b_counts = sum(ab_counts; dims=2)

    # frequency
    afreq = vec((a_counts) ./ (N))
    bfreq = vec((b_counts) ./ (N))
    abfreq =  vec((ab_counts) ./ (N))

    Δ_a = breaks_a[2] - breaks_a[1]
    Δ_b = breaks_b[2] - breaks_b[1]

    # approx entropy
    ha = entropy(afreq, base) + log(base, Δ_a)
    hb = entropy(bfreq, base) + log(base, Δ_b)
    hab = entropy(abfreq, base) + log(base, Δ_a * Δ_b)

    # mi
    mi = ha + hb - hab
    mi = normalize ? 2mi / (ha+hb) : mi
    # return (mi = mi, ha = ha, hb = hb, hab = hab)
    return mi
end

# value grouped by binary catagory
function empiricalMI(ab::AbstractVector{F}, mask::M; nbins=100, base=ℯ, normalize=false) where {F <: Number, M <: Union{BitVector, AbstractVector{<:Bool}}}
    length(ab) == length(mask) ||
        throw(ArgumentError("length of vals and meta must match; got vals=$(length(vals)), meta=$(length(meta))"))
    # num samples marginal then total
    N = length(ab)
    Na, Nb = sum(mask), (length(mask) - sum(mask))
    mask = Bool.(mask)
    # if mask is not grouping than there is no added information
    Na > 0 && Nb > 0 || return mi = 0.0
    
    ## otherwise ##

    # form edges
    edges = range(minimum(ab), maximum(ab), length=nbins)

    # fit hist and get counts
    a_counts = fit(Histogram, ab[mask], edges).weights
    b_counts = fit(Histogram, ab[.!mask], edges).weights
    ab_counts = a_counts .+ b_counts

    # get binwidth
    Δ = edges[2]-edges[1]

    # frequency
    afreq = vec((a_counts)./(Na)) #Δ/Δ
    bfreq = vec((b_counts)./(Nb)) #Δ/Δ
    abfreq =  vec((ab_counts)./(N)) #Δ/Δ

    # approx entropy
    ha = entropy(afreq, base) + log(base, Δ)
    hb = entropy(bfreq, base) + log(base, Δ)
    hab = entropy(abfreq, base) + log(base, Δ)

    # mi
    mi = hab - (Na/N*ha + Nb/N*hb) # original had flipped signs
    # return (mi = mi, ha = ha, hb = hb, hab = hab)
    mi = normalize ? 2mi / (ha+hb) : mi
    return mi
end

# value grouped by multiple categories
function empiricalMI(a::AbstractVector{V}, b::AbstractVector{M}; nbins=100, kwargs...) where 
    {V <: AbstractFloat,  M <: AbstractString}
    
    binned_a, a_bw = _bin_numbers(a, nbins)
    empiricalMI(binned_a, b; bw_a=a_bw, kwargs...)
end

# discrete catagories
empiricalMI(a::BitVector, b::BitVector; kwargs...) = empiricalMI(Int.(a), Int.(b); kwargs...)
function empiricalMI(a::AbstractVector{T}, b::AbstractVector{V}; bw_a=1., bw_b=1., base=ℯ, normalize=false) where 
    {T <: Union{Integer,AbstractString}, V <: Union{Integer,AbstractString}}
    
    counts = freqtable(a, b)
    N = sum(counts)
    Ha = entropy(sum(counts, dims=2)./N, base) + log(base, bw_a)
    Hb = entropy(sum(counts, dims=1)./N, base) + log(base, bw_b)
    Hab = entropy(counts./N, base) + log(base, bw_a * bw_b)
    mi = Ha + Hb - Hab
    mi = normalize ? 2mi / (Ha+Hb) : mi
end


_bin_numbers(v::AbstractVector{<:AbstractString}, nbins) = v, 1.
_bin_numbers(v::AbstractVector{<:Integer}, nbins) = string.(v), 1.
function _bin_numbers(v::AbstractVector{<:Number}, nbins)
    breaks = range(minimum(v), maximum(v), length=nbins)
    delta = breaks[2] - breaks[1]
    binned_v = cut(v, breaks; extend=true)
    return string.(binned_v), delta
end


"""
    vmeasure_homogeneity_completeness(labels_true, labels_pred; β=1.)

calculates and returns v-measure, homogeneity, completeness;
similar to f-score, precision, and recall respectively

Args:
* β, weighting term for v-measure, if β is greater than 1 completeness
 is weighted more strongly in the calculation, if β is less than 1,
 homogeneity is weighted more strongly

Citation:
* A. Rosenberg, J. Hirschberg, in Proceedings of the 2007 Joint Conference
 on Empirical Methods in Natural Language Processing and Computational Natural
 Language Learning (EMNLP-CoNLL) (Association for Computational Linguistics,
  Prague, Czech Republic, 2007; https://aclanthology.org/D07-1043), pp. 410–420.
"""
function vmeasure_homogeneity_completeness(labels_true, labels_pred; β=1.)
    if length(labels_true) == 0
        return 1.0, 1.0, 1.0
    end

    N = length(labels_true)
    entropy_C = entropy(freqtable(labels_true)./N)
    entropy_K = entropy(freqtable(labels_pred)./N)
    entropy_CK = entropy(freqtable(labels_true, labels_pred)./N)

    MI = entropy_C + entropy_K - entropy_CK

    homogeneity = entropy_C != 0. ? (MI / entropy_C) : 1.0
    completeness = entropy_K != 0. ? (MI / entropy_K) : 1.0

    if homogeneity + completeness == 0.
        v_measure_score = 0.0
    else
        v_measure_score = (
            (1 + β) * homogeneity * completeness / (β * homogeneity + completeness)
        )
    end
    return v_measure_score, homogeneity, completeness
end

"""
    adjustedrandindex(a::AbstractVector{<:Number}, b::AbstractVector{<:Number}; nbins=50)

Args:
* a, vector of numbers
* b, vector of numbers
* nbins, for continuous approximates discrete, for discrete choose nbins>max_number_of_classes

"""
function adjustedrandindex(a::AbstractVector{<:Number}, b::AbstractVector{<:Number}; nbins=50)

    # num samples marginal then total
    N = length(a)
    ab = vcat(a, b)
    breaks = range(minimum(ab), maximum(ab), length=nbins)

    # contingency and marginal counts
    ab_counts = fit(Histogram, (a, b), (breaks, breaks)).weights
    a_counts = sum(ab_counts; dims=1)
    b_counts = sum(ab_counts; dims=2)

    # count possible pairs
    contingentpairs = sum(binomial.(abcounts, 2))
    a_pairs = sum(binomial.(a_counts, 2))
    b_pairs = sum(binomial.(b_counts, 2))
    possiblepairs = binomial(N, 2)
    expectedpairs = a_pairs*b_pairs/possiblepairs

    # calculate ARI
    return (contingentpairs - expectedpairs) / ((a_pairs + b_pairs)/2 - expectedpairs) # ARI
end