var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SpectralInference","category":"page"},{"location":"#SpectralInference","page":"Home","title":"SpectralInference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SpectralInference.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SpectralInference]","category":"page"},{"location":"#SpectralInference.UPGMA_tree-Tuple{AbstractMatrix{<:Number}}","page":"Home","title":"SpectralInference.UPGMA_tree","text":"UPGMA_tree(Dij::AbstractMatrix{<:Number})\n\nshorthand for Clustering.hclust(Dij, linkage=:average, branchorder=:optimal)\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.adjustedrandindex-Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}}","page":"Home","title":"SpectralInference.adjustedrandindex","text":"adjustedrandindex(a::AbstractVector{<:Number}, b::AbstractVector{<:Number}; nbins=50)\n\nArgs:\n\na, vector of numbers\nb, vector of numbers\nnbins, for continuous approximates discrete, for discrete choose nbins>maxnumberof_classes\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.distancematrix_spaceneeded-Tuple{Any}","page":"Home","title":"SpectralInference.distancematrix_spaceneeded","text":"distancematrix_spaceneeded(n, p; bits=64) = Base.format_bytes(binomial(n,2) * p * bits)\n\nhow much memory is needed to store distance matrix Args:\n\nn: number of samples\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.distancetrace_spaceneeded-Tuple{Any, Any}","page":"Home","title":"SpectralInference.distancetrace_spaceneeded","text":"distancetrace_spaceneeded(n, p; bits=64) = Base.format_bytes(binomial(n,2) * p * bits)\n\nhow much memory is needed to store spectral residual trace\n\nArgs:\n\nn: number of samples\np: number of partitions/components\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.empiricalMI-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}}} where T<:AbstractFloat","page":"Home","title":"SpectralInference.empiricalMI","text":"empiricalMI(a::AbstractVector{<:Float}, b::AbstractVector{<:Float}[, nbins=50, normalize=false])\n\ncomputes empirical MI from identity of H(a) + H(b) - H(ab). where H = -sum(p(x)*log(p(x))) + log(Δ) the + log(Δ) corresponds to the log binwidth and unbiases the entropy estimate from binwidth choice. estimates are roughly stable from 32 (32^2  1000 total bins) to size of sample. going from a small undersestimate to a small overestimate across that range. We recommend choosing the sqrt(mean(1000, samplesize)) for nbins argument, or taking a few estimates across that range and averaging.\n\nArgs:\n\na, vecter of length N\nb, AbstractVector of length N\nnbins, number of bins per side, use 1000 < nbins^2 < length(a) for best results\nbase, base unit of MI (defaults to nats with base=ℯ)\nnormalize, bool, whether to normalize with mi / mean(ha, hb)\n\nReturns:\n\nMI\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.explainedvariance-Tuple{AbstractVector{<:Number}}","page":"Home","title":"SpectralInference.explainedvariance","text":"explainedvariance(s::AbstractVector{<:Number})\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.getintervals-Tuple{AbstractVector{<:Number}}","page":"Home","title":"SpectralInference.getintervals","text":"getintervals(S::AbstractVector{<:Number}; alpha=1.0, q=0.5)\n\nfinds spectral partitions. Computes log difference between each subsequent singular value and by default selects the differences that are larger than 1.0 * Q2(differences)\n\ni.e. finds breaks in the spectrum that explain smaller scales of variance\n\nArgs:\n\nS: singular values of a SVD factorization\nalpha: scalar multiple of q\nq: which quantile of log differences to use; by default Q2 \n\nReturns:\n\nAbstractVector{UnitRange} indices into S corresponding to the spectral partitions\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.getintervalsIQR-Tuple{AbstractVector{<:Number}}","page":"Home","title":"SpectralInference.getintervalsIQR","text":"getintervals_IQR(S::AbstractVector{<:Number}; alpha=1.5, ql=.25, qh=.75)\n\nfinds spectral partitions. Computes log difference between each subsequent singular value and by default selects the differences that are larger than 1.5 * IQR(differences)\n\ni.e. finds breaks in the spectrum that explain smaller scales of variance\n\nArgs:\n\nS: singular values of a SVD factorization\nalpha: scalar multiple of q\nq: which quantile of log differences to use; by default Q3 \n\nReturns:\n\nAbstractVector{UnitRange} indices into S corresponding to the spectral partitions\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.ij2k-Tuple{Any, Any, Any}","page":"Home","title":"SpectralInference.ij2k","text":"ij2k(i,j,n)\n\nwith pair (i,j) give index k to the pairs produced by combinations(vec), where vec is length n\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.k2ij-Tuple{Any, Any}","page":"Home","title":"SpectralInference.k2ij","text":"k2ij(k, n)\n\nwhich pair (i,j) produces the kth element of combinations(vec), where vec is length n\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.newickstring-Tuple{Clustering.Hclust}","page":"Home","title":"SpectralInference.newickstring","text":"newickstring(hc::Hclust[, tiplabels::AbstractVector[<:String]])\n\nconvert Hclust to newick tree string\n\nArgs:\n\nhc: Hclust object from Clustering package\ntiplabels: AbstractVector{<:String} names in same order as distance matrix\n\nReturns:\n\nnewick tree formated string\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.pairwise-Tuple{Function, AbstractMatrix}","page":"Home","title":"SpectralInference.pairwise","text":"pairwise(func::Function, m::AbstractMatrix)\n\nreturns the lower columnwise offdiagonal of result[k] = func(i, j)  where k is the kth pair and i and j are the ith and kth columns of m  calculated from enumerate(((i, j) for j in axes(m, 2) for i in (j+1):lastindex(m, 2)))\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.projectinLSV-Tuple{AbstractArray, SVD}","page":"Home","title":"SpectralInference.projectinLSV","text":"projectinLSV(data::AbstractArray{T}, usv::SVD{T}, [window])\n\nreturns estimated left singular vectors (aka: LSV or Û) for new data based on already calculated SVD factorization\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.projectinRSV-Tuple{AbstractArray, SVD}","page":"Home","title":"SpectralInference.projectinRSV","text":"projectinRSV(data::AbstractArray, usv::SVD, [window])\n\nreturns estimated transposed right singular vectors (RSV or V̂ᵗ) for new data based on already calculated SVD factorization\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.projectout-Tuple{SVD}","page":"Home","title":"SpectralInference.projectout","text":"projectout(usv::SVD, [window])\n\nrecreates original matrix i.e. calculates UΣV or if window is included  creates a spectrally filtered version of the original matrix off of the provided components in window.\n\ni.e., usv.U[:, window] * diagm(usv.S[window]) * usv.Vt[window, :]\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.readphylip-Tuple{AbstractString}","page":"Home","title":"SpectralInference.readphylip","text":"readphylip(fn::String)\n\nRead phylip alignment file, return dataframe of IDs and Sequences\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.scaledcumsum-Tuple{Any}","page":"Home","title":"SpectralInference.scaledcumsum","text":"scaledcumsum(c; dims=1)\n\ncumsum divided by maximum cumulative value\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.spectralcorrelations-Tuple{AbstractMatrix{<:Number}, Any}","page":"Home","title":"SpectralInference.spectralcorrelations","text":"calc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, window)\ncalc_spcorr_mtx(vecs::AbstractMatrix{<:Number}, vals::AbstractVector{<:Number}, window)\n\nCalculates pairwise spectral (pearson) correlations for a set of observations. \n\nArgs:\n\nvecs: set of left or right singular vectors with observations/features on rows and spectral components on columns\nvals: vector of singular values\nwindow: set of indices of vecs columns to compute correlations across\n\nReturns:\n\ncorrelation matrix where each pixel is the correlation between a pair of observations\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.spectraldistances-Tuple{AbstractMatrix{<:Number}}","page":"Home","title":"SpectralInference.spectraldistances","text":"spectraldistances(A::AbstractMatrix; [onrows=true, alpha=1.0, q=0.5])\nspectraldistances(usv::SVD; [onrows=true, alpha=1.0, q=0.5])\nspectraldistances(vecs::AbstactMatrix, vals::AbstractMatrix, intervals::AbstractVector{<:UnitRange})\n\ncomputes the cumulative spectral residual distance for spectral phylogenetic inference\n\n(∑_{p ∈ P} ||UₚΣₚ||₂)²\n\nwhere P are the spectral partitions found with getintervals. \n\nArgs:\n\nA, usv: AbstractMatrix or SVD factorization (AbstractMatrix is passed to svd() before calculation)\nonrows: if true will compute spectral distances on the left singular vectors (U matrix), if false will calculate on the right singular vectors or (V matrix)\nalpha, q: are passed to getintervals() see its documentation\n\nReturns:\n\ndistance matrix\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.spectraldistances_trace-Tuple{SVD}","page":"Home","title":"SpectralInference.spectraldistances_trace","text":"spectraldistances_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)\nspectraldistances_trace(vecs, vals, groups)\n\ncalculates spectral residual within each partition of spectrum and each pair of taxa\n\nreturns matrix where rows are spectral partitions and columns are taxa:taxa pairs ordered as the upper triangle in rowwise order, or lower triangle in colwise order.\n\nArgs:\n\nmethod: spectraldistances_trace(vecs, vals, groups)\nvecs: either usv.U or usv.V matrix\nvals: usv.S singular values vector\ngroups: usually calculated with getintervals(usv.S; alpha=alpha, q=q)\nmethod: spectraldistances_trace(usv::SVD; onrows=true, groups=nothing, alpha=1.0, q=0.5)     \nusv: SVD object\nonrows: true/false switch to calculate spectral distance on rows (U matrix) or columns (V matrix).\ngroups: if nothing groups are calculated with getintervals(usv.S; alpha=alpha, q=q),    otherwise they assume a vector of index ranges [1:1, 2:3, ...] to group usv.S with. \nalpha: passed to getintervals\nq: passed to getintervals\n\n\n\n\n\n","category":"method"},{"location":"#SpectralInference.squareform","page":"Home","title":"SpectralInference.squareform","text":"squareform(d::AbstractVector, fillvalue=zero(eltype(d)))\nsquareform(d::AbstractVector, fillvalue=zero(eltype(d)))\n\nIf d is a vector, squareform checks if it of n choose 2 length for integer n, then fills the values of a symetric square matrix with the values of d.\n\nIf d is a matrix, squareform checks if it is square then fills the values of vector with the lower offdiagonal of matrix d in column order form. \n\nfillvalue is the initial value of the produced vector or matrix. Only really apparant in a produced matrix where it will be the values on the diagonal.\n\n\n\n\n\n","category":"function"},{"location":"#SpectralInference.vmeasure_homogeneity_completeness-Tuple{Any, Any}","page":"Home","title":"SpectralInference.vmeasure_homogeneity_completeness","text":"vmeasure_homogeneity_completeness(labels_true, labels_pred; β=1.)\n\ncalculates and returns v-measure, homogeneity, completeness; similar to f-score, precision, and recall respectively\n\nArgs:\n\nβ, weighting term for v-measure, if β is greater than 1 completeness\n\nis weighted more strongly in the calculation, if β is less than 1,  homogeneity is weighted more strongly\n\nCitation:\n\nA. Rosenberg, J. Hirschberg, in Proceedings of the 2007 Joint Conference\n\non Empirical Methods in Natural Language Processing and Computational Natural  Language Learning (EMNLP-CoNLL) (Association for Computational Linguistics,   Prague, Czech Republic, 2007; https://aclanthology.org/D07-1043), pp. 410–420.\n\n\n\n\n\n","category":"method"}]
}
