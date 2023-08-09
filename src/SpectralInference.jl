module SpectralInference

using Reexport
using StatsBase: indicatormat, quantile, cor, fit, Histogram, entropy
using Distances: Distances, WeightedEuclidean
using Clustering: Clustering, hclust, Hclust
using LinearAlgebra: svd, SVD
using Printf: @sprintf
using FreqTables: freqtable
using CategoricalArrays: cut

@reexport using LinearAlgebra
@reexport using Clustering: hclust

const ALPHA = 1.0
const QUANT = 0.5

export spiresult, newickstring
export getintervals, getintervalsIQR, UPGMA_tree,
    spectraldistances, spectraldistances_trace, spectralcorrelations,
    projectinLSV, projectinRSV, projectout
include("core.jl")

export explainedvariance, scaledcumsum, 
    minspaceneeded, spimtx_spaceneeded,
    squareform
include("helpers.jl")

export empiricalMI, adjustedrandindex, vmeasure_homogeneity_completeness
include("empiricalMI.jl")

export readphylip, onehotencode
include("parsephylip.jl")

# extension functions
export getleafnames,
    network_distance, 
    network_distances, 
    patristic_distance, 
    patristic_distances,
    fscore_precision_recall, 
    cuttree,
    mapinternalnodes, 
    maplocalnodes, 
    collectiveLCA,
    as_polytomy, 
    as_polytomy!,
    pairedMI_across_treedepth, 
    clusters_per_cutlevel,
    ladderize!,
    spectral_lineage_encoding
include("treefunctions.jl")

include("precompile_workload.jl")

end

