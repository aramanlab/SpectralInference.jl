module SPI

using Reexport
using Printf
using StatsBase
using SparseArrays
using DataFrames
using Distances
using Clustering
using LinearAlgebra
using FreqTables: freqtable
using CategoricalArrays: cut
using Combinatorics: combinations

@reexport using LinearAlgebra
@reexport using Clustering: hclust

export spiresult
export getintervals, calc_spi_mtx, calc_spi_trace, 
    projectinLSV, projectinRSV, projectout
include("core.jl")

export zscore, scaledcumsum, nwstr
    minspaceneeded, spimtx_spaceneeded,
include("helpers.jl")

export empiricalMI, adjustedrandindex, vmeasure_homogeneity_completeness
include("empiricalMI.jl")

include("parsephylip.jl")

end
