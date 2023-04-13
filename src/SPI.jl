module SPI

using Reexport
using StatsBase
using SparseArrays
using DataFrames
using Distances
using Clustering
using LinearAlgebra
using Printf: @sprintf
using FreqTables: freqtable
using CategoricalArrays: cut

@reexport using LinearAlgebra
@reexport using Clustering: hclust

const ALPHA = 1.0
const QUANT = 0.5

export spiresult, nwstr
export getintervals, getintervalsIQR, 
    calc_spi_mtx, calc_spi_trace, calc_spcorr_mtx,
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

end
