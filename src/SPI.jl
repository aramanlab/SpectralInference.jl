module SPI

using Reexport
using Distances
using StatsBase
using Clustering
using LinearAlgebra

@reexport using LinearAlgebra
@reexport using Clustering: hclust

export spiresult
export getintervals, calc_spi_mtx, calc_spi_trace, 
    projectinLSV, projectinRSV, projectout
include("core.jl")

export zscore, scaledcumsum,
    minspaceneeded, spimtx_spaceneeded
include("helpers.jl")

end
