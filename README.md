# Spectral Inference

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aramanlab.github.io/SpectralInference.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aramanlab.github.io/SpectralInference.jl/dev)
[![Build Status](https://github.com/aramanlab/SpectralInference.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aramanlab/SpectralInference.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aramanlab/SpectralInference.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aramanlab/SpectralInference.jl)

## Install

type:

```julia-repl
] dev git@github.com:aramanlab/SpectralInference.jl.git
```
into the julia command line

## Examples

calculate Spectral distances
```julia
using SpectralInference

M = rand(100, 100)
usv = svd(M) # from linear algebra
spectralpartitions = getintervals(usv.S, alpha=1.0, q=0.5)
Dij = spectraldistances(usv.U, usv.S, spectralpartitions)
```
make a tree from the distances using basic hierarchical clusting
```julia
# shortcut of `hclust(spimtx; linkage=:average, branchorder=:optimal)` from Clustering package
spitree = UPGMA_tree(Dij)
```

or with neighbor joining
```julia
using NeighborJoining # https://github.com/BenjaminDoran/NeighborJoining.jl
njclusts = regNJ(Dij) # or fastNJ(Dij)
ids = string.(1:100);
nwstring = NeighborJoining.newickstring(njclusts, ids)
```

example script for running SpectralInference on a gene expression matrix in a `.csv` file to generate spectral tree.

```julia
# Example SpectralInference code for gene expression matrix in .csv format (genes = rows, cells = columns)
# Ben Doran / Noah Gamble
# 2023/08/07

# Load packages
using SpectralInference
using NeighborJoining
using CSV, DataFrames
using LinearAlgebra
using Clustering
using Muon
using JLD2

# Load matrix as data frame from .csv
df = CSV.read("/path/to/geneexpressionmatrix.csv", DataFrame)

# Copy to matrix
m = Matrix(df[:,2:end]) # 

# filter out low abundance genes
filt_thresh = 10
genecounts = mapslices(r -> sum(r.>0), m, dims=2) |> vec
idx_keep = genecounts .> filt_thresh
m = m[idx_keep,:]

# Organize into annotated data object
# m' is transpose of m (make samples rows)
adata = AnnData(X = m') 
adata.obs_names .= names(df)[2:end]
adata.var_names .= string.(df.Column1[idx_keep])
@info size(adata.X)

# Perform SVD
@info "Peforming SVD..."
@time usv = svd(adata.X)

# Perform SpectralInference
@info "Calculating SpectralInference matrix..."
@time dij = spectraldistances(usv.U, usv.S, getintervals(usv.S))

# Cluster to tree
@info "Clustering SpectralInference matrix..."
# @time spitree = UPGMA_tree(dij)
# @time nwtreestring = SpectralInference.newickstring(spitree, adata.obs_names.vals)
@time spitree = regNJ(dij)
@time nwtreestring = NeighborJoining.newickstring(spitree, adata.obs_names.vals)

# Add SpectralInference matrix to adata object
@info "Adding SpectralInference matrix to ADATA..."
adata.obsp["spectraldistances"] = dij
adata.obs[:,"order"] = spitree.order
adata.uns["cellmerges"] = spitree.merges
adata.uns["heights"] = spitree.heights
adata.uns["treelinkage"] = string(spitree.linkage)
adata.uns["newicktreestring"] = nwtreestring

# Write annotated data to and Newick tree to disk
@info "Writing results to disk..."
@time begin 
	writeh5ad("spi_test.h5ad", adata)
	open("test_tree.nw.txt", "w") do io
		println(io, nwtreestring)
	end
	jldsave("test_tree.jld2"; spitree)
end

```

example script for creating bootstrap spectral tree and calculating support values

```julia
using SpectralInference
using NeighborJoining
using Muon

spitst = readh5ad("spi_test.h5ad")
mtx = spitst.X[:, :]
NBOOT = 100
boottrees = map(1:NBOOT) do i
	ncols = size(mtx, 2)
	tmpmtx = mtx[:, rand(1:ncols, ncols)]
	usv = svd(tmpmtx)
	dij = spectraldistances(usv.U, usv.S, getintervals(usv.S))
	nwtreestring = NeighborJoining.newickstring(regNJ(dij))
end
open("test_tree_bootstraps.nw.txt", "w") do io
	for t in boottrees
		println(io, t)
	end
end

# calculate support values
using GoTree_jll
reftreefile = joinpath(pwd(), "test_tree.nw.txt")
boottreefile = joinpath(pwd(), "test_tree_bootstraps.nw.txt")
supporttreefile = joinpath(pwd(), "test_tree_bootstraps.nw.txt")
run(`$(gotree()) -i $reftreefile -b $boottreefile -o $supporttreefile`)
```

script for calculating MI between metavariables and tree depths

```julia
using SpectralInference
using NewickTree
using Muon, DataFrames

spitst = readh5ad("spi_test.h5ad")
spitree = readnw(readline(joinpath(pwd(), "test_tree_bootstraps.nw.txt")))
spitree_50pct = as_polytomy(n->NewickTree.support(n)<.5, spitree)
leafnames = getleafnames(spitree_50pct)
mi, treedepths = pairedMI_across_treedepth(eachcolumn(spitst.obs), leafnames, spitree_50pct; ncuts=100)


NBOOTS = 50
results = map(1:NBOOTS) do i
	pairedMI_across_treedepth(eachcolumn(spitst.obs), leafnames, spitree_50pct; bootstrap=true, ncuts=100)
end


rowmask = spitst.obs.count .> 15
results = map(1:NBOOTS) do i
	pairedMI_across_treedepth(eachcolumn(spitst.obs), leafnames, spitree_50pct; mask=rowmask, bootstrap=true, ncuts=100)
end
all_mi_results = first.(results)
all_treedepths_results = last.(results)

```



## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
