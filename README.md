# Spectral Phylogenetic Inference (SPI)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aramanlab.github.io/SPI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aramanlab.github.io/SPI.jl/dev)
[![Build Status](https://github.com/aramanlab/SPI.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aramanlab/SPI.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aramanlab/SPI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aramanlab/SPI.jl)

## Install

type:

```julia-repl
] dev git@github.com:aramanlab/SPI.jl.git
```
into the julia command line

## Example

```julia
using SPI

M = rand(100, 100)
usv = svd(M) # from linear algebra
spimtx = calc_spi_mtx(usv; alpha=1.5, q=0.75)
spitree = hclust(spimtx; linkage=:average, branchorder=:optimal) # hclust from Clustering

# or in one line
SPIresult = spiresult(M)
# with result having these fields
# SPIresult.usv, SPIresult.spimtx, SPIresult.spitree
```

example script for running SPI on a gene expression matrix in a `.csv` file

```julia
# Example SPI code for gene expression matrix in .csv format (genes = rows, cells = columns)
# Ben Doran / Noah Gamble
# 2022/04/13

# Load packages
using SPI
using CSV, DataFrames
using LinearAlgebra
using Clustering
using Muon
using JLD2

# Load matrix as data frame from .csv
df = CSV.read("/path/to/geneexpressionmatrix.csv", DataFrame)

# Copy to matrix
m = Matrix(df[:,2:end]) # 

filt_thresh = 10
b = (m .> 0)
idx_keep = vec(sum(b, dims = 2) .> filt_thresh)
m = m[idx_keep,:]

# Organize into annotated data object
adata = AnnData(X = m')
adata.obs_names .= names(df)[2:end]
adata.var_names .= string.(df.Column1[idx_keep])
@info size(adata.X)

# Perform SVD
@info "Peforming SVD..."
@time usv = svd(adata.X)

# Perform SPI
@info "Calculating SPI matrix..."
@time dij = calc_spi_mtx(usv)

# Cluster to tree
@info "Clustering SPI matrix..."
@time spitree = hclust(dij; linkage=:average, branchorder=:optimal)

# Converting to Newick string
@info "Converting to Newick string..."
@time nwtreestring = nwstr(spitree, adata.obs_names.vals)

# Add SPI matrix to adata object
@info "Adding SPI matrix to ADATA..."
adata.obsp["SPI_mtx"] = dij
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


## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
