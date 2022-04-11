# SPI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aramanlab.github.io/SPI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aramanlab.github.io/SPI.jl/dev)
[![Build Status](https://github.com/aramanlab/SPI.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aramanlab/SPI.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aramanlab/SPI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aramanlab/SPI.jl)

## Install

```
dev git@github.com:aramanlab/SPI.jl.git
```

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


## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
