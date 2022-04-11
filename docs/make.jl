using SPI
using Documenter

DocMeta.setdocmeta!(SPI, :DocTestSetup, :(using SPI); recursive=true)

makedocs(;
    modules=[SPI],
    authors="Benjamin Doran and collaborators",
    repo="https://github.com/aramanlab/SPI.jl/blob/{commit}{path}#{line}",
    sitename="SPI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aramanlab.github.io/SPI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aramanlab/SPI.jl",
    devbranch="main",
)
