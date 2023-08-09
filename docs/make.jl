using SpectralInference
using Documenter

DocMeta.setdocmeta!(SpectralInference, :DocTestSetup, :(using SpectralInference); recursive=true)

makedocs(;
    modules=[SpectralInference],
    authors="Benjamin Doran and collaborators",
    repo="https://github.com/aramanlab/SpectralInference.jl/blob/{commit}{path}#{line}",
    sitename="SpectralInference.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aramanlab.github.io/SpectralInference.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aramanlab/SpectralInference.jl",
    devbranch="main",
)
