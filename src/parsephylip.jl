"""
    readphylip(fn::String)

Read phylip alignment file, return dataframe of IDs and Sequences
"""
function readphylip(fn::AbstractString)
    smps, seqs = open(fn, "r") do reader
        nsmp, nfeat = parse.(Int, split(readline(reader)))
        smps = Vector()
        seqs = Vector()
        for l in 1:nsmp
            smp, seq = split(readline(reader))
            push!(smps, smp)
            push!(seqs, seq)
        end
        all(s->length(s)==nfeat, seqs) || throw(DimensionMismatch("Length of all sequences must be the equal"))

        ord = sortperm(smps)
        return string.(smps[ord]), seqs[ord]
    end
    return DataFrame("ids" => string.(smps), "seqs" => string.(seqs))
end # read phylip



onehotencode(seqs::AbstractVector{<:AbstractString}) = onehotencode(_stringcolumntocharmtx(seqs))
function onehotencode(chardf::AbstractMatrix{<:Char})
    ohemtx = Vector()
    for col in eachcol(chardf)
        push!(ohemtx, indicatormat(col)')
    end
    return sparse(hcat(ohemtx...))
end

function _stringcolumntocharmtx(seqs)
    Matrix(DataFrame(reduce(hcat, collect.(seqs)) |> permutedims, [string(i) for i in 1:length(first(seqs))]))
end