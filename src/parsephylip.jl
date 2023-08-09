"""
    readphylip(fn::String)

Read phylip alignment file, return dataframe of IDs and Sequences
"""
function readphylip(fn::AbstractString)
    smps, seqs = open(fn, "r") do reader
        nsmp, nfeat = parse.(Int, split(readline(reader)))
        smps = String[]
        seqs = String[]
        for l in 1:nsmp
            smp, seq = split(readline(reader))
            push!(smps, smp)
            push!(seqs, seq)
        end
        all(s->length(s)==nfeat, seqs) || throw(DimensionMismatch("Length of all sequences must be the equal"))
        return smps, seqs
    end
    return (;smps, seqs)
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
    reduce(vcat, permutedims.(collect.(seqs)))
end