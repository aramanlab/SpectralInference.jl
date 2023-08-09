using PrecompileTools

@setup_workload begin
    @compile_workload begin
        for T in [Float64, Float32, Int64, Int32]
            m = rand(T, 10, 10)
            usv = svd(m)
            dij = spectraldistances(usv.U, usv.S, getintervals(usv.S))
            dij = spectraldistances(usv, alpha=1.0, q=0.5)
            dij_trace = spectraldistances_trace(usv, alpha=1.0, q=0.5)
            dij_trace = spectraldistances_trace(usv.U, usv.S, getintervalsIQR(usv.S))
            spcorr = spectralcorrelations(usv.U, 1:4)
            spcorr = spectralcorrelations(usv.U, usv.S, 1:4)
            t = UPGMA_tree(dij)
            newickstring(t)
        end
    end
end