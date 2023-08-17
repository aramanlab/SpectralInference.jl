using SpectralInference
using Test
using Distributions

@testset "SpectralInference" begin
@testset "core" begin
    Dij_true = [
        0.0     2.0     5.3642  5.3642
        2.0     0.0     5.3642  5.3642
        5.3642  5.3642  0.0     2.0
        5.3642  5.3642  2.0     0.0
    ]
    
    M = Float64.([
        0 1 0 1 1 1
        0 1 1 0 1 1
        1 0 1 1 0 1
        1 0 1 1 1 0
    ])

    Mr2 = [
        -3.93055e-16   1.0         0.5  0.5  1.0  1.0
         9.01494e-16   1.0         0.5  0.5  1.0  1.0
         1.0          -1.5782e-16  1.0  1.0  0.5  0.5
         1.0          -1.5782e-16  1.0  1.0  0.5  0.5
    ]

    expS = [
        10.0, 3.0, 3.0, 1., 1., 1., 0.,
    ]   



    usv = svd(M)
    Dij_pred = spectraldistances(usv.U, usv.S, [1:1, 2:2, 3:4])
    @inferred spectraldistances(usv.U, usv.S, [1:1, 2:2, 3:4])
    # tests for getting partitions
    @test getintervals(expS) == [1:1, 2:3, 4:6, 7:7]
    @inferred getintervals(expS) == [1:1, 2:3, 4:6, 7:7]
    @test getintervalsIQR(expS) == [1:7]
    
    # test if SpectralInference distance works
    @test isapprox(Dij_pred, Dij_true, atol=1e-5)
    
    # test if we get the same SpectralInference tree
    hc = SpectralInference.UPGMA_tree(Dij_pred)
    nw = newickstring(hc, ["A", "B", "C", "D"], labelinternalnodes=false)
    @test replace(nw, r"(:)[^,)]*(?=[,);])"=>"") ∈ [
        "((A,B),(C,D));", "((A,B),(D,C));", "((B,A),(C,D));", "((B,A),(D,C));",
        "((C,D),(A,B));", "((D,C),(A,B));", "((C,D),(B,A));", "((D,C),(B,A));"
    ]
    # @test nw == "((B:2.000000e+00,A:2.000000e+00):3.364199e+00,(D:2.000000e+00,C:2.000000e+00):3.364199e+00):0.000000e+00;"

    # test projecting in and out of SVD space
    @test projectout(usv) ≈ M
    @test projectout(usv, 1:2) ≈ Mr2
    @test projectinLSV(M, usv) ≈ usv.U
    @test projectinLSV(M, usv, 1:3) ≈ usv.U[:, 1:3]
    @test projectinRSV(M, usv) ≈ usv.Vt
    @test projectinRSV(M, usv, 1:3) ≈ usv.Vt[1:3, :]

    # test partial spi_trace calculations result in summed SpectralInference distance matrix
    @test isapprox(vec(sum(spectraldistances_trace(usv,  q=.25), dims=1) .^ 2), Dij_true[triu(trues(4,4), 1)], atol=1e-5)
    # test spectral correlations 
    @inferred spectralcorrelations(usv.U, usv.S, 1:4)
end



@testset "helpers.jl" begin
    expS = [
        10.0, 3.0, 3.0, 1., 1., 1., 0.,
    ]
    v = explainedvariance(expS)
    @test v ≈ (expS .^ 2) ./ sum((expS .^ 2))
    @test scaledcumsum(v) ≈ cumsum(v) ./ maximum(cumsum(v))

    @test SpectralInference.distancetrace_spaceneeded(100,100; bits=64) == "30.212 MiB"
    @test SpectralInference.distancetrace_spaceneeded(1000,1000; bits=64) == "29.773 GiB"

    @test SpectralInference.distancematrix_spaceneeded(100; bits=64) == "625.000 KiB"
    @test SpectralInference.distancematrix_spaceneeded(1000; bits=64) == "61.035 MiB"
    @test SpectralInference.distancematrix_spaceneeded(1000; bits=32) == "30.518 MiB"

    M = Float64.([
        0 1 0 1 1 1
        0 1 1 0 1 1
        1 0 1 1 0 1
        1 0 1 1 1 0
    ])
    ps = SpectralInference.pairwise((i,j) -> i' * j , M')
    @test ps ≈ [3.0,2.0,2.0,2.0,2.0,3.0]

    @test SpectralInference.numpairs2N(6) == 4.0
    @test SpectralInference.numpairs2N(5) ≈ 3.7015621187164243

    @test squareform(ps) ≈ [
        0.0  3.0  2.0  2.0
        3.0  0.0  2.0  2.0
        2.0  2.0  0.0  3.0
        2.0  2.0  3.0  0.0
    ]

    @test squareform(ps, 4.) ≈ [
        4.0  3.0  2.0  2.0
        3.0  4.0  2.0  2.0
        2.0  2.0  4.0  3.0
        2.0  2.0  3.0  4.0
    ]

    
    @test SpectralInference.k2ij(4,4) == (3,2)
    @test SpectralInference.k2ij(6,4) == (4,3)
    @test SpectralInference.ij2k(3,2,4) == 4
    @test SpectralInference.ij2k(4,3,4) == 6
end

@testset "empiricalMI" begin
    
    x = [0, 0, 0, 1, 1, 1]
    for T in [Int64, Int32, Int16]
        t = T.(x)
        @test empiricalMI(t, t, base=2) == 1.0
        @test empiricalMI(t, t) ≈ 0.6931471805599453 atol=1e-10
    end

    # test masked MI
    x = randn(1000)
    for T in [Float64, Float32, Float16]
        t = T.(x)
        @test empiricalMI([t;t.+10.0], [trues(1000); falses(1000)], base=2) ≈ 1.0 atol=1e-3
    end

    
    ρ = 0.8
    mvN = MultivariateNormal([1 ρ; ρ 1])
    trueMI_e = (2*entropy(Normal(0,1))) - entropy(mvN) # entropy as marginals - joint
    trueMI_2 = (2*entropy(Normal(0,1), 2)) - entropy(mvN, 2) # entropy as marginals - joint
    x = rand(mvN, 10_000)'
    @test empiricalMI(x[:, 1], x[:,2], nbins=32) ≈ trueMI_e atol=1e-1
    @test empiricalMI(x[:, 1], x[:,2], nbins=32, base=2) ≈ trueMI_2 atol=1e-1

end


if VERSION >= v"1.9"
    @test isnothing(Base.get_extension(SpectralInference, :NewickTreeExt))
    using NewickTree
    @test !isnothing(Base.get_extension(SpectralInference, :NewickTreeExt))

    @testset "NewickTreeExt" begin
        m = rand(100, 40)
        usv = svd(m)
        dij = spectraldistances(usv.U, usv.S, getintervals(usv.S))
        hc = UPGMA_tree(dij)
        tree = readnw(newickstring(hc, string.(1:100)))
        leafnames = getleafnames(tree)
        @inferred network_distances(tree)
        @inferred patristic_distances(tree)
        @test fscore_precision_recall(tree, tree) == (1., 1., 1.)
        @inferred cuttree(network_distance, tree, 0.)
        @inferred Node collectiveLCA(getleaves(tree))
        @inferred Node as_polytomy(n->NewickTree.support(n)<0.5, tree)
        mi, treedepths = pairedMI_across_treedepth((;a=rand(Bool,100), b=rand([1,0], 100)), leafnames, tree; ncuts=50)
        @test length(treedepths) == length(mi.a) == length(mi.b)
        mi, treedepths = pairedMI_across_treedepth((;a=rand(UInt16,100), b=rand(Float64, 100)), leafnames, tree; ncuts=50, comparefun=(x,y)->abs(x-y))
        @test length(treedepths) == length(mi.a) == length(mi.b)
        mi, treedepths = pairedMI_across_treedepth((;a=rand(UInt16,100), b=rand(Float64, 100)), leafnames, tree; ncuts=50, treecut_distancefun=patristic_distance, comparefun=(x,y)->abs(x-y))
        @test length(treedepths) == length(mi.a) == length(mi.b)
        leafids = getleafids(tree)
        @inferred spectral_lineage_encoding(tree)
        @inferred spectral_lineage_encoding(tree, leafids)
        @inferred spectral_lineage_encoding(tree; filterfun=!isleaf)
        @inferred spectral_lineage_encoding(tree, sample(leafids, length(leafids), replace=false); filterfun=!isleaf)
    end
end # VERSION >= v"1.9"

end # SpectralInference testset