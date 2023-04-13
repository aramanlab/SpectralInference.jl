using SPI
using Test

2 + 2

@testset "SPI.jl core" begin
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
    Dij_pred = calc_spi_mtx(usv.U, usv.S, [1:1, 2:2, 3:4])
    @inferred calc_spi_mtx(usv.U, usv.S, [1:1, 2:2, 3:4])
    # tests for getting partitions
    @test getintervals(expS) == [1:1, 2:3, 4:6, 7:7]
    @inferred getintervals(expS) == [1:1, 2:3, 4:6, 7:7]
    @test getintervalsIQR(expS) == [1:7]
    
    # test if SPI distance works
    @test isapprox(Dij_pred, Dij_true, atol=1e-5)
    
    # test if we get the same SPI tree
    hc = SPI.UPGMA_tree(Dij_pred)
    nw = nwstr(hc, ["A", "B", "C", "D"], labelinternalnodes=false)
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

    # test partial spi_trace calculations result in summed SPI distance matrix
    @test isapprox(vec(sum(calc_spi_trace(usv,  q=.25), dims=1) .^ 2), Dij_true[triu(trues(4,4), 1)], atol=1e-5)
    @inferred calc_spcorr_mtx(usv.U, usv.S, 1:4)
end

@testset "SPI.jl helpers" begin
    expS = [
        10.0, 3.0, 3.0, 1., 1., 1., 0.,
    ]
    v = explainedvariance(expS)
    @test v ≈ (expS .^ 2) ./ sum((expS .^ 2))
    @test scaledcumsum(v) ≈ cumsum(v) ./ maximum(cumsum(v))

    @test SPI.minspaceneeded(100,100; bits=64) == "30.212 MiB"
    @test SPI.minspaceneeded(1000,1000; bits=64) == "29.773 GiB"

    @test SPI.spimtx_spaceneeded(100; bits=64) == "625.000 KiB"
    @test SPI.spimtx_spaceneeded(1000; bits=64) == "61.035 MiB"
    @test SPI.spimtx_spaceneeded(1000; bits=32) == "30.518 MiB"

    M = Float64.([
        0 1 0 1 1 1
        0 1 1 0 1 1
        1 0 1 1 0 1
        1 0 1 1 1 0
    ])
    ps = SPI.pairwise((i,j) -> i' * j , M')
    @test ps ≈ [3.0,2.0,2.0,2.0,2.0,3.0]

    @test SPI.numpairs2N(6) == 4.0
    @test SPI.numpairs2N(5) ≈ 3.7015621187164243

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

    
    @test SPI.k2ij(4,4) == (3,2)
    @test SPI.k2ij(6,4) == (4,3)
    @test SPI.ij2k(3,2,4) == 4
    @test SPI.ij2k(4,3,4) == 6
end
