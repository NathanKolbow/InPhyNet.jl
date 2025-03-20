using InPhyNet, PhyloNetworks, Test

gt1 = readnewick("((A:1,B:1):1,(C:1,D:1):1);")
gt2 = readnewick("((A:2,B:2):1,(C:2,D:2):1);")
gt3 = readnewick("(E:1,F:1):1;")

@testset "internodedistance" begin
    # Test single gene tree calls
    D, namelist = internodedistance(gt1)
    @test all(j -> namelist[j] == sort(tiplabels(gt1))[j], 1:length(namelist))

    @test D[1, 2] == 2
    @test D[3, 4] == 2
    @test D[1, 3] == 4
    @test D[1, 4] == 4
    @test D[2, 3] == 4
    @test D[2, 4] == 4


    D, namelist = internodedistance(gt2)
    @test all(j -> namelist[j] == sort(tiplabels(gt2))[j], 1:length(namelist))

    @test D[1, 2] == 4
    @test D[3, 4] == 4
    @test D[1, 3] == 6
    @test D[1, 4] == 6
    @test D[2, 3] == 6
    @test D[2, 4] == 6
end

@testset "AGID - missing pairs" begin
    D, namelist = calculateAGID([gt1, gt3], allow_missing_pairs=true, default_missing_value=Inf)
    @test all(j -> namelist[j] == ["A", "B", "C", "D", "E", "F"][j], 1:length(namelist))

    @test D[1, 5] == Inf
    @test D[5, 6] == 2
    @test D[1, 2] == 2
end


@testset "AGIC - missing pairs" begin
    # When some taxa do not appear together, make sure errors are thrown and default values are set
    try
        D, namelist = calculateAGIC([gt1, gt3])
        @test false
    catch e
        @test typeof(e) <: ArgumentError
    end

    
    D, namelist = calculateAGIC([gt1, gt3], allow_missing_pairs=true, default_missing_value=Inf)
    @test all(j -> namelist[j] == ["A", "B", "C", "D", "E", "F"][j], 1:length(namelist))

    @test D[1, 5] == Inf
    @test D[1, 2] == 1
    @test D[5, 6] == 0  # 0 instead of 1 b/c the tree is treated as unrooted
end

