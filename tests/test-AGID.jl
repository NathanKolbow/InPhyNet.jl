using InPhyNet, PhyloNetworks, NamedArrays, Test, Logging

Logging.disable_logging(Logging.Warn)

gt1 = readTopology("((A:1,B:1):1,(C:1,D:1):1);")
gt2 = readTopology("((A:2,B:2):1,(C:2,D:2):1);")
gt3 = readTopology("(E:1,F:1):1;")

@testset "internodedistance" begin
    # Test single gene tree calls
    D, namelist = internodedistance(gt1)
    n = NamedArray(D, (namelist, namelist))

    @test n["A", "B"] == 2
    @test n["C", "D"] == 2
    @test n["A", "C"] == 4
    @test n["A", "D"] == 4
    @test n["B", "C"] == 4
    @test n["B", "D"] == 4


    D, namelist = internodedistance(gt2)
    n = NamedArray(D, (namelist, namelist))

    @test n["A", "B"] == 4
    @test n["C", "D"] == 4
    @test n["A", "C"] == 6
    @test n["A", "D"] == 6
    @test n["B", "C"] == 6
    @test n["B", "D"] == 6
end

@testset "internodecount" begin

end

@testset "AGIC" begin
end

@testset "AGID" begin
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
    n = NamedArray(D, (namelist, namelist))

    @test n["A", "E"] == Inf
    @test n["A", "B"] == 1
    @test n["E", "F"] == 0  # 0 instead of 1 b/c the tree is treated as unrooted
end

@testset "AGID - missing pairs" begin
    D, namelist = calculateAGID([gt1, gt3], allow_missing_pairs=true, default_missing_value=Inf)
    n = NamedArray(D, (namelist, namelist))

    @test n["A", "E"] == Inf
    @test n["E", "F"] == 2
    @test n["A", "B"] == 2
end