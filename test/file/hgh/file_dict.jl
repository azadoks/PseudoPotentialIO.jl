using PseudoPotentialIO: hgh_to_dict, dict_to_hgh, _to_dict

@testset "HGH dictionaries" begin
    text = """
    Ni GTH-PADE-q18 GTH-LDA-q18
        4    6    8   zion
         0.35000000    2     3.61031072     0.44963832
        3
         0.24510489    2    12.16113071     3.51625420
                                           -9.07892931
         0.23474136    0
         0.21494950    1   -13.39506212
    """
    d = hgh_to_dict(text)
    @test d["title"] == "Ni GTH-PADE-q18 GTH-LDA-q18"
    @test d["zion"] == [4, 6, 8]  # the comment is dropped
    @test d["rloc"] === 0.35 && d["cloc"] == [3.61031072, 0.44963832]
    @test d["nonlocal"][1] == Dict("rp" => 0.24510489,
                                   "h" => [[12.16113071, 3.51625420], [-9.07892931]])
    @test d["nonlocal"][2]["h"] == []
    @test length(d["nonlocal"]) == 3

    io = IOBuffer()
    dict_to_hgh(io, d)
    @test hgh_to_dict(String(take!(io))) == d

    # the size of the data
    @test_throws ErrorException hgh_to_dict(replace(text, "    2     3.61" => "    3     3.61"))
    @test_throws ErrorException hgh_to_dict(replace(text, "-9.07892931" => ""))
    @test_throws ErrorException hgh_to_dict(replace(text, "\n    3\n" => "\n    4\n"))

    @testset "Round trip: $name" for (name, text) in HGH_TEXTS
        d = hgh_to_dict(text)
        io = IOBuffer()
        dict_to_hgh(io, d)
        @test hgh_to_dict(String(take!(io))) == d
        @test _to_dict(load_hgh(name)) == d
    end

    @testset "ABINIT" begin
        d = hgh_to_dict(HGH_TEXTS["gd-abinit.hgh"])
        @test d["zatom"] == 64 && d["zion"] == 18 && d["pspdat"] == 10605
        @test (d["pspcod"], d["pspxc"], d["lmax"]) == (3, 1, 3)
        # pspcod = 3: the diagonals; l = 3 has only h11 and k11
        @test d["nonlocal"][2] == Dict("rp" => 0.456953, "h" => [-0.058347, -1.697711, 0.0],
                                       "k" => [0.53554, 0.022024, -0.054301])
        @test d["nonlocal"][4] == Dict("rp" => 0.27339, "h" => [-21.92349], "k" => [0.014666])
        @test !haskey(d["nonlocal"][1], "k")

        # pspcod = 10, as ABINIT reads it: the numbers each line needs (nproj = 1 for l = 1,
        # so the k11 of l = 1 is the next line)
        text = """
        Goedecker pseudopotential for Si
           14   4  070301 zatom,zion,pspdat
        10 11  1 0 2001 0  pspcod,pspxc,lmax,lloc,mmax,r2well
             0.44000000    1    -6.26928833                    rloc nloc c1
            2                                                  nnonloc
             0.43563383    2     8.95174150    -2.70627082     rs ns hs11 hs12
                                                3.49378060     hs22
             0.49794218    1     2.43127673     0.00000000     rp np hp11 hp12
                                                0.00000000     hp22
        """
        d = hgh_to_dict(text)
        @test d["nonlocal"][2] == Dict("rp" => 0.49794218, "h" => [[2.43127673]], "k" => [[0.0]])
        io = IOBuffer()
        dict_to_hgh(io, d)
        @test hgh_to_dict(String(take!(io))) == d

        # the ABINIT formats are pspcod 3 and 10
        @test_throws ErrorException hgh_to_dict(replace(text, "10 11  1" => " 2 11  1"))
        # and a line must have the numbers it needs
        @test_throws ErrorException hgh_to_dict(replace(text, "    0.00000000     hp22" => ""))
    end
end
