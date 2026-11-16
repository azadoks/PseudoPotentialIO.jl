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
end
