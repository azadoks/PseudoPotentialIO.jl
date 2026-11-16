using PseudoPotentialIO: psp8_to_dict, dict_to_psp8, _to_dict

@testset "PSP8 dictionaries" begin
    text = PSP8_TEXTS["nlcc"]
    d = psp8_to_dict(text)
    h = d["header"]
    @test h["title"] == "X  test file" && h["zatom"] === 3.0 && h["pspd"] === 250101
    @test h["lloc"] === 4 && h["r2well"] === 0.0 && h["fchrg"] === 0.5
    @test !haskey(h, "mmax") && !haskey(h, "nproj")  # the sizes of the data
    @test d["rgrid"] == [0.0, 0.01, 0.02]
    @test only(d["projectors"]) == Dict("l" => 0, "ekb" => [2.5], "projectors" => [[1.0, 2.0, 3.0]])
    @test d["v_local"] == [-1.0, -2.0, -3.0]
    @test d["nlcc"]["d4_rhoc_dr4"] == [5.0, 5.0, 5.0]
    @test d["valence"] == Dict("rhov" => [7.0, 8.0, 9.0])
    @test d["trailing_text"] == "<INPUT>\n</INPUT>\n"

    io = IOBuffer()
    dict_to_psp8(io, d)
    @test psp8_to_dict(String(take!(io))) == d

    # the struct: an entry for every l
    file = PseudoPotentialIO.Psp8File(IOBuffer(text))
    @test file.header.nproj == [1, 0] && file.header.mmax == 3
    @test file.projectors == [[[1.0, 2.0, 3.0]], []] && file.ekb == [[2.5], []]
    @test isnothing(file.ae_rhov)

    # the local potential in place of l = 1
    text1 = replace(text, "8   11   1   4" => "8   11   1   1", "\n4\n" => "\n1\n")
    @test psp8_to_dict(text1)["v_local"] == d["v_local"]
    # blocks out of order, projectors at l = lloc, not psp8
    @test_throws ErrorException psp8_to_dict(replace(text, "\n4\n" => "\n3\n"))
    @test_throws ErrorException psp8_to_dict(replace(text1, "1  0  0  0  0    nproj" => "1  1    nproj"))
    @test_throws ErrorException psp8_to_dict(replace(text, "8   11" => "3   11"))
    @test_throws ErrorException psp8_to_dict(first(text, 300))  # cut short

    # spin-orbit (blocks for l >= 1), and the all-electron valence and core densities
    so = psp8_to_dict(PSP8_TEXTS["spin-orbit"])
    @test only(so["projectors_so"]) == Dict("l" => 1, "ekb" => [0.5],
                                            "projectors" => [[0.0, 0.1, 0.2]])
    @test so["valence"]["ae_rhoc"] == [70.0, 80.0, 90.0] && !haskey(so, "nlcc")
    file = load_psp8("spin-orbit")
    @test file.header.nproj == [1, 1] && file.header.nprojso == [0, 1]
    @test PseudoPotentialIO.has_spin_orbit(file) && functional(file) == [:lda_x, :lda_c_pw]

    @testset "Round trip: $name" for (name, text) in PSP8_TEXTS
        d = psp8_to_dict(text)
        io = IOBuffer()
        dict_to_psp8(io, d)
        @test psp8_to_dict(String(take!(io))) == d
        delete!(d, "trailing_text")  # not in the struct
        @test _to_dict(load_psp8(name)) == d
    end
end
