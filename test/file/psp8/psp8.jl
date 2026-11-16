@testset "PSP8" begin
    @testset "pspxc = 2 is Perdew-Zunger with the rounded constants" begin
        # ABINIT's and ONCVPSP's `xcpzca` use c = 0.0020, d = -0.0116: Libxc's LDA_C_PZ
        i = findfirst(e -> e["i"] == 2, PseudoPotentialIO.PSP8_EXCHANGE_CORRELATION)
        @test PseudoPotentialIO.PSP8_EXCHANGE_CORRELATION[i]["libxc"] == "LDA_X LDA_C_PZ"
        @test PseudoPotentialIO.libxc_to_abinit("LDA_X LDA_C_PZ") == 2
        @test PseudoPotentialIO.libxc_to_abinit("LDA_X LDA_C_PZ_MOD") == -1010
    end
end
