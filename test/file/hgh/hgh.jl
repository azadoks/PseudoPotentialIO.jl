@testset "HGH" begin
    @testset "c-q4.hgh" begin
        file = load_hgh("c-q4.hgh")

        @test occursin("c", lowercase(file.title))
        @test occursin("pade", lowercase(file.title))
        @test file.zion == [2, 2]
        @test file.rloc == 0.34883045
        @test file.cloc == [-8.51377110, 1.22843203]
        @test file.lmax == 1
        @test file.rp == [0.30455321, 0.2326773]
        @test file.h[1] == 9.52284179 * ones(1, 1)
        @test file.h[2] == zeros(0, 0)
    end

    @testset "ni-q18.hgh" begin
        file = load_hgh("ni-q18.hgh")

        @test occursin("ni", lowercase(file.title))
        @test occursin("pade", lowercase(file.title))
        @test file.zion == [4, 6, 8]
        @test file.rloc == 0.35000000
        @test file.cloc == [3.61031072, 0.44963832]
        @test file.lmax == 2
        @test file.rp == [0.24510489, 0.23474136, 0.21494950]
        @test file.h[1] == [[12.16113071, 3.51625420] [3.51625420, -9.07892931]]
        @test file.h[2] == [[-0.82062357, 2.54774737] [2.54774737, -6.02907069]]
        @test file.h[3] == -13.39506212 * ones(1, 1)
        @test isnothing(file.abinit) && all(isempty, file.k)
        @test !has_spin_orbit(file)
    end

    # ABINIT, pspcod = 3: the off-diagonal coefficients from the HGH relations
    @testset "gd-abinit.hgh" begin
        file = load_hgh("gd-abinit.hgh")
        @test file.abinit.pspcod == 3
        @test element(file).symbol == "Gd"
        @test valence_charge(file) == 18 && isempty(file.zion)
        @test functional(file) == [:lda_xc_teter93]  # pspxc = 1
        @test has_spin_orbit(file)
        @test file.lmax == 3
        @test file.rloc == 0.5175 && file.cloc == [17.512556, -0.719534, 0.0, 0.0]
        @test file.rp == [0.462014, 0.456953, 0.482368, 0.27339]
        @test size.(file.h) == [(2, 2), (2, 2), (1, 1), (1, 1)]  # up to the last nonzero
        @test size.(file.k) == [(0, 0), (3, 3), (1, 1), (1, 1)]
        h22, h33 = -1.697711, 0.0
        @test file.h[2] ≈ [-0.058347 -sqrt(5 / 7) / 2 * h22; -sqrt(5 / 7) / 2 * h22 h22]
        k = file.k[2]
        @test [k[1, 2], k[1, 3], k[2, 3]] ≈
              [-sqrt(5 / 7) / 2 * 0.022024, sqrt(35 / 11) / 6 * -0.054301,
               -14 / sqrt(11) / 6 * -0.054301]
        @test file.h[3] == fill(-0.562601, 1, 1) && file.k[4] == fill(0.014666, 1, 1)
    end

    # ABINIT, pspcod = 10: the coefficients as they are
    @testset "fe-abinit.hgh" begin
        file = load_hgh("fe-abinit.hgh")
        @test file.abinit.pspcod == 10
        @test element(file).symbol == "Fe"
        @test valence_charge(file) == 16
        @test functional(file) == [:gga_x_pbe, :gga_c_pbe]  # pspxc = 11
        @test file.cloc == [6.75678916, -0.22883251]
        @test file.h[1] == [0.62950570 7.91313242; 7.91313242 -10.21581002]
        @test file.k[1] == zeros(0, 0)
        @test file.k[2] == [0.09786820 0.08070002; 0.08070002 -0.09548555]
        @test file.k[3] == fill(0.01036288, 1, 1)
    end
end
