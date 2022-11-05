@testset "PSP v8" begin
    @testset "Internal data consistency" begin
        for (root, dirs, files) in walkdir("./psp8/"), file in files
            psp = Psp8PsP(joinpath(root, file))

            @test format(psp) == "PSP v8"

            @test psp.header.zatom >= psp.header.zion
            @test psp.header.pspcod == 8

            @test length(psp.projectors) == psp.header.lmax + 1
            if psp.header.lloc < psp.header.lmax
                @test isempty(psp.projectors[lloc + 1])
                @test isempty(psp.ekb[lloc + 1])
            end
            for projectors_l in psp.projectors, projector in projectors_l
                @test length(projector) == psp.header.mmax
            end
            for (l_plus_one, nproj_l) in enumerate(psp.header.nproj)
                @test length(psp.projectors[l_plus_one]) == nproj_l
                @test length(psp.ekb[l_plus_one]) == nproj_l
            end

            @test length(psp.v_local) == psp.header.mmax

            core_densities = [psp.rhoc, psp.d_rhoc_dr, psp.d2_rhoc_dr2, psp.d3_rhoc_dr3,
                              psp.d4_rhoc_dr4]
            if psp.header.fchrg > 0
                @test !any(isnothing, core_densities)
                for density in core_densities
                    @test length(density) == psp.header.mmax
                end
            else
                @test all(isnothing, core_densities)
            end

            if psp.header.extension_switch in (2, 3)
                @test length(psp.projectors_so) == psp.header.lmax + 1
                @test isempty(psp.projectors_so[0 + 1])
                for projectors_l in psp.projectors_so, projector in projectors_l
                    @test length(projector) == psp.header.mmax
                end
                @test length(psp.ekb_so) == psp.header.lmax + 1
                @test isempty(psp.projectors_so[0 + 1])
                for (l_plus_one, nprojso_l) in enumerate(psp.header.nprojso)
                    @test length(psp.projectors_so[l_plus_one]) == nprojso_l
                    @test length(psp.ekb_so[l_plus_one]) == nprojso_l
                end
            else
                @test isnothing(psp.projectors_so)
                @test isnothing(psp.ekb_so)
            end
        end
    end

    @testset "H_nc-sr-04_pbe_standard.psp8" begin
        psp = Psp8PsP("psp8/H_nc-sr-04_pbe_standard.psp8")

        @test psp.header.zatom == 1.0000
        @test psp.header.zion == 1.0000
        @test psp.header.pspd == 170501
        @test psp.header.pspcod == 8
        @test psp.header.pspxc == 11
        @test psp.header.lmax == 1
        @test psp.header.lloc == 4
        @test psp.header.mmax == 300
        @test psp.header.r2well == 0.0
        @test psp.header.rchrg == 2.99000000
        @test psp.header.fchrg == 0.00000000
        @test psp.header.qchrg == 0.00000000
        @test psp.header.nproj == [2, 1]
        @test psp.header.extension_switch == 1

        @test psp.ekb[0 + 1][1] == -1.7011234716222E+00
        @test psp.ekb[0 + 1][2] == -5.3448391737199E-01
        @test psp.projectors[0 + 1][1][5] == 1.6367055221032E-01
        @test psp.projectors[0 + 1][2][5] == 3.9125025705919E-01

        @test psp.v_local[5] == -3.1387511173374E+00
    end

    @testset "Ti_nc-fr-04_pbesol_standard.psp8" begin
        psp = Psp8PsP("psp8/Ti_nc-fr-04_pbesol_standard.psp8")

        @test psp.header.zatom == 22.0000
        @test psp.header.zion == 12.0000
        @test psp.header.pspd == 180429
        @test psp.header.pspcod == 8
        @test psp.header.pspxc == -116133
        @test psp.header.lmax == 2
        @test psp.header.lloc == 4
        @test psp.header.mmax == 600
        @test psp.header.r2well == 0.0
        @test psp.header.rchrg == 5.99000000
        @test psp.header.fchrg == 2.50000000
        @test psp.header.qchrg == 0.00000000
        @test psp.header.nproj == [2, 3, 2]
        @test psp.header.extension_switch == 3
        @test psp.header.nprojso == [0, 4, 3]

        @test psp.ekb[0 + 1][1] == 1.2539360126722E+00
        @test psp.ekb[0 + 1][2] == 1.1321100966516E+01
        @test psp.projectors[0 + 1][1][5] == -1.7889864901133E-01
        @test psp.projectors[0 + 1][2][5] == 3.5296270964743E-01

        @test psp.v_local[5] == -1.9130317133247E+01

        @test psp.ekb_so[1 + 1][1] == 6.9968670783987E-02
        @test psp.ekb_so[1 + 1][2] == -3.0732828917905E-02
        @test psp.ekb_so[1 + 1][3] == 4.6090926948408E-03
        @test psp.ekb_so[1 + 1][4] == -2.0693922436038E-04
        @test psp.projectors_so[1 + 1][1][5] == 1.7380194666161E-02
        @test psp.projectors_so[1 + 1][2][5] == -6.8327678530662E-02
        @test psp.projectors_so[1 + 1][3][5] == 2.3034649777788E-02
        @test psp.projectors_so[1 + 1][4][5] == -6.2709563805074E-03
    end

    @testset "Zn_nc-sr-04_pbesol_stringent.psp8" begin
        psp = Psp8PsP("psp8/Zn_nc-sr-04_pbesol_stringent.psp8")

        @test psp.header.zatom == 30.0000
        @test psp.header.zion == 20.0000
        @test psp.header.pspd == 180112
        @test psp.header.pspcod == 8
        @test psp.header.pspxc == -116133
        @test psp.header.lmax == 2
        @test psp.header.lloc == 4
        @test psp.header.mmax == 600
        @test psp.header.r2well == 0.0
        @test psp.header.rchrg == 5.99000000
        @test psp.header.fchrg == 2.50000000
        @test psp.header.qchrg == 0.00000000
        @test psp.header.nproj == [2, 2, 2]
        @test psp.header.extension_switch == 1

        @test psp.ekb[0 + 1][1] == 2.8308354242488E+00
        @test psp.ekb[0 + 1][2] == -1.5191074296972E+00
        @test psp.projectors[0 + 1][1][5] == 3.5070387789440E-01
        @test psp.projectors[0 + 1][2][5] == 1.2235313804112E-01

        @test psp.v_local[5] == -3.1208323049679E+01

        @test psp.rhoc[7] == 1.2735292737007E+02
        @test psp.d_rhoc_dr[7] == -1.3288155269898E+02
        @test psp.d2_rhoc_dr2[7] == -2.0848786910949E+03
        @test psp.d3_rhoc_dr3[7] == 6.3716538707879E+03
        @test psp.d4_rhoc_dr4[7] == 9.6377063748116E+04
    end
end