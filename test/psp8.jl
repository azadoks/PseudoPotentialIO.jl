@testset "[PSP v8] Loading and consistency" begin
    for (root, dirs, files) in walkdir("./psp8/"), file in files
        psp = Psp8PsP(joinpath(root, file))

        @test isa(psp, Psp8PsP)
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
            @test isempty(psp.projectors_so[0+1])
            for projectors_l in psp.projectors_so, projector in projectors_l
                @test length(projector) == psp.header.mmax
            end
            @test length(psp.ekb_so) == psp.header.lmax + 1
            @test isempty(psp.projectors_so[0+1])
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
