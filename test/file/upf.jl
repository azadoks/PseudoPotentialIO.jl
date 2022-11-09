@testset "UPF" begin
    dirs = [artifact"gbrv_pbe_1.5_upf", artifact"hgh_lda_upf",
            artifact"sg15_2022.02.06_upf", artifact"sssp_pbe_efficiency_1.1.2_upf"]

    filepaths = []
    for dir in dirs
        (_, _, files) = first(walkdir(dir))
        for file in files
            if file[1] != '.'  # Hack to avoid hidden files
                push!(filepaths, joinpath(dir, file))
            end
        end
    end

    @testset "Internal data consistency" begin
        @testset "$filepath" for filepath in filepaths
            psp = load_psp_file(filepath)

            @test haskey(PeriodicTable.elements, Symbol(psp.header.element))
            z_atom = PeriodicTable.elements[Symbol(psp.header.element)].number
            @test psp.header.z_valence <= z_atom

            if psp.header.pseudo_type == "PAW"
                @test psp.header.is_paw
                @test psp.header.is_ultrasoft
                @test !psp.header.is_coulomb
            elseif psp.header.pseudo_type in ("US", "USPP")
                @test psp.header.is_ultrasoft
                @test !any((psp.header.is_paw, psp.header.is_coulomb))
            elseif psp.header.pseudo_type == "1/r"
                @test psp.header.is_coulomb
                @test !any((psp.header.is_paw, psp.header.is_ultrasoft))
            else
                @test psp.header.pseudo_type in ("NC", "SL")
            end

            if psp.header.is_ultrasoft | psp.header.is_paw
                @test !isnothing(psp.nonlocal.augmentation)
                augmentation = psp.nonlocal.augmentation

                @test size(augmentation.q) == (psp.header.number_of_proj,
                                               psp.header.number_of_proj)

                if augmentation.q_with_l
                    @test isnothing(augmentation.qijs)
                    @test !isnothing(augmentation.qijls)
                    #TODO @test (length(augmentation.qijls) ==
                else
                    @test !isnothing(augmentation.qijs)
                    @test isnothing(augmentation.qijls)
                    n_proj = psp.header.number_of_proj
                    @test (length(augmentation.qijs) ==
                           psp.header.number_of_proj
                           *
                           (psp.header.number_of_proj + 1) / 2)
                end
            else
                @test isnothing(psp.nonlocal.augmentation)
            end
            if psp.header.is_coulomb
                #? What else should be missing in this case?
                @test isnothing(psp.local_)
            else
                @test !isnothing(psp.local_)
                @test length(psp.local_) == psp.header.mesh_size
            end

            if psp.header.core_correction
                @test !isnothing(psp.nlcc)
                @test length(psp.nlcc) == psp.header.mesh_size
            end

            if psp.header.has_so
                @test !isnothing(psp.spin_orb)
                @test length(psp.spin_orb.relbetas) == psp.header.number_of_proj
                @test length(psp.spin_orb.relwfcs) == psp.header.number_of_wfc
            end

            if psp.header.has_gipaw
                @test !isnothing(psp.gipaw)
                #TODO check the contents
            else
                @test isnothing(psp.gipaw)
            end

            if psp.header.number_of_wfc > 0
                @test length(psp.pswfc) == psp.header.number_of_wfc
            else
                @test isnothing(psp.pswfc)
            end

            if psp.header.has_wfc
                @test !isnothing(psp.full_wfc)
                #TODO check contents
            else
                @test isnothing(psp.full_wfc)
            end

            @test (length(psp.mesh.r) == length(psp.mesh.rab) == psp.header.mesh_size)
            if !isnothing(psp.mesh.mesh)
                @test psp.mesh.mesh == psp.header.mesh_size
            end

            @test length(psp.nonlocal.betas) == psp.header.number_of_proj
            for beta in psp.nonlocal.betas
                @test length(beta.beta) == beta.cutoff_radius_index
            end
        end
    end
end