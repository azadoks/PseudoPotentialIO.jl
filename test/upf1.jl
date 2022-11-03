@testset "[UPF v1.old] Loading and consistency" begin
    for (root, dirs, files) in walkdir("./upf1/"), file in files
        psp = UpfPsP(joinpath(root, file))

        @test isa(psp, UpfPsP)
        @test format(psp) == "UPF v1.old"

        @test haskey(PeriodicTable.elements, Symbol(psp.header.element))
        z_atom = PeriodicTable.elements[Symbol(psp.header.element)].number
        @test psp.header.z_valence <= z_atom
        
        @test psp.header.is_ultrasoft + psp.header.is_coulomb <= 1
        @test psp.header.is_paw + psp.header.is_coulomb <= 1
        @test psp.header.is_paw + psp.header.is_ultrasoft <= 1
        if psp.header.is_ultrasoft | psp.header.is_paw
            @test !isnothing(psp.nonlocal.augmentation)
            augmentation = psp.nonlocal.augmentation
            @test isnothing(augmentation.q)
            @test isnothing(augmentation.multipoles)
            @test length(augmentation.rinner) == 2psp.header.l_max + 1
            if augmentation.q_with_l
                @test isnothing(augmentation.qijs)
                @test !isnothing(augmentation.qijls)
                @test (length(augmentation.qijls) == psp.header.number_of_proj
                       * augmentation.nqlc)
            else
                @test !isnothing(augmentation.qijs)
                @test isnothing(augmentation.qijls)
                n_proj = psp.header.number_of_proj
                @test length(augmentation.qijs) == n_proj * (n_proj + 1) / 2
                @test length(augmentation.qfcoeff) == n_proj * (n_proj + 1) / 2
                @test all(length.(augmentation.qfcoeff) .== 
                          augmentation.nqf * (2psp.header.l_max + 1))
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
            @test len(psp.spin_orb.relbetas) == psp.header.number_of_proj
            @test len(psp.spin_orb.relwfcs) == psp.header.number_of_wfc
        end

        @test length(psp.nonlocal.betas) == psp.header.number_of_proj
        if psp.header.number_of_wfc > 0
            @test length(psp.pswfc) == psp.header.number_of_wfc
        else
            @test isnothing(psp.pswfc)
        end

        @test (length(psp.mesh.r) == length(psp.mesh.rab) == psp.mesh.mesh
               == psp.header.mesh_size)

        for beta in psp.nonlocal.betas
            @test length(beta.beta) == beta.cutoff_radius_index
        end
    end
end

@testset "[UPF v1.old] ag_lda_v1.4.uspp.F.upf" begin
    psp = UpfPsP("upf1/ag_lda_v1.4.uspp.F.upf")
    @test z_valence(psp) == 19.0
    @test element(psp).symbol == "Ag"
    @test is_ultrasoft(psp)
    @test l_max(psp) == 2
    @test psp.header.mesh_size == 887
    @test n_proj_radial(psp) == 6
    @test n_pseudo_wfc(psp) == 5
    @test psp.mesh.r[5] == 5.05789189118E-09
    @test psp.mesh.rab[5] == 1.32875017333E-09
    @test psp.nlcc[5] == 4.66361704568E-01
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 717
    @test psp.nonlocal.betas[1].beta[5] == -5.63186386318E-05
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == -1.76074960113E+00
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test length(psp.nonlocal.augmentation.rinner) == 5
    @test psp.nonlocal.augmentation.rinner[1] == 1.30000000000E+00
    @test psp.nonlocal.augmentation.qijs[1].qij[5] == -1.55354181933E-16
    @test length(psp.nonlocal.augmentation.qijs) == 21
    @test psp.nonlocal.augmentation.qfcoeff[1][5] == -3.73585933729E+00
    @test psp.pswfc[1].label == "4S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 1.53041946142E-08
end

@testset "[UPF v1.old] B_pbe_v1.01.uspp.F.upf" begin
    psp = UpfPsP("upf1/B_pbe_v1.01.uspp.F.upf")
    @test z_valence(psp) == 3.00
    @test element(psp).symbol == "B"
    @test is_ultrasoft(psp)
    @test has_nlcc(psp)
    @test l_max(psp) == 1
    @test psp.header.mesh_size == 781
    @test n_proj_radial(psp) == 4
    @test n_pseudo_wfc(psp) == 2
    @test psp.mesh.r[5] == 1.25728654505E-05
    @test psp.mesh.rab[5] == 3.24915430936E-06
    @test psp.nlcc[5] == 1.20782515483E+00
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 559
    @test psp.nonlocal.betas[1].beta[5] == -3.82483956865E-05
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == -1.22873815109E-01
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test length(psp.nonlocal.augmentation.rinner) == 3
    @test psp.nonlocal.augmentation.rinner[1] == 1.10000000000E+00
    @test psp.nonlocal.augmentation.qijs[1].qij[5] == -1.94859515081E-09
    @test length(psp.nonlocal.augmentation.qijs) == 10
    @test psp.nonlocal.augmentation.qfcoeff[1][5] == -7.21237350451E+01
    @test psp.pswfc[1].label == "2S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 5.96294641979E-06
end

@testset "[UPF v1.old] si_pbesol_v1.uspp.F.upf" begin
    psp = UpfPsP("upf1/si_pbesol_v1.uspp.F.upf")
    @test z_valence(psp) == 4.0
    @test element(psp).symbol == "Si"
    @test is_ultrasoft(psp)
    @test has_nlcc(psp)
    @test l_max(psp) == 2
    @test psp.header.mesh_size == 899
    @test n_proj_radial(psp) == 6
    @test n_pseudo_wfc(psp) == 2
    @test psp.mesh.r[5] == 4.49030908945E-06
    @test psp.mesh.rab[5] == 1.16041225334E-06
    @test psp.nlcc[5] == 6.79372951207E-01
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 627
    @test psp.nonlocal.betas[1].beta[5] == -8.50729690678E-06
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == 4.17504606973E-01
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test length(psp.nonlocal.augmentation.rinner) == 5
    @test psp.nonlocal.augmentation.rinner[1] == 9.00000000000E-01
    @test psp.nonlocal.augmentation.qijs[1].qij[5] == 2.30032305636E-12
    @test length(psp.nonlocal.augmentation.qijs) == 21
    @test psp.nonlocal.augmentation.qfcoeff[1][5] == -2.28802388683E+02
    @test psp.pswfc[1].label == "3S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 6.44327877084E-07
end

@testset "[UPF v1.old] mg_pbe_v1.4.uspp.F.upf" begin
    psp = UpfPsP("upf1/mg_pbe_v1.4.uspp.F.upf")
    @test z_valence(psp) == 10.0
    @test element(psp).symbol == "Mg"
    @test is_ultrasoft(psp)
    @test !has_nlcc(psp)
    @test l_max(psp) == 2
    @test psp.header.mesh_size == 821
    @test n_proj_radial(psp) == 7
    @test n_pseudo_wfc(psp) == 3
    @test psp.mesh.r[5] == 1.98100765738E-08
    @test psp.mesh.rab[5] == 5.20427151222E-09
    @test isnothing(psp.nlcc)
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 661
    @test psp.nonlocal.betas[1].beta[5] == -8.06372437750E-05
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == 2.22757816326E+01
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test length(psp.nonlocal.augmentation.rinner) == 5
    @test psp.nonlocal.augmentation.rinner[1] == 1.20000000000E+00
    @test psp.nonlocal.augmentation.qijs[1].qij[5] == 2.32704360102E-14
    @test length(psp.nonlocal.augmentation.qijs) == 28
    @test psp.nonlocal.augmentation.qfcoeff[1][5] == 3.43676888927E+02
    @test psp.pswfc[1].label == "2S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 4.74019527544E-08
end
