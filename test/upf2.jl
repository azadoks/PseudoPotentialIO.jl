# @testset "Load all UPF v2.0.1" begin
#     for (root, dirs, files) in walkdir("./upf2/"), file in files
#         psp = UpfPsP(joinpath(root, file))
#         @test isa(psp, UpfPsP)
#         @test format(psp) == "UPF v2.0.1"
#     end
# end

@testset "[UPF v2.0.1] Mg_nc-fr-04_pbesol_stringent.upf" begin
    psp = UpfPsP("upf2/Mg_nc-fr-04_pbesol_stringent.upf")
    @test z_valence(psp) == 10.0
    @test element(psp).symbol == "Mg"
    @test is_norm_conserving(psp)
    @test l_max(psp) == 1
    @test psp.header.mesh_size == 1510
    @test n_proj_radial(psp) == 6
    @test n_pseudo_wfc(psp) == 4
    @test has_spin_orbit(psp)
    @test psp.mesh.r[5] == 0.0400
    @test psp.mesh.rab[5] == 0.0100
    @test isnothing(psp.nlcc)
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 160
    @test psp.nonlocal.betas[1].cutoff_radius == 1.5900000000E+00
    @test psp.nonlocal.betas[1].beta[5] == 1.2872222840E-01
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == -1.6525827200E+00
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test isnothing(psp.nonlocal.augmentation)
    @test psp.pswfc[1].label == "2S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 1.3903075777E-01
    @test length(psp.spin_orb.relbetas) == 6
    @test psp.spin_orb.relbetas[3].lll == 1
    @test psp.spin_orb.relbetas[3].jjj == 0.5
    for i = eachindex(psp.spin_orb.relbetas)
        @test psp.spin_orb.relbetas[i].index == i
    end
    @test length(psp.spin_orb.relwfcs) == 4
    @test psp.spin_orb.relwfcs[3].lchi == 1
    @test psp.spin_orb.relwfcs[3].jchi == 0.5
    @test psp.spin_orb.relwfcs[3].nn == 2
    for i = eachindex(psp.spin_orb.relwfcs)
        @test psp.spin_orb.relwfcs[i].index == i
    end
end

@testset "[UPF v2.0.1] Si.pbe-n-rrkjus_psl.1.0.0.upf" begin
    psp = UpfPsP("upf2/Si.pbe-n-rrkjus_psl.1.0.0.upf")
    @test z_valence(psp) == 4.0
    @test element(psp).symbol == "Si"
    @test is_ultrasoft(psp)
    @test l_max(psp) == 2
    @test has_nlcc(psp)
    @test psp.header.l_max_rho == 4
    @test psp.header.mesh_size == 1141
    @test n_proj_radial(psp) == 6
    @test n_pseudo_wfc(psp) == 2
    @test psp.mesh.r[5] == 6.847393954957285E-005
    @test psp.mesh.rab[5] == 8.559242443696606E-007
    @test psp.nlcc[5] == 1.538616848020000E+000
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].label == "3S"
    @test psp.nonlocal.betas[1].cutoff_radius_index == 829
    @test psp.nonlocal.betas[1].ultrasoft_cutoff_radius == 1.800000000000000E+000
    @test psp.nonlocal.betas[1].beta[5] == -1.233266788373589E-003
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == 4.959233384252660E-001
    @test psp.nonlocal.augmentation.q_with_l
    @test length(psp.nonlocal.augmentation.q) == 36
    @test psp.nonlocal.augmentation.q[2] == -6.659970266131190E-002
    @test isnothing(psp.nonlocal.augmentation.qijs)
    @test length(psp.nonlocal.augmentation.qijls) == 34
    @test psp.nonlocal.augmentation.qijls[1].qijl[5] == 7.047762479021314E-010
    @test isnothing(psp.nonlocal.augmentation.rinner)
    @test isnothing(psp.nonlocal.augmentation.qfcoeff)
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test psp.pswfc[1].label == "3S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 3.542904978471305E-005
    @test psp.rhoatom[5] == 6.986543942450972E-009 
end

@testset "[UPF v2.0.1] Al.pbe-n-kjpaw_psl.1.0.0.upf" begin
    psp = UpfPsP("upf2/Al.pbe-n-kjpaw_psl.1.0.0.upf")
    @test z_valence(psp) == 3.0
    @test element(psp).symbol == "Al"
    @test is_ultrasoft(psp)
    @test is_paw(psp)
    @test psp.header.has_gipaw
    @test psp.header.paw_as_gipaw
    @test l_max(psp) == 2
    @test psp.header.l_max_rho == 4
    @test psp.header.mesh_size == 1135
    @test n_proj_radial(psp) == 6
    @test n_pseudo_wfc(psp) == 2
    @test psp.mesh.r[5] == 7.374116566877076E-005
    @test psp.mesh.rab[5] == 9.217645708596345E-007
    @test psp.nlcc[5] ==  2.656001002192721E-001
    @test psp.local_[5] == -6.296496837141098E+000
    @test psp.nonlocal.betas[1].label == "3S"
    @test psp.nonlocal.betas[1].angular_momentum == 0
    @test psp.nonlocal.betas[1].cutoff_radius_index == 827
    @test psp.nonlocal.betas[1].beta[5] == -1.792962988082909E-003
    for i = eachindex(psp.nonlocal.betas)
        @test psp.nonlocal.betas[i].index == i
    end
    @test e_kb(psp, 0, 1) == 4.162665025130763E-001
    @test psp.nonlocal.augmentation.q_with_l
    @test psp.nonlocal.augmentation.nqf == 0
    @test psp.nonlocal.augmentation.nqlc == 5
    @test psp.nonlocal.augmentation.shape == "PSQ"
    @test psp.nonlocal.augmentation.cutoff_r_index == 839
    @test psp.nonlocal.augmentation.l_max_aug == 4
    @test psp.nonlocal.augmentation.q[1] == -3.070968103524502E-002
    @test psp.nonlocal.augmentation.multipoles[7] == -3.491452538270029E-002
    @test isnothing(psp.nonlocal.augmentation.rinner)
    @test isnothing(psp.nonlocal.augmentation.qfcoeff)
    @test isnothing(psp.nonlocal.augmentation.qijs)
    @test length(psp.nonlocal.augmentation.qijls) == 34
    @test psp.nonlocal.augmentation.qijls[1].qijl[5] == 1.301443338375301E-009
    @test psp.pswfc[1].label == "3S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 2.099526933006344E-005
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
end

@testset "he-q2.upf" begin
    psp = UpfPsP("upf2/he-q2.upf")
    @test z_valence(psp) == 2.0
    @test element(psp).symbol == "He"
    @test l_max(psp) == -1
    @test psp.header.mesh_size == 985
    @test n_proj_radial(psp) == 0
    @test n_pseudo_wfc(psp) == 1
    @test psp.mesh.r[5] == 4.793175768470099E-004
    @test psp.mesh.rab[5] == 5.991469710587625E-006
    @test !has_nlcc(psp)
    @test isempty(psp.nonlocal.betas)
    @test isempty(psp.nonlocal.dij)
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test psp.pswfc[1].label == "1S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 1.831401620367775E-003
    @test psp.rhoatom[5] == 6.708063790171425E-006
end
