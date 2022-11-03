# @testset "Load all UPF v2.0.1" begin
#     for (root, dirs, files) in walkdir("./upf2/"), file in files
#         psp = UpfPsP(joinpath(root, file))
#         @test isa(psp, UpfPsP)
#         @test format(psp) == "UPF v2.0.1"
#     end
# end

@testset "Mg_nc-fr-04_pbesol_stringent.upf" begin
    psp = UpfPsP("./upf2/Mg_nc-fr-04_pbesol_stringent.upf")
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
    @test length(psp.spinorb.relbetas) == 6
    @test psp.spinorb.relbetas[3].lll == 1
    @test psp.spinorb.relbetas[3].jjj == 0.5
    for i = eachindex(psp.spinorb.relbetas)
        @test psp.spinorb.relbetas[i].index == i
    end
    @test length(psp.spinorb.relwfcs) == 4
    @test psp.spinorb.relwfcs[3].lchi == 1
    @test psp.spinorb.relwfcs[3].jchi == 0.5
    @test psp.spinorb.relwfcs[3].nn == 2
    for i = eachindex(psp.spinorb.relwfcs)
        @test psp.spinorb.relwfcs[i].index == i
    end
end

@testset "Si.pbe-n-rrkjus_psl.1.0.0.upf" begin
    psp = UpfPsP("./upf2/Si.pbe-n-rrkjus_psl.1.0.0.upf")
    @test z_valence(psp) == 8.0
    @test element(psp).symbol == "Al"
    @test is_ultrasoft(psp)
    @test l_max(psp) == 2
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
    @test length(psp.nonlocal.augmentation.q) == 36
    @test psp.nonlocal.augmentation.q[2] == -6.659970266131190E-002
    @test isnothing(psp.nonlocal.qijs)
    @test length(psp.nonlocal.qijls) == 68
    @test psp.nonlocal.qijls[1].qijl[5] == 7.047762479021314E-010
    @test isnothing(psp.nonlocal.augmentation.rinner)
    @test psp.nonlocal.augmentation.qfcoeff[1][5] == -3.73585933729E+00
    for i = eachindex(psp.pswfc)
        @test psp.pswfc[i].index == i
    end
    @test psp.pswfc[1].label == "4S"
    @test psp.pswfc[1].l == 0
    @test psp.pswfc[1].occupation == 2.0
    @test psp.pswfc[1].chi[5] == 1.53041946142E-08
end

# @testset "Al.pbe-n-kjpaw_psl.1.0.0.upf" begin
#     psp = UpfPsP("./upf2/Al.pbe-n-kjpaw_psl.1.0.0.upf")
#     @test z_valence(psp) == 19.0
#     @test element(psp).symbol == "Al"
#     @test is_ultrasoft(psp)
#     @test l_max(psp) == 2
#     @test psp.header.mesh_size == 887
#     @test n_proj_radial(psp) == 6
#     @test n_pseudo_wfc(psp) == 5
#     @test psp.mesh.r[5] == 5.05789189118E-09
#     @test psp.mesh.rab[5] == 1.32875017333E-09
#     @test psp.nlcc[5] == 4.66361704568E-01
#     @test psp.nonlocal.betas[1].angular_momentum == 0
#     @test psp.nonlocal.betas[1].cutoff_radius_index == 717
#     @test psp.nonlocal.betas[1].beta[5] == -5.63186386318E-05
#     for i = eachindex(psp.nonlocal.betas)
#         @test psp.nonlocal.betas[i].index == i
#     end
#     @test e_kb(psp, 0, 1) == -1.76074960113E+00
#     for i = eachindex(psp.pswfc)
#         @test psp.pswfc[i].index == i
#     end
#     @test length(psp.nonlocal.augmentation.rinner) == 5
#     @test psp.nonlocal.augmentation.rinner[1] == 1.30000000000E+00
#     @test psp.nonlocal.augmentation.qijs[1].qij[5] == -1.55354181933E-16
#     @test length(psp.nonlocal.augmentation.qijs) == 21
#     @test psp.nonlocal.augmentation.qfcoeff[1][5] == -3.73585933729E+00
#     @test psp.pswfc[1].label == "4S"
#     @test psp.pswfc[1].l == 0
#     @test psp.pswfc[1].occupation == 2.0
#     @test psp.pswfc[1].chi[5] == 1.53041946142E-08
# end

# @testset "he-q2.upf" begin
#     psp = UpfPsP("./upf2/he-q2.upf")
#     @test z_valence(psp) == 19.0
#     @test element(psp).symbol == "Al"
#     @test is_ultrasoft(psp)
#     @test l_max(psp) == 2
#     @test psp.header.mesh_size == 887
#     @test n_proj_radial(psp) == 6
#     @test n_pseudo_wfc(psp) == 5
#     @test psp.mesh.r[5] == 5.05789189118E-09
#     @test psp.mesh.rab[5] == 1.32875017333E-09
#     @test psp.nlcc[5] == 4.66361704568E-01
#     @test psp.nonlocal.betas[1].angular_momentum == 0
#     @test psp.nonlocal.betas[1].cutoff_radius_index == 717
#     @test psp.nonlocal.betas[1].beta[5] == -5.63186386318E-05
#     for i = eachindex(psp.nonlocal.betas)
#         @test psp.nonlocal.betas[i].index == i
#     end
#     @test e_kb(psp, 0, 1) == -1.76074960113E+00
#     for i = eachindex(psp.pswfc)
#         @test psp.pswfc[i].index == i
#     end
#     @test length(psp.nonlocal.augmentation.rinner) == 5
#     @test psp.nonlocal.augmentation.rinner[1] == 1.30000000000E+00
#     @test psp.nonlocal.augmentation.qijs[1].qij[5] == -1.55354181933E-16
#     @test length(psp.nonlocal.augmentation.qijs) == 21
#     @test psp.nonlocal.augmentation.qfcoeff[1][5] == -3.73585933729E+00
#     @test psp.pswfc[1].label == "4S"
#     @test psp.pswfc[1].l == 0
#     @test psp.pswfc[1].occupation == 2.0
#     @test psp.pswfc[1].chi[5] == 1.53041946142E-08
# end
