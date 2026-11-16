@testset "UPF v1.old" begin
    @testset "Internal data consistency" begin
        for filepath in values(UPF1_CASE_FILEPATHS)
            file = load_psp_file(filepath)

            @test format(file) == "UPF v1.old"

            # UPF v1.old has a different augmentation data format from UPF v2.0.1
            if file.header.is_ultrasoft | file.header.is_paw
                augmentation = file.nonlocal.augmentation

                @test isnothing(augmentation.multipoles)

                if augmentation.nqf > 0
                    @test !isnothing(augmentation.qfcoefs)
                    @test length(augmentation.qfcoefs) ==
                          (file.header.number_of_proj * (file.header.number_of_proj + 1) /
                           2)
                    for qfcoef in augmentation.qfcoefs
                        @test length(qfcoef.qfcoef) ==
                              augmentation.nqf * (2file.header.l_max + 1)
                    end
                else
                    @test isnothing(augmentation.qfcoefs)
                end

                @test !isnothing(augmentation.rinner)
                @test length(augmentation.rinner) == 2file.header.l_max + 1
            end
        end
    end

    @testset "B_pbe_v1.01.uspp.F.UPF" begin
        filename = "B_pbe_v1.01.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "B"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test header.core_correction
        @test header.functional == "SLA PW PBE PBE"  # "SLA PW PBE PBE PBE"
        @test header.z_valence == 3.00000000000
        @test header.total_psenergy == -6.00744913082
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 1
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 781
        @test header.number_of_wfc == 2
        @test header.number_of_proj == 4

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 1.25728654505E-05
        @test mesh.rab[5] == 3.24915430936E-06
        @test mesh.mesh == file.header.mesh_size
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test length(file.nlcc) == file.header.mesh_size
        @test file.nlcc[5] == 1.20782515483E+00

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -1.00137766482E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 559
        @test betas[1].beta[5] == -3.82483956865E-05

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 559
        @test betas[2].beta[5] == -8.67580981562E-05

        @test betas[3].angular_momentum == 1
        @test betas[3].cutoff_radius_index == 559
        @test betas[3].beta[5] == 1.99207032008E-07

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 559
        @test betas[4].beta[5] == 1.69150912265E-07

        @test size(nonlocal.dij) == (4, 4)
        @test nonlocal.dij[1, 1] == -1.22873815109E-01
        @test nonlocal.dij[1, 2] == -5.61912458144E+00
        @test nonlocal.dij[2, 2] == -5.87789890874E+00
        @test nonlocal.dij[3, 3] == 9.71393906034E+00
        @test nonlocal.dij[3, 4] == 2.46444027912E+01
        @test nonlocal.dij[4, 4] == 6.13899360806E+01
        @test count(d -> !iszero(d), nonlocal.dij) == 8

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 8

        @test size(augmentation.q) == (4, 4)
        @test augmentation.q[1, 1] == -5.06674546795E-01
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == -3.23289503777E-01
        @test augmentation.q[1, 3] == augmentation.q[3, 1] == 0.00000000000E+00
        @test augmentation.q[1, 4] == augmentation.q[4, 1] == 0.00000000000E+00
        @test augmentation.q[2, 2] == 1.89884871724E-01
        @test augmentation.q[2, 3] == augmentation.q[2, 3] == 0.00000000000E+00
        @test augmentation.q[2, 4] == augmentation.q[4, 2] == 0.00000000000E+00
        @test augmentation.q[3, 3] == 9.88682989485E-01
        @test augmentation.q[4, 3] == augmentation.q[3, 4] == 1.83449812909E+00
        @test augmentation.q[4, 4] == 3.27504485114E+00

        @test length(augmentation.rinner) == 3
        @test all(r -> r == 1.10000000000E+00, augmentation.rinner)

        @test length(augmentation.qijs) == 10
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test length(augmentation.qfcoefs) == 10

        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == -1.94859515081E-09
        @test augmentation.qfcoefs[1].qfcoef[5] == -7.21237350451E+01

        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == -1.53409489684E-09
        @test augmentation.qfcoefs[2].qfcoef[5] == -5.60931726792E+01

        @test augmentation.qijs[10].first_index == 4
        @test augmentation.qijs[10].second_index == 4
        @test augmentation.qijs[10].qij[5] == 1.57362794452E-08
        @test augmentation.qfcoefs[10].qfcoef[5] == 7.28024914160E+02

        @test isnothing(augmentation.qijls)
        @test augmentation.nqlc == 2file.header.l_max + 1
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "2S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 5.96294641979E-06

        @test pswfc[2].label == "2P"
        @test pswfc[2].l == 1
        @test pswfc[2].occupation == 1.00
        @test pswfc[2].chi[5] == 1.64367609052E-10

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 1.04916913718E-10

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end

    @testset "mg_pbe_v1.4.uspp.F.UPF" begin
        filename = "mg_pbe_v1.4.uspp.F.UPF"
        file = load_psp_file(UPF1_CASE_FILEPATHS[filename])

        header = file.header
        @test isnothing(header.generated)
        @test isnothing(header.author)
        @test isnothing(header.date)
        @test isnothing(header.comment)
        @test header.element == "Mg"
        @test header.pseudo_type == "US"
        @test header.relativistic == "scalar"
        @test header.is_ultrasoft
        @test !header.is_paw
        @test !header.is_coulomb
        @test !header.has_so
        @test !header.has_wfc
        @test !header.has_gipaw
        @test isnothing(header.paw_as_gipaw)
        @test !header.core_correction
        @test header.functional == "SLA PW PBX PBC"  # "SLA PW PBX PBC PBE"
        @test header.z_valence == 10.00000000000
        @test header.total_psenergy == -125.08981684000
        @test header.wfc_cutoff == 0.00000
        @test header.rho_cutoff == 0.00000
        @test header.l_max == 2
        @test isnothing(header.l_max_rho)
        @test isnothing(header.l_local)
        @test header.mesh_size == 821
        @test header.number_of_wfc == 3
        @test header.number_of_proj == 7

        mesh = file.mesh
        @test length(mesh.r) == file.header.mesh_size
        @test length(mesh.rab) == file.header.mesh_size
        @test mesh.r[5] == 1.98100765738E-08
        @test mesh.rab[5] == 5.20427151222E-09
        @test mesh.mesh == file.header.mesh_size
        @test isnothing(mesh.dx)
        @test isnothing(mesh.xmin)
        @test isnothing(mesh.rmax)
        @test isnothing(mesh.zmesh)

        @test isnothing(file.nlcc)

        @test length(file.local_) == file.header.mesh_size
        @test file.local_[5] == -3.96180531727E+01

        nonlocal = file.nonlocal

        betas = file.nonlocal.betas
        @test length(betas) == file.header.number_of_proj
        for (i, beta) in enumerate(betas)
            @test beta.index == i
            @test isnothing(beta.cutoff_radius)
            @test (length(beta.beta) == beta.cutoff_radius_index) | (length(beta.beta) == file.header.mesh_size)
            @test isnothing(beta.norm_conserving_radius)
            @test isnothing(beta.ultrasoft_cutoff_radius)
            @test isnothing(beta.label)
        end

        @test betas[1].angular_momentum == 0
        @test betas[1].cutoff_radius_index == 661
        @test betas[1].beta[5] == -8.06372437750E-05

        @test betas[2].angular_momentum == 0
        @test betas[2].cutoff_radius_index == 661
        @test betas[2].beta[5] == 1.45732756410E-05

        @test betas[3].angular_momentum == 0
        @test betas[3].cutoff_radius_index == 661
        @test betas[3].beta[5] == 6.78200641134E-04

        @test betas[4].angular_momentum == 1
        @test betas[4].cutoff_radius_index == 661
        @test betas[4].beta[5] == -4.68700375660E-09

        @test betas[5].angular_momentum == 1
        @test betas[5].cutoff_radius_index == 661
        @test betas[5].beta[5] == 1.69113364129E-07

        @test betas[6].angular_momentum == 2
        @test betas[6].cutoff_radius_index == 661
        @test betas[6].beta[5] == -6.58196838008E-14

        @test betas[7].angular_momentum == 2
        @test betas[7].cutoff_radius_index == 661
        @test betas[7].beta[5] == -3.36344971742E-13

        @test size(nonlocal.dij) == (7, 7)
        @test nonlocal.dij[1, 1] == 2.22757816326E+01
        @test nonlocal.dij[1, 2] == -5.12800750666E+00
        @test nonlocal.dij[1, 3] == 1.22489807665E+01
        @test nonlocal.dij[2, 2] == -5.11535718643E-01
        @test nonlocal.dij[2, 3] == -4.01674388929E+00
        @test nonlocal.dij[3, 3] == -1.72050231320E+00
        @test nonlocal.dij[4, 4] == 2.31491498927E+01
        @test nonlocal.dij[4, 5] == -3.23205785707E+01
        @test nonlocal.dij[5, 5] == 4.32591292604E+01
        @test nonlocal.dij[6, 6] == 5.06898497757E+00
        @test nonlocal.dij[6, 7] == -3.89652282790E+00
        @test nonlocal.dij[7, 7] == 3.00051462460E+00
        @test count(d -> !iszero(d), nonlocal.dij) == 17

        augmentation = nonlocal.augmentation
        @test isnothing(augmentation.multipoles)
        @test !augmentation.q_with_l
        @test augmentation.nqf == 8

        @test size(augmentation.q) == (7, 7)
        @test augmentation.q[1, 1] == 1.79874106443E+00
        @test augmentation.q[1, 2] == augmentation.q[2, 1] == -1.29615056877E+00
        @test augmentation.q[7, 7] == 9.57081971594E-02

        @test length(augmentation.rinner) == 5
        @test all(r -> r == 1.20000000000E+00, augmentation.rinner)

        @test length(augmentation.qijs) == 28
        for qij in augmentation.qijs
            @test length(qij.qij) == file.header.mesh_size
        end
        @test length(augmentation.qfcoefs) == 28

        @test augmentation.qijs[1].first_index == 1
        @test augmentation.qijs[1].second_index == 1
        @test augmentation.qijs[1].qij[5] == 2.32704360102E-14
        @test augmentation.qfcoefs[1].qfcoef[5] == 3.43676888927E+02

        @test augmentation.qijs[2].first_index == 1
        @test augmentation.qijs[2].second_index == 2
        @test augmentation.qijs[2].qij[5] == -1.08277556385E-14
        @test augmentation.qfcoefs[2].qfcoef[5] == -1.10516062036E+02

        @test augmentation.qijs[28].first_index == 7
        @test augmentation.qijs[28].second_index == 7
        @test augmentation.qijs[28].qij[5] == -4.33337772184E-16
        @test augmentation.qfcoefs[28].qfcoef[5] == -5.52962873366E+01

        @test isnothing(augmentation.qijls)
        @test augmentation.nqlc == 2file.header.l_max + 1
        @test isnothing(augmentation.shape)
        @test isnothing(augmentation.iraug)
        @test isnothing(augmentation.raug)
        @test isnothing(augmentation.l_max_aug)
        @test isnothing(augmentation.augmentation_epsilon)
        @test isnothing(augmentation.cutoff_r)
        @test isnothing(augmentation.cutoff_r_index)

        pswfc = file.pswfc
        @test length(pswfc) == file.header.number_of_wfc
        for (i, chi) in enumerate(file.pswfc)
            @test chi.index == i
            @test isnothing(chi.n)
            @test isnothing(chi.pseudo_energy)
            @test isnothing(chi.cutoff_radius)
            @test isnothing(chi.ultrasoft_cutoff_radius)
        end

        @test pswfc[1].label == "2S"
        @test pswfc[1].l == 0
        @test pswfc[1].occupation == 2.00
        @test pswfc[1].chi[5] == 4.74019527544E-08

        @test pswfc[3].label == "3S"
        @test pswfc[3].l == 0
        @test pswfc[3].occupation == 1.70
        @test pswfc[3].chi[5] == 1.02365544333E-08

        @test isnothing(file.full_wfc)

        @test length(file.rhoatom) == file.header.mesh_size
        @test file.rhoatom[5] == 3.70544359504E-14

        @test isnothing(file.spin_orb)
        @test isnothing(file.paw)
        @test isnothing(file.gipaw)
    end
end
