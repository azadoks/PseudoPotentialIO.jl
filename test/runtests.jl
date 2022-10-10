using UPF
using Test
using JSON
using PeriodicTable

Test.eval(quote
	function record(ts::DefaultTestSet, t::Union{Fail, Error})
		push!(ts.results, t)
	end
end)

PSP8_DIR = "./psp8/"
UPF_DIR = "./upf/"
JSON_DIR = "./json/"
HEADER_KEYS = [
    "number_of_proj", "core_correction", "element", "pseudo_type", "z_valence",
    "mesh_size", "number_of_wfc"
]
VECTOR_KEYS = ["radial_grid", "core_charge_density", "local_potential", "total_charge_density"]

pseudo_pairs = []
for (root, dirs, files) in walkdir(JSON_DIR)
    for file in files
        json_pseudo = JSON.parsefile(joinpath(root, file))["pseudo_potential"]
        upf_filename = replace(file, ".json" => ".upf")
        upf_pseudo = load_upf(joinpath(UPF_DIR, upf_filename))
        push!(pseudo_pairs, (json=json_pseudo, upf=upf_pseudo))
    end
end

@testset "UPF" begin
    @testset "mesh" begin
        for (_, upf) in pseudo_pairs
            @test length(upf["radial_grid"]) == upf["header"]["mesh_size"]
            @test length(upf["radial_grid_derivative"]) == upf["header"]["mesh_size"]
            r_params = upf["radial_grid_parameters"]
            if r_params["mesh_type"] in ["log_1", "log_2"]
                nr = r_params["mesh"]
                xmin = r_params["xmin"]
                dx = r_params["dx"]
                z = r_params["zmesh"]
                if r_params["mesh_type"] == "log_1"
                    r = map(i -> exp(xmin) * exp((i - 1) * dx) / z, 1:nr)
                    @test all(upf["radial_grid"] .≈ r)
                else r_params["mesh_type"] == "log_2"
                    r = map(i -> exp(xmin) * (exp((i - 1) * dx) - 1) / z, 1:nr)
                    @test all(upf["radial_grid"] .≈ r)
                end
            elseif r_params["mesh_type"] == "linear"
                @test all(isapprox.(
                    round.(diff(upf["radial_grid"]), digits=4),
                    upf["radial_grid_derivative"][1:end-1],
                ))
            end
        end
    end

    @testset "nlcc" begin
        for (_, upf) in pseudo_pairs
            if upf["header"]["core_correction"]
                @test length(upf["core_charge_density"]) == upf["header"]["mesh_size"]
            else
                @test !haskey(upf, "core_charge_density")
            end
        end
    end

    @testset "local" begin
        for (_, upf) in pseudo_pairs
            @test length(upf["local_potential"]) == upf["header"]["mesh_size"]
        end
    end

    @testset "beta" begin
        for (_, upf) in pseudo_pairs
            @test length(upf["beta_projectors"]) == upf["header"]["number_of_proj"]
            for proj in upf["beta_projectors"]
                @test length(proj["radial_function"]) == proj["cutoff_radius_index"]
                @test proj["index"] <= upf["header"]["number_of_proj"] 
                @test proj["angular_momentum"] <= upf["header"]["l_max"]
                if upf["header"]["has_so"]
                    @test haskey(proj, "total_angular_momentum")
                end
            end
        end
    end

    @testset "Dij" begin
        for (_, upf) in pseudo_pairs
            @test ndims(upf["D_ion"]) == 2
            @test size(upf["D_ion"], 1) == size(upf["D_ion"], 2) == upf["header"]["number_of_proj"]
        end
    end

    @testset "pswfc" begin
        for (_, upf) in pseudo_pairs
            @test length(upf["atomic_wave_functions"]) == upf["header"]["number_of_wfc"]
            @test (
                sum(wfc -> wfc["occupation"], upf["atomic_wave_functions"]) <=
                upf["header"]["z_valence"]
            )
            for wfc in upf["atomic_wave_functions"]
                @test length(wfc["radial_function"]) == upf["header"]["mesh_size"]
                @test wfc["angular_momentum"] <= upf["header"]["l_max"]
                if upf["header"]["has_so"]
                    @test haskey(wfc, "principal_quantum_number")
                    @test haskey(wfc, "total_angular_momentum")
                end
            end
        end
    end

    @testset "rhoatom" begin
        for (_, upf) in pseudo_pairs
            @test length(upf["total_charge_density"]) == upf["header"]["mesh_size"]
        end
    end
end

@testset "UPF--JSON" begin
    @testset "header" begin
        for pair in pseudo_pairs
            for key in HEADER_KEYS
                @test pair.json["header"][key] == pair.upf["header"][key]
            end
        end
    end

    @testset "radial_grid" begin
        for pair in pseudo_pairs
            @test all(pair.json["radial_grid"] .== pair.upf["radial_grid"])
        end
    end

    @testset "nlcc" begin
        for pair in pseudo_pairs
            if haskey(pair.json, "core_charge_density")
                @test all(pair.json["core_charge_density"] .== pair.upf["core_charge_density"])
            end
        end
    end

    @testset "local_potential" begin
        for pair in pseudo_pairs
            @test all(pair.json["local_potential"] .== pair.upf["local_potential"] ./ 2)
        end
    end

    @testset "beta_projectors" begin
        for pair in pseudo_pairs
            json_betas = pair.json["beta_projectors"]
            upf_betas = pair.upf["beta_projectors"]
            @test length(json_betas) == length(upf_betas)
            for i = eachindex(json_betas)
                if haskey(json_betas[i], "label")
                    @test json_betas[i]["label"] == upf_betas[i]["label"]
                end
                @test json_betas[i]["angular_momentum"] == upf_betas[i]["angular_momentum"]
                cutoff_idx = length(json_betas[i]["radial_function"])
                @test all(
                    json_betas[i]["radial_function"][begin:cutoff_idx] .==
                    upf_betas[i]["radial_function"][begin:cutoff_idx]
                )
            end
        end
    end

    # @testset "augmentation" begin
    #     aug_pairs = [pair for pair in pseudo_pairs if haskey(pair.json, "augmentation")]
    #     for pair in aug_pairs
    #         json_aug = pair.json["augmentation"]
    #         upf_aug = pair.upf["augmentation"]
    #         for i = eachindex(json_aug)
    #             @test json_aug[i]["i"] + 1 == upf_aug[i]["i"]
    #             @test json_aug[i]["j"] + 1 == upf_aug[i]["j"]
    #             @test json_aug[i]["angular_momentum"] == upf_aug[i]["angular_momentum"]
    #             @test isapprox(json_aug[i]["radial_function"], upf_aug[i]["radial_function"])
    #             # @test all(isapprox.(json_aug[i]["radial_function"], upf_aug[i]["radial_function"], atol=4e-10))  # rtol=5e-6))
    #         end
    #     end
    # end

    # @testset "paw_data" begin
    #     paw_pairs = [pair for pair in pseudo_pairs if haskey(pair.json, "paw_data")]
    #     @testset "aug_integrals" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["aug_integrals"]) == length(pair.upf["paw_data"]["aug_integrals"])
    #             @test all(pair.json["paw_data"]["aug_integrals"] .== pair.upf["paw_data"]["aug_integrals"])
    #         end
    #     end
    #     @testset "aug_multipoles" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["aug_multipoles"]) == length(pair.upf["paw_data"]["aug_multipoles"])
    #             @test all(pair.json["paw_data"]["aug_multipoles"] .== pair.upf["paw_data"]["aug_multipoles"])
    #         end
    #     end
    #     @testset "ae_wfc" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["ae_wfc"]) == length(pair.upf["paw_data"]["ae_wfc"])
    #             for i = 1:length(pair.json["paw_data"]["ae_wfc"])
    #                 @test pair.json["paw_data"]["ae_wfc"][i]["angular_momentum"] == pair.upf["paw_data"]["ae_wfc"][i]["angular_momentum"]
    #                 cutoff_idx = length(pair.json["paw_data"]["ae_wfc"][i]["radial_function"])
    #                 @test all(
    #                     pair.json["paw_data"]["ae_wfc"][i]["radial_function"][begin:cutoff_idx] .==
    #                     pair.upf["paw_data"]["ae_wfc"][i]["radial_function"][begin:cutoff_idx]
    #                 )
    #             end
    #         end
    #     end
    #     @testset "ps_wfc" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["ps_wfc"]) == length(pair.upf["paw_data"]["ps_wfc"])
    #             for i = 1:length(pair.json["paw_data"]["ps_wfc"])
    #                 @test pair.json["paw_data"]["ps_wfc"][i]["angular_momentum"] == pair.upf["paw_data"]["ps_wfc"][i]["angular_momentum"]
    #                 cutoff_idx = length(pair.json["paw_data"]["ps_wfc"][i]["radial_function"])
    #                 @test all(
    #                     pair.json["paw_data"]["ps_wfc"][i]["radial_function"][begin:cutoff_idx] .==
    #                     pair.upf["paw_data"]["ps_wfc"][i]["radial_function"][begin:cutoff_idx]
    #                 )
    #             end
    #         end
    #     end
    #     @testset "occupations" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["occupations"]) == length(pair.upf["paw_data"]["occupations"])
    #             @test all(pair.json["paw_data"]["occupations"] .== pair.upf["paw_data"]["occupations"])
    #         end
    #     end
    #     @testset "ae_core_charge_density" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["ae_core_charge_density"]) == length(pair.upf["paw_data"]["ae_core_charge_density"])
    #             @test all(pair.json["paw_data"]["ae_core_charge_density"] .== pair.upf["paw_data"]["ae_core_charge_density"])
    #         end
    #     end
    #     @testset "ae_local_potential" begin
    #         for pair in paw_pairs
    #             @test length(pair.json["paw_data"]["ae_local_potential"]) == length(pair.upf["paw_data"]["ae_local_potential"])
    #             @test all(pair.json["paw_data"]["ae_local_potential"] .== pair.upf["paw_data"]["ae_local_potential"] ./ 2)
    #         end
    #     end
    # end

    @testset "atomic_wave_functions" begin
        for pair in pseudo_pairs
            json_wfc = pair.json["atomic_wave_functions"]
            upf_wfc = pair.upf["atomic_wave_functions"]
            @test length(json_wfc) == length(upf_wfc)
            for i = eachindex(json_wfc)
                if haskey(json_wfc[i], "label")
                    @test json_wfc[i]["label"] == upf_wfc[i]["label"]
                end
                @test json_wfc[i]["occupation"] == upf_wfc[i]["occupation"]
                @test json_wfc[i]["angular_momentum"] == upf_wfc[i]["angular_momentum"]
                cutoff_idx = length(json_wfc[i]["radial_function"])
                @test all(
                    json_wfc[i]["radial_function"][begin:cutoff_idx] .==
                    upf_wfc[i]["radial_function"][begin:cutoff_idx]
                )
            end
        end
    end

    @testset "rhoatom" begin
        for pair in pseudo_pairs
            if haskey(pair.json, "total_charge_density")
                @test all(pair.json["total_charge_density"] .== pair.upf["total_charge_density"])
            end
        end
    end
end

pseudo_pairs = []
for (root, dirs, files) in walkdir(PSP8_DIR)
    for file in files
        upf_filename = replace(file, ".psp8" => ".upf")
        upf_path = joinpath(UPF_DIR, upf_filename)
        if isfile(upf_path)
            psp = load_psp(joinpath(root, file))
            upf = load_upf(upf_path)
            push!(pseudo_pairs, (psp=psp, upf=upf))
        end
    end
end

@testset "PSP8--UPF" begin
    @testset "header" begin
        for pair in pseudo_pairs
            psphead = pair.psp["header"]
            upfhead = pair.upf["header"]
            @test psphead["z_atom"] == elements[Symbol(upfhead["element"])].number
            @test psphead["z_valence"] == upfhead["z_valence"]
            @test psphead["l_max"] == upfhead["l_max"]
            @test psphead["core_correction"] == upfhead["core_correction"]
            @test psphead["has_so"] == upfhead["has_so"]
            if !psphead["has_so"]
                @test sum(psphead["number_of_proj"]) == upfhead["number_of_proj"]
            end
        end
    end

    @testset "radial_grid" begin
        for pair in pseudo_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            @test isapprox(
                pair.psp["radial_grid"],
                pair.upf["radial_grid"][begin:mesh_size]
            )
        end
    end

    @testset "local_potential" begin
        for pair in pseudo_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            @test isapprox(
                pair.psp["local_potential"],
                pair.upf["local_potential"][begin:mesh_size] ./ 2,
            )
        end
    end

    @testset "nlcc" begin
        for pair in pseudo_pairs
            if pair.psp["header"]["core_correction"]
                mesh_size = pair.psp["header"]["mesh_size"]
                @test isapprox(
                    pair.psp["nlcc"]["core_charge_density"],
                    pair.upf["core_charge_density"][begin:mesh_size],
                )
            end
        end
    end

    @testset "nonlocal_potential" begin
        non_so_pairs = filter(pair -> !pair.psp["header"]["has_so"], pseudo_pairs)
        # @testset "Vnl $(pair.psp["header"]["filename"])" for pair in non_so_pairs
        for pair in non_so_pairs
            mesh_size = pair.psp["header"]["mesh_size"]
            idx_upf = 1
            for l = 0:pair.psp["header"]["l_max"]
                for i = 1:length(pair.psp["beta_projectors"][l+1])
                    psp_β = pair.psp["beta_projectors"][l+1][i]
                    upf_β = pair.upf["beta_projectors"][idx_upf]["radial_function"]
                    ir_max = min(length(psp_β), length(upf_β))
                    
                    psp_ekb = pair.psp["ekb"][l+1][i]
                    upf_Dij = pair.upf["D_ion"][idx_upf,idx_upf]

                    psp_βekbβ = psp_β[1:ir_max] * psp_ekb * psp_β[1:ir_max]'
                    upf_βekbβ = upf_β[1:ir_max] * upf_Dij * upf_β[1:ir_max]'

                    @test psp_βekbβ ≈ upf_βekbβ ./ 2 atol=1e-6 rtol=1e-6
                    idx_upf += 1
                end
            end
        end
    end
end