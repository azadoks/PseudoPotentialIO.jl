using UPF
using Test
using JSON

UPF_DIR = "./upf/"
JSON_DIR = "./json/"
HEADER_KEYS = ["number_of_proj", "core_correction", "element", "pseudo_type", "z_valence", "mesh_size", "number_of_wfc"]
VECTOR_KEYS = ["radial_grid", "core_charge_density", "local_potential", "total_charge_density"]

pseudo_pairs = []
for (root, dirs, files) in walkdir(JSON_DIR)
    for file in files
        json_pseudo = JSON.parsefile(joinpath(root, file))["pseudo_potential"]
        upf_filename = splitpath(json_pseudo["header"]["original_upf_file"])[end]
        upf_pseudo = load_upf(joinpath(UPF_DIR, upf_filename))
        push!(pseudo_pairs, (json=json_pseudo, upf=upf_pseudo))
    end
end

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
        @test all(pair.json["local_potential"] .== pair.upf["local_potential"])
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

@testset "augmentation" begin
    aug_pairs = [pair for pair in pseudo_pairs if haskey(pair.json, "augmentation")]
    for pair in aug_pairs
        json_aug = pair.json["augmentation"]
        upf_aug = pair.upf["augmentation"]
        for i = eachindex(json_aug)
            @test json_aug[i]["i"] + 1 == upf_aug[i]["i"]
            @test json_aug[i]["j"] + 1 == upf_aug[i]["j"]
            @test json_aug[i]["angular_momentum"] == upf_aug[i]["angular_momentum"]
            @test all(isapprox.(json_aug[i]["radial_function"], upf_aug[i]["radial_function"], atol=4e-10))  # rtol=5e-6))
        end
    end
end

@testset "paw_data" begin
    paw_pairs = [pair for pair in pseudo_pairs if haskey(pair.json, "paw_data")]
    @testset "aug_integrals" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["aug_integrals"]) == length(pair.upf["paw_data"]["aug_integrals"])
            @test all(pair.json["paw_data"]["aug_integrals"] .== pair.upf["paw_data"]["aug_integrals"])
        end
    end
    @testset "aug_multipoles" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["aug_multipoles"]) == length(pair.upf["paw_data"]["aug_multipoles"])
            @test all(pair.json["paw_data"]["aug_multipoles"] .== pair.upf["paw_data"]["aug_multipoles"])
        end
    end
    @testset "ae_wfc" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["ae_wfc"]) == length(pair.upf["paw_data"]["ae_wfc"])
            for i = 1:length(pair.json["paw_data"]["ae_wfc"])
                @test pair.json["paw_data"]["ae_wfc"][i]["angular_momentum"] == pair.upf["paw_data"]["ae_wfc"][i]["angular_momentum"]
                cutoff_idx = length(pair.json["paw_data"]["ae_wfc"][i]["radial_function"])
                @test all(
                    pair.json["paw_data"]["ae_wfc"][i]["radial_function"][begin:cutoff_idx] .==
                    pair.upf["paw_data"]["ae_wfc"][i]["radial_function"][begin:cutoff_idx]
                )
            end
        end
    end
    @testset "ps_wfc" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["ps_wfc"]) == length(pair.upf["paw_data"]["ps_wfc"])
            for i = 1:length(pair.json["paw_data"]["ps_wfc"])
                @test pair.json["paw_data"]["ps_wfc"][i]["angular_momentum"] == pair.upf["paw_data"]["ps_wfc"][i]["angular_momentum"]
                cutoff_idx = length(pair.json["paw_data"]["ps_wfc"][i]["radial_function"])
                @test all(
                    pair.json["paw_data"]["ps_wfc"][i]["radial_function"][begin:cutoff_idx] .==
                    pair.upf["paw_data"]["ps_wfc"][i]["radial_function"][begin:cutoff_idx]
                )
            end
        end
    end
    @testset "occupations" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["occupations"]) == length(pair.upf["paw_data"]["occupations"])
            @test all(pair.json["paw_data"]["occupations"] .== pair.upf["paw_data"]["occupations"])
        end
    end
    @testset "ae_core_charge_density" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["ae_core_charge_density"]) == length(pair.upf["paw_data"]["ae_core_charge_density"])
            @test all(pair.json["paw_data"]["ae_core_charge_density"] .== pair.upf["paw_data"]["ae_core_charge_density"])
        end
    end
    @testset "ae_local_potential" begin
        for pair in paw_pairs
            @test length(pair.json["paw_data"]["ae_local_potential"]) == length(pair.upf["paw_data"]["ae_local_potential"])
            @test all(pair.json["paw_data"]["ae_local_potential"] .== pair.upf["paw_data"]["ae_local_potential"])
        end
    end
end

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