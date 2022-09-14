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

@testset "load_upf" begin
    for (root, dirs, files) in walkdir(UPF_DIR)
        for file in files
            upf = load_upf(joinpath(root, file))
        end
    end
end;

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
        @test length(pair.json["beta_projectors"]) == length(pair.upf["beta_projectors"])
        for i = 1:length(pair.json["beta_projectors"])
            if haskey(pair.json["beta_projectors"][i], "label")
                @test pair.json["beta_projectors"][i]["label"] == pair.upf["beta_projectors"][i]["label"]
            end
            @test pair.json["beta_projectors"][i]["angular_momentum"] == pair.upf["beta_projectors"][i]["angular_momentum"]
            cutoff_idx = length(pair.json["beta_projectors"][i]["radial_function"])
            @test all(
                pair.json["beta_projectors"][i]["radial_function"][begin:cutoff_idx] .==
                pair.upf["beta_projectors"][i]["radial_function"][begin:cutoff_idx]
            )
        end
    end
end

@testset "augmentation" begin
    aug_pairs = [pair for pair in pseudo_pairs if haskey(pair.json, "augmentation")]
    for pair in aug_pairs
        for i = 1:length(pair.json["augmentation"])
            @test pair.json["augmentation"][i]["i"] + 1 == pair.upf["augmentation"][i]["i"]
            @test pair.json["augmentation"][i]["j"] + 1 == pair.upf["augmentation"][i]["j"]
            @test pair.json["augmentation"][i]["angular_momentum"] == pair.upf["augmentation"][i]["angular_momentum"]
            @test all(pair.json["augmentation"][i]["radial_function"] .≈ pair.upf["augmentation"][i]["radial_function"])
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
        @test length(pair.json["atomic_wave_functions"]) == length(pair.upf["atomic_wave_functions"])
        for i = 1:length(pair.json["atomic_wave_functions"])
            if haskey(pair.json["atomic_wave_functions"][i], "label")
                @test pair.json["atomic_wave_functions"][i]["label"] == pair.upf["atomic_wave_functions"][i]["label"]
            end
            @test pair.json["atomic_wave_functions"][i]["occupation"] == pair.upf["atomic_wave_functions"][i]["occupation"]
            @test pair.json["atomic_wave_functions"][i]["angular_momentum"] == pair.upf["atomic_wave_functions"][i]["angular_momentum"]
            cutoff_idx = length(pair.json["atomic_wave_functions"][i]["radial_function"])
            @test all(
                pair.json["atomic_wave_functions"][i]["radial_function"][begin:cutoff_idx] .==
                pair.upf["atomic_wave_functions"][i]["radial_function"][begin:cutoff_idx]
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