upf1_filepaths = Dict(
    "ag_lda_v1.4.uspp.F.upf" => artifact"gbrv_lda_1.5_upf",
    "B_pbe_v1.01.uspp.F.upf" => artifact"sssp_pbe_precision_1.1.2_upf",
    "si_pbesol_v1.uspp.F.upf" => artifact"gbrv_pbesol_1.5_upf",
    "mg_pbe_v1.4.uspp.F.upf" => artifact"gbrv_pbe_1.5_upf"
)
upf1_filepaths = Dict(key => joinpath(value, key) for (key, value) in upf1_filepaths)

upf2_filepaths = Dict(
    "Mg.upf" => artifact"pd_nc_fr_pbesol_stringent_0.4_upf",
    "Si.pbe-n-rrkjus_psl.1.0.0.upf" => artifact"sssp_pbe_efficiency_1.1.2_upf",
    "Al.pbe-n-kjpaw_psl.1.0.0.upf" => artifact"sssp_pbe_efficiency_1.1.2_upf",
    "Dy.GGA-PBE-paw-v1.0.UPF" => artifact"sssp_pbe_efficiency_1.1.2_upf",
    "He.pbe-hgh.UPF" => artifact"hgh_pbe_upf"
)
upf2_filepaths = Dict(key => joinpath(value, key) for (key, value) in upf2_filepaths)

psp8_filepaths = Dict(
    "H.psp8" => artifact"pd_nc_sr_pbe_standard_0.4.1_psp8",
    "Ti.psp8" => artifact"pd_nc_fr_pbesol_standard_0.4_psp8",
    "Zn.psp8" => artifact"pd_nc_sr_pbesol_stringent_0.4.1_psp8"
)
psp8_filepaths = Dict(key => joinpath(value, key) for (key, value) in psp8_filepaths)
