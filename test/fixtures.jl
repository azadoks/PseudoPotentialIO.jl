families = (artifact"gbrv_lda_1.5_upf", artifact"sssp_pbe_precision_1.1.2_upf",
            artifact"hgh_lda_hgh", artifact"hgh_pbe_upf",
            artifact"pd_nc_sr_pbe_standard_0.4.1_psp8",
            artifact"pd_nc_fr_pbesol_standard_0.4_psp8", artifact"sg15_2022.02.06_upf",
            artifact"pd_nc_sr_pbe_standard_0.4.1_upf")

upf1_filepaths = Dict("ag_lda_v1.4.uspp.F.UPF" => artifact"gbrv_lda_1.5_upf",
                      "B_pbe_v1.01.uspp.F.UPF" => artifact"sssp_pbe_precision_1.1.2_upf",
                      "si_pbesol_v1.uspp.F.UPF" => artifact"gbrv_pbesol_1.5_upf",
                      "mg_pbe_v1.4.uspp.F.UPF" => artifact"gbrv_pbe_1.5_upf")
upf1_filepaths = Dict(key => joinpath(value, key) for (key, value) in upf1_filepaths)

upf2_filepaths = Dict("Mg.upf" => artifact"pd_nc_fr_pbesol_stringent_0.4_upf",
                      "Si.pbe-n-rrkjus_psl.1.0.0.UPF" => artifact"sssp_pbe_efficiency_1.1.2_upf",
                      "Al.pbe-n-kjpaw_psl.1.0.0.UPF" => artifact"sssp_pbe_efficiency_1.1.2_upf",
                      "Dy.GGA-PBE-paw-v1.0.UPF" => artifact"sssp_pbe_efficiency_1.1.2_upf",
                      "He.pbe-hgh.UPF" => artifact"hgh_pbe_upf",
                      "Al.pbe-hgh.UPF" => artifact"hgh_pbe_upf",
                      "H.upf" => artifact"pd_nc_sr_pbe_standard_0.4.1_upf",
                      "Fe.upf" => artifact"pd_nc_sr_pbe_standard_0.4.1_upf")
upf2_filepaths = Dict(key => joinpath(value, key) for (key, value) in upf2_filepaths)

psp8_filepaths = Dict("H.psp8" => artifact"pd_nc_sr_pbe_standard_0.4.1_psp8",
                      "Fe.psp8" => artifact"pd_nc_sr_pbe_standard_0.4.1_psp8",
                      "Ti.psp8" => artifact"pd_nc_fr_pbesol_standard_0.4_psp8",
                      "Zn.psp8" => artifact"pd_nc_sr_pbesol_stringent_0.4.1_psp8")
psp8_filepaths = Dict(key => joinpath(value, key) for (key, value) in psp8_filepaths)

hgh_filepaths = Dict("c-q4.hgh" => artifact"hgh_lda_hgh",
                     "ni-q18.hgh" => artifact"hgh_lda_hgh")
hgh_filepaths = Dict(key => joinpath(value, key) for (key, value) in hgh_filepaths)

upf2_hgh_pairs = [("H.pbe-hgh.UPF", "h-q1.hgh"),
                  ("He.pbe-hgh.UPF", "he-q2.hgh"),
                  ("B.pbe-hgh.UPF", "b-q3.hgh"),
                  ("C.pbe-hgh.UPF", "c-q4.hgh"),
                  ("N.pbe-hgh.UPF", "n-q5.hgh"),
                  ("O.pbe-hgh.UPF", "o-q6.hgh"),
                  ("F.pbe-hgh.UPF", "f-q7.hgh"),
                  ("Ne.pbe-hgh.UPF", "ne-q8.hgh"),
                  ("Mg.pbe-hgh.UPF", "mg-q2.hgh"),
                  ("Al.pbe-hgh.UPF", "al-q3.hgh"),
                  ("Si.pbe-hgh.UPF", "si-q4.hgh"),
                  ("P.pbe-hgh.UPF", "p-q5.hgh"),
                  ("S.pbe-hgh.UPF", "s-q6.hgh"),
                  ("Cl.pbe-hgh.UPF", "cl-q7.hgh"),
                  ("Ar.pbe-hgh.UPF", "ar-q8.hgh"),
                  ("Ga.pbe-hgh.UPF", "ga-q3.hgh"),
                  ("Ge.pbe-hgh.UPF", "ge-q4.hgh"),
                  ("As.pbe-hgh.UPF", "as-q5.hgh"),
                  ("Se.pbe-hgh.UPF", "se-q6.hgh"),
                  ("Br.pbe-hgh.UPF", "br-q7.hgh"),
                  ("Kr.pbe-hgh.UPF", "kr-q8.hgh"),
                  ("Ru.pbe-hgh.UPF", "ru-q8.hgh"),
                  ("Rh.pbe-hgh.UPF", "rh-q9.hgh"),
                  ("Pd.pbe-hgh.UPF", "pd-q10.hgh"),
                  ("In.pbe-hgh.UPF", "in-q3.hgh"),
                  ("Sn.pbe-hgh.UPF", "sn-q4.hgh"),
                  ("Sb.pbe-hgh.UPF", "sb-q5.hgh"),
                  ("Te.pbe-hgh.UPF", "te-q6.hgh"),
                  ("I.pbe-hgh.UPF", "i-q7.hgh"),
                  ("Xe.pbe-hgh.UPF", "xe-q8.hgh"),
                  ("La.pbe-hgh.UPF", "la-q11.hgh"),
                  ("Ta.pbe-hgh.UPF", "ta-q5.hgh"),
                  ("W.pbe-hgh.UPF", "w-q6.hgh"),
                  ("Re.pbe-hgh.UPF", "re-q7.hgh"),
                  ("Os.pbe-hgh.UPF", "os-q8.hgh"),
                  ("Ir.pbe-hgh.UPF", "ir-q9.hgh"),
                  ("Pt.pbe-hgh.UPF", "pt-q10.hgh"),
                  ("Tl.pbe-hgh.UPF", "tl-q3.hgh"),
                  ("Pb.pbe-hgh.UPF", "pb-q4.hgh"),
                  ("Bi.pbe-hgh.UPF", "bi-q5.hgh"),
                  ("Po.pbe-hgh.UPF", "po-q6.hgh"),
                  ("At.pbe-hgh.UPF", "at-q7.hgh"),
                  ("Rn.pbe-hgh.UPF", "rn-q8.hgh")]
upf2_hgh_pairs = [(joinpath(artifact"hgh_pbe_upf", pair[1]),
                   joinpath(artifact"hgh_pbe_hgh", pair[2])) for pair in upf2_hgh_pairs]

upf2_psp8_elements = ("Ag", "Al", "Ar", "As", "Au", "B", "Ba", "Be", "Bi", "Br", "C", "Ca",
                      # "Cr", "Cu", # TODO Chromium and Copper agreement is terrible for some reason
                      "Cd", "Cl", "Co", "Cs", "F", "Fe", "Ga", "Ge", "H", "He",
                      "Hf", "Hg", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lu", "Mg", "Mn",
                      "Mo", "N", "Na", "Nb", "Ne", "Ni", "O", "Os", "P", "Pb", "Pd", "Po",
                      "Pt", "Rb", "Re", "Rh", "Rn", "Ru", "S", "Sb", "Sc", "Se", "Si", "Sn",
                      "Sr", "Ta", "Tc", "Te", "Ti", "Tl", "V", "W", "Xe", "Y", "Zn", "Zr")
upf2_psp8_pairs = [(joinpath(artifact"pd_nc_sr_pbe_standard_0.4.1_upf", "$(element).upf"),
                    joinpath(artifact"pd_nc_sr_pbe_standard_0.4.1_psp8", "$(element).psp8"))
                   for element in upf2_psp8_elements]
