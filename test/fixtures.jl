# Test files: pseudopotential families of PseudoPotentialData (downloaded on first use),
# and the files of `data/`.

# Every file of a family, in the order of the elements
function family_filepaths(identifier::AbstractString)
    family = PseudoFamily(identifier)
    return [family[element] for element in sort(collect(keys(family)))]
end

const DATA_DIR = joinpath(@__DIR__, "data")

# PSML: small files of each version, from ABINIT's tests and PseudoDojo (see data/psml/README.md)
PSML_FILEPATHS = filter(endswith(".psml"), readdir(joinpath(DATA_DIR, "psml"); join=true))

# For the tests of the `PsPFile` interface
TEST_FAMILIES = (
    "sg15.nc.sr.pbe.v1_2.canonical.upf",
    "sssp.mixed.sr.pbe.v1_3_0.precision.upf",
    "dojo.nc.sr.pbe.v0_4_1.standard.upf",
    "dojo.nc.fr.pbesol.v0_4.standard.upf",
)
TEST_FILEPATHS = vcat(family_filepaths.(collect(TEST_FAMILIES))..., PSML_FILEPATHS)

# For the consistency tests of UPF (v1 and v2; NC, US and PAW; with and without spin-orbit)
UPF_FAMILIES = (
    "sg15.nc.sr.pbe.v1_2.canonical.upf",
    "sg15.nc.fr.pbe.v1_1.canonical.upf",
    "sssp.mixed.sr.pbe.v1_3_0.efficiency.upf",
    "sssp.mixed.sr.pbesol.v1_3_0.precision.upf",
)
UPF_FILEPATHS = vcat(family_filepaths.(collect(UPF_FAMILIES))...,
                     joinpath(DATA_DIR, "56_Ba_m.upf"))

# Files whose values are tested, by their original name (PseudoPotentialData names the files
# by their element)
UPF1_CASES = [
    ("sssp.mixed.sr.pbe.v1_3_0.precision.upf", :B, "B_pbe_v1.01.uspp.F.UPF"),
    ("sssp.mixed.sr.pbe.v1_3_0.precision.upf", :Mg, "mg_pbe_v1.4.uspp.F.UPF"),
]
UPF1_CASE_FILEPATHS = Dict(name => PseudoFamily(family)[element]
                           for (family, element, name) in UPF1_CASES)

UPF2_CASES = [
    ("dojo.nc.fr.pbesol.v0_4.stringent.upf", :Mg, "Mg.upf"),
    ("sssp.mixed.sr.pbe.v1_3_0.efficiency.upf", :Si, "Si.pbe-n-rrkjus_psl.1.0.0.UPF"),
    ("sssp.mixed.sr.pbe.v1_3_0.efficiency.upf", :Al, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ("dojo.nc.sr.pbe.v0_4_1.standard.upf", :H, "H.upf"),
    ("dojo.nc.sr.pbe.v0_4_1.standard.upf", :Fe, "Fe.upf"),
]
UPF2_CASE_FILEPATHS = Dict(name => PseudoFamily(family)[element]
                           for (family, element, name) in UPF2_CASES)
UPF2_CASE_FILEPATHS["56_Ba_m.upf"] = joinpath(DATA_DIR, "56_Ba_m.upf")
UPF2_CASE_FILEPATHS["al.blps.lda.upf"] = joinpath(DATA_DIR, "al.blps.lda.upf")

# Files that break the specification of their format (METAPSP writes relativistic="F"): the
# tests read them without validation, but that of the file itself (in upf2.jl)
const OFF_SPEC_FILES = ("56_Ba_m.upf",)
test_validation(name::AbstractString) = basename(name) in OFF_SPEC_FILES ? :none : :error
load_test_file(path) = load_psp_file(path; validation=test_validation(path))

# Small PSP8 files
const PSP8_TEXTS = Dict(
    # l = 0 with projectors, l = 1 without (no block), the local potential after lmax, a
    # core charge, a valence density with only its first column, and trailing text
    "nlcc" => """
X  test file
3.0000  2.0000  250101    zatom,zion,pspd
8   11   1   4   3   0   pspcod,pspxc,lmax,lloc,mmax,r2well
1.0  0.5  0.0   rchrg fchrg qchrg
1  0  0  0  0    nproj
1  1   extension_switch
0   2.5D+00
1  0.0D+00  1.0D+00
2  1.0D-02  2.0D+00
3  2.0D-02  3.0D+00
4
1  0.0D+00 -1.0D+00
2  1.0D-02 -2.0D+00
3  2.0D-02 -3.0D+00
1  0.0D+00  1.0  2.0  3.0  4.0  5.0
2  1.0D-02  1.0  2.0  3.0  4.0  5.0
3  2.0D-02  1.0  2.0  3.0  4.0  5.0
1  0.0D+00  7.0
2  1.0D-02  8.0
3  2.0D-02  9.0
<INPUT>
</INPUT>
""",
    # spin-orbit blocks, a Libxc functional, and all the valence density columns
    "spin-orbit" => """
X  spin-orbit test file
3.0000  2.0000  250101    zatom,zion,pspd
8   -001012   1   4   3   0   pspcod,pspxc,lmax,lloc,mmax,r2well
1.0  0.0  0.0   rchrg fchrg qchrg
1  1    nproj
3   extension_switch
1   nprojso
0   2.5D+00
1  0.0D+00  1.0D+00
2  1.0D-02  2.0D+00
3  2.0D-02  3.0D+00
1   1.5D+00
1  0.0D+00  0.0D+00
2  1.0D-02  4.0D-01
3  2.0D-02  6.0D-01
4
1  0.0D+00 -1.0D+00
2  1.0D-02 -2.0D+00
3  2.0D-02 -3.0D+00
1   5.0D-01
1  0.0D+00  0.0D+00
2  1.0D-02  1.0D-01
3  2.0D-02  2.0D-01
1  0.0D+00  7.0  7.5  70.0
2  1.0D-02  8.0  8.5  80.0
3  2.0D-02  9.0  9.5  90.0
""")

load_psp8(name) = Psp8File(IOBuffer(PSP8_TEXTS[name]); identifier=name)

# Two files of the HGH LDA family (Hartwigsen, Goedecker, Hutter), as they are
const HGH_TEXTS = Dict(
    "c-q4.hgh" => """
C GTH-PADE-q4 GTH-LDA-q4
    2    2
     0.34883045    2    -8.51377110     1.22843203
    2
     0.30455321    1     9.52284179
     0.23267730    0
""",
    "ni-q18.hgh" => """
Ni GTH-PADE-q18 GTH-LDA-q18
    4    6    8
     0.35000000    2     3.61031072     0.44963832
    3
     0.24510489    2    12.16113071     3.51625420
                                       -9.07892931
     0.23474136    2    -0.82062357     2.54774737
                                       -6.02907069
     0.21494950    1   -13.39506212
""")

load_hgh(name) = HghFile(IOBuffer(HGH_TEXTS[name]); identifier=name)
