using PseudoPotentialIO: psml_to_dict, dict_to_psml, _to_dict

@testset "PSML dictionaries" begin
    text = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <psml version="1.1" energy_unit="hartree" length_unit="bohr" uuid="ecf3"
          xmlns="http://esl.cecam.org/PSML/ns/1.1">
      <provenance creator="ATM" date="2017-08-01">
        <annotation action="generation" scale="0.5E-04" />
        <input-file name="INP"><![CDATA[ n=F c=ca & <x> ]]></input-file>
      </provenance>
      <pseudo-atom-spec atomic-label="F" atomic-number="9" z-pseudo="7."
                        relativity="no" spin-dft="no" core-corrections="yes"
                        flavor="Hamann&apos;s oncvpsp">
        <valence-configuration total-valence-charge="7">
          <shell n="2" l="s" occupation="2" />
        </valence-configuration>
      </pseudo-atom-spec>
      <grid npts="3"><grid-data> 0.0 1.0D-1 2.0 </grid-data></grid>
      <valence-charge total-charge="7"><radfunc><data> 1 2 3 </data></radfunc></valence-charge>
      <local-potential type="l=0">
        <grid npts="2"><grid-data> 0 1 </grid-data></grid>
        <radfunc><data> 1 2 </data></radfunc>
      </local-potential>
      <nonlocal-projectors set="lj">
        <proj l="p" j="1.5" seq="1" ekb="2.5" type="oncv">
          <radfunc><data npts="2"> 1 2 </data></radfunc>
        </proj>
      </nonlocal-projectors>
    </psml>
    """
    d = psml_to_dict(text)
    @test d["version"] == "1.1" && d["xmlns"] == "http://esl.cecam.org/PSML/ns/1.1"
    p = only(d["provenance"])
    @test p["date"] == "2017-08-01"
    @test p["annotation"]["scale"] == "0.5E-04"  # annotations are free text
    @test only(p["input-file"])["content"] == " n=F c=ca & <x> "
    a = d["pseudo-atom-spec"]
    @test a["atomic-label"] == "F" && a["atomic-number"] === 9 && a["z-pseudo"] === 7.0
    @test a["relativity"] == "no" && a["spin-dft"] === false && a["core-corrections"] === true
    @test a["flavor"] == "Hamann's oncvpsp"
    @test only(a["valence-configuration"]["shell"])["l"] == "s"
    @test d["grid"] == Dict("grid-data" => [0.0, 0.1, 2.0])
    @test d["valence-charge"]["radfunc"]["data"] == [1.0, 2.0, 3.0]
    @test d["local-potential"]["radfunc"]["data"] == [1.0, 2.0]  # the local grid
    proj = only(only(d["nonlocal-projectors"])["proj"])
    @test proj["j"] === 1.5 && proj["seq"] === 1
    @test proj["radfunc"]["data"] == [1.0, 2.0]  # shorter than the grid, with npts

    io = IOBuffer()
    dict_to_psml(io, d)
    @test psml_to_dict(String(take!(io))) == d

    # data must have the size of its grid, or npts no larger than it
    @test_throws ErrorException psml_to_dict(replace(text, "<data> 1 2 3 </data>" => "<data> 1 2 </data>"))
    @test_throws ErrorException psml_to_dict(replace(text, """<data npts="2">""" => """<data npts="4">"""))
    @test_throws ErrorException psml_to_dict(replace(text, """<grid npts="3">""" => """<grid npts="4">"""))
    @test_throws ErrorException psml_to_dict(replace(text, """<grid npts="3">""" => "<grid>"))
    @test_throws ErrorException psml_to_dict("""<UPF version="2.0.1"></UPF>""")

    @testset "Round trip: $(basename(path))" for path in PSML_FILEPATHS
        d = psml_to_dict(read(path, String))
        io = IOBuffer()
        dict_to_psml(io, d)
        @test psml_to_dict(String(take!(io))) == d
        # the layout of the struct is that of PSML 1.1
        d["version"] in ("0.8", "1.0") || @test _to_dict(load_psp_file(path)) == d
    end
end
