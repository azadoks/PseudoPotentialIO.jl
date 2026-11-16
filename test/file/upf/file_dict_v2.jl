using PseudoPotentialIO: upf_to_dict, dict_to_upf

@testset "UPF v2 dictionaries" begin
    @testset "Values" begin
        text = """
        <UPF version="2.0.1">
          <PP_INFO>a &amp; b, R &lt; 2 & c<PP_INPUTFILE>&input /</PP_INPUTFILE></PP_INFO>
          <PP_HEADER element="F" date="180509" relativistic="full" is_paw=".false."
                     has_so="T" core_correction="true" z_valence="7.00" l_max="2" />
          <PP_R size="4"> 1.0D-3 2.0-319 2*3.5 </PP_R>
          <PP_NONLOCAL>
            <PP_BETA.1 index="*" angular_momentum="1.0"> 1 2 </PP_BETA.1>
            <PP_BETA.2> 3 4 </PP_BETA.2>
          </PP_NONLOCAL>
          <PP_QIJ.1.2 composite_index="2"> 0 </PP_QIJ.1.2>
        </UPF>
        trailing junk"""
        d = upf_to_dict(text)
        @test d["version"] == "2.0.1"
        @test d["info"]["content"] == "a & b, R < 2 & c"
        @test d["info"]["inputfile"]["content"] == "&input /"
        h = d["header"]
        @test h["element"] == "F" && h["date"] == "180509" && h["relativistic"] == "full"
        @test h["is_paw"] === false && h["has_so"] === true && h["core_correction"] === true
        @test h["z_valence"] === 7.0 && h["l_max"] === 2
        @test d["r"] == [1e-3, 2.0e-319, 3.5, 3.5]
        # the suffix wins over the index attribute; elements without attributes are dicts
        @test [b["index"] for b in d["nonlocal"]["beta"]] == [1, 2]
        @test d["nonlocal"]["beta"][1]["angular_momentum"] === 1.0
        @test d["nonlocal"]["beta"][2]["content"] == [3.0, 4.0]
        @test d["qij"][1]["first_index"] == 1 && d["qij"][1]["second_index"] == 2

        # written and read back: the same dictionaries
        io = IOBuffer()
        dict_to_upf(io, d)
        @test upf_to_dict(String(take!(io))) == d

        @test_throws ErrorException upf_to_dict("""<UPF version="2.0.1"><PP_R size="3"> 1 2 </PP_R></UPF>""")
        @test_throws ErrorException upf_to_dict("""<qe_pp:pseudo></qe_pp:pseudo>""")
    end

    @testset "Round trip: $name" for (name, filepath) in UPF2_CASE_FILEPATHS
        d = upf_to_dict(read(filepath, String))
        io = IOBuffer()
        dict_to_upf(io, d)
        @test upf_to_dict(String(take!(io))) == d
    end
end
