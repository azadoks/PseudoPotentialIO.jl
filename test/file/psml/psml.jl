using PseudoPotentialIO: psml_to_dict, save_psp, psml_grid, _to_dict

@testset "PsmlFile" begin
    text = """
    <?xml version="1.0" encoding="UTF-8" ?>
    <psml version="1.1" energy_unit="hartree" length_unit="bohr" uuid="ecf3"
          xmlns="http://esl.cecam.org/PSML/ns/1.1">
      <provenance creator="ONCVPSP" date="2017-08-01">
        <annotation action="generation" />
        <input-file name="INP"><![CDATA[Fe 26.00 3 4 4 psp8]]></input-file>
      </provenance>
      <pseudo-atom-spec atomic-label="Fe" atomic-number="26" z-pseudo="16"
                        relativity="dirac" core-corrections="yes" flavor="oncvpsp">
        <exchange-correlation>
          <libxc-info number-of-functionals="2">
            <functional name="PBE" type="exchange" id="101" />
            <functional name="PBE" type="correlation" id="130" />
          </libxc-info>
        </exchange-correlation>
        <valence-configuration total-valence-charge="16">
          <shell n="3" l="s" occupation="2" />
          <shell n="3" l="d" occupation="6" occupation-up="5" occupation-down="1" />
          <shell n="3" l="p" occupation="6" />
          <shell n="4" l="s" occupation="2" />
        </valence-configuration>
        <core-configuration total-core-charge="10">
          <shell n="1" l="s" occupation="2" />
          <shell n="2" l="s" occupation="2" />
          <shell n="2" l="p" occupation="6" />
        </core-configuration>
      </pseudo-atom-spec>
      <grid npts="4"><grid-data> 0 1 2 3 </grid-data></grid>
      <valence-charge total-charge="16" is-unscreening-charge="yes">
        <radfunc><data> 1 2 3 4 </data></radfunc>
      </valence-charge>
      <pseudocore-charge matching-radius="0.5" number-of-continuous-derivatives="4">
        <radfunc><data npts="2"> 5 6 </data></radfunc>
      </pseudocore-charge>
      <semilocal-potentials set="lj">
        <slps n="3" l="p" j="1.5" rc="1.1" eref="-2.2">
          <radfunc><data> 1 2 3 4 </data></radfunc>
        </slps>
      </semilocal-potentials>
      <local-potential type="oncv-fit">
        <radfunc>
          <grid npts="2"><grid-data> 0 0.5 </grid-data></grid>
          <data> 7 8 </data>
        </radfunc>
        <local-charge><radfunc><data> 1 1 1 1 </data></radfunc></local-charge>
      </local-potential>
      <nonlocal-projectors set="lj">
        <grid npts="3"><grid-data> 0 0.1 0.2 </grid-data></grid>
        <proj l="d" j="2.5" seq="1" ekb="12.7" type="oncv">
          <radfunc><data> 1 2 3 </data></radfunc>
        </proj>
      </nonlocal-projectors>
      <pseudo-wave-functions set="lj">
        <pswf n="3" l="s" j="0.5" energy_level="-3.5"><radfunc><data> 1 2 3 4 </data></radfunc></pswf>
      </pseudo-wave-functions>
    </psml>
    """
    file = PsmlFile(IOBuffer(text); identifier="Fe.psml")
    @test file isa PsmlFile && file.version == "1.1" && file.uuid == "ecf3"
    @test identifier(file) == "Fe.psml" && format(file) == "PSML"
    @test element(file) == elements[:Fe]
    @test functional(file) == [:gga_x_pbe, :gga_c_pbe]
    @test valence_charge(file) === 16.0
    @test is_norm_conserving(file) && !is_ultrasoft(file) && !is_paw(file)
    @test has_spin_orbit(file) && has_model_core_charge_density(file)

    p = only(file.provenance)
    @test p.annotation == Dict("action" => "generation") && isnothing(p.record_number)
    @test only(p.input_files).content == "Fe 26.00 3 4 4 psp8"
    spec = file.pseudo_atom_spec
    @test spec.core_corrections && isnothing(spec.spin_dft) && spec.flavor == "oncvpsp"
    @test [s.l for s in spec.valence_configuration.shells] == [0, 2, 1, 0]
    @test spec.valence_configuration.shells[2].occupation_up === 5.0
    @test spec.core_configuration.total_charge === 10.0
    @test file.valence_charge.is_unscreening_charge === true
    @test isnothing(file.valence_charge.rescaled_to_z_pseudo)
    @test file.pseudocore_charge.number_of_continuous_derivatives === 4
    slps = only(only(file.semilocal_potentials).slps)
    @test (slps.l, slps.j, slps.n, slps.eref) == (1, 1.5, 3, -2.2)
    pswf = only(only(file.pseudo_wave_functions).pswfs)
    @test pswf.energy_level === -3.5 && pswf.j === 0.5

    # the grid of a radial function: its own, its group's, or the top-level one
    group = only(file.nonlocal_projectors)
    proj = only(group.projectors)
    @test (proj.l, proj.j, proj.seq, proj.ekb, proj.type) == (2, 2.5, 1, 12.7, "oncv")
    @test psml_grid(file, proj.radfunc, group) == [0.0, 0.1, 0.2]
    @test psml_grid(file, file.local_potential.radfunc) == [0.0, 0.5]
    @test psml_grid(file, file.pseudocore_charge.radfunc) == [0.0, 1.0]  # npts = 2
    @test file.local_potential.local_charge.data == ones(4)

    # struct -> dictionaries -> file -> struct
    @test _to_dict(file) == psml_to_dict(text)
    io = IOBuffer()
    save_psp(io, file)
    @test _to_dict(PsmlFile(IOBuffer(take!(io)))) == _to_dict(file)

    @test_throws ErrorException PsmlFile(IOBuffer(replace(text, """ l="d" j""" => """ l="x" j""")))
    @test_throws ArgumentError PsmlFile(IOBuffer(replace(text, "number-of-functionals=\"2\"" =>
                                                                 "number-of-functionals=\"3\""));
                                              validation=:error)
    @test_throws ErrorException PsmlFile(IOBuffer(replace(text, " ekb=\"12.7\"" => "")))

    # PSML 0.8: the header, the operator, and a text annotation
    old = """
    <psml version="0.8" energy_unit="hartree" length_unit="bohr">
      <provenance creator="oncvpsp" date="2017" />
      <header atomic-label="Au" atomic-number="79." z-pseudo="19." relativity="dirac"
              polarized="no" core-corrections="yes">
        <exchange-correlation>
          <libxc-info number-of-functionals="1">
            <functional name="XC_NOT_IMPL" type="XC_EXCHANGE_CORRELATION" id="-1" />
          </libxc-info>
        </exchange-correlation>
        <valence-configuration total-valence-charge="19.">
          <shell n="5" l="s" occupation="2." />
          <shell n="5" l="p" occupation="6." />
          <shell n="5" l="d" occupation="10." />
          <shell n="6" l="s" occupation="1." />
        </valence-configuration>
      </header>
      <grid npts="2"><grid-data> 0 1 </grid-data></grid>
      <valence-charge total-charge="19."><radfunc><data> 1 2 </data></radfunc></valence-charge>
      <pseudocore-charge annotation="Monotonic polynomial">
        <radfunc><data> 3 4 </data></radfunc>
      </pseudocore-charge>
      <pseudopotential-operator version="0.1" energy_unit="hartree" length_unit="bohr">
        <provenance creator="oncvpsp 2.1.2+psml" />
        <local-potential type="oncv-fit"><radfunc><data> 5 6 </data></radfunc></local-potential>
        <projectors set="lj">
          <proj l="s" j="0.0" seq="1" ekb="9.2" type="oncv"><radfunc><data> 7 8 </data></radfunc></proj>
        </projectors>
      </pseudopotential-operator>
    </psml>
    """
    file = PsmlFile(IOBuffer(old))
    @test file.version == "0.8" && element(file) == elements[:Au] && isempty(functional(file))
    @test file.pseudo_atom_spec.spin_dft === false
    @test [p.creator for p in file.provenance] == ["oncvpsp", "oncvpsp 2.1.2+psml"]
    @test file.pseudocore_charge.annotation == Dict("text" => "Monotonic polynomial")
    @test file.local_potential.radfunc.data == [5.0, 6.0]
    @test only(only(file.nonlocal_projectors).projectors).ekb === 9.2
    d = _to_dict(file)  # written as PSML 1.1
    @test d["version"] == "1.1" && haskey(d, "pseudo-atom-spec") && haskey(d, "nonlocal-projectors")
    @test !haskey(d, "header") && !haskey(d, "pseudopotential-operator")
    io = IOBuffer()
    save_psp(io, file)
    @test _to_dict(PsmlFile(IOBuffer(take!(io)))) == d
end
