# # Tutorial
#
# Here, we'll take a look at the basic usage of PseudoPotentialIO.jl: searching for, loading, and investigating the contents of pseudopotentials.

using PseudoPotentialIO
using CairoMakie
using Colors

# ## 1. Searching for pseudopotentials
# PseudoPotentialIO provides easy access to a variety of pseudopotential families with zero setup using Julia [Artifacts](https://docs.julialang.org/en/v1/stdlib/Artifacts/), [LazyArtifacts.jl](https://github.com/JuliaPackaging/LazyArtifacts.jl), and [PseudoLibrary](https://github.com/JuliaMolSim/PseudoLibrary).
# In order to list the available pre-bundled families, use `list_families`

list_families(with_info=true)

# Because families are downloaded lazily, detailed information on some of the families is missing (they need to be downloaded first).
# In order to download a pseudopotential family, you can load it using `load_family`.
#
# !!! note "Loading pseudopotential families"
#     `load_family` can also load all the pseudopotentials in a local directory!
#
# You can then see an overview of which elements the family supports with `show_family_periodic_table`.

family = load_family("hgh_lda_upf");
show_family_periodic_table(family)

# For more detailed information, use `show_family_list`.
# You can restrict the output by providing a list of elements that you're interested in.

show_family_list(family)  # Show all the pseudopotentials
show_family_list(family; elements=["Ba", "Ti", "O"])  # Only show the pseudos for Ba, Ti, and O

# ## 2. Loading pseudopotential files
# To load an individual pseudopotential _file_, use `load_psp_file`, specifying the family name or directory and the filename of the pseudopotential

Ba_psp_file = load_psp_file("hgh_lda_upf", "Ba.pz-sp-hgh.UPF")

# PseudoPotentialIO distinguishes between pseudopotential _files_ and the pseudopotentials themselves.
# Structures like `HghFile`, `UpfFile`, and `Psp8File` correspond to pseudopotential file formats and make the quantities that these files contain directly available, with no unit conversion or processing.
# For example, we can take a look at the author field from the header in the barium UPF pseudopotential we just loaded

Ba_psp_file.header.author

# We can also check that the properties of Ba_psp_file match up with the sections of a UPF file

propertynames(Ba_psp_file)

# ## 3. Loading pseudopotentials
# Once you've decided that you would like to use a given pseudopotential for a calculation, either convert its `PsPFile` struct a corresponding pseudopotential structure

Ba_psp_from_File = load_psp(Ba_psp_file)

# , or load the pseudopotential structure directly from the file using `load_psp`

Ba_psp_from_disk = load_psp("hgh_lda_upf", "Ba.pz-sp-hgh.UPF")

# We can confirm that these pseudopotentials are identical
Ba_psp_from_File == Ba_psp_from_disk

# This procedure has processed the contents of the UPF file (on disk) or the `UpfFile` struct into a common and consistent data representation that PseudoPotentialIO uses for calculations. We can see that the contents have changed by looking at the property names of our new `UpfPsP` struct

propertynames(Ba_psp_from_disk)

# ## 4. Inspecting pseudopotential quantities
# One thing that we might want to do with a processed pseudopotential is to visualize some of the quantities it contains.
# Let's plot the Kleinman-Bylander projectors from a PseudoDojo barium pseudopotential

Ba_psp = load_psp("pd_nc_sr_pbesol_standard_0.4.1_upf", "Ba.upf");
let
    linestyles = [:solid, :dash, :dot]
    colors = Colors.JULIA_LOGO_COLORS
    fig = Figure(); ax = Axis(fig[1,1], xlabel="r [a₀]", ylabel="β(r)")
    for l in angular_momenta(Ba_psp)                  # Iterate over each angular momentum 0:lmax
        color = colors[l+1]
        for n in beta_projector_radial_indices(Ba_psp, l)  # Iterate over each projector at l 1:nmax
            linestyle = linestyles[n]
            r²βln = Ba_psp.β[l][n]                    # Projector multiplied by r²
            i_rc_ln = lastindex(r²βln)                # Cutoff radius index
            βln = r²βln ./ Ba_psp.r[1:i_rc_ln].^2     # Remove the r² prefactor
            lines!(ax, Ba_psp.r[1:i_rc_ln], βln, label="|β[$l][$n]⟩",
                   linestyle=linestyle, color=color)
        end
    end
    axislegend()
    fig
end
