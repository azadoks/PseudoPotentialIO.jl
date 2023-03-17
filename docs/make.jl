using PseudoPotentialIO
using Documenter
using Literate

DEBUG = false  # Set to true to disable some checks and cleanup
SRCPATH   = joinpath(@__DIR__, "src")
BUILDPATH = joinpath(@__DIR__, "build")
ROOTPATH  = joinpath(@__DIR__, "..")
CONTINUOUS_INTEGRATION = get(ENV, "CI", nothing) == "true"
PAGES = ["Home" => "index.md",
         "Tutorial" => "tutorial.jl",
         "Quantities" => "quantities.md",
         "Formats" => "formats.md",
         "api.md"]

# Get list of files from PAGES
extract_paths(pages::AbstractArray) = collect(Iterators.flatten(extract_paths.(pages)))
extract_paths(file::AbstractString) = [file]
extract_paths(pair::Pair) = extract_paths(pair.second)

# Transform files to *.md
transform_to_md(pages::AbstractArray) = transform_to_md.(pages)
transform_to_md(file::AbstractString) = first(splitext(file)) * ".md"
transform_to_md(pair::Pair) = (pair.first => transform_to_md(pair.second))

# Collect files to treat with Literate (i.e. the examples and the .jl files in the docs)
# The examples go to docs/literate_build/examples, the .jl files stay where they are
literate_files = map(filter!(endswith(".jl"), extract_paths(PAGES))) do file
       if startswith(file, "examples/")
           (src=joinpath(ROOTPATH, file), dest=joinpath(SRCPATH, "examples"), example=true)
       else
           (src=joinpath(SRCPATH, file), dest=joinpath(SRCPATH, dirname(file)), example=false)
       end
   end

# Run Literate on them all
for file in literate_files
       preprocess = file.example ? add_badges : identity
       Literate.markdown(file.src, file.dest;
                         flavor=Literate.DocumenterFlavor(),
                         credit=false, preprocess)
       Literate.notebook(file.src, file.dest; credit=false,
                         execute=CONTINUOUS_INTEGRATION || DEBUG)
   end

DocMeta.setdocmeta!(PseudoPotentialIO, :DocTestSetup, :(using PseudoPotentialIO);
                    recursive=true)

makedocs(;
         modules=[PseudoPotentialIO],
         authors="Austin Zadoks",
         repo="https://github.com/azadoks/PseudoPotentialIO.jl/blob/{commit}{path}#{line}",
         sitename="PseudoPotentialIO.jl",
         format=Documenter.HTML(;
                                prettyurls=CONTINUOUS_INTEGRATION,
                                canonical="https://azadoks.github.io/PseudoPotentialIO.jl",
                                assets=String[]),
         pages=transform_to_md(PAGES),
         checkdocs=:exports,
         strict=true)

deploydocs(;
           repo="github.com/azadoks/PseudoPotentialIO.jl",
           devbranch="master")
