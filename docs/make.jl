using PseudoPotentialIO
using Documenter

DocMeta.setdocmeta!(PseudoPotentialIO, :DocTestSetup, :(using PseudoPotentialIO);
                    recursive=true)

makedocs(;
         modules=[PseudoPotentialIO],
         authors="Austin Zadoks",
         repo="https://github.com/azadoks/PseudoPotentialIO.jl/blob/{commit}{path}#{line}",
         sitename="PseudoPotentialIO.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://azadoks.github.io/PseudoPotentialIO.jl",
                                assets=String[]),
         pages=["Home" => "index.md",
                "api.md"],
         checkdocs=:exports,
         strict=true)

deploydocs(;
           repo="github.com/azadoks/PseudoPotentialIO.jl",
           devbranch="master")
