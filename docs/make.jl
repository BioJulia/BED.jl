using Pkg

Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, BED

makedocs(
    format = Documenter.HTML(
        edit_branch = "develop"
    ),
    sitename = "BED.jl",
    modules = [BED],
    pages = [
        "Home" => "index.md",
        "API Reference" => "lib/public.md"
    ],
    authors = "Kenta Sato, D. C. Jones, Ben J. Ward, Ciarán O'Mara, The BioJulia Organisation and other contributors.",
)
deploydocs(
    repo = "github.com/BioJulia/BED.jl.git",
    devbranch = "develop"
)
