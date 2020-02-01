using Pkg
using Documenter, BED

makedocs(
    format = Documenter.HTML(
        edit_link = "develop"
    ),
    modules = [BED],
    sitename = "BED.jl",
    pages = [
        "Home" => "index.md",
        "BED" => "man/bed.md",
        "API Reference" => "man/api.md"
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."

)

deploydocs(
    repo = "github.com/BioJulia/BED.jl.git",
    devbranch = "develop",
    push_preview = true
)
