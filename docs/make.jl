using Pkg
using Documenter, BED

format = Documenter.HTML(
    edit_link = "develop"
)

makedocs(
    format = format,
    checkdocs = :all,
    linkcheck = true,
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
