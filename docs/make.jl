using TOML
using Documenter, BED

format = Documenter.HTML(
    edit_link = "develop"
)

authors = let
    project_path = joinpath(Base.pkgdir(BED), "Project.toml")
    toml = TOML.parsefile(project_path)
    authors_with_email = join(toml["authors"], ", ")
    authors = replace(authors_with_email, r" <.*?>" => "")
    authors * ", The BioJulia Organisation, and other contributors."
end
    

makedocs(
    format = format,
    checkdocs = :all,
    #linkcheck = true,
    modules = [BED],
    sitename = "BED.jl",
    pages = [
        "Home" => "index.md",
        "BED" => "man/bed.md",
        "API Reference" => "man/api.md"
    ],
    authors = authors,
)

deploydocs(
    repo = "github.com/BioJulia/BED.jl.git",
    devbranch = "develop",
    push_preview = true
)
