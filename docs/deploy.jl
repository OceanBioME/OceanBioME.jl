using Documenter

deploydocs(
    repo = "github.com/OceanBioME/OceanBioME.jl",
    versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
    forcepush = true,
    push_preview = true,
    devbranch = "main"
)
