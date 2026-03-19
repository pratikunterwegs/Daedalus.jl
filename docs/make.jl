using Daedalus
using Documenter

DocMeta.setdocmeta!(Daedalus, :DocTestSetup, :(using Daedalus); recursive = true)

makedocs(;
    modules = [Daedalus],
    authors = "Pratik Gupte <pratikgupte16@gmail.com> and contributors",
    sitename = "Daedalus.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://jameel-institute.github.io/Daedalus.jl",
        edit_link = "main",
        assets = String[],
        size_threshold_ignore = [
            "ensemble.md",
            "country_data.md"  # prevent HTML size errors during docs build
        ]
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Modelling parameter uncertainty" => "ensemble.md",
        "Implementing interventions" => "npis.md",
        "Multiple contact settings" => "settings.md",
        "Benchmarking" => "benchmarking.md",
        "Parallelisation" => "parallelisation.md",
        "Country and pathogen data" => "country_data.md",
        "Implementing reactive events" => "musings.md",
        "Index" => "pkg_index.md",
        "Function Reference" => "reference.md"
    ]
)

deploydocs(;
    repo = "github.com/jameel-institute/Daedalus.jl",
    devbranch = "main"
)
