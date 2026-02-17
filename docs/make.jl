using Daedalus
using Documenter

DocMeta.setdocmeta!(Daedalus, :DocTestSetup, :(using Daedalus); recursive=true)

makedocs(;
    modules=[Daedalus],
    authors="Pratik Gupte <pratikgupte16@gmail.com> and contributors",
    sitename="Daedalus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pratikunterwegs.github.io/Daedalus.jl",
        edit_link="main",
        assets=String[],
    ),
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Index" => "pkg_index.md",
        "Function Reference" => "reference.md",
        "Reactive NPIs" => "npis.md",
        "Implementing reactive events" => "musings.md",
        "Benchmarking" => "benchmarking.md"
    ],
)

deploydocs(;
    repo="github.com/pratikunterwegs/Daedalus.jl",
    devbranch="main",
)
