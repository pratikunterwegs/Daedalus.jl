using Daedalus
using Documenter

DocMeta.setdocmeta!(Daedalus, :DocTestSetup, :(using Daedalus); recursive=true)

makedocs(;
    modules=[Daedalus],
    authors="Pratik Gupte <pratikgupte16@gmail.com> and contributors",
    sitename="Daedalus.jl",
    format=Documenter.HTML(;
        canonical="https://pratikgupte16@gmail.com.github.io/Daedalus.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pratikgupte16@gmail.com/Daedalus.jl",
    devbranch="main",
)
