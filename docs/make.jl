using Documenter, RRIFT

makedocs(;
    modules=[RRIFT],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/notZaki/RRIFT.jl/blob/{commit}{path}#L{line}",
    sitename="RRIFT.jl",
    authors="Zaki A",
    assets=String[],
)

deploydocs(;
    repo="github.com/notZaki/RRIFT.jl",
)
