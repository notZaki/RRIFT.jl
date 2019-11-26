using Documenter, RRIFT

makedocs(;
    modules=[RRIFT],
    format=Documenter.HTML(),
    pages=[
        "Introduction" => "index.md",
        "Guide" => Any[
            "Downloading" => "guide/downloading.md",
            "Pre-processing" => "guide/preprocessing.md",
            "Fitting" => "guide/fitting.md"
         ]
    ],
    repo="https://github.com/notZaki/RRIFT.jl/blob/{commit}{path}#L{line}",
    sitename="RRIFT.jl",
    authors="Zaki A",
)

deploydocs(;
    repo="github.com/notZaki/RRIFT.jl",
)
