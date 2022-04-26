using Documenter
using AMRVW

DocMeta.setdocmeta!(AMRVW, :DocTestSetup, :(using AMRVW); recursive=true)

makedocs(
    sitename = "AMRVW",
    format = Documenter.HTML(ansicolor=true),
    modules = [AMRVW]
)

deploydocs(
    repo = "github.com/jverzani/AMRVW.jl"
)
