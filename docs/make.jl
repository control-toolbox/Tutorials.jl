using Documenter
using OptimalControl

# For reproducibility
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(joinpath(@__DIR__, "Manifest.toml"), joinpath(@__DIR__, "src", "assets", "Manifest.toml"), force = true)
cp(joinpath(@__DIR__, "Project.toml"), joinpath(@__DIR__, "src", "assets", "Project.toml"), force = true)

repo_url = "github.com/control-toolbox/Tutorials.jl"

makedocs(;
    draft=true, # if draft is true, then the julia code from .md is not executed
    # to disable the draft mode in a specific markdown file, use the following:
    # ```@meta
    # Draft = false
    # ```
    warnonly=[:cross_references, :autodocs_block],
    sitename="Tutorials",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=[""],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Getting Started" => "index.md",
        "Tutorials and Advanced Features" => [
            "Discrete continuation" => "tutorial-continuation.md",
            "Discretisation methods" => "tutorial-discretisation.md",
            "NLP manipulations" => "tutorial-nlp.md",
            "Indirect simple shooting" => "tutorial-iss.md",
            "Goddard: direct, indirect" => "tutorial-goddard.md",
            "Linear–quadratic regulator" => "tutorial-lqr.md",
            "Minimal action" => "tutorial-mam.md",
            "Model Predictive Control" => "tutorial-mpc.md",
        ],
    ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
