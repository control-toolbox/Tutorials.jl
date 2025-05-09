using Documenter
using OptimalControl

repo_url = "github.com/control-toolbox/Tutorials.jl"

makedocs(;
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
            "Linear–quadratic regulator" => "tutorial-lqr-basic.md",
            "Minimal action" => "tutorial-mam.md",
        ],
    ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
