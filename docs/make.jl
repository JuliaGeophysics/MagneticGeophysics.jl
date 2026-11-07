using Documenter, DocumenterVitepress

makedocs(;
    sitename = "MagneticGeophysics.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = DocumenterVitepress.MarkdownVitepress(;
        repo       = "github.com/JuliaGeophysics/MagneticGeophysics.jl",
        devbranch  = "main",
        devurl     = "dev",
        # full URL with https://, otherwise the host is taken as part of the base path
        deploy_url = "https://juliageophysics.com/MagneticGeophysics.jl",
        description = "A software repository for magnetic geophysics research and applications",
        # dev is the only published version, so let search engines index it
        noindex_non_stable = false,
    ),
    # nested lists become the dropdown menus of the top navigation bar
    pages = [
        "Home" => "index.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/MagneticGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)
