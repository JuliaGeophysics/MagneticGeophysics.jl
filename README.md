<h1 align="center">MagneticGeophysics.jl</h1>

<p align="center"><em>A software repository for magnetic research and applications.</em></p>

<p align="center">
	<a href="https://github.com/JuliaGeophysics/MagneticGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

3D magnetic forward modelling and LSQR inversion on the tensor meshes of MTGeophysics.jl, with topography and lakes.
This is an early version.

## Installation

> [!TIP]
> **New to Julia?** Install it first from [julialang.org/install](https://julialang.org/install/), which sets up Julia and its version manager, [juliaup](https://github.com/JuliaLang/juliaup), in one step. Then open a terminal, type `julia`, and follow the steps below.

Install MagneticGeophysics.jl into its own environment (press `]` to enter the package manager):

```julia
julia> ]
pkg> activate @magneticgeophysics
pkg> add https://github.com/JuliaGeophysics/MagneticGeophysics.jl
```

Then load it in a Julia session or script:

```julia
julia> using MagneticGeophysics
```

In a new session, run `pkg> activate @magneticgeophysics` again before `using MagneticGeophysics`.

## Examples

The KDISynthetic example: a synthetic airborne survey over the Koillismaa Deep Intrusion, Finland.

```
julia --project=. --threads=8 examples/KDISynthetic/prepare.jl
julia --project=. --threads=8 examples/KDISynthetic/run_KDISynthetic.jl
julia --project=. examples/KDISynthetic/plot_KDISynthetic.jl
```

Contributions are welcome, see [CONTRIBUTING.md](CONTRIBUTING.md).

> [!NOTE]
> **Also check out** other JuliaGeophysics packages: [MTGeophysics.jl](https://github.com/JuliaGeophysics/MTGeophysics.jl) for magnetotellurics, [Timekeepers.jl](https://github.com/JuliaGeophysics/Timekeepers.jl) for time-series analysis in electromagnetic geophysics, and [TerraScope.jl](https://github.com/JuliaGeophysics/TerraScope.jl) for looking at diverse geoscientific data together.
