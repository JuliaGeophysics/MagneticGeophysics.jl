# 3D magnetic forward modelling with prisms
# Author: @pankajkmishra
# Closed-form field of uniformly magnetised right rectangular prisms (Nagy et al. 2000; Fukushima
# 2020), with the kernels and their safe logarithm/arctangent as in choclo (Soler & Uieda 2021,
# BSD-3): B = C_m Σ_vertices (−1)^(i+j+k) K(vertex − station) M, K the second derivatives of ∫1/r
# (third derivatives with one upward index for the vertical gradients). The kernels are written
# east-north-up; stations and results are x north, y east, z DOWN. Each kernel depends on the
# vertex only, so it is evaluated once per mesh node and summed into the cells sharing that node.
# Induced magnetisation along the inducing field, M = κ F/μ₀, so C_m M = κ F/4π: F in nT gives nT.
# One row routine serves both solvers, so the dense G and the matrix-free products are identical.

#---------- kernels (east, north, up; vertex − station) ----------

@inline function _safe_atan(y, x)
    x == 0 && return y > 0 ? π / 2 : y < 0 ? -π / 2 : 0.0
    atan(y / x)
end

# log(x + r) without cancellation for x < 0
@inline function _safe_log(x, y, z, r)
    r == 0 && return 0.0
    x >= 0 && return log(x + r)
    y == 0 && z == 0 && return -log(-2x)
    log((y^2 + z^2) / (r - x))
end

@inline function _k_iii(a, b, c, r)
    ((a == 0 && b == 0) || (a == 0 && c == 0)) && return 0.0
    a2, b2, c2 = a^2, b^2, c^2
    -b * c * (2a2 + b2 + c2) / ((a2 + b2) * (a2 + c2) * r)
end

@inline _k_iij(a, b, c, r) = a == 0 && b == 0 ? 0.0 : a * c / ((a^2 + b^2) * r)

# second derivatives: symmetric 3 × 3 (ee, nn, uu, en, eu, nu)
@inline function _kernel2(e, n, u, r)
    (-_safe_atan(n * u, e * r), -_safe_atan(e * u, n * r), -_safe_atan(e * n, u * r),
     _safe_log(u, e, n, r), _safe_log(n, e, u, r), _safe_log(e, n, u, r))
end

# third derivatives with one upward index: (eeu, nnu, uuu, enu, euu, nuu)
@inline function _kernel3u(e, n, u, r)
    (_k_iij(e, u, n, r), _k_iij(n, u, e, r), _k_iii(u, n, e, r),
     -1 / r, _k_iij(u, e, n, r), _k_iij(u, n, e, r))
end

#---------- linear functionals ----------

# A datum is σ wᵀ K f per node (w, f east-north-up), from second (order 2) or third (order 3)
# derivative kernels. NED → ENU: (n, e, d) ↦ (e, n, −d).
_enu(v) = (v[2], v[1], -v[3])

struct Functional
    w::NTuple{3,Float64}
    f::NTuple{3,Float64}
    order::Int
    sign::Float64
end

# the functionals of a data set: one, or three (Bx, By, Bz along the frame) for the amplitude
function functionals(d::Data)
    f = _enu(field_direction(d))
    c = d.component
    c == :tmi && return [Functional(f, f, 2, 1.0)]
    c in (:bx, :by, :bz) && return [Functional(_enu(sensor_direction(d)), f, 2, 1.0)]
    c == :rtp && return [Functional((0.0, 0.0, -1.0), (0.0, 0.0, -1.0), 2, 1.0)]
    # ∂/∂z_down = −∂/∂u: dBz/dz = ∂b_u/∂u, dTMI/dz = −f·∂b/∂u
    c == :dbz_dz && return [Functional((0.0, 0.0, 1.0), f, 3, 1.0)]
    c == :dtmi_dz && return [Functional(f, f, 3, -1.0)]
    c == :amplitude && return [Functional(_enu(e), f, 2, 1.0) for e in ((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0))]
    throw(ArgumentError("no functional for component $c"))
end

#---------- rows ----------

# S[node] = σ wᵀ K(node − station) f for station (x0, y0, z0)
function _node_values!(S::Array{Float64,3}, m::Model, x0, y0, z0, q::Functional)
    we, wn, wu = q.w; fe, fn, fu = q.f
    @inbounds for k in 1:m.nz+1
        u = z0 - m.z[k]                       # upward shift: vertex up − station up
        for j in 1:m.ny+1
            e = m.y[j] - y0
            for i in 1:m.nx+1
                n = m.x[i] - x0
                r = sqrt(e^2 + n^2 + u^2)
                kee, knn, kuu, ken, keu, knu = q.order == 2 ? _kernel2(e, n, u, r) : _kernel3u(e, n, u, r)
                S[i, j, k] = q.sign * (we * (kee * fe + ken * fn + keu * fu) +
                                       wn * (ken * fe + knn * fn + knu * fu) +
                                       wu * (keu * fe + knu * fn + kuu * fu))
            end
        end
    end
    S
end

# row of the active cells: Σ over the 8 vertices of (−1)^(a+b+c) S[i+a, j+b, k+c], times F/4π
function _row!(row::AbstractVector, S::Array{Float64,3}, cells::Vector{CartesianIndex{3}}, F::Float64)
    c = F / 4π
    @inbounds for (n, I) in enumerate(cells)
        i, j, k = I[1], I[2], I[3]
        row[n] = c * (S[i+1, j+1, k] - S[i, j+1, k] - S[i+1, j, k] + S[i, j, k] -
                      S[i+1, j+1, k+1] + S[i, j+1, k+1] + S[i+1, j, k+1] - S[i, j, k+1])
    end
    row
end

#---------- operator ----------

"""
    MagForward(m, d; solver=:dense)

Forward operator of the data `d` on the active (non-air) cells of `m`. `solver = :dense` builds and
stores the sensitivity G (nT per SI, or nT/m per SI); `:matrixfree` recomputes the rows for every
product and needs no storage. Both give the same numbers. Checks the data frame against the mesh.
"""
struct MagForward
    m::Model
    d::Data
    cells::Vector{CartesianIndex{3}}
    zs::Vector{Float64}                   # station z on the mesh (d.z + m.datum)
    q::Vector{Functional}
    solver::Symbol
    G::Vector{Matrix{Float64}}            # one per functional; empty for :matrixfree
end

function MagForward(m::Model, d::Data; solver::Symbol = :dense)
    solver in (:dense, :matrixfree) || throw(ArgumentError("solver must be :dense or :matrixfree"))
    check_frame(d, m)
    cells = findall(active_cells(m))
    q = functionals(d)
    op = MagForward(m, d, cells, _mesh_z(d, m), q, solver, Matrix{Float64}[])
    if solver == :dense
        bytes = 8.0 * ndata(d) * length(cells) * length(q)
        bytes > 0.8 * Sys.total_memory() &&
            throw(ArgumentError(@sprintf("dense G needs %.1f GB; use solver = :matrixfree", bytes / 2^30)))
        for k in eachindex(q)
            G = Matrix{Float64}(undef, ndata(d), length(cells))
            _foreach_row(op, k) do n, row
                G[n, :] .= row
            end
            push!(op.G, G)
        end
    end
    op
end

ncells(op::MagForward) = length(op.cells)

# run f(n, row) for every station n with row = functional k, threaded over chunks of stations
function _foreach_row(f, op::MagForward, k::Int)
    m, d = op.m, op.d
    N = ndata(d)
    chunks = collect(Iterators.partition(1:N, cld(N, max(1, Threads.nthreads()))))
    @sync for chunk in chunks
        Threads.@spawn begin
            S = Array{Float64,3}(undef, m.nx + 1, m.ny + 1, m.nz + 1)
            row = Vector{Float64}(undef, length(op.cells))
            for n in chunk
                _node_values!(S, m, d.x[n], d.y[n], op.zs[n], op.q[k])
                f(n, _row!(row, S, op.cells, d.F))
            end
        end
    end
    nothing
end

# G_k v for every functional k
function _gv(op::MagForward, k::Int, v::AbstractVector)
    op.solver == :dense && return op.G[k] * v
    out = zeros(ndata(op.d))
    _foreach_row(op, k) do n, row
        out[n] = dot(row, v)
    end
    out
end

# G_kᵀ w
function _gtw(op::MagForward, k::Int, w::AbstractVector)
    op.solver == :dense && return transpose(op.G[k]) * w
    acc = [zeros(ncells(op)) for _ in 1:max(1, Threads.nthreads())]
    N = ndata(op.d)
    chunks = collect(Iterators.partition(1:N, cld(N, length(acc))))
    m, d = op.m, op.d
    @sync for (t, chunk) in enumerate(chunks)
        Threads.@spawn begin
            S = Array{Float64,3}(undef, m.nx + 1, m.ny + 1, m.nz + 1)
            row = Vector{Float64}(undef, ncells(op))
            a = acc[t]
            for n in chunk
                _node_values!(S, m, d.x[n], d.y[n], op.zs[n], op.q[k])
                _row!(row, S, op.cells, d.F)
                a .+= w[n] .* row
            end
        end
    end
    sum(acc)
end

_check_len(op, v) = length(v) == ncells(op) ||
    throw(ArgumentError("model vector has $(length(v)) values, the operator $(ncells(op)) active cells"))

"""
    forward(op, κ) -> predicted data

`κ` holds the susceptibility (SI) of the active cells, in `findall(active_cells(m))` order.
"""
function forward(op::MagForward, κ::AbstractVector)
    _check_len(op, κ)
    length(op.q) == 1 && return _gv(op, 1, κ)
    sqrt.(sum(_gv(op, k, κ) .^ 2 for k in 1:3))
end

"""
    jvec(op, κ, v) / jtvec(op, κ, w)

Jacobian products at κ. Linear data: J = G. Amplitude a = |B|: J = Σ_c (B_c / a) G_c.
"""
function jvec(op::MagForward, κ::AbstractVector, v::AbstractVector)
    _check_len(op, v)
    length(op.q) == 1 && return _gv(op, 1, v)
    B = [_gv(op, k, κ) for k in 1:3]; a = sqrt.(sum(b .^ 2 for b in B))
    sum((B[k] ./ a) .* _gv(op, k, v) for k in 1:3)
end

function jtvec(op::MagForward, κ::AbstractVector, w::AbstractVector)
    length(w) == ndata(op.d) || throw(ArgumentError("data vector length differs from the data"))
    length(op.q) == 1 && return _gtw(op, 1, w)
    B = [_gv(op, k, κ) for k in 1:3]; a = sqrt.(sum(b .^ 2 for b in B))
    sum(_gtw(op, k, w .* B[k] ./ a) for k in 1:3)
end

"""
    model_vector(op, m) / model_array(op, κ)

Active-cell vector of a model on the operator's mesh, and the full array (NaN in air) of a vector.
"""
model_vector(op::MagForward, m::Model) = (same_mesh(op.m, m) || throw(ArgumentError("model mesh differs")); m.A[op.cells])
function model_array(op::MagForward, κ::AbstractVector)
    _check_len(op, κ)
    A = fill(NaN, size(op.m.A)); A[op.cells] .= κ; A
end

#---------- file workflow ----------

"""
    ForwardSolve3D(model_path, data_path, fwd_path; output_path, add_noise=false, rng_seed=20261017) -> path

Forward run from files: a susceptibility model, a native data file giving the stations, data type,
inducing field and errors, and `FwdCtrl` giving the solver. Writes the predicted data, by default
`data.pred` next to the data file, with the errors of the data file; `add_noise` adds Gaussian
noise of those errors.
"""
function ForwardSolve3D(model_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString;
                        output_path::AbstractString = joinpath(dirname(abspath(data_path)), "data.pred"),
                        add_noise::Bool = false, rng_seed::Integer = 20261017)
    fwd = ReadFwdCtrl3D(fwd_path)
    m, d = load_model(model_path), load_data_mag(data_path)
    op = MagForward(m, d; solver = fwd.solver)
    p = forward(op, model_vector(op, m))
    add_noise && (p .+= d.error .* randn(Xoshiro(rng_seed), length(p)))
    write_data_mag(output_path, d; value = p, comment = "predicted by $(basename(model_path))" * (add_noise ? ", with noise" : ""))
end
