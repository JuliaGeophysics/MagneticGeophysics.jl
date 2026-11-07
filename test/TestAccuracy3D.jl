# Solver and gradient accuracy. Each check states its tolerance and prints the measured error.

using HCubature

relmax(a, b) = maximum(abs, a .- b) / maximum(abs, b)

# induced dipole field (NED) of moment κ V along f at R = station − source, and its z derivative
function dipole(F, f, κV, R)
    r = norm(R)
    κV * F / 4π * (3 * dot(f, R) .* R ./ r^5 .- f ./ r^3)
end
function dipole_dz(F, f, κV, R)
    r = norm(R)
    T = [3(i == 3) * R[j] + 3R[i] * (j == 3) for i in 1:3, j in 1:3] ./ r^5 .- 15 .* (R * R') .* R[3] ./ r^7 .+
        3 .* [i == j for i in 1:3, j in 1:3] .* R[3] ./ r^5
    κV * F / 4π * (T * f)
end

# value of each data type from a field vector B and its z derivative G (NED)
function reference(d, B, G)
    f = collect(field_direction(d))
    c = d.component
    c == :tmi ? dot(f, B) : c == :amplitude ? norm(B) : c == :bz ? B[3] :
    c in (:bx, :by) ? dot(collect(sensor_direction(d)), B) : c == :dbz_dz ? G[3] : c == :dtmi_dz ? dot(f, G) :
    error("no reference for $c")
end

# 10 m cube, κ = 0.1, at 100–110 m depth, in a mesh padded with zero-κ cells
cube_mesh() = Model([1000.0, 10, 1000], [1000.0, 10, 1000], [100.0, 10, 100]; origin = [-1005.0, -1005.0, 0.0])
cube_κ() = (κ = zeros(27); κ[LinearIndices((3, 3, 3))[2, 2, 2]] = 0.1; κ)

@testset "4.1 cube vs dipole (far field)" begin
    m, κ = cube_mesh(), cube_κ()
    x = [400.0, -300, 50, 0]; y = [100.0, 250, -500, 0]; z = [-50.0, -20, -80, -350]
    for c in (:tmi, :bx, :by, :bz, :amplitude, :dtmi_dz, :dbz_dz)
        d = station_data(c, x, y, z; saz = 7.0)
        f = collect(field_direction(d))
        ref = [reference(d, dipole(d.F, f, 1000 * 0.1, [x[k], y[k], z[k] - 105]),
                         dipole_dz(d.F, f, 1000 * 0.1, [x[k], y[k], z[k] - 105])) for k in eachindex(x)]
        # a cube has no quadrupole: the difference is the l = 4 term, (a/r)⁴ ≈ 1e-7 here
        report("cube vs dipole, $c", relmax(forward(MagForward(m, d), κ), ref), 1e-6)
    end
    d = station_data(:rtp, x, y, z)
    ref = [dipole(d.F, [0, 0, 1.0], 100.0, [x[k], y[k], z[k] - 105])[3] for k in eachindex(x)]
    report("cube vs dipole, rtp", relmax(forward(MagForward(m, d), κ), ref), 1e-6)
end

@testset "4.2 near field vs numerical integration" begin
    # prism column sticks up out of the topography: air down to 110 m beside it
    m = Model([100.0, 10, 100], [100.0, 10, 100], [100.0, 10, 100]; origin = [-100.0, -100.0, 0.0])
    m.A[:, :, 1:2] .= NaN; m.A[2, 2, 2] = 0.0
    m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    op_cells = findall(active_cells(m)); κ = [I == CartesianIndex(2, 2, 2) ? 0.1 : 0.0 for I in op_cells]
    stations = [(5.0, 5.0, 99.0, "1 m above the centre"), (0.0, 5.0, 99.0, "1 m above an edge"),
                (0.0, 0.0, 99.0, "1 m above a corner"), (-1.0, 5.0, 105.0, "1 m beside, below the top"),
                (12.0, 13.0, 95.0, "off a corner, above")]
    for c in (:tmi, :bx, :bz, :dtmi_dz, :dbz_dz)
        for (x, y, z, label) in stations
            d = station_data(c, [x], [y], [z])
            f = collect(field_direction(d))
            B = hcubature(p -> dipole(d.F, f, 0.1, [x, y, z] .- p), [0.0, 0, 100], [10.0, 10, 110]; rtol = 1e-11, maxevals = 10^7)[1]
            G = c in (:dtmi_dz, :dbz_dz) ?
                hcubature(p -> dipole_dz(d.F, f, 0.1, [x, y, z] .- p), [0.0, 0, 100], [10.0, 10, 110]; rtol = 1e-11, maxevals = 10^7)[1] : B
            ref = reference(d, B, G)
            got = forward(MagForward(m, d), κ)[1]
            report("$c, $label", abs(got - ref) / abs(ref), 1e-8)
        end
    end
end

@testset "4.3 superposition and limits" begin
    one = Model([1000.0, 10, 1000], [1000.0, 10, 1000], [100.0, 10, 100]; origin = [-1005.0, -1005.0, 0.0])
    split = Model([1000.0, 5, 5, 1000], [1000.0, 5, 5, 1000], [100.0, 5, 5, 100]; origin = [-1005.0, -1005.0, 0.0])
    κs = [I[1] in 2:3 && I[2] in 2:3 && I[3] in 2:3 ? 0.1 : 0.0 for I in CartesianIndices((4, 4, 4))]
    x = [0.0, 30, -7]; y = [0.0, -12, 9]; z = [-1.0, -40, -3]
    for c in (:tmi, :by, :dtmi_dz)
        d = station_data(c, x, y, z; saz = 3.0)
        report("one prism = 8 sub-prisms, $c", relmax(forward(MagForward(split, d), vec(κs)), forward(MagForward(one, d), cube_κ())), 1e-12)
    end
    m, κ = cube_mesh(), cube_κ()
    a = forward(MagForward(m, station_data(:rtp, x, y, z; I = 30.0, D = 40.0)), κ)
    b = forward(MagForward(m, station_data(:tmi, x, y, z; I = 90.0, D = 0.0)), κ)
    report("RTP = TMI at the pole", relmax(a, b), 1e-14)
    report("zero model gives zero data", maximum(abs, forward(MagForward(m, station_data(:tmi, x, y, z)), zeros(27))), 0.0)
end

@testset "4.4 round-off: far stations and thin cells vs BigFloat" begin
    setprecision(BigFloat, 256) do
        function prism(T, order, st, lo, hi, f)
            x0, y0, z0 = T.(st); s = zero(T)
            for a in 0:1, b in 0:1, c in 0:1
                n = T(a == 1 ? hi[1] : lo[1]) - x0; e = T(b == 1 ? hi[2] : lo[2]) - y0; u = z0 - T(c == 0 ? lo[3] : hi[3])
                r = sqrt(e^2 + n^2 + u^2)
                K = order == 2 ? MG._kernel2(e, n, u, r) : MG._kernel3u(e, n, u, r)
                kee, knn, kuu, ken, keu, knu = K
                s += (-1)^(a + b + c) * (f[1] * (kee * f[1] + ken * f[2] + keu * f[3]) + f[2] * (ken * f[1] + knn * f[2] + knu * f[3]) +
                                         f[3] * (keu * f[1] + knu * f[2] + kuu * f[3]))
            end
            s
        end
        f = MG._enu((cosd(65) * cosd(20), cosd(65) * sind(20), sind(65)))
        # cancellation between the two z-corners grows as cells thin: about 1e-10 at 1 m
        # thickness, 1.6e-8 for a 1 cm stress case; both far below any data error
        for order in (2, 3), (lo, hi, label, tol) in (((-5.0, -5.0, 100.0), (5.0, 5.0, 110.0), "10 m cell", 1e-10),
                                                       ((-5.0, -5.0, 100.0), (5.0, 5.0, 101.0), "1 m thick cell", 5e-10),
                                                       ((-5.0, -5.0, 100.0), (5.0, 5.0, 100.01), "1 cm thick cell (stress)", 1e-7))
            near = abs(prism(BigFloat, order, (0.0, 0.0, -100.0), lo, hi, f))
            for R in (5_000.0, 50_000.0, 200_000.0)
                st = (0.6R, 0.8R, -100.0)
                err = Float64(abs(prism(Float64, order, st, lo, hi, f) - prism(BigFloat, order, st, lo, hi, f)) / near)
                report("order $order, $label at $(Int(R ÷ 1000)) km, err / response above", err, tol)
            end
        end
    end
end

@testset "4.5 vertical gradients vs central differences" begin
    m, κ = cube_mesh(), cube_κ()
    x = [3.0, 40, -25]; y = [2.0, -30, 15]; z = [-2.0, -20, -60]
    for (g, base) in ((:dbz_dz, :bz), (:dtmi_dz, :tmi))
        exact = forward(MagForward(m, station_data(g, x, y, z)), κ)
        fd(h) = (forward(MagForward(m, station_data(base, x, y, z .+ h)), κ) .-
                 forward(MagForward(m, station_data(base, x, y, z .- h)), κ)) ./ 2h
        e = [relmax(fd(h), exact) for h in (1.0, 0.1, 0.01)]
        report("$g: central difference h = 1 m", e[1], 1e-1)
        report("$g: h² convergence, error(1 m) / error(0.1 m) ≥ 50 → 1/ratio", e[2] / e[1], 2e-2)
        report("$g: h = 1 cm", e[3], 1e-5)
    end
end

@testset "4.9 rotation: the same body on a 90°-turned frame" begin
    # mesh B is mesh A seen from a frame turned 90° clockwise: x' = y, y' = −x
    A = Model([100.0, 20, 30, 20, 100], [100.0, 10, 20, 40, 100], [50.0, 20, 40]; origin = [-135.0, -135.0, 0.0])
    κA = [0.01 * (I[1] + 2I[2] + 3I[3]) for I in CartesianIndices(size(A.A))]
    B = Model(A.dy, reverse(A.dx), A.dz, [κA[A.nx + 1 - j, i, k] for i in 1:A.ny, j in 1:A.nx, k in 1:A.nz];
              origin = [A.origin[2], -(A.origin[1] + sum(A.dx)), A.origin[3]], rotation = 90.0)
    x = [10.0, -40, 33]; y = [5.0, 20, -27]; z = [-5.0, -15, -30]
    for c in (:tmi, :bx, :by, :amplitude, :dtmi_dz)
        d = station_data(c, x, y, z; saz = 15.0)
        a = forward(MagForward(A, d), vec(κA))
        b = forward(MagForward(B, rotate_data(d, 90.0)), vec(B.A))
        report("rotated frame, $c", relmax(b, a), 1e-12)
    end
end

@testset "4.10 air cells above the earth" begin
    flat = Model([200.0, 50, 50, 200], [200.0, 50, 50, 200], [30.0, 60]; origin = [-250.0, -250.0, 0.0])
    κ = [0.01 * (I[1] + I[2] * I[3]) for I in CartesianIndices((4, 4, 2))]
    air = Model(flat.dx, flat.dy, [40.0, 30, 60], cat(fill(NaN, 4, 4, 1), κ; dims = 3); origin = [-250.0, -250.0, -40.0])
    d = station_data(:tmi, [10.0, -60], [5.0, 70], [-45.0, -50])
    report("air layer on top changes nothing", relmax(forward(MagForward(air, d), vec(κ)), forward(MagForward(flat, d), vec(κ))), 1e-13)
end

@testset "4.11 SimPEG benchmark (cross-gradient tutorial TMI)" begin
    pad = [5.0 * 1.3^k for k in 1:5]
    h = vcat(reverse(pad), fill(5.0, 40), pad)
    v = parse.(Float64, split(read(joinpath(DATA, "simpeg_true_susc.txt"), String)))
    A = permutedims(reshape(v, 50, 50, 20), (2, 1, 3))[:, :, end:-1:1]     # east fastest, z up → [north, east, depth]
    nair = permutedims(reduce(hcat, [parse.(Int, split(l)) for l in eachline(joinpath(DATA, "simpeg_air_layers.txt")) if !startswith(l, "#")]))
    for j in 1:50, i in 1:50; A[i, j, 1:nair[i, j]] .= NaN; end            # nair[north, east]
    m = Model(h, h, vcat(fill(5.0, 15), pad), A; origin = [-sum(h) / 2, -sum(h) / 2, 0.0])
    obs = reduce(hcat, [parse.(Float64, split(l)) for l in eachline(joinpath(DATA, "simpeg_magnetic_data.obs"))])
    d = station_data(:tmi, obs[2, :], obs[1, :], -obs[3, :]; I = 90.0, D = 0.0)
    op = MagForward(m, d)
    pred = forward(op, model_vector(op, m))
    c = sum(obs[4, :] .- pred) / ndata(d)
    @printf("  SimPEG data − ours: constant %.4f nT (%.2f %% of max |d|)\n", c, 100c / maximum(abs, obs[4, :]))
    report("SimPEG TMI after the constant, max |r| / max |d|", relmax(pred .+ c, obs[4, :]), 1e-3)
end
