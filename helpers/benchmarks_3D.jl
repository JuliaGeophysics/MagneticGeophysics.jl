# 3D magnetic benchmark case for the examples: a regional synthetic at the Koillismaa Deep Hole, Kuusamo, Finland
# Author: @pankajkmishra
# Coordinates: ETRS-TM35FIN (EPSG:3067). The frame is that grid (x northing, y easting, m; quadratic fit to the
# projection, rms 2 cm), so models, maps and sections are in EPSG:3067; grid north is 1.6° west of true north there.
# The survey is centred on the hole's collar, N 7 290 550, E 579 780 (Contrib. Mineral. Petrol. 2025, doi:10.1007/s00410-025-02226-4).
# Stations every 1 km on a 21 × 21 grid over 20 × 20 km, 500 m above flat ground at 250 m a.s.l. (a regional
# compilation). Inducing field IGRF-14 at the collar, 2025-01-01 (BGS): F = 53 861 nT, I = 76.87°, D = 14.09°.
# The mesh comes from MakeMesh3D: 500 m cells, 250 m layers, 10 km deep core.
#
# Two blocks, positions relative to the collar (northing n, easting e, depth z below ground):
#   deep     4 × 4 km in plan, centred under the collar, 1.4–3.4 km deep (the hole met ultramafic rock at ~1.4 km),
#            κ 0.04 SI
#   shallow  2 × 2 km, centred at n = +5, e = +5 km, 0.5–1.5 km deep, κ 0.03 SI
# Gaussian noise of 2 % of the largest |anomaly| (seeded). Writes examples/synthetic/: model.start and model.prior
# (κ = 0), model.true, model.mask (earth cells) and data.dat. Generated, not tracked.
# Usage: julia --project=. helpers/benchmarks_3D.jl [output_dir]

using MagneticGeophysics
using Printf
using Random
import Proj

const OUT = mkpath(length(ARGS) >= 1 ? ARGS[1] : joinpath(dirname(@__DIR__), "examples", "synthetic"))
const N0, E0 = 7_290_550.0, 579_780.0                                   # Koillismaa Deep Hole collar, EPSG:3067
const GROUND, CLEARANCE = 250.0, 500.0                                  # m a.s.l., m above ground
const FIELD = (53_861.0, 76.868, 14.093)                                # IGRF-14, 2025-01-01, F (nT), I, D (°)

elevation(x, y) = GROUND                                                 # flat ground (m a.s.l.)

# the blocks at northing n, easting e (m from the collar) and depth z below the ground (m)
deep(n, e, z) = abs(n) <= 2000 && abs(e) <= 2000 && 1400 <= z <= 3400
shallow(n, e, z) = abs(n - 5000) <= 1000 && abs(e - 5000) <= 1000 && 500 <= z <= 1500

#---------- stations: a regular EPSG:3067 grid, imported as survey data ----------

g = -10_000:1000:10_000
E = [E0 + e for _ in g for e in g]; N = [N0 + n for n in g for _ in g]
csv = joinpath(OUT, "stations.csv")
open(csv, "w") do io
    println(io, "line,easting,northing,elevation,value,error")
    for (k, (e, n)) in enumerate(zip(E, N))
        @printf(io, "L%02d,%.1f,%.1f,%.1f,0,1\n", cld(k, length(g)), e, n, GROUND + CLEARANCE)
    end
end
tr = Proj.Transformation("EPSG:3067", "EPSG:4326"; always_xy = true)
ll = [tr((e, n)) for (e, n) in zip(E, N)]
pad(v) = (minimum(v) - 0.05, maximum(v) + 0.05)
frame = Frame("EPSG:3067", pad(last.(ll)), pad(first.(ll)))
d = load_data_mag(csv; crs = "EPSG:3067", frame, field = FIELD,
                  columns = (easting = "easting", northing = "northing", elevation = "elevation", value = "value",
                             error = "error", line = "line"))
rm(csv)
d.name = "synthetic"
data_path = write_data_mag(joinpath(OUT, "data.dat"), d; comment = "stations of the Koillismaa synthetic")

#---------- mesh and true model ----------

mesh = MakeMesh3D(data_path; mode = :nogui, out_model = joinpath(OUT, "model.start"), out_mask = joinpath(OUT, "model.mask"),
                  surface = elevation, first_layer_div = 4.0)                # 250 m layers under 1 km stations
start = load_model(joinpath(OUT, "model.start"))
truth = deepcopy(start)
for i in 1:truth.nx, j in 1:truth.ny, k in 1:truth.nz
    isnan(truth.A[i, j, k]) && continue
    n, e = truth.cx[i] - N0, truth.cy[j] - E0
    z = truth.cz[k] + elevation(truth.cx[i], truth.cy[j])               # depth below the local ground
    truth.A[i, j, k] = deep(n, e, z) ? 0.04 : shallow(n, e, z) ? 0.03 : 0.0
end
check_frame(d, start)
op = MagForward(truth, d)
clean = forward(op, model_vector(op, truth))
σ = 0.02 * maximum(abs, clean)
d.value .= clean .+ σ .* randn(Xoshiro(1), ndata(d))
d.error .= σ

write_model(joinpath(OUT, "model.prior"), start; comment = "prior (reference) model, κ = 0")
write_model(joinpath(OUT, "model.true"), truth; comment = "true model of the Koillismaa synthetic, two blocks")
write_data_mag(joinpath(OUT, "data.dat"), d; comment = "synthetic TMI at the Koillismaa Deep Hole, EPSG:3067, 2 % noise")
@printf("mesh %d × %d × %d (%d regular layers), %d stations, anomaly %.1f to %.1f nT, σ = %.2f nT, frame %s → %s\n",
        start.nx, start.ny, start.nz, mesh.nz_regular, ndata(d), extrema(clean)..., σ, frame.source, OUT)
