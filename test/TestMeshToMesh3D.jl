# MeshToMesh: a susceptibility model onto a finer mesh with its own topography; mask and data files.

@testset "MeshToMesh" begin
    dir = mktempdir()
    src = Model(fill(100.0, 8), fill(100.0, 8), fill(50.0, 6); origin = [-400.0, -400.0, 0.0])
    src.A[5:8, :, :] .= 0.02                                     # northern half magnetic
    sp = write_model(joinpath(dir, "coarse.sus"), src)
    d = station_data(:tmi, [-100.0, 100.0], [0.0, 0.0], [-30.0, -30.0])
    dp = write_data_mag(joinpath(dir, "obs.dat"), d)
    # target: a resistivity mesh, 50 m cells, top layer air over the west half and sea water east of y = 200 m
    ρ = fill(100.0, 16, 16, 12); ρ[:, 1:8, 1] .= 1e17; ρ[:, 13:16, 1] .= 0.3
    open(joinpath(dir, "fine.rho"), "w") do io
        println(io, "16 16 12 0 LOGE"); println(io, join(fill(50.0, 16), " ")); println(io, join(fill(50.0, 16), " "))
        println(io, join(fill(25.0, 12), " "))
        for k in 1:12, j in 1:16
            println(io, join((log(ρ[i, j, k]) for i in 16:-1:1), " "))
        end
        println(io, "-400 -400 0"); println(io, "0")
    end
    r = MeshToMesh(sp, dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "out"), water_log10 = 0.0)
    m = load_model(r.model)
    @test (m.nx, m.ny, m.nz) == (16, 16, 12) && r.dst_grid == (16, 16, 12)
    @test all(isnan, m.A[:, 1:8, 1]) && all(isnan, m.A[:, 13:16, 1])       # target air and sea water stay out
    @test all(m.A[9:16, :, 2:end] .== 0.02) && all(m.A[1:8, :, 2:end] .== 0.0)
    @test r.n_air == 16 * 12 && r.cells_outside_source == 0
    @test load_mask(r.mask, m) == .!isnan.(m.A)
    @test load_data_mag(r.data).z == d.z                                     # heights kept
    lin = MeshToMesh(sp, dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "lin"), method = :linear, snap_sites = true,
                     clearance = 5.0)
    ml = load_model(lin.model)
    @test 0 < ml.A[8, 9, 6] < 0.02 && ml.A[12, 9, 6] == 0.02                  # interpolated across the contact
    dl = load_data_mag(lin.data)
    @test all(dl.z .≈ 25.0 - 5.0) || all(isapprox.(dl.z, ml.Z[MG._column_index(ml, d.x[1], d.y[1])...] - 5.0))
    @test_throws ErrorException MeshToMesh(sp, dp, joinpath(dir, "missing.rho"))
    # a source with a datum keeps each value at its elevation
    lay = Model(src.dx, src.dy, src.dz, src.A .+ reshape(0.001 .* (1:6), 1, 1, 6); origin = src.origin)
    up = Model(lay.dx, lay.dy, lay.dz, lay.A; origin = lay.origin .+ [0, 0, 100.0], datum = 100.0)
    a = MeshToMesh(write_model(joinpath(dir, "lay.sus"), lay), dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "a"))
    b = MeshToMesh(write_model(joinpath(dir, "up.sus"), up), dp, joinpath(dir, "fine.rho"); outdir = joinpath(dir, "b"))
    @test isequal(load_model(b.model).A, load_model(a.model).A) && b.cells_outside_source == 0
    # a target turned 90°: a body 300 m north of the centre lies at y′ = −300 m in its frame, and so do the stations
    body = Model(src.dx, src.dy, src.dz; origin = src.origin)
    body.A[7:8, 4:5, :] .= 0.02                                  # x 200–400 m (north), y −100–100 m
    turned = Model(src.dx, src.dy, src.dz; origin = src.origin, rotation = 90.0)
    rt = MeshToMesh(write_model(joinpath(dir, "body.sus"), body), dp, write_model(joinpath(dir, "turned.sus"), turned);
                    outdir = joinpath(dir, "turned"))
    mt = load_model(rt.model)
    @test all(mt.A[4:5, 1:2, :] .== 0.02) && count(>(0), mt.A) == 2 * 2 * 6   # x′ −100–100 m, y′ −400–−200 m
    dt = load_data_mag(rt.data)
    @test dt.azimuth == 90.0 && dt.x ≈ d.y && dt.y ≈ -d.x
end
