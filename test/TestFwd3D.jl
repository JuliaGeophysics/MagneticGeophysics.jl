# Forward operator: dense = matrix-free, adjoint, amplitude Jacobian, model vectors, frame check.

# 6 × 5 × 4 mesh with topography (two air cells) and a random model
function small_problem(c; seed = 1)
    m = Model([200.0, 50, 50, 50, 50, 200], [200.0, 50, 50, 50, 200], [25.0, 25, 50, 100])
    m.A[1:2, 1, 1] .= NaN
    m = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin)
    x = range(-150, 150; length = 7); y = range(-100, 100; length = 5)
    d = station_data(c, [a for a in x for _ in y], [b for _ in x for b in y], fill(-30.0, length(x) * length(y)); saz = 12.0)
    κ = 0.01 .+ 0.02 .* abs.(sin.(seed .* (1:count(active_cells(m)))))
    m, d, κ
end

@testset "datum moves the mesh, not the stations" begin
    m, d, κ = small_problem(:tmi)
    s = Model(m.dx, m.dy, m.dz, m.A; origin = m.origin .+ [0, 0, 120.0], datum = 120.0)   # the same earth
    a, b = MagForward(m, d), MagForward(s, d)
    @test b.zs ≈ d.z .+ 120
    @test forward(b, κ) ≈ forward(a, κ) rtol = 1e-12
    opt = Inv3DOptions(weighting = :depth)
    var = 1:length(a.cells)
    @test MG._weights(b, nothing, var, opt) ≈ MG._weights(a, nothing, var, opt)
end

@testset "dense = matrix-free, adjoint" begin
    for c in COMPONENTS
        m, d, κ = small_problem(c)
        D, F = MagForward(m, d), MagForward(m, d; solver = :matrixfree)
        v = cos.(1:length(κ)); w = sin.(1:ndata(d))
        @test isapprox(forward(D, κ), forward(F, κ); rtol = 1e-12)
        @test isapprox(jvec(D, κ, v), jvec(F, κ, v); rtol = 1e-12)
        @test isapprox(jtvec(D, κ, w), jtvec(F, κ, w); rtol = 1e-12)
        for op in (D, F)
            Jv = jvec(op, κ, v)
            @test abs(dot(Jv, w) - dot(v, jtvec(op, κ, w))) / (norm(Jv) * norm(w)) < 1e-12
        end
    end
end

@testset "amplitude Jacobian (Taylor test)" begin
    m, d, κ = small_problem(:amplitude)
    op = MagForward(m, d)
    v = cos.(1:length(κ)) .* 0.01
    J = jvec(op, κ, v); d0 = forward(op, κ)
    r(h) = norm(forward(op, κ .+ h .* v) .- d0 .- h .* J)
    @test r(1e-2) / r(1e-3) > 90          # second order: ×100 per decade
    @test r(1e-3) < 1e-6 * norm(d0)
end

@testset "linear data and model vectors" begin
    m, d, κ = small_problem(:tmi)
    op = MagForward(m, d)
    @test forward(op, 2κ) ≈ 2forward(op, κ)
    @test isequal(isnan.(model_array(op, κ)), air_cells(m))
    mm = Model(m.dx, m.dy, m.dz, model_array(op, κ); origin = m.origin)
    @test model_vector(op, mm) == κ
    @test_throws ArgumentError forward(op, κ[1:end-1])
    @test_throws ArgumentError MagForward(m, rotate_data(d, 5.0))
end
