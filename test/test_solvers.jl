using Statistics

using OceanBioME: NewtonRaphsonSolver, 
                  NewtonRaphsonSafeSolver, 
                  DampedNewtonRaphsonSolver, 
                  Bisection

f1(x, p) = tanh(x)
f1′(x, p) = sech(x)^2
X1 = 0

f2(x, p) = x^2-2
f2′(x, p) = 2*x
X2 = (-√2, √2)

f3(x, p) = exp(x) - 1000
f3′(x, p) = exp(x)
X3 = 6.907755278982137

f4(x, p) = 1/(1+x^2)-0.5
f4′(x, p) = -2*x/(1+x^2)^2
X4 = (-1, 1)

f5(x, p) = p.a * x^2 + p.b
f5′(x, p) = 2 * p.a * x

@testset "Testing solvers" begin
    @info "Benchmark times are given for each solver, but should be taken with a pinch of salt as they vary a lot depending on the problem and start point"
    # Plain newton raphson
    solver = NewtonRaphsonSolver()

    # make all this parallel in a kernel to check it will run on a GPU correctly

    correct1 = map(x0 -> solver(f1, f1′, x0, nothing) == X1, (-10, -1, 0, 1, 10, 100))
    @test all(correct1 .== (false, true, true, true, false, false))

    correct2 = (map(x0 -> solver(f2, f2′, x0, nothing) ≈ X2[1], (-5, -1, -0.1, 0))...,
                map(x0 -> solver(f2, f2′, x0, nothing) ≈ X2[2], (0.1, 1, 5))...)
    @test all(correct2 .== (true, true, true, false, true, true, true))

    correct3 = map(x0 -> solver(f3, f3′, x0, nothing) == X3, (-1, 0, 1, 5, 7, 10))
    @test all(correct3 .== (false, false, true, true, true, true))

    correct4 = (map(x0 -> solver(f4, f4′, x0, nothing) ≈ X4[1], (-10, -1.5, -1, 0))...,
                map(x0 -> solver(f4, f4′, x0, nothing) ≈ X4[2], (1, 1.5, 5))...)
    @test all(correct4 .== (false, true, true, false, true, true, false))

    @test solver(f5, f5′, -5, (a = 1, b = -2)) ≈ -√2

    rt = map(n -> (@timed solver(f1, f1′, -1, nothing)).time, 1:100)

    @info "Tested NewtonRaphsonSolver with benchmark time: $(round(mean(rt)*10^9, digits = 1))ns per solve"

    # Newton raphson safe (Newt safe) - "garanteed" convergence, but expensive
    solver = NewtonRaphsonSafeSolver()
    
    correct1 = map(x0 -> isapprox(solver(f1, f1′, x0, nothing), X1, atol = solver.atol), ((-1000, 1000), (-100, 1), (-10, 0)))
    @test all(correct1)

    correct2 = (map(x0 -> isapprox(solver(f2, f2′, x0, nothing), X2[1], atol = solver.atol), ((-1000, 0), (-10, √2 - 0.1)))...,
                isapprox(solver(f2, f2′, (0, 1000), nothing), X2[2], atol = solver.atol))
    @test all(correct2)

    correct3 = solver(f3, f3′, (-10000, 10000), nothing) == X3
    @test correct3

    correct4 = (map(x0 -> isapprox(solver(f4, f4′, x0, nothing), X4[1], atol = solver.atol), ((-1000, 0), (-10, 0.9)))...,
                isapprox(solver(f4, f4′, (0, 1000), nothing), X4[2], atol = solver.atol))
    @test all(correct4)

    @test solver(f5, f5′, (-100, 0), (a = 1, b = -2)) ≈ -√2

    rt = map(n -> (@timed solver(f1, f1′, (-1, 1), nothing)).time, 1:100)

    @info "Tested NewtonRaphsonSafeSolver with benchmark time: $(round(mean(rt)*10^9, digits = 1))ns per solve"

    # DampedNewtonRaphsonSolver
    solver = DampedNewtonRaphsonSolver(; damping = 0.1)

    correct1 = map(x0 -> solver(f1, f1′, x0, nothing) == X1, (-10, -1, 0, 1, 10))
    @test all(correct1)

    correct2 = (map(x0 -> solver(f2, f2′, x0, nothing) ≈ X2[1], (-5, -1, -0.1, 0))...,
                map(x0 -> solver(f2, f2′, x0, nothing) ≈ X2[2], (0.1, 1, 5))...)
    @test all(correct2 .== (true, true, true, false, true, true, true))

    correct3 = map(x0 -> solver(f3, f3′, x0, nothing) == X3, (-1, 0, 1, 5, 7, 10))
    @test all(correct3)

    correct4 = (map(x0 -> solver(f4, f4′, x0, nothing) ≈ X4[1], (-10, -1.5, -1, 0))...,
                map(x0 -> solver(f4, f4′, x0, nothing) ≈ X4[2], (1, 1.5, 5))...)
    @test all(correct4 .== (true, true, true, false, true, true, true))

    @test solver(f5, f5′, -5, (a = 1, b = -2)) ≈ -√2

    solver = DampedNewtonRaphsonSolver()

    rt = map(n -> (@timed solver(f1, f1′, -1, nothing)).time, 1:100)

    @info "Tested DampedNewtonRaphsonSolver with benchmark time: $(round(mean(rt)*10^9, digits = 1))ns per solve"

    # Bisection
    solver = Bisection()
    
    correct1 = map(x0 -> isapprox(solver(f1, f1′, x0, nothing), X1, atol = solver.atol), ((-1000, 1000), (-100, 1), (-10, 0)))
    @test all(correct1)

    correct2 = (map(x0 -> isapprox(solver(f2, f2′, x0, nothing), X2[1], atol = solver.atol), ((-1000, 0), (-10, √2 - 0.1)))...,
                isapprox(solver(f2, f2′, (0, 1000), nothing), X2[2], atol = solver.atol))
    @test all(correct2)

    correct3 = solver(f3, f3′, (-10000, 10000), nothing) ≈ X3
    @test correct3

    correct4 = (map(x0 -> isapprox(f4(solver(f4, f4′, x0, nothing), nothing), 0, atol = solver.atol), ((-1000, 0), (-10, 0.9)))...,
                isapprox(f4(solver(f4, f4′, (0, 1000), nothing), nothing), 0, atol = solver.atol))
    @test all(correct4)

    @test solver(f5, f5′, (-100, 0), (a = 1, b = -2)) ≈ -√2

    rt = map(n -> (@timed solver(f1, f1′, (-10, 1), nothing)).time, 1:100)

    @info "Tested Bisection with benchmark time: $(round(mean(rt)*10^9, digits = 1))ns per solve"
end