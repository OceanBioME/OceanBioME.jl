@kwdef struct NewtonRaphsonSolver{FT, IT}
    max_iters :: IT = 100
         atol :: FT = 10^-10
end

@inline function (nrs::NewtonRaphsonSolver)(f, f′, x0, params)
    x = x0
    x⁻ = Inf
    N = 0

    while (abs(x - x⁻) > nrs.atol) & (N < nrs.max_iters)
        x⁻ = x
        x -= step(f, f′, x, params)
        N += 1
    end

    return x
end

@inline step(f, f′, x, params) = f(x, params) / f′(x, params)
