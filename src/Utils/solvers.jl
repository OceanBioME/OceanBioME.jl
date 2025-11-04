@kwdef struct NewtonRaphsonSolver{FT, IT}
    max_iters :: IT = 1000
         atol :: FT = 10^-10
end

@inline function (nrs::NewtonRaphsonSolver)(f, f′, x0, params; warn_fail = false)
    x = x0
    x⁻ = Inf
    N = 0
    fx = f(x, params)

    while (abs(fx) > nrs.atol) & (N < nrs.max_iters)
        fx = f(x, params)
        x⁻ = x
        x -= fx / f′(x, params)
        N += 1
    end

    (warn_fail && (abs(fx) > nrs.atol)) && @warn "Failed to converge"

    return x
end

@kwdef struct NewtonRaphsonSafeSolver{FT, IT}
    max_iters :: IT = 100
         atol :: FT = 10^-10
end

@inline function (nrss::NewtonRaphsonSafeSolver)(f, f′, x0, params; warn_fail = false)
    x⁻ = minimum(x0)
    x⁺ = maximum(x0)
    xⁿ = (x⁻ + x⁺)/2
    N = 0

    while (abs(f((x⁻ + x⁺)/2, params)) > nrss.atol) & (N < nrss.max_iters)
        xⁿ = (x⁻ + x⁺)/2
        xⁿ⁺¹ = xⁿ - f(xⁿ, params) / f′(xⁿ, params)
          
        x⁻, x⁺ = bisection(x⁺, x⁻, ifelse(x⁻ < xⁿ⁺¹ < x⁺, xⁿ⁺¹, xⁿ), f, params)

        xⁿ = (x⁻ + x⁺)/2
        N += 1
    end

    (warn_fail && (abs(f(xⁿ, params)) > nrss.atol)) && @warn "Failed to converge"

    return xⁿ
end

@kwdef struct Bisection{FT, IT}
    max_iters :: IT = 100
         atol :: FT = 10^-10
      damping :: FT = 1.0
end

@inline (bs::Bisection)(f, x0, params; kwargs...) = bs(f, nothing, x0, params; kwargs...)
@inline function (bs::Bisection)(f, f′, x0, params; warn_fail = false)
    x⁻, x⁺ = x0
    xⁿ = (x⁻ + x⁺)/2
    N = 0

    while (abs(f(xⁿ, params)) > bs.atol) & (N < bs.max_iters)
        x⁻, x⁺ = bisection(x⁺, x⁻, xⁿ, f, params)

        xⁿ = (x⁻ + x⁺)/2
        N += 1
    end

    (warn_fail && (abs(f(xⁿ, params)) > bs.atol)) && @warn "Failed to converge"

    return xⁿ
end

@inline function bisection(x⁺, x⁻, x, f, params)
    fx = f(x, params)

    return (ifelse((sign(f(x⁻, params)) == sign(fx)) | (fx == 0), x, x⁻),
            ifelse((sign(f(x⁺, params)) == sign(fx)) | (fx == 0), x, x⁺))
end

@kwdef struct DampedNewtonRaphsonSolver{FT, IT, BO}
    max_iters :: IT = 100
         atol :: FT = 10^-20
      damping :: FT = 0.5
       bounds :: BO = (lower = -Inf, upper = Inf)
end

@inline function (dnrs::DampedNewtonRaphsonSolver)(f, f′, x0, params; warn_fail = false)
    x = x0
    x⁻ = Inf
    N = 0
    fx = f(x, params)

    while (abs(fx) > dnrs.atol) & (N < dnrs.max_iters)
        fx = f(x, params)
        x⁻ = x
        δ = fx / f′(x, params)
        λ = 1
        while ((abs(f(x - λ * δ, params)) >= abs(fx)) | !(dnrs.bounds.lower <= x - λ * δ <= dnrs.bounds.upper)) & (λ > dnrs.damping^10)
            λ *= dnrs.damping
        end
        x -= λ * δ
        N += 1
    end

    (warn_fail && (abs(fx) > dnrs.atol)) && @warn "Failed to converge"

    return x
end