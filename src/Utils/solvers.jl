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

"""
    DRTSafeSolver(; max_iters = 100, bracket_grows = 3, atol = 10^-10,
                    H_lower = 10^-9, H_upper = 10^-6)

A bracketed Newton-Raphson solver which falls back to bisection whenever a Newton step would leave the
bracket, or would not at least halve the step. This is a port of MARBL's `drtsafe`
(`marbl_co2calc_mod.F90`, itself a modified Numerical Recipes `rtsafe`) and exists so that carbonate
chemistry can be compared with MARBL's.

Unlike the other solvers here it stops on the *step* size (`|δ| < atol`) rather than on the residual,
and returns the iterate without re-evaluating it, both as MARBL does. At MARBL's tolerance (`xacc =
10^-10`) this deliberately under converges, leaving a residual of order `5×10⁻¹⁰`; that is the point,
since reproducing MARBL means reproducing its convergence rather than bettering it.

`H_lower` and `H_upper` bracket the hydrogen ion concentration and default to MARBL's cold start
bracket of pH 6 to 9, which is used in preference to the initial guess `CarbonChemistry` supplies. The
bracket width matters: too narrow and deep, low pH cells fall outside it. Where the bracket does not
contain a sign change it is widened geometrically up to `bracket_grows` times.

Keyword Arguments
=================
- `max_iters`: maximum number of iterations (MARBL `maxit`)
- `bracket_grows`: maximum number of times the bracket may be widened (MARBL `max_bracket_grow_it`)
- `atol`: absolute tolerance on the step (MARBL `xacc`)
- `H_lower`, `H_upper`: the initial bracket on `[H⁺]`
"""
@kwdef struct DRTSafeSolver{FT, IT}
        max_iters :: IT = 100
    bracket_grows :: IT = 3
             atol :: FT = 10^-10
          H_lower :: FT = 10^-9
          H_upper :: FT = 10^-6
end

@inline initial_bracket(drts::DRTSafeSolver, x0) = (drts.H_lower, drts.H_upper)
@inline initial_bracket(::DRTSafeSolver, x0::Tuple) = x0

@inline function (drts::DRTSafeSolver)(f, f′, x0, params; warn_fail = false)
    x₁, x₂ = initial_bracket(drts, x0)

    f₁ = f(x₁, params)
    f₂ = f(x₂, params)

    # widen the bracket geometrically until the ends straddle the root
    N = 0
    same_sign = ((f₁ > 0) & (f₂ > 0)) | ((f₁ < 0) & (f₂ < 0))

    while same_sign & (N < drts.bracket_grows)
        δ = sqrt(x₂ / x₁)

        x₂ *= δ
        x₁ /= δ

        f₁ = f(x₁, params)
        f₂ = f(x₂, params)

        same_sign = ((f₁ > 0) & (f₂ > 0)) | ((f₁ < 0) & (f₂ < 0))

        N += 1
    end

    (warn_fail && same_sign) && @warn "Failed to bracket the root"

    # `x⁻` is the end where the residual is negative. Alkalinity falls as `[H⁺]` rises, so that is the
    # *larger* concentration and `x⁻ > x⁺`; the bisection step below is negative in consequence
    x⁻, x⁺ = ifelse(f₁ < 0, (x₁, x₂), (x₂, x₁))

    x  = (x⁻ + x⁺) / 2
    δ⁻ = abs(x⁻ - x⁺)
    δ  = δ⁻

    fx  = f(x, params)
    f′x = f′(x, params)

    for _ in 1:drts.max_iters
        leaves_bracket  = ((x - x⁺) * f′x - fx) * ((x - x⁻) * f′x - fx) >= 0
        step_decreases  = abs(2fx) <= abs(δ⁻ * f′x)

        xᵖ = x

        if leaves_bracket | !step_decreases
            δ⁻, δ = δ, (x⁺ - x⁻) / 2

            x = x⁻ + δ

            stalled = x⁻ == x
        else
            δ⁻, δ = δ, -fx / f′x

            x += δ

            stalled = xᵖ == x
        end

        # a step which does not move the iterate, or which is smaller than the tolerance, ends the
        # solve. The iterate is returned *without* re-evaluating it, so the residual is left at the
        # order of `atol * |∂f|` rather than at zero
        (stalled | (abs(δ) < drts.atol)) && return x

        fx  = f(x, params)
        f′x = f′(x, params)

        # keep `x⁻` on the negative side of the root
        x⁻, x⁺ = ifelse(fx < 0, (x, x⁺), (x⁻, x))
    end

    warn_fail && @warn "Failed to converge"

    return x
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
  armijo_constant :: FT = 0.5
      min_damping :: FT = damping^10
           bounds :: BO = (lower = nothing, upper = nothing)
end

@inline function (dnrs::DampedNewtonRaphsonSolver)(f, f′, x0, params; warn_fail = false)
    x = x0
    N = 0
    U = dnrs.bounds.upper
    L = dnrs.bounds.lower
    λₘ = dnrs.min_damping
    c = dnrs.armijo_constant

    fx = f(x, params)

    while (abs(fx) > dnrs.atol) & (N < dnrs.max_iters)
        δ = fx / f′(x, params)

        λ = bounded_λ(x, L, U, δ)

        fnew = f(x - λ * δ, params)

        while ((abs(fnew) >= abs(fx) - c * λ * fx * sign(fx))) & (λ > λₘ)
            λ *= dnrs.damping
            fnew = f(x - λ * δ, params)
        end

        x -= λ * δ
        N += 1

        fx = fnew
    end

    (warn_fail && (abs(fx) > dnrs.atol)) && @warn "Failed to converge"

    return x
end

@inline bounded_λ(x, ::Nothing, ::Nothing, δ) = one(x)
@inline bounded_λ(x, ::Nothing, U, δ) = bounded_λ(x, U, nothing, δ)

@inline function bounded_λ(x, B, ::Nothing, δ)
    λb = (x - B)/δ
    λb = ifelse(λb < 0, 1, λb)

    return min(1, λb)
end

@inline function bounded_λ(x, L, U, δ)
    λl = (x - L)/δ
    λu = (x - U)/δ

    λl = ifelse(λl < 0, 1, λl)
    λu = ifelse(λu < 0, 1, λu)

    return min(1, λl, λu)
end
