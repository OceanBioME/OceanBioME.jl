struct PolynomialParameterisation{N, C}
    coefficients :: C
end

function PolynomialParameterisation{N}(coefficients) where N
    length(coefficients) == N + 1 || 
        throw(ArgumentError("You must provide N+1 coefficients for an order N polynomial"))

    return PolynomialParameterisation{N, typeof(coefficients)}(coefficients)
end

# generic implementation
@inline function (p::PolynomialParameterisation{N})(x) where N
    y = 0
    for n in 0:N
        y += @inbounds p.coefficients[n+1] * x^N
    end
    return y
end

# fast implementations for low order
@inline (p::PolynomialParameterisation{0})(x) = @inbounds p.coefficients[1]
@inline (p::PolynomialParameterisation{1})(x) = @inbounds p.coefficients[1] + p.coefficients[2] * x
@inline (p::PolynomialParameterisation{2})(x) = @inbounds p.coefficients[1] + p.coefficients[2] * x + p.coefficients[3] * x^2

@inline (p::PolynomialParameterisation{3})(x) = @inbounds (p.coefficients[1] + p.coefficients[2] * x + p.coefficients[3] * x^2
                                                           + p.coefficients[4] * x^3)

@inline (p::PolynomialParameterisation{4})(x) = @inbounds (p.coefficients[1] + p.coefficients[2] * x + p.coefficients[3] * x^2
                                                           + p.coefficients[4] * x^3 + p.coefficients[5] * x^4)

@inline (p::PolynomialParameterisation{5})(x) = @inbounds (p.coefficients[1] + p.coefficients[2] * x + p.coefficients[3] * x^2
                                                           + p.coefficients[4] * x^3 + p.coefficients[5] * x^4 + p.coefficients[6] * x^5)


summary(::PolynomialParameterisation{N}) where N = "Order $N polynomial parameterisation"
show(io::IO, p::PolynomialParameterisation{N}) where N = 
    println(io, summary(p), "\n",
                "    p(x) = Σ{n ∈ Z : [0, $N]}(cₙ xⁿ⁻¹) where c = $(p.coefficients)")