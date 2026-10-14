# The integral function H_α(x) = ∫_{x₀}^{x} W(g(α φ(ξ))) dξ.

"""
    HFunction(problem, kernel, α; modulation = NoModulation(), x0 = basepoint(problem))

The function `H(x) = ∫_{x₀}^{x} W(g(α φ(ξ))) dξ` for a problem, a kernel `g`, a scale
`α` and a modulation `W`. An `HFunction` describes *what* to compute; an
[`AbstractEvaluator`](@ref) decides *how*.
"""
struct HFunction{P<:AbstractProblem,K,M}
    problem::P
    kernel::K
    α::Rational{BigInt}
    modulation::M
    x0::Rational{BigInt}
end

function HFunction(problem::AbstractProblem, kernel, α::Real;
                   modulation = NoModulation(), x0 = basepoint(problem))
    return HFunction(problem, kernel, Rational{BigInt}(α), modulation, Rational{BigInt}(x0))
end

"""
    HFunction(problem::IntegerZeroProblem; kernel = GaussKernel(), modulation = NoModulation())

`HFunction` with the certified scale `α` of [`mass_test_parameters`](@ref).
"""
function HFunction(problem::IntegerZeroProblem; kernel = GaussKernel(),
                   modulation = NoModulation())
    P = mass_test_parameters(problem; modulation)
    return HFunction(problem, kernel, P.α; modulation)
end

"""
    integrand(H, ξ)

The integrand `W(g(α φ(ξ)))` at `ξ`, in the number type of `ξ`.
"""
@inline function integrand(H::HFunction, ξ)
    s = lift(ξ, H.α) * phi(H.problem, ξ)
    return modulate(H.modulation, kernelvalue(H.kernel, s))
end

"""
    required_precision(H, x, atol; guard = 64)

Working precision in bits that suffices, with margin, to evaluate `H(x)` to absolute
accuracy `atol`: room for the reduction of the arguments of `φ` (`scalebits`), for the
scale `α`, for the length of the integration interval and for `atol` itself.
"""
function required_precision(H::HFunction, x, atol; guard::Integer = 64)
    L = abs(Rational{BigInt}(x) - H.x0) + 1
    bits = scalebits(H.problem) + ceillog2(max(H.α, 1 // 1)) + ceillog2(L) +
           ceillog2(1 / Rational{BigInt}(atol))
    return bits + guard
end
