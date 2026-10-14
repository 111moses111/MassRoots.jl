# The evaluator interface: how H(x) is computed. This is the extension point of the
# package; every approach to computing H is a subtype of AbstractEvaluator.

"""
    AbstractEvaluator

A way of computing `H(x)`. A new approach is a subtype with one method

    MassRoots.evaluate(ev::MyEvaluator, H::HFunction, x::Rational{BigInt},
                       atol::Rational{BigInt}) -> HValue

which returns `H(x)` with absolute error at most `atol` (guaranteed if the method is
rigorous, aimed for otherwise). The method may dispatch on the problem type
(`H::HFunction{<:FactorProblem}`) to use structure specific to that problem. The search
code sees evaluators only through [`Oracle`](@ref), which counts and times the queries.
"""
abstract type AbstractEvaluator end

"""
    HValue(value, radius; rigorous = false, info = (;))

A value of `H` returned by an evaluator: the approximation, an error bound (`radius`; a
guaranteed bound when `rigorous`, the method's own estimate otherwise, `NaN` when the
method has none) and method-specific diagnostics in `info` (number of integrand
evaluations, subintervals, working precision, ...).
"""
struct HValue
    value::BigFloat
    radius::BigFloat
    rigorous::Bool
    info::NamedTuple
end

function HValue(value, radius; rigorous::Bool = false, info::NamedTuple = (;))
    return HValue(tobigfloat(value), tobigfloat(radius), rigorous, info)
end

"""
    tobigfloat(v)

Convert a value returned by an evaluator to `BigFloat` without losing precision.
Extensions add methods (e.g. the midpoint of an Arb ball).
"""
tobigfloat(v::BigFloat) = v
tobigfloat(v::Real) = BigFloat(v; precision = max(precision(BigFloat), 256))

"""
    evaluate(ev::AbstractEvaluator, H::HFunction, x, atol) -> HValue

Compute `H(x)` to absolute accuracy `atol` with the approach `ev`.
"""
function evaluate end

evaluate(ev::AbstractEvaluator, H::HFunction, x::Real, atol::Real) =
    evaluate(ev, H, Rational{BigInt}(x), Rational{BigInt}(atol))

# Fallback that names the missing piece instead of a bare MethodError.
function evaluate(ev::AbstractEvaluator, H::HFunction, x::Rational{BigInt},
                  atol::Rational{BigInt})
    throw(ArgumentError("no method `MassRoots.evaluate(::$(typeof(ev)), ::HFunction, x, atol)`; " *
                        "implement it to use this evaluator"))
end

# ------------------------------------------------------------------------------------
# Helpers shared by quadrature-based evaluators.

"""
    default_segment_length(H)

Length of the uniform initial subdivision used by quadrature evaluators when none is
given. For a certified problem (`α = S²`) it is `2w = 1/(Sν)`: by the proof of Theorem 5
the kernel is at least 1/2 on `[d - w, d + w]` around every zero `d`, so every peak
covers at least half of a segment and cannot fall between the quadrature nodes. The
choice uses only `α`, not the location of the zeros. For other problems there is no
default subdivision (`nothing`).
"""
function default_segment_length(H::HFunction)
    H.problem isa IntegerZeroProblem || return nothing
    isone(denominator(H.α)) || return nothing
    S = isqrt(numerator(H.α))
    S * S == numerator(H.α) || return nothing
    return 1 // (S * ceilsqrt(S))
end

"Uniform subdivision of `[a, b]` (exact rationals) into pieces of length at most `h`."
function uniform_breaks(a::Rational{BigInt}, b::Rational{BigInt}, h)
    h === nothing && return [a, b]
    n = max(1, cld(numerator((b - a) / h), denominator((b - a) / h)))
    return [a + (b - a) * k // n for k in 0:n]
end
