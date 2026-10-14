"""
    IdealEvaluator(; noise = :random, seed = 0)

**Test oracle that knows the zeros.** It finds the zeros of an
[`IntegerZeroProblem`](@ref) by checking every candidate with `isroot`, and returns a
value of `H` that is consistent with the guarantees of Theorem 5 and otherwise
unfavourable:

* every zero contributes exactly the guaranteed minimum mass `B₁`;
* each zero-free stretch carries background mass up to `B₀`, growing linearly across it;
* with `noise = :random` the value is shifted by `±atol` (sign from a hash of `x` and
  `seed`); `noise = :none` returns the model value.

It exists to test the search logic and the query counts independently of any
integration method. It is brute force in disguise and must never be used as a solver.
"""
struct IdealEvaluator <: AbstractEvaluator
    noise::Symbol
    seed::UInt
    function IdealEvaluator(; noise::Symbol = :random, seed::Integer = 0)
        noise in (:random, :none) || throw(ArgumentError("noise must be :random or :none"))
        return new(noise, UInt(seed))
    end
end

function evaluate(ev::IdealEvaluator, H::HFunction{<:IntegerZeroProblem},
                  x::Rational{BigInt}, atol::Rational{BigInt})
    prob = H.problem
    cands = candidates(prob)
    Q = isqrt(numerator(H.α))
    ν = ceilsqrt(Q)
    V = modulate_exact(H.modulation, big(1) // 2)
    B1 = V / (Q * ν)
    B0 = V / (2Q * ν)
    roots = [m for m in cands if isroot(prob, m)]
    model(y) = ideal_model(y, roots, first(cands), last(cands), B0, B1)
    v = model(x) - model(H.x0)
    if ev.noise === :random
        s = isodd(hash((x, ev.seed))) ? 1 : -1
        v += s * atol
    end
    return HValue(BigFloat(v; precision = 64 + ceillog2(1 / atol) + 64), atol;
                  rigorous = false, info = (model = :ideal,))
end

# Monotone model of H on the real line: a jump of B1 at every zero, plus background that
# grows linearly from 0 to B0 across each zero-free stretch between consecutive zeros
# (and between the ends of the candidate range and the first/last zero).
function ideal_model(y::Rational{BigInt}, roots, lo, hi, B0, B1)
    edges = vcat(big(lo) - 1 // 2, [big(r) // 1 for r in roots], big(hi) + 1 // 2)
    total = zero(B1)
    for k in 1:length(edges)-1
        a, b = edges[k], edges[k+1]
        frac = y <= a ? 0 // 1 : (y >= b ? 1 // 1 : (y - a) / (b - a))
        total += B0 * frac
    end
    total += B1 * count(r -> r < y, roots)
    return total
end
