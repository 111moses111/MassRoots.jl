"""
    ArbEvaluator(; precision = nothing, depth_limit = nothing, eval_limit = nothing,
                 deg_limit = nothing, max_attempts = 5, verbose = 0)

Rigorous integration with Arb (`acb_calc_integrate` through Arblib.jl): the integrand is
evaluated on complex balls and the result is a ball that contains `H(x)`. Requires
`using Arblib` (the method lives in a package extension).

The integrand must be holomorphic near the integration path; for the built-in problems
and kernels it is entire (quadratic problem) or holomorphic away from 0 (factorization).
The kernel grows quickly off the real axis, so long intervals need many subintervals. Arb
treats the tolerance as a goal, so the evaluator asks for `atol/64`; if the returned ball
is still wider than `atol`, it repeats with a tolerance `2⁻⁸` times tighter, 16 more bits
and the evaluation and depth limits multiplied by 4, up to `max_attempts` times
(`nothing` limits start from Arb's defaults).
"""
struct ArbEvaluator <: AbstractEvaluator
    precision::Union{Nothing,Int}
    depth_limit::Union{Nothing,Int}
    eval_limit::Union{Nothing,Int}
    deg_limit::Union{Nothing,Int}
    max_attempts::Int
    verbose::Int
end

ArbEvaluator(; precision = nothing, depth_limit = nothing, eval_limit = nothing,
             deg_limit = nothing, max_attempts::Integer = 5, verbose::Integer = 0) =
    ArbEvaluator(precision, depth_limit, eval_limit, deg_limit, Int(max_attempts), Int(verbose))

function evaluate(ev::ArbEvaluator, H::HFunction, x::Rational{BigInt}, atol::Rational{BigInt})
    return arb_evaluate(ev, H, x, atol)
end

# Implemented in ext/MassRootsArblibExt.jl.
function arb_evaluate(args...)
    throw(ArgumentError("ArbEvaluator needs Arblib.jl: run `using Arblib` first"))
end
