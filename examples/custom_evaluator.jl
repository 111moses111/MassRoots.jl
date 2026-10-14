# How to add a new way of computing H: a composite midpoint rule in Float64.
#
#   julia --project=examples examples/custom_evaluator.jl
#
# An evaluator is any subtype of AbstractEvaluator with one `evaluate` method. The search,
# the query counting and the logging come for free.

using MassRoots
import MassRoots: evaluate

struct MidpointEvaluator <: AbstractEvaluator
    points_per_segment::Int
end

function evaluate(ev::MidpointEvaluator, H::HFunction, x::Rational{BigInt},
                  atol::Rational{BigInt})
    x == H.x0 && return HValue(0, 0)
    lo, hi, sgn = x < H.x0 ? (x, H.x0, -1) : (H.x0, x, 1)
    # Same uniform subdivision as the built-in quadratures (2w for certified problems).
    breaks = MassRoots.uniform_breaks(lo, hi, MassRoots.default_segment_length(H))
    total = 0.0
    evals = 0
    for k in 1:length(breaks)-1
        a, b = Float64(breaks[k]), Float64(breaks[k+1])
        h = (b - a) / ev.points_per_segment
        for i in 0:ev.points_per_segment-1
            total += h * integrand(H, a + (i + 0.5) * h)
            evals += 1
        end
    end
    return HValue(sgn * total, NaN; info = (evals = evals,))   # NaN: no error estimate
end

for t in (15, 21, 13)
    r = find_divisor(t, MidpointEvaluator(4))
    evals = sum(q.info.evals for q in r.oracle.log)
    println("t = $t: ", r, ", integrand evaluations: ", evals)
end

# An evaluator can also specialise on the problem, e.g. to exploit structure of f_t:
#   evaluate(ev::MyEvaluator, H::HFunction{FactorProblem}, x, atol) = ...
