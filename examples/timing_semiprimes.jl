# Divisor search for semiprimes t = p·q with known prime factors of similar size:
# number of values of H and running time (rigorous Arb evaluator, erfc kernel).
#
#   julia --project=examples examples/timing_semiprimes.jl
#
# The number of queries grows like log₂ t, the time roughly like t^1.15 (see the table in
# the README). The last case, t ≈ 10⁶, takes about 45 minutes on a 2-core machine and is
# commented out.

using MassRoots, Arblib, Printf

find_divisor(21, ArbEvaluator(); kernel = ErfcKernel())          # compilation, not timed

cases = [(31, 37), (101, 103), (317, 331)]                       # add (1009, 1013) for t ≈ 10⁶

@printf("%-10s %-6s %-11s %-7s %-8s %-10s %s\n",
        "t", "digits", "factors", "found", "queries", "total, s", "longest query, s")
for (p, q) in cases
    t = p * q
    el = @elapsed r = find_divisor(t, ArbEvaluator(); kernel = ErfcKernel())
    qmax = maximum(e.seconds for e in r.oracle.log)
    @printf("%-10d %-6d %-11s %-7s %-8d %-10.2f %.2f\n",
            t, ndigits(t), "$p·$q", string(r.root), r.queries, el, qmax)
end
