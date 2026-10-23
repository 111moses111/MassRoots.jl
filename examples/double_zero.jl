# A zero without a sign change.
#
#   julia --project=examples examples/double_zero.jl
#
# f(x) = (x - 1)²(x - 2) has a double zero at x = 1, where f does not change sign, and a
# simple zero at x = 2. Methods that bracket a sign change cannot see x = 1; the integral
# mass test finds both zeros. There is no certified gap for a general f: α, the threshold
# and the resolution are chosen by hand.

using MassRoots

f(x) = (x - 1)^2 * (x - 2)
prob = GeneralProblem(f, 0, 3)          # φ = f² (the default)

println("f(0.9) = ", f(0.9), ",  f(1.1) = ", f(1.1), "   (no sign change at x = 1)")

for α in (100, 1000)
    ivals, oracle = localize(prob, QuadGKEvaluator(Float64; segment_length = 1 // 64);
                             α = α, threshold = 1e-3, resolution = 1 // 64)
    println("α = ", α, ": ", length(ivals), " intervals, ", nqueries(oracle), " queries")
    for (a, b) in ivals
        println("    [", Float64(a), ", ", Float64(b), "]")
    end
end
