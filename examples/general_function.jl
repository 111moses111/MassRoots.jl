# Zeros that are not integers: the example (2.6) of the paper,
# f(x) = sin(πx/2) + sin²(4πx) + 0.1 on [0, 4], with φ = f².
#
#   julia --project=examples examples/general_function.jl
#
# There is no certified gap here: α, the threshold and the resolution are chosen by hand.
# With α = 1000 the positive minima of f² near x = 0, 2, 4 (f ≈ 0.096) are suppressed
# (α f² ≈ 9.2), while every genuine zero keeps a peak of mass ~10⁻³.

using MassRoots

f(x) = sinpi(x / 2) + sinpi(4x)^2 + 0.1
prob = GeneralProblem(f, 0, 4; square = true)

ivals, oracle = localize(prob, QuadGKEvaluator(Float64; segment_length = 1 // 256);
                         α = 1000, threshold = 1e-4, resolution = 1 // 1024)

println(length(ivals), " intervals, ", nqueries(oracle), " queries")
for (a, b) in ivals
    println("  [", Float64(a), ", ", Float64(b), "]")
end
