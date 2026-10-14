# Quick start: divisor search, factorization and the quadratic problem.
#
#   julia --project=examples examples/quickstart.jl
#
# (examples/Project.toml lists MassRoots, Arblib and Primes.)

using MassRoots
using Arblib                     # loads the rigorous ArbEvaluator
using Printf

fmt(x) = @sprintf("%.6e", x)     # BigFloat values can be far below Float64 range

# --- one divisor ------------------------------------------------------------------------
r = find_divisor(1001, ArbEvaluator())
println(r)                                   # LocateResult(:found, root = 7, queries = 12)
println(r.params)                            # S = 2t², α = S², p ≤ 3n + 7
for q in r.oracle.log                        # every value of H that was requested
    println("  H(", q.x, ") = ", fmt(q.value), " ± ", fmt(q.radius), "   ",
            round(q.seconds; digits = 3), " s")
end

# --- full factorization -----------------------------------------------------------------
factors, queries = factorize(2^2 * 3 * 7^2, ArbEvaluator())   # 588
println("factors = ", factors, ", queries = ", queries)

# --- quadratic feasibility a m² + b k = c ------------------------------------------------
Q = QuadraticProblem(3, 7, 103)              # 3·m² + 7·k = 103
ans, d = decide(Q, ArbEvaluator())           # exactly two queries
println("solvable: ", ans, " (", d.queries, " queries)")
s = locate(Q, ArbEvaluator())
println("solution m = ", s.root, ", k = ", (103 - 3 * s.root^2) ÷ 7)
println("3·m² + 7·k = 100 solvable: ", first(decide(QuadraticProblem(3, 7, 100), ArbEvaluator())))

# --- H itself ---------------------------------------------------------------------------
H = HFunction(FactorProblem(91))             # Gaussian kernel, certified α
atol = mass_test_parameters(FactorProblem(91)).atol
v = evaluate(ArbEvaluator(), H, 89 // 2, atol)
println("H(89/2) = ", fmt(v.value), " ± ", fmt(v.radius))
