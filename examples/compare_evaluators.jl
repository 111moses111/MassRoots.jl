# Compare approaches to computing H on divisor search: correctness, queries, time and
# integrand evaluations as t grows. Writes examples/results.csv.
#
#   julia --project=examples examples/compare_evaluators.jl
#
# Each evaluator has a cap on t, because the cost of computing H grows with t (the
# number of queries is O(log t), but each query integrates over an interval of length ~t
# with peaks of width ~t⁻² … t⁻³).

using MassRoots
using Arblib
using Primes

evaluators = [
    # name                       evaluator                                max t
    ("QuadGK Float64",           QuadGKEvaluator(Float64),                25),
    ("Gauss-Legendre Float64",   GaussLegendreEvaluator(Float64),         25),
    ("Arb (rigorous)",           ArbEvaluator(),                          3000),
    ("Arb, erfc kernel",         ArbEvaluator(),                          3000),
]
kernels = Dict("Arb, erfc kernel" => ErfcKernel())

ts = [15, 21, 23, 91, 97, 143, 323, 1001, 2021, 2999]

rows = String["t,evaluator,status,root,correct,queries,seconds,evals,max_precision"]
for t in ts, (name, ev, maxt) in evaluators
    t <= maxt || continue
    kernel = get(kernels, name, GaussKernel())
    stats = @timed find_divisor(t, ev; kernel)
    r = stats.value
    correct = isprime(t) ? r.status === :none : (r.status === :found && t % r.root == 0)
    evals = sum((get(q.info, :evals, 0) for q in r.oracle.log); init = 0)
    prec = maximum((get(q.info, :precision, 0) for q in r.oracle.log); init = 0)
    push!(rows, join((t, name, r.status, something(r.root, ""), correct, r.queries,
                      round(stats.time; digits = 3), evals, prec), ","))
    println(rpad("t = $t", 10), rpad(name, 26), rpad(string(r.status), 12),
            "queries = ", lpad(r.queries, 2), "   ", lpad(round(stats.time; digits = 2), 8), " s",
            correct ? "" : "   WRONG")
end

out = joinpath(@__DIR__, "results.csv")
write(out, join(rows, "\n") * "\n")
println("written ", out)
