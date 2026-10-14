# Search procedures. They see H only through an Oracle.

"""
    LocateResult

Outcome of [`locate`](@ref) or [`decide`](@ref):

* `status` – `:found` (a verified zero), `:none` (the initial block tested negative),
  `:unverified` (bisection ended at a candidate that fails the exact check) or
  `:inconsistent` (a positive block with two negative children; impossible with an
  oracle that meets its accuracy, so it signals an evaluator that missed it);
* `root` – the zero found, or `nothing`;
* `queries` – number of distinct evaluations of `H`;
* `params`, `oracle` – the certified parameters and the oracle with its query log
  (`nothing` when no query was needed).
"""
struct LocateResult{O}
    status::Symbol
    root::Union{Nothing,BigInt}
    queries::Int
    params::MassTestParameters
    oracle::O
end

function Base.show(io::IO, r::LocateResult)
    print(io, "LocateResult(:", r.status, ", root = ", something(r.root, "nothing"),
          ", queries = ", r.queries, ")")
end

positive(o::Oracle, ℓ, r, θ) = mass(o, half(ℓ - 1), half(r)) > θ

"""
    locate(prob::IntegerZeroProblem, evaluator; kernel = GaussKernel(),
           modulation = NoModulation(), params = mass_test_parameters(prob; modulation))

Find one zero of `prob` by bisection over candidate blocks `[ℓ - 1/2, r + 1/2]`
(§5.4 and Theorem 6 of the paper). Uses at most `2 + ⌈log₂ M⌉` queries for `M`
candidates; the zero is checked exactly with `isroot`.
"""
function locate(prob::IntegerZeroProblem, ev::AbstractEvaluator; kernel = GaussKernel(),
                modulation = NoModulation(),
                params::MassTestParameters = mass_test_parameters(prob; modulation))
    H = HFunction(prob, kernel, params.α; modulation)
    o = Oracle(ev, H, params.atol)
    θ = params.threshold
    c = candidates(prob)
    ℓ, r = first(c), last(c)
    positive(o, ℓ, r, θ) || return LocateResult(:none, nothing, nqueries(o), params, o)
    while ℓ < r
        j = fld(ℓ + r, 2)
        if positive(o, ℓ, j, θ)
            r = j
        elseif positive(o, j + 1, r, θ)          # no new query: both ends are cached
            ℓ = j + 1
        else
            return LocateResult(:inconsistent, nothing, nqueries(o), params, o)
        end
    end
    status = isroot(prob, ℓ) ? :found : :unverified
    return LocateResult(status, ℓ, nqueries(o), params, o)
end

"""
    decide(prob::IntegerZeroProblem, evaluator; kwargs...) -> (Bool, LocateResult)

Does `prob` have a zero? One mass test on the whole candidate range: exactly two
nonadaptive queries, `H(first - 1/2)` and `H(last + 1/2)` (Theorem 6.2).
"""
function decide(prob::IntegerZeroProblem, ev::AbstractEvaluator; kernel = GaussKernel(),
                modulation = NoModulation(),
                params::MassTestParameters = mass_test_parameters(prob; modulation))
    H = HFunction(prob, kernel, params.α; modulation)
    o = Oracle(ev, H, params.atol)
    c = candidates(prob)
    ans = positive(o, first(c), last(c), params.threshold)
    return ans, LocateResult(ans ? :positive : :none, nothing, nqueries(o), params, o)
end

"""
    find_divisor(t, evaluator; kwargs...) -> LocateResult

A nontrivial divisor of `t ≥ 2`, or `:none` if `t` is prime. `t = 2, 3` are answered
without queries.
"""
function find_divisor(t::Integer, ev::AbstractEvaluator; kwargs...)
    t >= 2 || throw(ArgumentError("t must be ≥ 2"))
    if t < 4
        return LocateResult(:none, nothing, 0, mass_test_parameters(8), nothing)
    end
    return locate(FactorProblem(t), ev; kwargs...)
end

"""
    factorize(t, evaluator; kwargs...) -> (factors, queries)

Prime factorization of `t ≥ 2` by repeated [`find_divisor`](@ref); `factors` is sorted
with multiplicity, `queries` is the total number of evaluations of `H`.
"""
function factorize(t::Integer, ev::AbstractEvaluator; kwargs...)
    t >= 2 || throw(ArgumentError("t must be ≥ 2"))
    stack = BigInt[big(t)]
    primes = BigInt[]
    total = 0
    while !isempty(stack)
        n = pop!(stack)
        r = find_divisor(n, ev; kwargs...)
        total += r.queries
        if r.status === :none
            push!(primes, n)
        elseif r.status === :found
            push!(stack, r.root, div(n, r.root))
        else
            error("factorization of $n stopped with status :$(r.status); " *
                  "the evaluator did not reach the required accuracy")
        end
    end
    return sort!(primes), total
end

"""
    localize(prob::AbstractProblem, evaluator; α, threshold, resolution,
             kernel = GaussKernel(), modulation = NoModulation(), atol = threshold / 8,
             interval = domain(prob), firstonly = false)

Bisection with mass tests for any problem, including zeros that are not integers.
Returns the intervals of width at most `resolution` whose computed mass exceeds
`threshold` (with `merge = true`, adjacent ones joined, so a peak wider than the
resolution gives one interval), and the oracle. No guarantee: `α`, `threshold` and `resolution` are the
caller's choice, and a zero-free interval may pass the test near a small positive
minimum of `φ`, while a zero may be missed if its peak carries less mass than the
threshold. For problems with integer zeros and a certified scale use [`locate`](@ref).
"""
function localize(prob::AbstractProblem, ev::AbstractEvaluator; α::Real, threshold::Real,
                  resolution::Real, kernel = GaussKernel(), modulation = NoModulation(),
                  atol::Real = Rational{BigInt}(threshold) / 8,
                  interval = domain(prob), firstonly::Bool = false, merge::Bool = true)
    H = HFunction(prob, kernel, α; modulation)
    o = Oracle(ev, H, atol)
    θ = Rational{BigInt}(threshold)
    δ = Rational{BigInt}(resolution)
    found = Tuple{Rational{BigInt},Rational{BigInt}}[]
    work = [(Rational{BigInt}(interval[1]), Rational{BigInt}(interval[2]))]
    while !isempty(work)
        a, b = popfirst!(work)
        mass(o, a, b) > θ || continue
        if b - a <= δ
            push!(found, (a, b))
            firstonly && break
        else
            m = (a + b) / 2
            pushfirst!(work, (a, m), (m, b))
        end
    end
    sort!(found)
    merge && (found = merge_adjacent(found))
    return found, o
end

"Join intervals that share an endpoint: one cluster per peak."
function merge_adjacent(iv::Vector{Tuple{Rational{BigInt},Rational{BigInt}}})
    out = similar(iv, 0)
    for (a, b) in iv
        if !isempty(out) && out[end][2] == a
            out[end] = (out[end][1], b)
        else
            push!(out, (a, b))
        end
    end
    return out
end
