"""
    Oracle(evaluator, H, atol)

Answers queries `H(x)` with absolute accuracy `atol` through `evaluator`, caches them,
and records every query (point, value, error bound, time, evaluator diagnostics). The
search code talks only to an `Oracle`, so the number of queries reported is exactly the
number of distinct evaluations of `H`.
"""
mutable struct Oracle{E<:AbstractEvaluator,HF<:HFunction}
    evaluator::E
    H::HF
    atol::Rational{BigInt}
    cache::Dict{Rational{BigInt},HValue}
    log::Vector{NamedTuple}
end

Oracle(ev::AbstractEvaluator, H::HFunction, atol::Real) =
    Oracle(ev, H, Rational{BigInt}(atol), Dict{Rational{BigInt},HValue}(), NamedTuple[])

"""
    query(oracle, x) -> HValue

Value of `H(x)`; evaluated once per distinct `x`. A rigorous evaluator whose error bound
exceeds the requested accuracy raises an error, because the test would no longer be
certified.
"""
function query(o::Oracle, x::Real)
    xr = Rational{BigInt}(x)
    haskey(o.cache, xr) && return o.cache[xr]
    t = @elapsed v = evaluate(o.evaluator, o.H, xr, o.atol)
    if v.rigorous && !(v.radius <= o.atol)
        error("rigorous evaluator returned radius $(Float64(v.radius)) > atol = " *
              "$(Float64(o.atol)) at x = $xr; raise its precision or subdivision limits")
    end
    o.cache[xr] = v
    push!(o.log, (x = xr, value = v.value, radius = v.radius, rigorous = v.rigorous,
                  seconds = t, info = v.info))
    return v
end

"Number of distinct queries answered so far."
nqueries(o::Oracle) = length(o.cache)

"Computed mass of `[u, v]`."
mass(o::Oracle, u::Real, v::Real) = query(o, v).value - query(o, u).value
