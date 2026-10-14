"""
    GaussLegendreEvaluator(T = Float64; degree = 8, precision = nothing,
                           segment_length = :auto)

Composite Gauss–Legendre rule of fixed `degree` on a uniform subdivision (non-adaptive).
Its cost is known in advance (`segments × degree` integrand evaluations), which makes it
convenient for cost experiments. Not rigorous; the radius is the difference from the
rule of degree `degree - 1` on the same subdivision, a heuristic estimate.
"""
struct GaussLegendreEvaluator{T<:Number} <: AbstractEvaluator
    degree::Int
    precision::Union{Nothing,Int}
    segment_length::Any
end

function GaussLegendreEvaluator(::Type{T} = Float64; degree::Integer = 8, precision = nothing,
                                segment_length = :auto) where {T<:Number}
    degree >= 2 || throw(ArgumentError("degree must be ≥ 2"))
    return GaussLegendreEvaluator{T}(Int(degree), precision, segment_length)
end

function evaluate(ev::GaussLegendreEvaluator{T}, H::HFunction, x::Rational{BigInt},
                  atol::Rational{BigInt}) where {T}
    x == H.x0 && return HValue(0, 0; info = (evals = 0, segments = 0))
    lo, hi, sgn = x < H.x0 ? (x, H.x0, -1) : (H.x0, x, 1)
    h = ev.segment_length === :auto ? default_segment_length(H) : ev.segment_length
    breaks = uniform_breaks(lo, hi, h)
    nseg = length(breaks) - 1
    prec = something(ev.precision, required_precision(H, x, atol))
    return with_precision(T, prec) do
        xs, ws = gauss(T, ev.degree)          # nodes and weights on [-1, 1]
        xl, wl = gauss(T, ev.degree - 1)
        hi_sum = zero(liftto(T, 0, prec))
        lo_sum = zero(hi_sum)
        for k in 1:nseg
            a = liftto(T, breaks[k], prec)
            b = liftto(T, breaks[k+1], prec)
            c, r = (a + b) / 2, (b - a) / 2
            s = zero(hi_sum)
            for i in eachindex(xs)
                s += ws[i] * integrand(H, c + r * xs[i])
            end
            hi_sum += r * s
            s = zero(hi_sum)
            for i in eachindex(xl)
                s += wl[i] * integrand(H, c + r * xl[i])
            end
            lo_sum += r * s
        end
        HValue(sgn * hi_sum, abs(hi_sum - lo_sum); rigorous = false,
               info = (evals = nseg * (2ev.degree - 1), segments = nseg, precision = prec))
    end
end
