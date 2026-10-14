"""
    QuadGKEvaluator(T = Float64; precision = nothing, order = 7, maxevals = 10^8,
                    segment_length = :auto)

Adaptive Gauss–Kronrod quadrature (QuadGK.jl) in the number type `T` (`Float64`,
`BigFloat`, `Double64`, ...). Not rigorous: the returned radius is QuadGK's error
estimate.

* `precision` – working precision in bits for `BigFloat`; by default
  [`required_precision`](@ref).
* `segment_length` – the integration interval is first split uniformly into pieces of
  this length, each integrated adaptively. `:auto` uses [`default_segment_length`](@ref);
  `nothing` disables the split. Without a split, an adaptive rule can step over a peak
  much narrower than the interval and return a value near zero.
"""
struct QuadGKEvaluator{T<:Number} <: AbstractEvaluator
    precision::Union{Nothing,Int}
    order::Int
    maxevals::Int
    segment_length::Any
end

function QuadGKEvaluator(::Type{T} = Float64; precision = nothing, order::Integer = 7,
                         maxevals::Integer = 10^8, segment_length = :auto) where {T<:Number}
    return QuadGKEvaluator{T}(precision, Int(order), Int(maxevals), segment_length)
end

function evaluate(ev::QuadGKEvaluator{T}, H::HFunction, x::Rational{BigInt},
                  atol::Rational{BigInt}) where {T}
    x == H.x0 && return HValue(0, 0; info = (evals = 0, segments = 0))
    lo, hi, sgn = x < H.x0 ? (x, H.x0, -1) : (H.x0, x, 1)
    h = ev.segment_length === :auto ? default_segment_length(H) : ev.segment_length
    breaks = uniform_breaks(lo, hi, h)
    nseg = length(breaks) - 1
    prec = something(ev.precision, required_precision(H, x, atol))
    return with_precision(T, prec) do
        f = ξ -> integrand(H, ξ)
        segatol = liftto(T, atol / (2 * nseg), prec)
        total = zero(liftto(T, 0, prec))
        err = zero(total)
        evals = 0
        for k in 1:nseg
            a = liftto(T, breaks[k], prec)
            b = liftto(T, breaks[k+1], prec)
            I, E, n = quadgk_count(f, a, b; atol = segatol, rtol = zero(segatol),
                                   order = ev.order, maxevals = ev.maxevals)
            total += I
            err += E
            evals += n
        end
        HValue(sgn * total, err; rigorous = false,
               info = (evals = evals, segments = nseg, precision = prec))
    end
end
