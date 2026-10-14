# Number-type plumbing. All formulas in the package are written for a generic number
# type `T` (Float64, BigFloat, Double64, Arb, Acb, ...). Exact constants (BigInt and
# Rational{BigInt}) are converted to the type of the current argument with `lift`, which
# keeps the working precision of that argument. Extensions add methods for their types.

"""
    lift(x, v)

Convert the exact constant `v` (an `Integer` or `Rational`) to the number type of `x`,
using the precision of `x` where the type has one.
"""
lift(x::T, v) where {T<:Number} = convert(T, v)
lift(x::BigFloat, v) = BigFloat(v; precision = precision(x))

"""
    liftto(T, v, prec)

Convert the exact value `v` to the number type `T` at `prec` bits (ignored by
fixed-precision types).
"""
liftto(::Type{T}, v, prec::Integer) where {T<:Number} = convert(T, v)
liftto(::Type{BigFloat}, v, prec::Integer) = BigFloat(v; precision = prec)

"""
    sinpi2(y)

`sin(π y)^2`, computed with `sinpi` so that the reduction of a large argument is done by
the number type itself.
"""
@inline function sinpi2(y)
    s = sinpi(y)
    return s * s
end

"""
    with_precision(f, T, prec)

Run `f()` with working precision `prec` bits for types that have a global precision
setting (`BigFloat`); fixed-precision types run `f()` unchanged.
"""
with_precision(f, ::Type{BigFloat}, prec::Integer) = setprecision(f, BigFloat, prec)
with_precision(f, ::Type, prec::Integer) = f()

"Smallest integer `p ≥ 0` with `2^p ≥ r`."
function ceillog2(r::Rational)
    r <= 1 && return 0
    m = cld(numerator(r), denominator(r))
    return ndigits(m - 1; base = 2)
end
ceillog2(m::Integer) = ceillog2(m // 1)

"Smallest integer `ν` with `ν^2 ≥ S`."
function ceilsqrt(S::Integer)
    ν = isqrt(S)
    return ν * ν < S ? ν + 1 : ν
end

"Number of binary digits of a positive integer."
bitlength(n::Integer) = n == 0 ? 1 : ndigits(abs(n); base = 2)

"Half-integer `k + 1/2` as an exact rational."
half(k::Integer) = (2 * big(k) + 1) // 2
